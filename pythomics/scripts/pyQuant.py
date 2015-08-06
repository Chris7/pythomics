#!/usr/bin/env python
from __future__ import division

# TODO: Threading for single files. Since much time is spent in fetching m/z records, it may be pointless
# because I/O is limiting. If so, create a thread just for I/O for processing requests the other threads
# interface with

description = """
This will quantify labeled peaks (such as SILAC) in ms1 spectra. It relies solely on the distance between peaks,
 which can correct for errors due to amino acid conversions.
"""

import sys
import os
import copy
import operator
from scipy.stats import mode
import traceback
import pandas as pd
import numpy as np
import re
from collections import OrderedDict, defaultdict
from itertools import izip_longest
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Process, Queue, Manager
try:
    from profilestats import profile
except ImportError:
    pass
from Queue import Empty
from scipy import integrate, optimize
from scipy.stats import signaltonoise, linregress

from matplotlib.backends.backend_pdf import PdfPages
from cStringIO import StringIO
import argparse
import urllib
import base64
import scipy.stats

from pythomics.templates import CustomParser
from pythomics.proteomics.parsers import GuessIterator
from pythomics.proteomics import peaks, config

RESULT_ORDER = [('peptide', 'Peptide'), ('modifications', 'Modifications'),
                ('charge', 'Charge'), ('ms1', 'MS1 Spectrum ID'), ('scan', 'MS2 Spectrum ID')]


parser = CustomParser(description=description)
parser.add_processed_ms()
parser.add_argument('--mzml', help="The mzML files for the raw data or directory containing them. Defaults to the directory of the processed file.", type=argparse.FileType('r'), nargs='*')
parser.add_argument('--filemap', help="By default, the mzML file is assumed to be a derivation of the scan filename. " \
                                     "If this is not true, this file provides a correct mapping.")
parser.add_argument('--rt-width', help="The width of the retention time window to search for. Default: 1.0 minute", type=float, default=1.0)
parser.add_argument('--precision', help="The precision for storing m/z values. Defaults to 6 decimal places.", type=int, default=6)
parser.add_argument('--no-spread', help="Assume there is no spread of the label and only compare pairwise peaks.", action='store_true')
parser.add_argument('--no-sparse', help="Don't store data as sparse entities. Major memory savings at the expense of speed.", action='store_false')
parser.add_argument('--no-temp', help="This will disable saving of the intermediate files for SILAC analysis.", action='store_false')
parser.add_argument('--debug', help="This will output debug information and graphs.", action='store_true')
parser.add_argument('--skip', help="If true, skip scans with missing files in the mapping.", action='store_true')
parser.add_argument('--peptide', help="The peptide to limit quantification to.", type=str)
parser.add_argument('--html', help="Output a HTML table summary.", action='store_true')
parser.add_argument('--resume', help="Will resume from the last run. Only works if not directing output to stdout.", action='store_true')
parser.add_argument('--sample', help="How much of the data to sample. Enter as a decimal (ie 1.0 for everything, 0.1 for 10%%)", type=float, default=1.0)
parser.add_argument('-o', '--out', nargs='?', help='Stuff', type=str)

import random

class Reader(Process):
    def __init__(self, incoming, outgoing, scan_dict=None, raw_file=None, ms2_ms1_map=None, scan_rt_map=None, rt_index_map=None):
        super(Reader, self).__init__()
        self.incoming = incoming
        self.outgoing = outgoing
        self.ms2_ms1_map = ms2_ms1_map
        self.scan_rt_map = scan_rt_map
        self.rt_index_map = rt_index_map
        self.scan_dict = scan_dict
        self.raw = GuessIterator(raw_file, full=False, store=False, ms_filter=1)
        sys.stderr.write('Thread processing peptide file.\n')
        for index, i in enumerate(self.raw):
            if index % 100 == 0:
                sys.stderr.write('.')
            if i is None:
                break
            if i.ms_level != 1:
                self.ms2_ms1_map[int(i.id)] = scan_id
                continue
            scan_id = int(i.id)
            rt = i.rt
            self.scan_rt_map[scan_id] = rt
            self.rt_index_map[rt] = scan_id
        sys.stderr.write('Scans Loaded\n')

    def run(self):
        for scan_id in iter(self.incoming.get, None):
            if scan_id in self.scan_dict:
                continue
            scan = self.raw.getScan(scan_id)
            self.outgoing.put(scan)
        sys.stderr.write('reader done\n')


class Worker(Process):
    def __init__(self, queue=None, results=None, precision=6, raw_name=None, silac_labels=None,
                 temp=False, sparse=False, debug=False, html=False, mono=False,
                 scan_dict=None, reader_in=None, reader_out=None, ms2_ms1_map=None, scan_rt_map=None, rt_index_map=None):
        super(Worker, self).__init__()
        self.precision = precision
        self.sparse = sparse
        self.queue=queue
        self.reader_in, self.reader_out = reader_in, reader_out
        self.scan_dict, self.ms2_ms1_map, self.scan_rt_map, self.rt_index_map = scan_dict, ms2_ms1_map, scan_rt_map, rt_index_map
        self.results = results
        self.silac_labels = {} if silac_labels is None else silac_labels
        self.shifts = {0: "Light"}
        self.shifts.update({sum(silac_masses.keys()): silac_label for silac_label, silac_masses in self.silac_labels.iteritems()})
        self.raw_name = raw_name
        self.filename = os.path.split(self.raw_name)[1]
        self.temp = temp
        self.rt_width = 1
        self.rt_tol = 0.2 # for fitting
        self.debug = debug
        self.html = html
        self.mono = mono
        self.last_clean = 0 # for cleaning up scans before our current RT window
        sys.stderr.write('Thread spawned: Sparse: {0}\n'.format(self.sparse))

    def convertScan(self, scan):
        import numpy as np
        scan_vals = np.array(scan.scans)
        res = pd.Series(scan_vals[:, 1].astype(np.uint64), index=np.round(scan_vals[:, 0], self.precision), name=scan.rt, dtype='uint64')
        # due to precision, we have multiple m/z values at the same place. We can eliminate this by grouping them and summing them.
        # Summation is the correct choice here because we are combining values of a precision higher than we care about.
        return res.groupby(level=0).sum()

    def getScan(self, ms1):
        ms1 = int(ms1)
        try:
            scan = self.reader_out.get(timeout=0.1)
            while scan:
                self.scan_dict[int(scan.title)] = self.convertScan(scan)
                scan = self.reader_out.get(timeout=0.1)
        except Empty:
            pass
        if ms1 not in self.scan_dict:
            self.reader_in.put(ms1)
            try:
                scan = self.reader_out.get(timeout=0.1)
            except Empty:
                return self.getScan(ms1)
            while int(scan.title) != ms1:
                self.scan_dict[int(scan.title)] = self.convertScan(scan)
                try:
                    scan = self.reader_out.get(timeout=0.1)
                except Empty:
                    return self.getScan(ms1)
            self.scan_dict[int(scan.title)] = self.convertScan(scan)
            # put scans +- 50 into the queue just to keep it going
            scans = sorted(self.rt_index_map.values())
            current_index = scans.index(ms1)
            left, right = current_index-25, current_index+25
            if left < 0:
                left = 0
            if right > len(scans):
                right = len(scans)
            for i in xrange(left, right):
                self.reader_in.put(scans[i])
            # clean up scan_dict to kick out old entries to make serialization not so burdensome
            stored_scans = sorted(self.scan_dict.keys())
            far_left = left - 50
            if far_left > 0:
                to_remove = [i for i in stored_scans if i < scans[far_left]]
                for i in to_remove:
                    del self.scan_dict[i]
        return self.scan_dict[ms1]

    def run(self):
        import time
        scan_rt_map = dict(self.scan_rt_map.items())
        rt_index_map = dict(self.rt_index_map.items())
        ms2_ms1_map = dict(self.ms2_ms1_map.items())

        rt_index_map = pd.Series(dict(rt_index_map.items()))

        silac_shifts = {}
        for silac_label, silac_masses in self.silac_labels.items():
            for mass, aas in silac_masses.iteritems():
                pmass = float('{0:0.5f}'.format(float(mass)))
                try:
                    silac_shifts[mass] |= aas
                except:
                    silac_shifts[mass] = aas

        # @profile(print_stats=10, dump_stats=True)
        def run_thing(params):
            try:
                import operator
                from pythomics.proteomics import peaks, config
                html_images = {}
                start_time = time.time()
                scan_info = params.get('scan_info')
                scan = params.get('scan')
                scanId = scan_info['id']
                if not scan:
                    return
                try:
                    ms1 = int(scan.ms1_scan.title)
                except AttributeError:
                    ms1 = ms2_ms1_map.get(int(scan.rawId))
                rt = scan_rt_map.get(ms1)
                if not rt:
                    return
                mods = scan_info.get('modifications')
                charge = int(scan.charge)
                precursor = float(scan.mass)+config.HYDROGEN*(charge-1)
                shift = 0
                for mod in filter(lambda x: x, mods.split('|')):
                    aa, pos, mass, type = mod.split(',', 3)
                    mass = float('{0:0.5f}'.format(float(mass)))
                    if aa in silac_shifts.get(mass, {}):
                        shift += mass
                if shift > 0:
                    precursor -= shift
                peptide = scan.peptide
                if self.debug:
                    sys.stderr.write('on ms {0} {1} {2} {3}\n'.format(ms1, rt, precursor, scan_info))
                precursors = {'Light': 0.0}
                silac_dict = {'data': None, 'df': pd.DataFrame(), 'precursor': 'NA',
                              'isotopes': {}, 'peaks': OrderedDict(), 'intensity': 'NA'}
                data = OrderedDict()
                data['Light'] = copy.deepcopy(silac_dict)
                combined_data = pd.DataFrame()
                for silac_label, silac_masses in self.silac_labels.items():
                    silac_shift=0
                    for label_mass, label_masses in silac_masses.items():
                        labels = [label_mass for mod_aa in peptide if mod_aa in label_masses]
                        silac_shift += sum(labels)
                    precursors[silac_label] = silac_shift
                    data[silac_label] = copy.deepcopy(silac_dict)
                precursors = OrderedDict(sorted(precursors.items(), key=operator.itemgetter(1)))
                shift_maxes = {i: j for i,j in zip(precursors.keys(), precursors.values()[1:])}
                base_rt = rt_index_map.index.searchsorted(rt)
                finished = set([])
                finished_isotopes = {i: set([]) for i in precursors.keys()}
                result_dict = {'peptide': scan_info.get('mod_peptide', peptide),
                               'scan': scanId, 'ms1': ms1,
                               'charge': charge, 'modifications': mods}
                rt_window = []
                ms_index = 0
                delta = -1
                theo_dist = peaks.calculate_theoretical_distribution(peptide.upper())
                spacing = config.NEUTRON/float(charge)
                isotope_labels = {}
                chosen_window = defaultdict(dict)
                while True:
                    if len(finished) == len(precursors.keys()):
                        break
                    rt_window.append(rt_index_map.index[base_rt+ms_index])
                    df = self.getScan(rt_index_map.iloc[base_rt+ms_index])
                    found = False
                    for precursor_label, precursor_shift in precursors.items():
                        if precursor_label in finished:
                            continue
                        precursor_mass = precursor+precursor_shift
                        precursor_mz = precursor_mass/float(charge)
                        data[precursor_label]['precursor'] = precursor_mz
                        shift_max = shift_maxes.get(precursor_label)
                        shift_max = precursor+shift_max if shift_max is not None else None
                        envelope = peaks.findEnvelope(df, start_mz=precursor_mass, max_mz=shift_max,
                                                      charge=charge, ppm=5, ppm2=2, heavy=True,
                                                      theo_dist=theo_dist, label=precursor_label, skip_isotopes=finished_isotopes[precursor_label])
                        peaks_found = data[precursor_label]['peaks']

                        if not envelope['envelope']:
                            finished.add(precursor_label)
                            continue
                        selected = {}
                        for isotope, vals in envelope['micro_envelopes'].iteritems():
                            if isotope in finished_isotopes[precursor_label]:
                                continue
                            if vals.get('int') == 0:
                                finished_isotopes[precursor_label].add(isotope)
                                continue
                            else:
                                found = True
                            selected[precursor_mz+isotope*spacing] = vals.get('int')#df.iloc[range(*vals)].sum()
                            chosen_window[df.name][precursor_mz+isotope*spacing] = vals
                        selected = pd.Series(selected, name=df.name).to_frame()
                        for i in ['envelope', 'micro_envelopes', 'ppms']:
                            for isotope, vals in envelope.get(i, {}).iteritems():
                                if isotope in finished_isotopes[precursor_label]:
                                    continue
                                val_dict = {'info': vals, 'df': df}
                                try:
                                    peaks_found[isotope][i].append(val_dict)
                                except KeyError:
                                    try:
                                        peaks_found[isotope][i] = [val_dict]
                                    except KeyError:
                                        peaks_found[isotope] = {i: [val_dict]}

                        for i, isotope_index in zip(selected.index, envelope['envelope'].keys()):
                            isotope_labels[i] = {'label': precursor_label, 'isotope_index': isotope_index}
                        if df.name in combined_data.columns:
                            combined_data = combined_data.add(selected, axis='index', fill_value=0)
                        else:
                            combined_data = pd.concat([combined_data, selected], axis=1).fillna(0)

                    if found is False:
                        if delta == -1:
                            delta = 1
                            ms_index=0
                            finished = set([])
                            finished_isotopes = {i: set([]) for i in precursors.keys()}
                        elif ms_index >= 2:
                            break
                    ms_index += delta
                # bookmark with zeros, do the right end first because pandas will by default append there
                combined_data[rt_index_map.index[rt_index_map.index.searchsorted(combined_data.columns[-1])+1]] = 0
                combined_data[rt_index_map.index[rt_index_map.index.searchsorted(combined_data.columns[0])-1]] = 0
                combined_data = combined_data[sorted(combined_data.columns)]
                # shared_isotopes = set([])
                rt_window.sort()
                combined_data = combined_data.sort(axis='index').sort(axis='columns')
                from scipy.signal import argrelmax, argrelmin
                from scipy import integrate
                from scipy.optimize import minimize
                start_rt = rt_index_map.index[base_rt]
                quant_vals = defaultdict(dict)
                isotope_labels = pd.DataFrame(isotope_labels).T
                fig_map = {}

                if self.html:
                    fname = '{2}_{0}_{1}_{3}_clusters.png'.format(peptide, ms1, self.filename, scanId)
                    fig = plt.figure(figsize=(10,10))
                    subplot_rows = len(precursors.keys())+1
                    subplot_columns = pd.Series(isotope_labels['label']).value_counts().iloc[0]+1
                    subplot_count = int(len(combined_data)/4+1)
                    ax = fig.add_subplot(subplot_rows, subplot_columns, 1, projection='3d')
                    X=combined_data.columns.astype(float).values
                    Y=combined_data.index.astype(float).values
                    Z=combined_data.fillna(0).values
                    # Z/=Z.max()
                    Xi,Yi = np.meshgrid(X, Y)
                    ax.plot_wireframe(Yi, Xi, Z, cmap=plt.cm.coolwarm)

                combined_peaks = defaultdict(dict)
                plot_index = {}
                quan_start = None
                fig_nums = defaultdict(list)
                for mz, label in isotope_labels['label'].iteritems():
                    fig_nums[label].append(mz)
                labelx = False
                labely = False
                ymax = combined_data.max().max()
                label_fig_row = {v: i+1 for i,v in enumerate(precursors.keys())}
                for row_num, (index, values) in enumerate(combined_data.iterrows()):
                    quant_label = isotope_labels.loc[index, 'label']
                    fig_index = label_fig_row.get(quant_label)*subplot_columns+fig_nums[quant_label].index(index)+2
                    current_row = int(fig_index/subplot_columns+1)
                    if (fig_index-2)%subplot_columns == 0:
                        labely = True
                    if current_row == subplot_rows:
                        labelx = True
                    # print quant_label, subplot_rows, subplot_columns, fig_index
                    fig_map[index] = fig_index
                    plot_index[index, quant_label] = row_num
                    res, all_peaks = peaks.findAllPeaks(values)
                    xdata = values.index.values.astype(float)
                    ydata = values.fillna(0).values
                    mval = ydata.max()
                    rt_means = res.x[1::3]
                    rt_amps = res.x[::3]
                    rt_vars = res.x[2::3]
                    combined_peaks[quant_label][index] = [{'mean': i, 'amp': j*mval, 'var': k, 'peak': l, 'total': values.sum()}
                                                          for i,j,k,l in zip(rt_means, rt_amps, rt_vars, all_peaks)]
                    if self.html:
                        ax = fig.add_subplot(subplot_rows, subplot_columns, fig_index)
                        ax.plot(xdata, ydata, 'bo-', alpha=0.7)
                        ax.plot(xdata, peaks.gauss_ndim(xdata, *res.x)*mval, color='r')
                        if labely:
                            labely = False
                        else:
                            ax.set_yticklabels([])
                        if labelx:
                            labelx = False
                            ax.set_xticklabels(['{0:.2f}'.format(i) for i in xdata], rotation=45)
                        else:
                            ax.set_xticklabels([])
                        ax.set_ylim(0,combined_data.max().max())
                # get two most common peak, pick the closest to our RT
                # we may need to add a check for a minimal # of in for max distance from the RT as well here.
                common_peaks = pd.Series([peak['peak'] for i, values in combined_peaks.items() for index, value_peaks in values.iteritems() for peak in value_peaks]).value_counts()
                common_peaks = common_peaks.sort_index()
                tcommon_peaks = common_peaks[common_peaks>=4]

                # combine peaks that are separated by a single scan
                spillover = {}
                spillover_peaks = tcommon_peaks.index.to_series().apply(lambda x: np.where(xdata==x)[0][0])
                spillover_peaks = spillover_peaks.sort_index()
                spillover = defaultdict(list)
                for index, value in spillover_peaks.iteritems():
                    spillover_matches = spillover_peaks==(value+1)
                    if spillover_matches.any():
                        spillover[spillover_peaks[spillover_matches].index[0]].extend(spillover.get(index, [index]))
                    else:
                        spillover[index].extend([index])
                new_common = pd.Series(0, index=spillover.keys())
                for i,v in spillover.iteritems():
                    new_common[i] += sum([tcommon_peaks[val] for val in v])

                common_peaks = new_common if new_common.any() else common_peaks
                common_peaks_deltas = sorted([(i, np.abs(i-start_rt)) for i in common_peaks.index], key=operator.itemgetter(1))
                common_peak = common_peaks_deltas[0][0]
                common_loc = np.where(xdata==common_peak)[0][0]
                common_var = pd.Series([peak['var'] for i, values in combined_peaks.items() for index, value_peaks in values.iteritems() for peak in value_peaks if peak['peak'] == common_peak]).median()
                for quant_label, quan_values in combined_peaks.items():
                    for index, values in quan_values.items():
                        if not values:
                            continue
                        rt_values = combined_data.loc[index]
                        xdata = rt_values.index.values.astype(float)
                        ydata = rt_values.fillna(0)
                        closest_rts = sorted([(i, np.abs(i['peak']-common_peak)) for i in values], key=operator.itemgetter(1))[0][0]
                        # if we move more than a # of ms1 to the dominant peak, update to our known peak
                        gc = 'k'
                        peak_loc = np.where(xdata == closest_rts['peak'])[0][0]
                        mean = closest_rts['mean']
                        amp = closest_rts['amp']
                        mean_diff = np.abs(mean-xdata[common_loc])
                        if len(xdata) >= 3 and (np.abs(peak_loc-common_loc) > 2 or mean_diff > 0.3):
                            mean = common_peak
                            amp = ydata[common_peak]
                            gc = 'g'
                        var_rat = closest_rts['var']/common_var
                        var = closest_rts['var']
                        peak_params = (amp,  mean, var)
                        # int_args = (res.x[rt_index]*mval, res.x[rt_index+1], res.x[rt_index+2])
                        # gauss beats simps/sumtraps/quadrature/fixed_quad
                        int_val = integrate.quad(peaks.gauss, xdata[0], xdata[-1], args=peak_params)[0]
                        isotope_index = isotope_labels.loc[index, 'isotope_index']
                        if int_val and not pd.isnull(int_val) and gc != 'c':
                            try:
                                quant_vals[quant_label][isotope_index] += int_val
                            except KeyError:
                                quant_vals[quant_label][isotope_index] =  int_val
                    #     quant_vals[quant_label] += integrate.simps(ydata*mval, xdata)
                        if self.html:
                            ax = fig.add_subplot(subplot_rows, subplot_columns, fig_map.get(index))
                            ax.plot(xdata, peaks.gauss(xdata, *peak_params), '{}o-'.format(gc), alpha=0.7)
                            ax.plot([start_rt, start_rt], ax.get_ylim(),'k-')
                            # ax.set_ylim(0, amp)

                for silac_label1 in data.keys():
                    qv1 = quant_vals.get(silac_label1)
                    for silac_label2 in data.keys():
                        if silac_label1 == silac_label2:
                            continue
                        qv2 = quant_vals.get(silac_label2)
                        ratio = 'NA'
                        if qv1 is not None and qv2 is not None:
                            #common_isotopes = set(qv1.keys()).intersection(qv2.keys())
                            common_isotopes = set(qv1.keys()).union(qv2.keys())
                            quant1 = sum([qv1.get(i, 0) for i in common_isotopes])
                            quant2 = sum([qv2.get(i, 0) for i in common_isotopes])
                            ratio = quant1/quant2 if quant1 and quant2 else 'NA'
                        result_dict.update({'{}_{}_ratio'.format(silac_label1, silac_label2): ratio})

                if self.html:
                    plt.tight_layout()
                    ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=100)
                    html_images['clusters'] = os.path.join(self.html['rel'], fname)
                # if self.debug:
                #     pdf.close()
                if self.debug or self.html:
                    plt.close('all')
                result_dict.update({'html_info': html_images})
                for silac_label, silac_data in data.iteritems():
                    result_dict.update({
                        '{}_intensity'.format(silac_label): silac_data['intensity'],
                        '{}_precursor'.format(silac_label): silac_data['precursor']
                    })
                self.results.put(result_dict)
            except:
                sys.stderr.write('ERROR ON {0}\n{1}\n{2}\n'.format(ms1, scan_info, traceback.format_exc()))
                return
        for params in iter(self.queue.get, None):
            run_thing(params)
        self.results.put(None)



def main():
    args = parser.parse_args()
    source = args.processed.name
    threads = args.p
    skip = args.skip
    out = args.out
    hdf = args.mzml
    save_files = args.no_temp
    sparse = args.no_sparse
    html = args.html
    resume = args.resume
    manager = Manager()

    hdf_filemap = {}

    if isinstance(hdf, list):
        nfunc = lambda i: (os.path.splitext(os.path.split(i.name)[1])[0], os.path.abspath(i.name)) if hasattr(i, 'name') else (os.path.splitext(os.path.split(i)[1])[0], os.path.abspath(i))
        hdf_filemap = dict([nfunc(i) for i in hdf])
    else:
        hdf = hdf.name if hdf else hdf
        if not hdf:
            hdf = os.path.abspath(os.path.split(source)[0])
        if os.path.isdir(hdf):
            hdf_filemap = dict([(os.path.splitext(i)[0], os.path.abspath(os.path.join(hdf, i))) for i in os.listdir(hdf) if i.lower().endswith('mzml')])
        else:
            hdf_filemap[os.path.splitext(os.path.split(hdf)[1])[0]] = os.path.abspath(hdf)

    results = GuessIterator(source, full=True, store=False, peptide=args.peptide)

    mass_scans = {}
    msf_scan_index = {}

    silac_labels = results.getSILACLabels()
    labels = silac_labels.keys()
    for silac_label in labels:
        RESULT_ORDER.extend([('{}_intensity'.format(silac_label), '{} Intensity'.format(silac_label)),
                            ('{}_precursor'.format(silac_label), '{} Precursor'.format(silac_label))])
        for silac_label2 in labels:
            if silac_label != silac_label2:
                RESULT_ORDER.extend([('{}_{}_ratio'.format(silac_label, silac_label2),
                                     '{}/{}'.format(silac_label, silac_label2))])
                RESULT_ORDER.extend([('{}_{}_residual'.format(silac_label, silac_label2),
                                     '{}/{} residual'.format(silac_label, silac_label2))])
    # ('light_clusters', 'Light Peaks Identified'), ('heavy_clusters', 'Heavy Peaks Identified'), ('light_rt_dev', 'Light Retention Time Deviation'),
    # ('heavy_rt_dev', 'Heavy Retention Time Deviation'), ('light_residuals', 'Light Intensity Residuals of Fit'),
    # ('heavy_residuals', 'Heavy Intensity Residuals of Fit'), ('cluster_ratio_deviation', 'Deviation of Heavy Clusters From Light'),
    # ('light_snr', 'Light SNR'), ('heavy_snr', 'Heavy SNR'), ('heavy_ks', 'Heavy KS'), ('light_ks', 'Light KS')

    sys.stderr.write('Loading Scans:\n')
    sample = args.sample
    for index, i in enumerate(results):
        if index%1000 == 0:
            sys.stderr.write('.')
        if args.peptide and args.peptide.lower() != i.peptide.lower():
            continue
        if not args.peptide and (sample != 1.0 and random.random() > sample):
            continue
        specId = str(i.rawId)
        if specId in mass_scans:
            continue
        mass_scans[specId] = {'id': specId, 'peptide': i.peptide, 'mod_peptide': i.modifiedPeptide, 'rt': i.rt, 'charge': i.charge, 'modifications': i.getModifications()}

    msf_scans = {}

    for index, scan in enumerate(results.getScans(modifications=False)):
        if index%1000 == 0:
            sys.stderr.write('.')
        if scan is None:
            break
        msf_scans[(str(scan.rawId), scan.peptide)] = scan
    sys.stderr.write('\nScans loaded.\n')

    raw_files = {}

    sys.stderr.write('Processing retention times.\n')
    for specId, i in mass_scans.iteritems():
        if index%1000 == 0:
            sys.stderr.write('.')
        scan = msf_scans.get((specId, i.get('peptide')))
        fname = os.path.splitext(scan.file)[0]
        if fname not in hdf_filemap:
            fname = os.path.split(fname)[1]
            if fname not in hdf_filemap:
                if skip:
                    continue
                sys.stderr.write('{0} not found in filemap. Filemap is {1}.'.format(fname, hdf_filemap))
                return 1
        i['file'] = fname
        try:
            raw_files[fname].append((i.get('rt'), i))
        except KeyError:
            raw_files[fname] = [(i.get('rt'), i)]

    # sort by RT so we can minimize our memory footprint by throwing away scans we no longer need
    for i in raw_files.keys():
        raw_files[i].sort(key=operator.itemgetter(0))
        raw_files[i] = [j[1] for j in raw_files[i]]

    sys.stderr.write('\nRetention times loaded.\n')

    workers = []
    completed = 0
    sys.stderr.write('Beginning SILAC quantification.\n')
    scan_count = len(mass_scans)
    headers = ['Raw File']+[i[1] for i in RESULT_ORDER]
    if resume:
        if not out:
            sys.stderr.write('You may only resume runs with a file output.\n')
            return -1
        out = open(out, 'ab')
        out_path = out.name
    else:
        if out:
            out = open(out, 'wb')
            out_path = out.name
        else:
            out = sys.stdout
            out_path = source
        out.write('{0}\n'.format('\t'.join(headers)))
    html_list = []
    if html:
        import unicodedata
        value = unicodedata.normalize('NFKD', unicode(os.path.splitext(os.path.split(out_path)[1])[0])).encode('ascii', 'ignore')
        value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
        value = unicode(re.sub('[-\s]+', '-', value))
        html = os.path.join(os.path.split(out_path)[0], os.path.normpath(value)+'_images')
        try:
            os.mkdir(html)
        except OSError:
            pass
        html = {'full': os.path.join(os.path.split(out_path)[0], os.path.normpath(value)+'_images'),
                'rel':os.path.normpath(value)+'_images' }
        def table_rows(html_list, res=None):
            # each item is a string like a\tb\tc
            if html_list:
                d = html_list.pop(0)
            else:
                return res
            l = d['table']
            images = d['images']
            if res is None:
                res = '<tr>'
            out = []
            for i,v in zip(l.split('\t'), headers):
                rt_format = '{}_rt'.format(v.split(' ')[0])
                if v.endswith('Intensity') and rt_format in images:
                    out.append("""<td data-toggle="popover" data-trigger="click" data-content='<img src="{0}">'>{1}</td>""".format(images.get(rt_format), i))
                elif v == 'Peptide' and 'clusters' in images:
                    out.append("""<td data-toggle="popover" data-trigger="click" data-content='<img src="{0}">'>{1}</td>""".format(images.get('clusters'), i))
                else:
                    out.append('<td>{0}</td>'.format(i))
            res += '\n'.join(out)+'</tr>'
            return table_rows(html_list, res=res)
        if resume:
            html_out = open('{0}.html'.format(out_path), 'ab')
        else:
            html_out = open('{0}.html'.format(out_path), 'wb')
            html_out.write(
                    """<!DOCTYPE html>
                    <html>
                    <head lang="en">
                        <meta charset="UTF-8">
                        <title>{0}</title>
                        <link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css" type="text/css">
                        <link rel="stylesheet" href="http://cdn.datatables.net/1.10.5/css/jquery.dataTables.css" type="text/css">
                    </head>
                    <body>
                        <table id="raw-table" class="table table-striped table-bordered table-hover">
                            <thead>
                                <tr>
                                {1}
                                </tr>
                            </thead>
                            <tbody>
                    """.format(
                        source,
                        '\n'.join(['<th>{0}</th>'.format(i) for i in ['Raw File']+[i[1] for i in RESULT_ORDER]])
                    )
            )

    skip_map = set([])
    if resume:
        import csv
        key = None
        for entry in csv.reader(open(out.name, 'rb'), delimiter='\t'):
            # key is filename, peptide, modifications, charge, ms2 id
            key = (entry[0], entry[1], entry[5], entry[6], entry[8])
            skip_map.add(key)
        skip_map.discard(key)

    for filename, mass_scans in raw_files.iteritems():
        filepath = hdf_filemap[filename]
        sys.stderr.write('Processing {0}.\n'.format(filename))
        in_queue = Queue()
        result_queue = Queue()
        reader_in, reader_out = Queue(), Queue()
        scan_dict = manager.dict()

        ms2_ms1_map = manager.dict()
        scan_rt_map = manager.dict()
        rt_index_map = manager.dict()
        reader = Reader(reader_in, reader_out, scan_dict=scan_dict, raw_file=filepath, ms2_ms1_map=ms2_ms1_map, scan_rt_map=scan_rt_map,
                        rt_index_map=rt_index_map)
        reader.start()
        for i in xrange(threads):
            worker = Worker(queue=in_queue, results=result_queue, raw_name=filepath, silac_labels=silac_labels,
                            temp=save_files, sparse=sparse, debug=args.debug, html=html, mono=args.no_spread,
                            reader_in=reader_in, reader_out=reader_out, scan_dict=scan_dict, ms2_ms1_map=ms2_ms1_map,
                            scan_rt_map=scan_rt_map, rt_index_map=rt_index_map)
            workers.append(worker)
            worker.start()

        for scan_index, v in enumerate(mass_scans):
            if resume:
                key = (v.get('file'), v.get('peptide'), v.get('modifications'), v.get('charge'), v.get('id'))
                if key in skip_map:
                    completed += 1
                    continue
            # if scan_index % 100 != 0:
            #     continue
            params = {}
            params['scan_info'] = v
            scan = msf_scans.get((v.get('id'), v.get('peptide')))
            params['scan'] = scan
            in_queue.put(params)

        sys.stderr.write('{0} processed and placed into queue.\n'.format(filename))

        in_queue.put(None)
        none_count = 0
        while workers:
            try:
                result = result_queue.get(timeout=0.1)
            except Empty:
                continue
            if result is None:
                none_count += 1
            if none_count == threads:
                [i.terminate() for i in workers]
                workers = []
            elif result is not None:
                completed += 1
                if completed % 10 == 0:
                    sys.stderr.write('\r{0:2.2f}% Completed'.format(completed/scan_count*100))
                    sys.stderr.flush()
                    # sys.stderr.write('Processed {0} of {1}...\n'.format(completed, scan_count))
                res = '{0}\t{1}\n'.format(filename, '\t'.join([str(result.get(i[0], 'NA')) for i in RESULT_ORDER]))
                out.write(res)
                out.flush()
                if html:
                    html_out.write(table_rows([{'table': res.strip(), 'images': result['html_info']}]))
                    html_out.flush()
        reader_in.put(None)

    while workers:
        try:
            result = result_queue.get(timeout=0.1)
        except Empty:
            continue
        if result is None:
            none_count += 1
        if none_count == threads:
            [i.terminate() for i in workers]
            workers = []
        else:
            completed += 1
            if completed % 10 == 0:
                sys.stderr.write('\r{0:2.2f}% Completed'.format(completed/scan_count*100))
                sys.stderr.flush()
                # sys.stderr.write('Processed {0} of {1}...\n'.format(completed, scan_count))
            res = '{0}\t{1}\n'.format(filename, '\t'.join([str(result.get(i[0], 'NA')) for i in RESULT_ORDER]))
            out.write(res)
            out.flush()
            if html:
                html_out.write(table_rows([{'table': res.strip(), 'images': result['html_info']}]))
                html_out.flush()



    out.flush()
    out.close()

    if html:
        html_out.write(
                """
                        </tbody>
                    </table>
                </body>
                <footer></footer>
                <script type="text/javascript" src="http://code.jquery.com/jquery-1.11.1.min.js"></script>
                <script type="text/javascript" src="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/js/bootstrap.min.js"></script>
                <script type="text/javascript" src="http://cdn.datatables.net/1.10.5/js/jquery.dataTables.min.js"></script>
                <script>
                    $(document).ready(function() {
                        $('#raw-table').DataTable({
                            "iDisplayLength": 100,
                        });
                        var reload = false;
                        $('[data-toggle="popover"]').popover({ html : true, container: 'footer', placement: 'bottom' });
                        $('#raw-table').on( 'page.dt search.dt init.dt order.dt length.dt', function () {
                        	reload = true;
                        });
                        $('#raw-table').on( 'draw.dt', function () {
                        	if(reload){
                        		$('[data-toggle="popover"]').popover({ html : true, container: 'footer', placement: 'bottom' });
                        		reload = false;
                        	}
                        });
                        // from http://stackoverflow.com/questions/20466903/bootstrap-popover-hide-on-click-outside
                        $(document).on('click', function (e) {
                            $('[data-toggle=popover]').each(function () {
                                // hide any open popovers when the anywhere else in the body is clicked
                                if (!$(this).is(e.target) && $(this).has(e.target).length === 0 && $('.popover').has(e.target).length === 0) {
                                    $(this).popover('hide');
                                }
                            });
                        });
                    });
                </script>
                </html>
                """
            )


if __name__ == "__main__":
    sys.exit(main())
