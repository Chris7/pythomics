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
import traceback
import pandas as pd
import numpy as np
import re
from collections import OrderedDict
from itertools import izip_longest
from matplotlib import pyplot as plt
from multiprocessing import Process, Queue
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

class Worker(Process):
    def __init__(self, queue=None, results=None, precision=6, raw_name=None, silac_labels=None,
                 temp=False, sparse=False, debug=False, html=False, mono=False):
        super(Worker, self).__init__()
        self.precision = precision
        self.sparse = sparse
        self.queue=queue
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

    def getScan(self, ms1, dense=True):
        import numpy as np
        if ms1 not in self.data:
            scan = self.raw.getScan(ms1)
            scan_vals = np.array(scan.scans)
            res = pd.Series(scan_vals[:, 1].astype(np.uint64), index=np.round(scan_vals[:, 0], self.precision), name=scan.rt, dtype='uint64')
            # due to precision, we have multiple m/z values at the same place. We can eliminate this by grouping them and summing them.
            # Summation is the correct choice here because we are combining values of a precision higher than we care about.
            res = res.groupby(level=0).sum()
            if self.sparse:
                res = res.to_sparse(fill_value=0)
            self.data[ms1] = res
        return self.data[ms1].to_dense() if dense else self.data[ms1]

    @staticmethod
    def sign(a):
        return -1 if a < 0 else 1

    def cleanScans(self, rt_map, current_rt):
        min_rt = current_rt - self.rt_width
        if min_rt > self.last_clean:
            old_scans = rt_map[self.last_clean:min_rt].values
            self.last_clean = current_rt - self.rt_width
            for i in old_scans:
                if i in self.data:
                    del self.data[i]
            import gc
            gc.collect()

    def gfit(self, params, X, Y, amp, rt):
        #mu, std, offset = params[0], params[1], params[2]
        mu, std = params[0], params[1]
        if abs(mu-rt) > self.rt_tol:
            return np.inf
        fit = amp*np.exp(-(X-mu)**2/(2*std**2))
        lb, rb = mu-2*abs(std), mu+2*abs(std)
        fit = np.array([j for i,j in zip(X, fit) if i>=lb or i<=rb])
        nfit = fit/max(fit)
        ny = np.array([j for i,j in zip(X, Y) if i>=lb or i<=rb])
        ny = ny/max(ny)
        #residuals = (nfit-ny)**2
        #nfit = fit/max(fit)
        #ny = Y/max(Y)
        residuals = sum((nfit-ny)**2)
        #residuals = sum((fit-Y)**2)
        return residuals

    def run(self):
        import time
        args = parser.parse_args()
        source = self.raw_name
        self.raw = GuessIterator(source, full=False, store=False, ms_filter=1)
        self.data = {}
        self.data2 = pd.DataFrame(dtype='uint16')
        scan_rt_map = {}
        rt_index_map = {}
        ms2_ms1_map = {}
        sys.stderr.write('Thread processing peptide file.\n')
        for index, i in enumerate(self.raw):
            if index % 100 == 0:
                sys.stderr.write('.')
            if i is None:
                break
            if i.ms_level != 1:
                ms2_ms1_map[int(i.id)] = scan_id
                continue
            scan_id = int(i.id)
            rt = i.rt
            scan_rt_map[scan_id] = rt
            rt_index_map[rt] = scan_id
        sys.stderr.write('Scans Loaded\n')
        rt_index_map = pd.Series(rt_index_map)

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
                        # count += len(labels)
                        silac_shift += sum(labels)
                    precursors[silac_label] = silac_shift
                    data[silac_label] = copy.deepcopy(silac_dict)
                precursors = OrderedDict(sorted(precursors.items(), key=operator.itemgetter(1)))
                shift_maxes = {i: j for i,j in zip(precursors.keys(), precursors.values()[1:])}
                rt_width = int(len(rt_index_map[rt-self.rt_width:rt+self.rt_width])/2)
                # sys.stderr.write('{}\n'.format(rt_width))
                base_rt = rt_index_map.index.searchsorted(rt)
                finished = set([])
                result_dict = {'peptide': scan_info.get('mod_peptide', peptide),
                               'scan': scanId, 'ms1': ms1,
                               'charge': charge, 'modifications': mods}
                rt_window = []
                ms_index = 0
                delta = -1
                theo_dist = peaks.calculate_theoretical_distribution(peptide.upper())
                spacing = config.NEUTRON/float(charge)
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
                                                      charge=charge, ppm=5, heavy=False, theo_dist=theo_dist, label=precursor_label)
                        peaks_found = data[precursor_label]['peaks']
                        # look for Proline/Glutamate/Glutamines
                        # if 'P' in peptide or 'E' in peptide or 'Q' in peptide:
                        #     heavy2 = peaks.findEnvelope(df, start_mz=heavy_precursor, isotope_offset=6, charge=charge, ppm=10, heavy=True)
                        #     print rt_index_map.iloc[base_rt+ms_index], heavy_precursor
                        #     print heavy
                        #     print heavy2
                        if not envelope['envelope']:
                            if not peaks_found:
                                finished.add(precursor_label)
                            continue
                        found = True
                        for i in ['envelope', 'micro_envelopes', 'ppms']:
                            for isotope, vals in envelope.get(i, {}).iteritems():
                                val_dict = {'info': vals, 'df': df}
                                try:
                                    peaks_found[isotope][i].append(val_dict)
                                except KeyError:
                                    try:
                                        peaks_found[isotope][i] = [val_dict]
                                    except KeyError:
                                        peaks_found[isotope] = {i: [val_dict]}
                        selected = {}
                        for isotope, vals in envelope['micro_envelopes'].iteritems():
                            selected[precursor_mz+isotope*spacing] = df.iloc[range(*vals)].sum()
                        selected = pd.Series(selected, name=df.name).to_frame()
                        if df.name in combined_data.columns:
                            combined_data = combined_data.add(selected, axis='index', fill_value=0)
                        else:
                            combined_data = pd.concat([combined_data, selected], axis=1).fillna(0)
                    if found is False:
                        if delta == -1:
                            delta = 1
                            ms_index=1
                        elif ms_index >= 2:
                            break
                    ms_index += delta
                # shared_isotopes = set([])
                rt_window.sort()
                # if self.mono:
                #     for silac_label, silac_data in data.iteritems():
                #         peaks_found = silac_data.get('peaks')
                #         found_isotopes = set(peaks_found.keys())
                #         if found_isotopes:
                #             shared_isotopes = shared_isotopes.intersection(found_isotopes) if shared_isotopes else found_isotopes
                # combine all our peaks and get our RT boundaries. Do it this way fo ra
                combined_data = combined_data.sort(axis='index').sort(axis='columns')
                for silac_label, silac_data in data.iteritems():
                    peaks_found = silac_data.get('peaks')
                    peak_data = peaks.buildEnvelope(peaks_found=peaks_found, rt_window=rt_window, start_rt=rt, silac_label=silac_label)
                    silac_data['df'] = peak_data.get('data')

                    if self.debug or self.html:
                        if self.debug:
                            pdf = PdfPages('{2}_{0}_{1}_fix.pdf'.format(peptide, ms1, self.raw_name))
                            plt.figure()
                    light_rt_data = peak_data.get('rt')
                    light_int = peak_data.get('integration')
                    silac_data['intensity'] = light_int
                    if self.debug or self.html:
                        plt.figure()
                        if not light_rt_data.empty:
                            ax = light_rt_data.plot(color='b', title=str(light_int))
                            if self.debug:
                                pdf.savefig(figure=ax.get_figure())
                            if self.html:
                                fname = '{2}_{0}_{1}_{3}_{4}_rt.png'.format(peptide, ms1, self.filename, scanId, silac_label)
                                ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                                html_images['{}_rt'.format(silac_label)] = os.path.join(self.html['rel'], fname)
                data_keys = data.keys()
                for silac_label in data_keys:
                    label1_keys = set(data[silac_label]['peaks'].keys())
                    label1_data = data[silac_label]
                    if not self.mono:
                        #Chris_Ecoli_1-2-4_quant_ns - pearson medium 0.82, heavy, 0.93
                        label1_int = label1_data.get('intensity')
                        label1_array = label1_data.get('df').fillna(0).sum(axis=1)
                        if not label1_array.empty:
                            array1_fit = peaks.fit_data(label1_array, charge=float(charge), peptide=peptide.upper())
                            #label1_int = array1_fit.get('fit').sum()
                            label1_res = array1_fit.get('residual')
                            label1_data['nb'] = array1_fit.get('fit')
                        else:
                            #label1_int = label1_array.sum()
                            label1_res = np.inf
                            #Chris_Ecoli_1-2-4_quant_ns2 -  pearson medium 0.71, heavy 0.92
                            #label1_int = label1_data['df'].fillna(0).sum(axis=1).sum()
                    for silac_label2 in data_keys:
                        if silac_label == silac_label2:
                            continue
                        label2_data = data[silac_label2]
                        if self.mono:
                            # compare the ratios of monoisotopic peaks
                            shared_isotopes = label1_keys.intersection(set(label2_data['peaks'].keys()))
                            if not shared_isotopes:
                                continue
                            # force the peak patterns to be in the same rough shape
                            label1_diff = label1_data['df'].fillna(0).sum(axis=1).diff().fillna(0)
                            label2_diff = label2_data['df'].fillna(0).sum(axis=1).diff().fillna(0)
                            # si2 = set([])
                            # # print label1_diff
                            # # print label2_diff
                            # # print shared_isotopes
                            # for i in sorted(list(shared_isotopes)):
                            #     if self.sign(label1_diff.iloc[i]) == self.sign(label2_diff.iloc[i]):
                            #         si2.add(i)
                            #     else:
                            #         break
                            # shared_isotopes = si2

                            # This method is inferior
                            # peaks_found = label1_data.get('peaks')
                            # peak_data = peaks.buildEnvelope(peaks_found=peaks_found, isotopes=shared_isotopes,
                            #                                 rt_window=rt_window)
                            ##Chris_Ecoli_1-2-4_quant_wide3 pearson: medium 0.55, heavy 0.73
                            # label1_int = peak_data.get('integration')
                            ##Chris_Ecoli_1-2-4_quant_wide4 pearson: medium 0.72, heavy 0.96
                            # label1_int = peak_data.get('rt').sum()
                            ##Chris_Ecoli_1-2-4_quant_wide2 pearson: medium 0.72, heavy 0.96
                            i1_data = label1_data['df'].fillna(0).sum(axis=1).iloc[list(shared_isotopes)]
                            if not i1_data.empty:
                                array1_fit = peaks.fit_data(i1_data , charge=float(charge), peptide=peptide.upper())
                                label1_array = array1_fit.get('fit')
                                label1_int = array1_fit.get('fit').sum()
                                label1_res = array1_fit.get('residual')
                                label1_data['nb'] = array1_fit.get('fit')
                            else:
                                label1_int = i1_data.sum()
                                label1_res = np.inf

                            # peaks_found = label2_data.get('peaks')
                            # peak_data = peaks.buildEnvelope(peaks_found=peaks_found, isotopes=shared_isotopes,
                            #                                 rt_window=rt_window)
                            # label2_int = peak_data.get('integration')
                            # label2_int = peak_data.get('rt').sum()
                            i2_data = label2_data['df'].fillna(0).sum(axis=1).iloc[list(shared_isotopes)]
                            if not i2_data.empty:
                                array2_fit = peaks.fit_data(i2_data , charge=float(charge), peptide=peptide.upper())
                                label2_array = array2_fit.get('fit')
                                label2_int = array2_fit.get('fit').sum()
                                label2_res = array2_fit.get('residual')
                                label2_data['nb'] = array2_fit.get('fit')
                            else:
                                label2_int = i2_data.sum()
                                label2_res = np.inf
                        else:
                            label2_int = label2_data.get('intensity')
                            label2_array = label2_data.get('df').fillna(0).sum(axis=1)
                            if not label2_array.empty:
                                array2_fit = peaks.fit_data(label2_array, charge=float(charge), peptide=peptide.upper())
                                #label2_int = array2_fit.get('fit').sum()
                                label2_res = array2_fit.get('residual')
                                label2_data['nb'] = array2_fit.get('fit')
                            else:
                                #label2_int = label2_array.sum()
                                label2_res = np.inf

                        if self.mono and not label1_array.empty and not label2_array.empty:
                            # we want the intercept for the ratio
                            slope, intercept, r_value, p_value, std_err = linregress(range(len(shared_isotopes)), label1_array.iloc[list(shared_isotopes)].values/label2_array.iloc[list(shared_isotopes)].values)
                            result_dict.update({'{}_{}_ratio'.format(silac_label, silac_label2):
                                                    intercept})
                            result_dict.update({'{}_{}_residual'.format(silac_label, silac_label2):
                                                    '{}'.format(abs(slope))})
                        else:
                            #return
                            result_dict.update({'{}_{}_ratio'.format(silac_label, silac_label2):
                                                    label1_int/label2_int if label1_int and label2_int else 'NA'})
                            result_dict.update({'{}_{}_residual'.format(silac_label, silac_label2):
                                                    '{}'.format(label1_res+label2_res)})

                if self.html:
                    fname = '{2}_{0}_{1}_{3}_clusters.png'.format(peptide, ms1, self.filename, scanId)
                    plt.figure()
                    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
                    legends = []
                    for silac_label, _ in precursors.iteritems():
                        color = colors.pop(0) if colors else 'black'
                        silac_data = data.get(silac_label)
                        x, y = [],[]
                        ax = plt.subplot(2, 1, 1)
                        cseries = silac_data.get('nb')
                        if cseries is not None:
                            for index, val in zip(cseries.index, cseries):
                                x += [index,index,index]
                                y += [0,val,0]
                            ax.bar(cseries.index, cseries.values, width=1.0/float(charge)/2,
                                   color='{}'.format(color), alpha=0.7, label=silac_label)
                        x, y = [],[]
                        light_df = silac_data.get('df')
                        ax = plt.subplot(2, 1, 2)
                        cseries = light_df.fillna(0).sum(axis=1)
                        for index, val in zip(cseries.index, cseries):
                            x += [index,index,index]
                            y += [0,val,0]
                        line = ax.plot(x,y, '{}'.format(color), alpha=0.7)
                        legends.append((line[0], silac_label))

                    # Put a legend above our top plot
                    # shift our plots down
                    for i in xrange(1,3):
                        ax = plt.subplot(2, 1, i)
                        box = ax.get_position()
                        # set_position is [left, bottom, width, height]
                        ax.set_position([box.x0, box.y0 - box.height * 0.1,
                                         box.width, box.height * 0.9])
                    ax.get_figure().legend([i[0] for i in legends], [i[1] for i in legends], loc='upper center',
                                           fancybox=True, shadow=True, ncol=3)
                    ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                    html_images['clusters'] = os.path.join(self.html['rel'], fname)
                if self.debug:
                    pdf.close()
                if self.debug or self.html:
                    plt.close('all')
                result_dict.update({'html_info': html_images})
                for silac_label, silac_data in data.iteritems():
                    result_dict.update({
                        '{}_intensity'.format(silac_label): silac_data['intensity'],
                        '{}_precursor'.format(silac_label): silac_data['precursor']
                    })
                # print result_dict
                self.results.put(result_dict)
                self.cleanScans(rt_index_map, rt)
                # print 'run tim eis', time.time()-start_time
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
        if not args.peptide and (sample != 1.0 and random.random() < sample):
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

    workers = {}
    in_queues = {}
    result_queues = {}
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
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images.get(rt_format), i))
                elif v == 'Peptide' and 'clusters' in images:
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images.get('clusters'), i))
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
        for entry in csv.reader(open(out.name, 'rb'), delimiter='\t'):
            # key is filename, peptide, modifications, charge, ms2 id
            key = (entry[0], entry[1], entry[5], entry[6], entry[8])
            skip_map.add(key)

    for filename, mass_scans in raw_files.iteritems():
        filepath = hdf_filemap[filename]
        sys.stderr.write('Processing {0}.\n'.format(filename))
        in_queue = Queue()
        result_queue = Queue()
        in_queues[filename] = in_queue
        result_queues[filename] = result_queue
        worker = Worker(queue=in_queue, results=result_queue, raw_name=filepath, silac_labels=silac_labels,
                        temp=save_files, sparse=sparse, debug=args.debug, html=html, mono=args.no_spread)
        workers[filename] = worker
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
        if len(workers) == threads:
            # clear the queues
            while workers:
                empty_queue = set([])
                for rawstr, queue in result_queues.iteritems():
                    if rawstr in empty_queue:
                        continue
                    try:
                        result = queue.get(timeout=0.1)
                    except Empty:
                        continue
                    if result is None:
                        worker = workers[rawstr]
                        if worker.is_alive():
                            worker.terminate()
                        del workers[rawstr]
                        empty_queue.add(rawstr)
                    else:
                        completed += 1
                        if completed % 10 == 0:
                            sys.stderr.write('\r{0:2.2f}% Completed'.format(completed/scan_count*100))
                            sys.stderr.flush()
                            # sys.stderr.write('Processed {0} of {1}...\n'.format(completed, scan_count))
                        res = '{0}\t{1}\n'.format(rawstr, '\t'.join([str(result.get(i[0], 'NA')) for i in RESULT_ORDER]))
                        out.write(res)
                        out.flush()
                        if html:
                            html_out.write(table_rows([{'table': res.strip(), 'images': result['html_info']}]))
                            html_out.flush()
            workers = {}
            result_queues = {}

    while workers:
        empty_queue = set([])
        for rawstr, queue in result_queues.iteritems():
            if rawstr in empty_queue:
                continue
            try:
                result = queue.get(timeout=0.1)
            except Empty:
                continue
            if result is None:
                worker = workers[rawstr]
                if worker.is_alive():
                    worker.terminate()
                del workers[rawstr]
                empty_queue.add(rawstr)
            else:
                completed += 1
                if completed % 10 == 0:
                    sys.stderr.write('\r{0:2.2f}% Completed'.format(completed/scan_count*100))
                    sys.stderr.flush()
                    # sys.stderr.write('Processed {0} of {1}...\n'.format(completed, scan_count))
                res = '{0}\t{1}\n'.format(rawstr, '\t'.join([str(result.get(i[0],'NA')) for i in RESULT_ORDER]))
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
                    });
                </script>
                </html>
                """
            )


if __name__ == "__main__":
    sys.exit(main())
