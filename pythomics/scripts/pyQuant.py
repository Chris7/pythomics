#!/usr/bin/env python
from __future__ import division, unicode_literals

# TODO: Threading for single files. Since much time is spent in fetching m/z records, it may be pointless
# because I/O is limiting. If so, create a thread just for I/O for processing requests the other threads
# interface with

description = """
This will quantify labeled peaks (such as SILAC) in ms1 spectra. It relies solely on the distance between peaks,
 which can correct for errors due to amino acid conversions.
"""
from guppy import hpy
import sys
import os
import copy
import operator
from scipy.stats import mode
import operator
from pythomics.proteomics import peaks, config
import traceback
import pandas as pd
import numpy as np
import re
from collections import OrderedDict, defaultdict
from itertools import izip_longest
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Process, Queue, Manager, Array
import ctypes
import sqlite3
from scipy.stats import linregress
try:
    #from profilestats import profile
    from memory_profiler import profile
except ImportError:
    pass
from Queue import Empty
from scipy import integrate, optimize
from scipy.stats import signaltonoise, linregress

from matplotlib.backends.backend_pdf import PdfPages
import argparse
from cStringIO import StringIO
from scipy import integrate

from pythomics.templates import CustomParser
from pythomics.proteomics.parsers import GuessIterator
from pythomics.proteomics import peaks, config

RESULT_ORDER = [('peptide', 'Peptide'), ('modifications', 'Modifications'),
                ('charge', 'Charge'), ('ms1', 'MS1 Spectrum ID'), ('scan', 'MS2 Spectrum ID'), ('rt', 'Retention Time')]


parser = CustomParser(description=description)
parser.add_processed_ms()
parser.add_argument('--mzml', help="The mzML files for the raw data. If not provided, raw files assumed to be in the directory of the processed file.", type=argparse.FileType('r'), nargs='*')
parser.add_argument('--raw-dir', help="The directory containing raw data.", type=str)

parser.add_argument('--filemap', help="By default, the mzML file is assumed to be a derivation of the scan filename. " \
                                     "If this is not true, this file provides a correct mapping.")
parser.add_argument('--rt-width', help="The width of the retention time window to search for. Default: 1.0 minute", type=float, default=1.0)
parser.add_argument('--precision', help="The precision for storing m/z values. Defaults to 6 decimal places.", type=int, default=6)
parser.add_argument('--precursor-ppm', help="The mass accuracy for the first monoisotopic peak in ppm.", type=float, default=5)
parser.add_argument('--isotope-ppm', help="The mass accuracy for the isotopic cluster.", type=float, default=2.5)
parser.add_argument('--spread', help="Assume there is spread of the isotopic label.", action='store_true')
parser.add_argument('--debug', help="This will output debug information and graphs.", action='store_true')
parser.add_argument('--skip', help="If true, skip scans with missing files in the mapping.", action='store_true')
parser.add_argument('--peptide', help="The peptide to limit quantification to.", type=str)
parser.add_argument('--html', help="Output a HTML table summary.", action='store_true')
parser.add_argument('--resume', help="Will resume from the last run. Only works if not directing output to stdout.", action='store_true')
parser.add_argument('--sample', help="How much of the data to sample. Enter as a decimal (ie 1.0 for everything, 0.1 for 10%%)", type=float, default=1.0)
parser.add_argument('--disable-stats', help="Disable confidence statistics on data.", action='store_true')
parser.add_argument('-o', '--out', nargs='?', help='The prefix for the file output', type=str)
tsv_group = parser.add_argument_group('Tabbed File Input')
tsv_group.add_argument('--tsv', help='Indicate the procesed argument is a delimited file.', action='store_true')
tsv_group.add_argument('--label', help='The column indicating the label state of the peptide. If not found, entry assumed to be light variant.', default='Labeling State')
tsv_group.add_argument('--peptide-col', help='The column indicating the peptide.', default='Sequence')
tsv_group.add_argument('--rt', help='The column indicating the retention time.', default='Retention time')
tsv_group.add_argument('--mz', help='The column indicating the MZ value of the precursor ion. This is not the MH+.', default='m/z')
tsv_group.add_argument('--ms2', help='The column indicating the ms2 scan corresponding to the ion.', default='MS/MS Scan Number')
tsv_group.add_argument('--charge', help='The column indicating the charge state of the ion.', default='Charge')
tsv_group.add_argument('--source', help='The column indicating the raw file the scan is contained in.', default='Raw file')
tsv_group.add_argument('--label-scheme', help='The file corresponding to the labeling scheme utilized.', type=argparse.FileType('r'))

import random

class Reader(Process):
    def __init__(self, incoming, outgoing, raw_file=None, ms2_ms1_map=None, scan_rt_map=None, rt_index_map=None, scans_loaded=None):
        super(Reader, self).__init__()
        self.incoming = incoming
        self.outgoing = outgoing
        self.ms2_ms1_map = ms2_ms1_map
        self.scan_rt_map = scan_rt_map
        self.rt_index_map = rt_index_map
        self.scan_dict = {}
        self.scans_loaded = scans_loaded
        self.raw = raw_file

    def run(self):
        for scan_request in iter(self.incoming.get, None):
            # print scan_request
            thread, scan_id = scan_request
            d = self.scan_dict.get(scan_id)
            if not d:
                scan = self.raw.getScan(scan_id)
                scan_vals = np.array(scan.scans)
                # add to our database
                d = {'vals': scan_vals, 'rt': scan.rt, 'title': int(scan.title)}
                self.scan_dict[scan_id] = d
                # the scan has been stored, delete it
                del scan
            self.outgoing[thread].put(d)
        sys.stderr.write('reader done\n')

class Worker(Process):
    def __init__(self, queue=None, results=None, precision=6, raw_name=None, silac_labels=None, isotope_ppms=None,
                 debug=False, html=False, mono=False, scans_loaded=None, precursor_ppm=5.0, isotope_ppm=2.5,
                 reader_in=None, reader_out=None, ms2_ms1_map=None, scan_rt_map=None, rt_index_map=None,
                 last_scan_accessed=None, thread=None, fitting_run=False):
        super(Worker, self).__init__()
        self.scans_loaded = scans_loaded
        self.precision = precision
        self.precursor_ppm = precursor_ppm
        self.isotope_ppm = isotope_ppm
        self.queue=queue
        self.reader_in, self.reader_out = reader_in, reader_out
        self.ms2_ms1_map, self.scan_rt_map, self.rt_index_map = ms2_ms1_map, scan_rt_map, rt_index_map
        self.results = results
        self.silac_labels = {'Light': {}} if silac_labels is None else silac_labels
        self.shifts = {0: "Light"}
        self.shifts.update({sum(silac_masses.keys()): silac_label for silac_label, silac_masses in self.silac_labels.iteritems()})
        self.raw_name = raw_name
        self.filename = os.path.split(self.raw_name)[1]
        self.rt_width = 1
        self.rt_tol = 0.2 # for fitting
        self.debug = debug
        self.html = html
        self.mono = mono
        self.last_scan_accessed = last_scan_accessed
        self.thread = thread
        self.fitting_run = fitting_run
        self.isotope_ppms = isotope_ppms

    def convertScan(self, scan):
        import numpy as np
        scan_vals = scan['vals']
        res = pd.Series(scan_vals[:, 1].astype(np.uint64), index=np.round(scan_vals[:, 0], self.precision), name=scan['rt'], dtype='uint64')
        del scan_vals
        # due to precision, we have multiple m/z values at the same place. We can eliminate this by grouping them and summing them.
        # Summation is the correct choice here because we are combining values of a precision higher than we care about.
        return res.groupby(level=0).sum()

    def getScan(self, ms1):
        self.last_scan_accessed[self.thread] = ms1
        ms1 = int(ms1)
        self.reader_in.put((self.thread, ms1))
        scan = self.reader_out.get()
        return self.convertScan(scan)

    def run_thing(self, params, scan_rt_map, rt_index_map, ms2_ms1_map, silac_shifts):
        try:
            html_images = {}
            scan_info = params.get('scan_info')
            ms2_info = params.get('ms2_info')
            scanId = scan_info['id']
            if not ms2_info:
                return
            ms1 = ms2_info['ms1']
            infer = ms1 is None
            if infer:
                ms1 = ms2_ms1_map.get(int(ms2_info['rawId']))
            else:
                try:
                    ms1 = int(ms1)
                except AttributeError:
                    ms1 = ms2_ms1_map.get(int(ms2_info['rawId']))
            rt = scan_info.get('rt', scan_rt_map.get(ms1))
            if not rt:
                return
            mods = scan_info.get('modifications')
            charge = int(ms2_info['charge'])
            precursor = (float(ms2_info['mass'])-config.HYDROGEN)*charge if infer else float(ms2_info['mass'])+config.HYDROGEN*(charge-1)
            shift = 0
            if mods is None:
                shift += self.silac_labels.get('label', 0)
            else:
                for mod in filter(lambda x: x, mods.split('|')):
                    aa, pos, mass, type = mod.split(',', 3)
                    mass = float('{0:0.5f}'.format(float(mass)))
                    if aa in silac_shifts.get(mass, {}):
                        shift += mass
            if shift > 0:
                precursor -= shift
            theor_mass = scan_info['theor_mass']-shift/float(charge)
            peptide = ms2_info['peptide']
            if self.debug:
                sys.stderr.write('thread {4} on ms {0} {1} {2} {3}\n'.format(ms1, rt, precursor, scan_info, id(self)))
            precursors = {'Light': 0.0}
            silac_dict = {'data': None, 'df': pd.DataFrame(), 'precursor': 'NA',
                          'isotopes': {}, 'peaks': OrderedDict(), 'intensity': 'NA'}
            data = OrderedDict()
            data['Light'] = copy.deepcopy(silac_dict)
            combined_data = pd.DataFrame()
            for silac_label, silac_masses in self.silac_labels.items():
                silac_shift=0
                global_mass = None
                added_residues = set([])
                cterm_mass = 0
                nterm_mass = 0
                for label_mass, label_masses in silac_masses.items():
                    if 'X' in label_masses:
                        global_mass = label_mass
                    if ']' in label_masses:
                        cterm_mass = label_mass
                    if '[' in label_masses:
                        nterm_mass = label_mass
                    added_residues = added_residues.union(label_masses)
                    labels = [label_mass for mod_aa in peptide if mod_aa in label_masses]
                    silac_shift += sum(labels)
                if global_mass is not None:
                    silac_shift += sum([global_mass for mod_aa in peptide if mod_aa not in added_residues])
                silac_shift += cterm_mass+nterm_mass
                # get the non-specific ones
                precursors[silac_label] = silac_shift
                data[silac_label] = copy.deepcopy(silac_dict)
            precursors = OrderedDict(sorted(precursors.items(), key=operator.itemgetter(1)))
            shift_maxes = {i: j for i,j in zip(precursors.keys(), precursors.values()[1:])}
            base_rt = rt_index_map.index.searchsorted(rt)
            finished = set([])
            finished_isotopes = {i: set([]) for i in precursors.keys()}
            result_dict = {'peptide': scan_info.get('mod_peptide', peptide),
                           'scan': scanId, 'ms1': ms1, 'charge': charge,
                           'modifications': mods, 'rt': rt}
            rt_window = []
            ms_index = 0
            delta = -1
            theo_dist = peaks.calculate_theoretical_distribution(peptide.upper())
            spacing = config.NEUTRON/float(charge)
            isotope_labels = {}
            isotopes_chosen = {}
            # chosen_window = defaultdict(dict)
            last_precursors = {-1: {}, 1: {}}
            while True:
                if len(finished) == len(precursors.keys()):
                    break
                if base_rt+ms_index == len(rt_index_map) or base_rt+ms_index < 0:
                    break
                rt_window.append(rt_index_map.index[base_rt+ms_index])
                df = self.getScan(rt_index_map.iloc[base_rt+ms_index])
                found = False
                if df is not None:
                    for precursor_label, precursor_shift in precursors.items():
                        if precursor_label in finished:
                            continue
                        measured_precursor = (precursor+precursor_shift)/float(charge)
                        precursor_mass = theor_mass*float(charge)+precursor_shift
                        precursor_mz = precursor_mass/float(charge)
                        theo_mass = theor_mass+precursor_shift/float(charge)
                        # print precursor_label, precursor_mz, theo_mass
                        data[precursor_label]['precursor'] = precursor_mz
                        shift_max = shift_maxes.get(precursor_label)
                        shift_max = precursor+shift_max if shift_max is not None else None
                        envelope = peaks.findEnvelope(df, start_mz=precursor_mass, max_mz=shift_max,
                                                      charge=charge, precursor_ppm=self.precursor_ppm, isotope_ppm=self.isotope_ppm, heavy=True,
                                                      isotope_ppms=self.isotope_ppms if self.fitting_run else None,
                                                      theo_dist=theo_dist, label=precursor_label, skip_isotopes=finished_isotopes[precursor_label],
                                                      last_precursor=last_precursors[delta].get(precursor_label, measured_precursor))
                        # print df.name, precursor_mass, envelope['envelope'], envelope['micro_envelopes'].get(0, {}).get('int'), last_precursors[delta].get(precursor_label, measured_precursor)
                        if not envelope['envelope']:
                            finished.add(precursor_label)
                            continue
                        if 0 in envelope['micro_envelopes'] and envelope['micro_envelopes'][0].get('int'):
                            if ms_index == 0:
                                last_precursors[delta*-1][precursor_label] = envelope['micro_envelopes'][0]['params'][1]
                            last_precursors[delta][precursor_label] = envelope['micro_envelopes'][0]['params'][1]
                        selected = {}

                        for isotope, vals in envelope['micro_envelopes'].iteritems():
                            if isotope in finished_isotopes[precursor_label]:
                                continue
                            if vals.get('int') == 0:
                                finished_isotopes[precursor_label].add(isotope)
                                continue
                            else:
                                found = True
                            selected[precursor_mz+isotope*spacing] = vals.get('int')
                            vals['isotope'] = isotope
                            isotope_labels[precursor_mz+isotope*spacing] = {'label': precursor_label, 'isotope_index': isotope}
                            key = (df.name, precursor_mz+isotope*spacing)
                            isotopes_chosen[key] = {'label': precursor_label, 'isotope_index': isotope, 'amplitude': vals.get('int')}

                        selected = pd.Series(selected, name=df.name).to_frame()
                        if df.name in combined_data.columns:
                            combined_data = combined_data.add(selected, axis='index', fill_value=0)
                        else:
                            combined_data = pd.concat([combined_data, selected], axis=1).fillna(0)
                        del envelope
                        del selected
                if found is False or np.abs(ms_index) > 75:
                    # the 75 check is in case we're in something crazy. We should already have the elution profile of the ion
                    # of interest, else we're in an LC contaminant that will never end.
                    if delta == -1:
                        delta = 1
                        ms_index=0
                        finished = set([])
                        finished_isotopes = {i: set([]) for i in precursors.keys()}
                    elif ms_index >= 2:
                        break
                del df
                ms_index += delta

            if isotope_labels:
                # bookend with zeros, do the right end first because pandas will by default append there
                combined_data[rt_index_map.index[rt_index_map.index.searchsorted(combined_data.columns[-1])+1]] = 0
                combined_data[rt_index_map.index[rt_index_map.index.searchsorted(combined_data.columns[0])-1]] = 0
                combined_data = combined_data[sorted(combined_data.columns)]
                rt_window.sort()
                combined_data = combined_data.sort(axis='index').sort(axis='columns')
                start_rt = rt
                quant_vals = defaultdict(dict)
                isotope_labels = pd.DataFrame(isotope_labels).T

                fig_map = {}

                isotopes_chosen = pd.DataFrame(isotopes_chosen).T
                isotopes_chosen.index.names = ['RT', 'MZ']

                if self.html:
                    # make the figure of our isotopes selected
                    isotope_figure = '{2}_{0}_{1}_{3}_isotopes.png'.format(peptide, ms1, self.filename, scanId)
                    isotopes_chosen['RT'] = isotopes_chosen.index.get_level_values('RT')
                    subplots = len(isotopes_chosen.index.get_level_values('RT').drop_duplicates())
                    fig = plt.figure()
                    all_x = sorted(isotopes_chosen.index.get_level_values('MZ').drop_duplicates())
                    for counter, (index, row) in enumerate(isotopes_chosen.groupby('RT')):
                        ax = fig.add_subplot(subplots, 1, counter+1)
                        colors = 'bmrk'
                        for group, color in zip(precursors.keys(), colors):
                            label_df = row[row['label'] == group]
                            x = label_df['amplitude'].index.get_level_values('MZ')
                            ax.bar(x, label_df['amplitude'], width=spacing/4.0, facecolor=color, align='center')
                        ax.set_xticks(all_x)
                        ax.set_xticklabels([])
                        ax.set_xlim([all_x[0]-0.5, all_x[-1]+0.5])
                    ax.set_xticklabels(['{0:.1f}'.format(i) for i in all_x], rotation=45, ha='right')
                    html_images['isotopes'] = os.path.join(self.html['full'], isotope_figure)
                    ax.get_figure().savefig(html_images['isotopes'], format='png', dpi=100)

                    fname = '{2}_{0}_{1}_{3}_clusters.png'.format(peptide, ms1, self.filename, scanId)
                    fig = plt.figure(figsize=(10, 10))
                    subplot_rows = len(precursors.keys())+1
                    subplot_columns = pd.Series(isotope_labels['label']).value_counts().iloc[0]+1
                    ax = fig.add_subplot(subplot_rows, subplot_columns, 1, projection='3d')
                    X=combined_data.columns.astype(float).values
                    Y=combined_data.index.astype(float).values
                    Z=combined_data.fillna(0).values
                    Xi,Yi = np.meshgrid(X, Y)
                    ax.plot_wireframe(Yi, Xi, Z, cmap=plt.cm.coolwarm)

                combined_peaks = defaultdict(dict)
                plot_index = {}
                fig_nums = defaultdict(list)
                for mz, label in isotope_labels['label'].iteritems():
                    fig_nums[label].append(mz)
                labelx = False
                labely = False
                label_fig_row = {v: i+1 for i,v in enumerate(precursors.keys())}
                for row_num, (index, values) in enumerate(combined_data.iterrows()):
                    quant_label = isotope_labels.loc[index, 'label']
                    if self.html:
                        fig_index = label_fig_row.get(quant_label)*subplot_columns+fig_nums[quant_label].index(index)+2
                        current_row = int(fig_index/subplot_columns+1)
                        if (fig_index-2)%subplot_columns == 0:
                            labely = True
                        if current_row == subplot_rows:
                            labelx = True
                        fig_map[index] = fig_index
                        plot_index[index, quant_label] = row_num

                    xdata = values.index.values.astype(float)
                    ydata = values.fillna(0).values.astype(float)
                    if ydata.any():
                        res, all_peaks = peaks.findAllPeaks(values, filter=True, bigauss_fit=True)
                        # res2, all_peaks2 = peaks.findAllPeaks2(values, filter=True)
                        # if len(res.x) > 4:
                        #     print res
                        # res2, all_peaks2 = peaks.findAllPeaks2(values, filter=True)
                        #print res1.success, res2.success, res1.fun, res2.fun, res1.bic, res2.bic
                        # res, all_peaks = (res, all_peaks) if res.bic < res2.bic else (res2, all_peaks2)
                        mval = ydata.max()
                        rt_means = res.x[1::4]
                        rt_amps = res.x[::4]
                        rt_std = res.x[2::4]
                        rt_std2 = res.x[3::4]
                        combined_peaks[quant_label][index] = [{'mean': i, 'amp': j*mval, 'std': l, 'std2': k, 'peak': m, 'total': values.sum()}
                                                              for i, j, l, k, m in zip(rt_means, rt_amps, rt_std, rt_std2, all_peaks)]
                    if self.html:
                        ax = fig.add_subplot(subplot_rows, subplot_columns, fig_index)
                        ax.plot(xdata, ydata, 'bo-', alpha=0.7)
                        if ydata.any():
                            ax.plot(xdata, peaks.bigauss_ndim(xdata, *res.x)*mval, color='r')
                        # if labely:
                        #     labely = False
                        # else:
                        #     ax.set_yticklabels([])
                        # if labelx:
                        #     labelx = False
                        ax.set_xticks(xdata[::int(len(xdata)/5)+1])
                        ax.set_xticklabels(['{0:.2f} '.format(i) for i in xdata[::int(len(xdata)/5)+1]], rotation=45, ha='right')
                        # else:
                        #     ax.set_xticklabels([])
                # get two most common peak, pick the closest to our RT
                # we may need to add a check for a minimal # of in for max distance from the RT as well here.
                common_peaks = pd.Series([peak['peak'] for i, values in combined_peaks.items() for index, value_peaks in values.iteritems() for peak in value_peaks]).value_counts()
                common_peaks = common_peaks.sort_index()
                tcommon_peaks = common_peaks[common_peaks>=4]

                # combine peaks that are separated by a single scan
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
                peak_info = {i: {'amp': -1, 'var': 0} for i in data.keys()}
                for quant_label, quan_values in combined_peaks.items():
                    for index, values in quan_values.items():
                        if not values:
                            continue
                        rt_values = combined_data.loc[index]
                        xdata = rt_values.index.values.astype(float)
                        ydata = rt_values.fillna(0)
                        # pick the biggest within a rt cutoff of 0.2, otherwise pick closest
                        # closest_rts = sorted([(i, i['amp']) for i in values if np.abs(i['peak']-common_peak) < 0.2], key=operator.itemgetter(1), reverse=True)
                        # if not closest_rts:
                        closest_rts = sorted([(i, np.abs(i['peak']-common_peak)) for i in values], key=operator.itemgetter(1))
                        closest_rts = closest_rts[0][0]
                        # if we move more than a # of ms1 to the dominant peak, update to our known peak
                        gc = 'k'
                        peak_loc = np.where(xdata == closest_rts['peak'])[0][0]
                        mean = closest_rts['mean']
                        amp = closest_rts['amp']
                        mean_diff = np.abs(mean-xdata[common_loc])
                        if quant_label == 'Light':
                            pass
                        if len(xdata) >= 3 and (mean_diff > 0.3 or (np.abs(peak_loc-common_loc) > 2 and mean_diff > 0.2)):
                            mean = common_peak
                            amp = ydata[common_peak]
                            gc = 'g'
                        #var_rat = closest_rts['var']/common_var
                        std = closest_rts['std']
                        std2 = closest_rts['std2']
                        peak_params = (amp,  mean, std, std2)
                        # int_args = (res.x[rt_index]*mval, res.x[rt_index+1], res.x[rt_index+2])
                        left, right = xdata[0]-4*std, xdata[-1]+4*std2
                        xr = np.linspace(left, right, 1000)
                        int_val = integrate.simps(peaks.bigauss(xr, *peak_params), x=xr)
                        isotope_index = isotope_labels.loc[index, 'isotope_index']

                        if int_val and not pd.isnull(int_val) and gc != 'c':
                            try:
                                quant_vals[quant_label][isotope_index] += int_val
                            except KeyError:
                                quant_vals[quant_label][isotope_index] =  int_val
                        if peak_info.get(quant_label, {}).get('amp', -1) < amp:
                            peak_info[quant_label].update({'amp': amp, 'std': std, 'std2': std2, 'mean_diff': mean_diff})
                        if self.html:
                            ax = fig.add_subplot(subplot_rows, subplot_columns, fig_map.get(index))
                            ax.plot(xdata, peaks.bigauss(xdata, *peak_params), '{}o-'.format(gc), alpha=0.7)
                            ax.plot([start_rt, start_rt], ax.get_ylim(),'k-')
                            ax.set_ylim(0,combined_data.max().max())
                            # ax.set_ylim(0, amp)

                for silac_label1 in data.keys():
                    qv1 = quant_vals.get(silac_label1)
                    for silac_label2 in data.keys():
                        if silac_label1 == silac_label2:
                            continue
                        qv2 = quant_vals.get(silac_label2)
                        ratio = 'NA'
                        if qv1 is not None and qv2 is not None:
                            if self.mono:
                                common_isotopes = set(qv1.keys()).intersection(qv2.keys())
                            else:
                                common_isotopes = set(qv1.keys()).union(qv2.keys())
                            quant1 = sum([qv1.get(i, 0) for i in common_isotopes])
                            quant2 = sum([qv2.get(i, 0) for i in common_isotopes])
                            ratio = quant1/quant2 if quant1 and quant2 else 'NA'
                        result_dict.update({'{}_{}_ratio'.format(silac_label1, silac_label2): ratio})

                if self.html:
                    plt.tight_layout()
                    ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=100)
                    html_images['clusters'] = os.path.join(self.html['rel'], fname)
                if self.debug or self.html:
                    plt.close('all')
                result_dict.update({'html_info': html_images})
                for silac_label, silac_data in data.iteritems():
                    result_dict.update({
                        '{}_intensity'.format(silac_label): sum(quant_vals[silac_label].values()),
                        '{}_isotopes'.format(silac_label): sum(isotopes_chosen['label'] == silac_label),
                        '{}_rt_width'.format(silac_label): peak_info.get(silac_label, {}).get('var', 'NA'),
                        '{}_mean_diff'.format(silac_label): peak_info.get(silac_label, {}).get('mean_diff', 'NA'),
                    })
                del combined_peaks
            for silac_label, silac_data in data.iteritems():
                result_dict.update({
                    '{}_precursor'.format(silac_label): silac_data['precursor']
                })
            self.results.put(result_dict)
            del result_dict
            del data
            del combined_data
            del isotopes_chosen
        except:
            sys.stderr.write('ERROR ON {0}\n{1}\n{2}\n'.format(ms1, scan_info, traceback.format_exc()))
            return

    def run(self):
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

        for params in iter(self.queue.get, None):
            self.run_thing(params, scan_rt_map, rt_index_map, ms2_ms1_map, silac_shifts)
        self.results.put(None)



def main():
    args = parser.parse_args()
    source = args.processed.name
    threads = args.p
    skip = args.skip
    out = args.out
    raw_file = args.mzml
    html = args.html
    resume = args.resume
    manager = Manager()
    calc_stats = not args.disable_stats

    scan_filemap = {}

    if isinstance(raw_file, list):
        nfunc = lambda i: (os.path.splitext(os.path.split(i.name)[1])[0], os.path.abspath(i.name)) if hasattr(i, 'name') else (os.path.splitext(os.path.split(i)[1])[0], os.path.abspath(i))
        scan_filemap = dict([nfunc(i) for i in raw_file])
    else:
        if args.raw_dir:
            raw_file = args.raw_dir
        else:
            raw_file = raw_file.name if raw_file else raw_file
            if not raw_file:
                raw_file = os.path.abspath(os.path.split(source)[0])
        if os.path.isdir(raw_file):
            scan_filemap = dict([(os.path.splitext(i)[0], os.path.abspath(os.path.join(raw_file, i))) for i in os.listdir(raw_file) if i.lower().endswith('mzml')])
        else:
            scan_filemap[os.path.splitext(os.path.split(raw_file)[1])[0]] = os.path.abspath(raw_file)

    mass_scans = manager.dict()
    msf_scans = manager.dict()
    raw_files = {}

    if args.label_scheme:
        label_info = pd.read_table(args.label_scheme.name, sep='\t', header=None, dtype='str')
        try:
            label_info.columns = ['Label', 'AA', 'Mass', 'UserName']
            name_mapping = dict([(v['Label'],v['UserName']) for i,v in label_info.iterrows()])
        except ValueError:
            label_info.columns = ['Label', 'AA', 'Mass']
            name_mapping = dict([(v['Label'],v['Label']) for i,v in label_info.iterrows()])
        silac_labels = {'Light': {0: set([])}}
        for group_name, group_info in label_info.groupby('Label'):
            masses = {}
            label_name = name_mapping.get(group_name, group_name)
            for mass, mass_aas in group_info.groupby('Mass'):
                mass_val = float(mass)
                mass_list = mass_aas['AA'].values.tolist()
                try:
                    masses[mass_val].add(mass_list)
                except KeyError:
                    masses[mass_val] = set(mass_list)
            silac_labels.update({label_name: masses})

    if args.tsv:
        dat = pd.read_table(source, sep='\t')

        sys.stderr.write('Loading Scans:\n')
        sample = args.sample

        peptide_col = args.peptide_col
        ms2_col = args.ms2
        mz_col = args.mz
        rt_col = args.rt
        charge_col = args.charge
        file_col = args.source
        label_col = args.label

        for index, row in enumerate(dat.iterrows()):
            if index%1000 == 0:
                sys.stderr.write('.')
            row_index, i = row
            peptide = i[peptide_col].strip()
            if args.peptide and args.peptide.lower() != peptide.lower():
                continue
            if not args.peptide and (sample != 1.0 and random.random() > sample):
                continue
            specId = str(i[ms2_col])
            if specId in mass_scans:
                continue
            fname = i[file_col]
            if fname not in scan_filemap:
                fname = os.path.split(fname)[1]
                if fname not in scan_filemap:
                    if skip:
                        continue
                    sys.stderr.write('{0} not found in filemap. Filemap is {1}.'.format(fname, scan_filemap))
                    return 1
            charge = float(i[charge_col])
            if not charge or not peptide:
                continue
            d = {'id': specId, 'theor_mass': i[mz_col], 'peptide': peptide, 'mod_peptide': peptide, 'rt': i[rt_col], 'charge': charge, 'modifications': None, 'label': name_mapping.get(i[label_col]) if label_col in i else None}
            mass_scans[specId] = d
            msf_scans[(specId, peptide)] = {'file': fname, 'ms1': None, 'rawId': specId, 'charge': charge, 'mass': i[mz_col], 'peptide': i[peptide_col]}
            try:
                raw_files[i[file_col]].append((float(i[rt_col]), d))
            except:
                raw_files[i[file_col]] = [(float(i[rt_col]), d)]

        # sort by RT so we can minimize our memory footprint by throwing away scans we no longer need
        for i in raw_files.keys():
            raw_files[i].sort(key=operator.itemgetter(0))
            raw_files[i] = [j[1] for j in raw_files[i]]

    else:
        results = GuessIterator(source, full=True, store=False, peptide=args.peptide)

        if not args.label_scheme:
            silac_labels = {'Light': {0: set([])}}
            silac_labels.update(results.getSILACLabels())

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
            mass_scans[specId] = {'id': specId, 'theor_mass': i.getTheorMass(), 'peptide': i.peptide, 'mod_peptide': i.modifiedPeptide, 'rt': i.rt, 'charge': i.charge, 'modifications': i.getModifications()}
            del i

        for index, scan in enumerate(results.getScans(modifications=False)):
            if index%1000 == 0:
                sys.stderr.write('.')
            if scan is None:
                break
            msf_scans[(str(scan.rawId), scan.peptide)] = {'file': scan.file, 'ms1': scan.ms1_scan.title, 'rawId': scan.rawId, 'charge': scan.charge, 'mass': scan.mass, 'peptide': scan.peptide}
            del scan
        sys.stderr.write('\nScans loaded.\n')

        raw_files = {}

        sys.stderr.write('Processing retention times.\n')
        for specId in mass_scans.keys():
            i = mass_scans[specId]
            if index%1000 == 0:
                sys.stderr.write('.')
            scan = msf_scans.get((specId, i.get('peptide')))
            if scan is None:
                continue
            fname = os.path.splitext(scan.get('file'))[0]
            if fname not in scan_filemap:
                fname = os.path.split(fname)[1]
                if fname not in scan_filemap:
                    if skip:
                        continue
                    sys.stderr.write('{0} not found in filemap. Filemap is {1}.'.format(fname, scan_filemap))
                    return 1
            i['file'] = fname
            try:
                raw_files[fname].append((float(i.get('rt')), i))
            except KeyError:
                raw_files[fname] = [(float(i.get('rt')), i)]


        # sort by RT so we can minimize our memory footprint by throwing away scans we no longer need
        for i in raw_files.keys():
            raw_files[i].sort(key=operator.itemgetter(0))
            raw_files[i] = [j[1] for j in raw_files[i]]

        raw_files = manager.dict(raw_files)

        sys.stderr.write('\nRetention times loaded.\n')

    raw_files = manager.dict(raw_files)

    labels = silac_labels.keys()
    for silac_label in labels:
        RESULT_ORDER.extend([('{}_intensity'.format(silac_label), '{} Intensity'.format(silac_label)),
                             ('{}_precursor'.format(silac_label), '{} Precursor'.format(silac_label)),
                             ('{}_rt_width'.format(silac_label), '{} RT Width'.format(silac_label)),
                             ('{}_isotopes'.format(silac_label), '{} Isotopes Found'.format(silac_label)),
                             ('{}_mean_diff'.format(silac_label), '{} Mean Offset'.format(silac_label))
                             ])
        for silac_label2 in labels:
            if silac_label != silac_label2:
                RESULT_ORDER.extend([('{}_{}_ratio'.format(silac_label, silac_label2), '{}/{}'.format(silac_label, silac_label2)),
                                     ])
                if calc_stats:
                    RESULT_ORDER.extend([('{}_{}_confidence'.format(silac_label, silac_label2), '{}/{} Confidence'.format(silac_label, silac_label2)),
                                         ])

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
                elif v == 'Retention Time' and 'isotopes' in images:
                    out.append("""<td data-toggle="popover" data-trigger="click" data-content='<img src="{0}">'>{1}</td>""".format(images.get('isotopes'), i))
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

    for filename in raw_files.keys():
        raw_scans = raw_files[filename]
        filepath = scan_filemap[filename]
        if not len(raw_scans):
            continue
        in_queue = Queue()
        result_queue = Queue()
        reader_in = Queue()
        reader_outs = {}
        for i in xrange(threads):
            reader_outs[i] = Queue()

        ms2_ms1_map = manager.dict()
        scan_rt_map = manager.dict()
        rt_index_map = manager.dict()

        raw = GuessIterator(filepath, full=False, store=False, ms_filter=None if args.tsv else 1)
        sys.stderr.write('Processing raw file.\n')
        for index, i in enumerate(raw):
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

        scans_loaded = Array(ctypes.c_bool, scan_id)
        last_scan_accessed = Array(ctypes.c_uint32, threads)

        reader = Reader(reader_in, reader_outs, scans_loaded=scans_loaded, raw_file=raw, ms2_ms1_map=ms2_ms1_map, scan_rt_map=scan_rt_map,
                        rt_index_map=rt_index_map)
        reader.start()

        # # establish our test parameters first
        #
        # for i in xrange(threads):
        #     worker = Worker(queue=in_queue, results=result_queue, raw_name=filepath, silac_labels=silac_labels,
        #                     debug=args.debug, html=html, mono=not args.spread, precursor_ppm=args.precursor_ppm,
        #                     isotope_ppm=args.isotope_ppm*5,
        #                     reader_in=reader_in, reader_out=reader_outs[i], ms2_ms1_map=ms2_ms1_map,
        #                     scan_rt_map=scan_rt_map, rt_index_map=rt_index_map, scans_loaded=scans_loaded,
        #                     last_scan_accessed=last_scan_accessed, thread=i, fitting_run=True)
        #     workers.append(worker)
        #     worker.start()
        #
        # try:
        #     sampling = random.sample(range(len(raw_scans)), 100)
        # except ValueError:
        #     sampling = range(len(raw_scans))
        #
        # for scan_index in sampling:
        #     v = raw_scans[scan_index]
        #     if resume:
        #         key = (v.get('file'), v.get('peptide'), v.get('modifications'), v.get('charge'), v.get('id'))
        #         if key in skip_map:
        #             completed += 1
        #             continue
        #     # if scan_index % 100 != 0:
        #     #     continue
        #     params = {}
        #     params['scan_info'] = v
        #     ms2_scan = msf_scans.get((v.get('id'), v.get('peptide')))
        #     params['ms2_info'] = ms2_scan
        #     in_queue.put(params)
        # # worker.run()
        #
        # sys.stderr.write('{0} processed and placed into queue.\n'.format(filename))
        #
        # # kill the workers
        # [in_queue.put(None) for i in xrange(threads)]
        # fit_results = []
        # while workers or result is not None:
        #     try:
        #         result = result_queue.get(timeout=0.1)
        #     except Empty:
        #         # kill expired workers
        #         result = None
        #         to_del = []
        #         for i, v in enumerate(workers):
        #             if not v.is_alive():
        #                 v.terminate()
        #                 to_del.append(i)
        #         for i in sorted(to_del, reverse=True):
        #             del workers[i]
        #     if result is not None:
        #         fit_results.append(result)
        #
        # # do it for real now
        # isotope_ppms = {}
        # for res in fit_results:
        #     for isotope, isotope_params in res.items():
        #         error = isotope_params.get('error')
        #         if isotope_params.get('int') and error and not pd.isnull(error):
        #             try:
        #                 isotope_ppms[isotope].append(error)
        #             except KeyError:
        #                 isotope_ppms[isotope] = [error]
        # for i in isotope_ppms.keys():
        #     ppm_info = pd.Series(isotope_ppms[i])
        #     med_devs = np.abs(np.log2(ppm_info/ppm_info.median()))
        #     # remove entries that deviate from the median by > 2
        #     ppm_info = ppm_info[med_devs<1]
        #     isotope_ppms[i] = args.isotope_ppm if i != 0 else args.precursor_ppm#(ppm_info.mean()+ppm_info.std()*2)*1e6
        # print isotope_ppms
        for i in xrange(threads):
            worker = Worker(queue=in_queue, results=result_queue, raw_name=filepath, silac_labels=silac_labels,
                            debug=args.debug, html=html, mono=not args.spread, precursor_ppm=args.precursor_ppm,
                            isotope_ppm=args.isotope_ppm, isotope_ppms=None,
                            reader_in=reader_in, reader_out=reader_outs[i], ms2_ms1_map=ms2_ms1_map,
                            scan_rt_map=scan_rt_map, rt_index_map=rt_index_map, scans_loaded=scans_loaded,
                            last_scan_accessed=last_scan_accessed, thread=i)
            workers.append(worker)
            worker.start()

        for scan_index, v in enumerate(raw_scans):
            if resume:
                key = (v.get('file'), v.get('peptide'), v.get('modifications'), v.get('charge'), v.get('id'))
                if key in skip_map:
                    completed += 1
                    continue
            params = {}
            params['scan_info'] = v
            ms2_scan = msf_scans.get((v.get('id'), v.get('peptide')))
            params['ms2_info'] = ms2_scan
            in_queue.put(params)

        sys.stderr.write('{0} processed and placed into queue.\n'.format(filename))

        # kill the workers
        [in_queue.put(None) for i in xrange(threads)]
        all_results = []
        html_results = []
        while workers or result is not None:
            try:
                result = result_queue.get(timeout=0.1)
            except Empty:
                # kill expired workers
                result = None
                to_del = []
                for i, v in enumerate(workers):
                    if not v.is_alive():
                        v.terminate()
                        to_del.append(i)
                for i in sorted(to_del, reverse=True):
                    del workers[i]
            if result is not None:
                completed += 1
                if completed % 10 == 0:
                    sys.stderr.write('\r{0:2.2f}% Completed'.format(completed/scan_count*100))
                    sys.stderr.flush()
                    # sys.stderr.write('Processed {0} of {1}...\n'.format(completed, scan_count))
                res_list = [filename]+[result.get(i[0], 'NA') for i in RESULT_ORDER]
                if calc_stats:
                    all_results.append(res_list)
                    html_results.append(result.get('html_info', {}))
                res = '{0}\n'.format('\t'.join(map(str, res_list)))
                out.write(res)
                out.flush()
                if html and not calc_stats:
                    html_out.write(table_rows([{'table': res.strip(), 'images': result.get('html_info', {})}]))
                    html_out.flush()
        reader_in.put(None)
        del ms2_ms1_map
        del scan_rt_map
        del rt_index_map

    if calc_stats:
        from scipy import stats
        data = pd.DataFrame.from_records(all_results, columns=[i for i in headers if i != 'Confidence'])
        data = data.replace('NA', np.nan)
        for silac_label1 in silac_labels.keys():
            label1_log = 'L{}'.format(silac_label1)
            label1_logp = 'L{}_p'.format(silac_label1)
            label1_int = '{} Intensity'.format(silac_label1)
            label1_hif = '{} Isotopes Found'.format(silac_label1)
            label1_hifp = '{} Isotopes Found p'.format(silac_label1)
            for silac_label2 in silac_labels.keys():
                if silac_label1 == silac_label2:
                    continue
                label2_log = 'L{}'.format(silac_label2)
                label2_logp = 'L{}_p'.format(silac_label2)
                label2_int = '{} Intensity'.format(silac_label2)
                label2_hif = '{} Isotopes Found'.format(silac_label2)
                label2_hifp = '{} Isotopes Found p'.format(silac_label2)

                mixed = '{}/{}'.format(silac_label1, silac_label2)
                mixed_p = '{}/{}_p'.format(silac_label1, silac_label2)
                mixed_mean = '{}_Mean_Diff'.format(mixed)
                mixed_mean_p = '{}_Mean_Diff_p'.format(mixed)
                mixed_rt_diff = '{}_RT_Diff'.format(mixed)
                mixed_rt_diff_p = '{}_p'.format(mixed_rt_diff)
                mixed_isotope_diff = '{}_Isotope_Diff'.format(mixed)
                mixed_isotope_diff_p = '{}_Isotope_Diff_p'.format(mixed)

                data[label1_log] = np.log(data[label1_int].astype(float)+1)
                data[label1_logp] = stats.norm.cdf((data[label1_log] - data[data[label1_log]>0][label1_log].mean())/data[data[label1_log]>0][label1_log].std(ddof=0))
                data[label2_log] = np.log(data[label2_int].astype(float)+1)
                data[label2_logp] = stats.norm.cdf((data[label2_log] - data[data[label2_log]>0][label2_log].mean())/data[data[label2_log]>0][label2_log].std(ddof=0))

                nz = data[(data[label2_log] > 0) & (data[label1_log] > 0)]
                mu = pd.Series(np.ravel(nz.loc[:,(label2_log, label1_log)])).mean()
                std = pd.Series(np.ravel(nz.loc[:,(label2_log, label1_log)])).std()

                data[mixed_p] = stats.norm.cdf((data.loc[:,(label2_log, label1_log)].mean(axis=1)-mu)/std)
                data[mixed_rt_diff] = np.log2(np.abs(data['{} RT Width'.format(silac_label2)].astype(float)-data['{} RT Width'.format(silac_label1)].astype(float)))
                data[mixed_mean] = np.abs(data['{} Mean Offset'.format(silac_label1)].astype(float)-data['{} Mean Offset'.format(silac_label1)].astype(float))
                data[mixed_rt_diff] = data[mixed_rt_diff].replace([np.inf, -np.inf], np.nan)
                data[mixed_rt_diff_p] = stats.norm.cdf((data[mixed_rt_diff] - data[mixed_rt_diff].mean())/data[mixed_rt_diff].std(ddof=0))
                data[mixed_mean_p] = stats.norm.cdf((data[mixed_mean] - data[mixed_mean].mean())/data[mixed_mean].std(ddof=0))

                data[label2_hifp] = np.log2(data[label2_hif]).replace([np.inf, -np.inf], np.nan)
                data[label2_hifp] = stats.norm.cdf((data[label2_hifp]-data[label2_hifp].median())/data[label2_hifp].std())

                data[label1_hifp] = np.log2(data[label1_hif]).replace([np.inf, -np.inf], np.nan)
                data[label1_hifp] = stats.norm.cdf((data[label1_hifp]-data[label1_hifp].median())/data[label2_hifp].std())

                data[mixed_isotope_diff] = np.log2(data[label2_hif]/data[label1_hif]).replace([np.inf, -np.inf], np.nan)
                data[mixed_isotope_diff_p] = stats.norm.cdf((data[mixed_isotope_diff] - data[mixed_isotope_diff].median())/data[mixed_isotope_diff].std(ddof=0))

                # confidence assessment
                mixed_confidence = '{}/{} Confidence'.format(silac_label1, silac_label2)
                data[mixed_confidence] = 10
                data.loc[(data[mixed_mean_p] > 0.90), mixed_confidence] -= 1
                data.loc[(data[mixed_rt_diff_p] > 0.90), mixed_confidence] -= 1
                data.loc[(data[label2_logp] < 0.10), mixed_confidence] -= 0.5
                data.loc[(data[label1_logp] < 0.10), mixed_confidence] -= 0.5
                data.loc[(data[mixed_p] < 0.10), mixed_confidence] -= 0.5
                data.loc[((data[mixed_isotope_diff_p] > 0.90) | (data[mixed_isotope_diff_p] < 0.10)), mixed_confidence] -= 1
                data.loc[(data[label2_hifp] < 0.10), mixed_confidence] -= 1
                data.loc[(data[label2_hifp] < 0.10), mixed_confidence] -= 1

        data.to_csv(out.name, sep=str('\t'), index=None)
        if html:
            for index, (row_index, row) in enumerate(data.iterrows()):
                res = '\t'.join(row.astype(str))
                html_out.write(table_rows([{'table': res.strip(), 'images': html_results[index]}]))
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
