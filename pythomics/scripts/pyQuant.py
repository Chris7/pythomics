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
from scipy.stats import signaltonoise

from matplotlib.backends.backend_pdf import PdfPages
from cStringIO import StringIO
import argparse
import urllib
import base64
import scipy.stats

from pythomics.templates import CustomParser
from pythomics.proteomics.parsers import GuessIterator
from pythomics.proteomics import peaks, config

RESULT_ORDER = (('peptide', 'Peptide'), ('heavy', 'Heavy'), ('light', 'Light'), ('h_l', 'Heavy/Light'), ('modifications', 'Modifications'),
                ('charge', 'Charge'), ('ms1', 'MS1 Spectrum ID'), ('scan', 'MS2 Spectrum ID'), ('light_precursor', 'Light Precursor M/Z'),
                ('heavy_precursor', 'Heavy Precursor M/Z'),
                ('light_clusters', 'Light Peaks Identified'), ('heavy_clusters', 'Heavy Peaks Identified'), ('light_rt_dev', 'Light Retention Time Deviation'),
                ('heavy_rt_dev', 'Heavy Retention Time Deviation'), ('light_residuals', 'Light Intensity Residuals of Fit'),
                ('heavy_residuals', 'Heavy Intensity Residuals of Fit'), ('cluster_ratio_deviation', 'Deviation of Heavy Clusters From Light'),
                ('light_snr', 'Light SNR'), ('heavy_snr', 'Heavy SNR'), ('heavy_ks', 'Heavy KS'), ('light_ks', 'Light KS'))


parser = CustomParser(description=description)
parser.add_processed_ms()
parser.add_argument('--mzml', help="The mzML files for the raw data or directory containing them. Defaults to the directory of the processed file.", type=argparse.FileType('r'), nargs='?')
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
parser.add_argument('-o', '--out', nargs='?', help='Stuff', type=str)

class Worker(Process):
    def __init__(self, queue=None, results=None, precision=6, raw_name=None, heavy_masses=None,
                 temp=False, sparse=False, debug=False, html=False, mono=False):
        super(Worker, self).__init__()
        self.precision = precision
        self.sparse = sparse
        self.queue=queue
        self.results = results
        if heavy_masses is None:
            self.heavy_masses = {'R': 10.008269, 'K': 8.0142}
        else:
            self.heavy_masses = heavy_masses
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
        sys.stderr.write('Thread processing peptide file.\n')
        for index, i in enumerate(self.raw):
            if index % 100 == 0:
                sys.stderr.write('.')
            if i is None:
                break
            if i.ms_level != 1:
                continue
            scan_title = int(i.title)
            rt = i.rt
            scan_rt_map[scan_title] = rt
            rt_index_map[rt] = scan_title
        sys.stderr.write('Scans Loaded\n')
        rt_index_map = pd.Series(rt_index_map)

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
                ms1 = int(scan.get('ms1'))
                rt = scan_rt_map.get(ms1)
                if not rt:
                    return
                mods = scan_info.get('modifications')
                charge = int(scan.get('charge'))
                precursor = scan.get('mass')+config.HYDROGEN*(charge-1)

                shift = scan_info['mass_shift']
                peptide = scan_info['peptide']
                mod_pep = scan_info['mod_peptide']
                if self.debug:
                    sys.stderr.write('on ms {0} {1} {2} {3}\n'.format(ms1, rt, precursor, scan_info))
                if shift == 0:
                    count = 0
                    for label_mass, label_masses in self.heavy_masses.iteritems():
                        labels = [label_mass for mod_aa in mod_pep if mod_aa in label_masses]
                        count += len(labels)
                        shift += sum(labels)
                    if count != 1:
                        return
                    light_precursor = precursor
                    heavy_precursor = precursor+shift
                else:
                    light_precursor = precursor-shift
                    heavy_precursor = precursor

                quant = light_precursor != heavy_precursor

                light_data, heavy_data = None, None
                light_rt_data = pd.Series()
                heavy_rt_data = pd.Series()
                rt_width = int(len(rt_index_map[rt-self.rt_width:rt+self.rt_width])/2)
                # sys.stderr.write('{}\n'.format(rt_width))
                base_rt = rt_index_map.index.searchsorted(rt)
                light_df = pd.DataFrame()
                heavy_df = pd.DataFrame()
                light_isotopes = {}
                heavy_isotopes = {}
                peaks_found = {'light': OrderedDict(), 'heavy': OrderedDict()}
                for ms_index in xrange(-3,rt_width+1):
                    df = self.getScan(rt_index_map.iloc[base_rt+ms_index])
                    light = peaks.findEnvelope(df, start_mz=light_precursor, max_mz=heavy_precursor if quant else None, charge=charge, ppm=5, heavy=False)
                    if quant:
                        heavy = peaks.findEnvelope(df, start_mz=heavy_precursor, charge=charge, ppm=10, heavy=True)
                        # look for Proline/Glutamate/Glutamines
                        # if 'P' in peptide or 'E' in peptide or 'Q' in peptide:
                        #     heavy2 = peaks.findEnvelope(df, start_mz=heavy_precursor, isotope_offset=6, charge=charge, ppm=10, heavy=True)
                        #     print rt_index_map.iloc[base_rt+ms_index], heavy_precursor
                        #     print heavy
                        #     print heavy2
                    if not light['envelope'] and (quant and not heavy['envelope']):
                        if ms_index < 0:
                            continue
                        if not peaks_found['light'] and peaks_found['heavy']:
                            return
                        # sys.stderr.write('break on {}\n'.format(ms_index))
                        break
                    for i in ['envelope', 'micro_envelopes', 'ppms']:
                        for isotope, vals in light.get(i, {}).iteritems():
                            val_dict = {'info': vals, 'df': df}
                            try:
                                peaks_found['light'][isotope][i].append(val_dict)
                            except KeyError:
                                try:
                                    peaks_found['light'][isotope][i] = [val_dict]
                                except KeyError:
                                    try:
                                        peaks_found['light'][isotope] = {i: [val_dict] }
                                    except KeyError:
                                        peaks_found['light'] = {isotope: {i: [val_dict]}}
                        for isotope, vals in heavy.get(i, {}).iteritems():
                            val_dict = {'info': vals, 'df': df}
                            try:
                                peaks_found['heavy'][isotope][i].append(val_dict)
                            except KeyError:
                                try:
                                    peaks_found['heavy'][isotope][i] = [val_dict]
                                except KeyError:
                                    try:
                                        peaks_found['heavy'][isotope] = {i: [val_dict] }
                                    except KeyError:
                                        peaks_found['heavy'] = {isotope: {i: [val_dict]}}
                zip_func = zip if self.mono else izip_longest
                light_keys = peaks_found.get('light', {}).keys()
                heavy_keys = peaks_found.get('heavy', {}).keys()
                if self.mono:
                    light_keys = sorted(list(set(light_keys).intersection(heavy_keys)))
                    heavy_keys = light_keys
                light_clusters, heavy_clusters = len(light_keys), len(heavy_keys)
                for light_isotope, heavy_isotope in zip_func(light_keys, heavy_keys):
                    light_info = peaks_found.get('light', {}).get(light_isotope, {})
                    heavy_info = peaks_found.get('heavy', {}).get(heavy_isotope, {})
                    for light_envelope, light_micro, heavy_envelope, heavy_micro in izip_longest(
                            light_info.get('envelope', []),
                            light_info.get('micro_envelopes', []),
                            heavy_info.get('envelope', []),
                            heavy_info.get('micro_envelopes', [])):
                        if light_envelope:
                            iso_light_df = light_envelope['df']
                            light_rt = iso_light_df.name
                        if heavy_envelope:
                            iso_heavy_df = heavy_envelope['df'] if heavy_envelope else None
                            heavy_rt = iso_heavy_df.name

                        if light_micro and light_micro.get('info'):
                            start, end = light_micro['info'][0], light_micro['info'][1]
                            selector = iso_light_df.index[range(start,end+1)]
                            temp_series = pd.Series(iso_light_df[selector], index=selector)
                            try:
                                series_index = light_isotopes[light_isotope]
                            except KeyError:
                                series_index = temp_series.idxmax()
                                light_isotopes[light_isotope] = series_index
                            # temp_int = temp_series.sum()#integrate.simps(temp_series.values, temp_series.index.values)
                            temp_int = integrate.simps(temp_series.values)#, temp_series.index.values)
                            if series_index in light_df.index:
                                if light_rt in light_df and not pd.isnull(light_df.loc[series_index, light_rt]):
                                    light_df.loc[series_index, light_rt] += temp_int
                                else:
                                    light_df.loc[series_index, light_rt] = temp_int
                            else:
                                light_df = light_df.append(pd.DataFrame(temp_int, columns=[light_rt], index=[series_index]))

                        if heavy_micro and heavy_micro.get('info', ()):
                            start, end = heavy_micro['info'][0], heavy_micro['info'][1]
                            selector = iso_heavy_df.index[range(start, end+1)]
                            temp_series = pd.Series(iso_heavy_df[selector], index=selector)
                            try:
                                series_index = heavy_isotopes[heavy_isotope]
                            except KeyError:
                                series_index = temp_series.idxmax()
                                heavy_isotopes[heavy_isotope] = series_index
                            # temp_int = temp_series.sum()#integrate.simps(temp_series.values, temp_series.index.values)
                            temp_int = integrate.simps(temp_series.values)#, temp_series.index.values)
                            if series_index in heavy_df.index:
                                if heavy_rt in heavy_df and not pd.isnull(heavy_df.loc[series_index, heavy_rt]):
                                    heavy_df.loc[series_index, heavy_rt] += temp_int
                                else:
                                    heavy_df.loc[series_index, heavy_rt] = temp_int
                            else:
                                heavy_df = heavy_df.append(pd.DataFrame(temp_int, columns=[heavy_rt], index=[series_index]))

                if self.debug or self.html:
                    if self.debug:
                        pdf = PdfPages('{2}_{0}_{1}_fix.pdf'.format(peptide, ms1, self.raw_name))
                        plt.figure()
                light_rt_data = light_df.fillna(0).sum(axis=0).sort_index()
                light_int = integrate.trapz(light_rt_data.values, light_rt_data.index.values) if not light_rt_data.empty else 0#light_df.fillna(0).sum(axis=1).sum()
                # if light_int < 0:
                #     1
                # light_int = light_df.fillna(0).sum(axis=1).sum()
                # light_rt_data.plot(color='b', title=str('abc')).get_figure().savefig('/home/chris/test.png', format='png', dpi=50)
                if self.debug or self.html:
                    plt.figure()
                    if not light_rt_data.empty:
                        ax = light_rt_data.plot(color='b', title=str(light_int))
                        # X = light_rt_data.index.values
                        # pd.Series(dict(zip(X, fit(X)+opt_offset))).plot(color='r')
                        # plt.plot([rt, rt], [0, (amp-opt_offset)], color='k')
                        # plt.plot([mu, mu], [0, (amp-opt_offset)], color='r')
                        if self.debug:
                            pdf.savefig(figure=ax.get_figure())
                        if self.html:
                            fname = '{2}_{0}_{1}_{3}_light_rt.png'.format(peptide, ms1, self.filename, scanId)
                            ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                            html_images['light_rt'] = os.path.join(self.html['rel'], fname)
                        ax.get_figure().clf()

                if not heavy_df.empty:
                    heavy_rt_data = heavy_df.fillna(0).sum(axis=0).sort_index()
                    heavy_int = integrate.trapz(heavy_rt_data.values, heavy_rt_data.index.values) if not heavy_rt_data.empty else 0##integrate.simps(heavy_rt_data.values, heavy_rt_data.index.values)
                    # heavy_int = heavy_df.fillna(0).sum(axis=1).sum()
                    light_rats = []#[df.iloc[v]/df.iloc[light['envelope'][i]] for i,v in enumerate(light['envelope'][1:])]
                    heavy_rats = []#[df.iloc[v]/df.iloc[heavy['envelope'][i]] for i,v in enumerate(heavy['envelope'][1:])]
                    cluster_rat = sum([(i-v) for i, v in zip(light_rats, heavy_rats[:-1])])
                    if self.debug or self.html:
                        plt.figure()
                        if not heavy_rt_data.empty:
                            ax = heavy_rt_data.plot(color='b', title=str(heavy_int))
                            # X = heavy_rt_data.index.values
                            # pd.Series(dict(zip(X, fit(X)+opt_offset))).plot(color='r')
                            # plt.plot([rt, rt], [0, (amp-opt_offset)], color='k')
                            # plt.plot([mu, mu], [0, (amp-opt_offset)], color='r')
                            if self.debug:
                                pdf.savefig(figure=ax.get_figure())
                            if self.html:
                                fname = '{2}_{0}_{1}_{3}_heavy_rt.png'.format(peptide, ms1, self.filename, scanId)
                                ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                                html_images['heavy_rt'] = os.path.join(self.html['rel'], fname)
                            ax.get_figure().clf()
                            plt.close('all')

                if self.html:
                    fname = '{2}_{0}_{1}_{3}_clusters.png'.format(peptide, ms1, self.filename, scanId)
                    x, y = [],[]
                    cseries = light_df.fillna(0).sum(axis=1)
                    for index, val in zip(cseries.index, cseries):
                        x += [index,index,index]
                        y += [0,val,0]
                    ax = plt.plot(x,y, 'b')
                    x, y = [],[]
                    cseries = heavy_df.fillna(0).sum(axis=1)
                    for index, val in zip(cseries.index, cseries):
                        x += [index,index,index]
                        y += [0,val,0]
                    ax = plt.plot(x,y, 'r')
                    ax[0].get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                    html_images['clusters'] = os.path.join(self.html['rel'], fname)
                if self.debug:
                    pdf.close()

                quant = not heavy_df.empty

                self.results.put({'peptide': scan_info.get('mod_peptide', peptide), 'heavy': heavy_int if quant else 'NA', 'light': light_int,
                                  'scan': scanId, 'light_clusters': light_clusters,
                                  'h_l': heavy_int/light_int if quant and light_int else 'NA',
                                  'heavy_clusters': heavy_clusters if quant else 'NA', 'ms1': ms1,
                                  'charge': charge, 'modifications': mods, 'light_precursor': light_precursor,
                                  'heavy_precursor': heavy_precursor if quant else 'NA',
                                  'cluster_ratio_deviation': cluster_rat if quant else 'NA',
                                  'html_info': html_images,
                                  'light_snr': signaltonoise(light_data.values) if light_data is not None else 'NA',
                                  'heavy_snr': signaltonoise(heavy_data.values) if quant and heavy_data is not None else 'NA',})
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
    hdf = args.mzml.name if args.mzml else args.mzml
    save_files = args.no_temp
    sparse = args.no_sparse
    html = args.html
    resume = args.resume


    hdf_filemap = {}
    if not hdf:
        hdf = os.path.abspath(os.path.split(source)[0])
    if os.path.isdir(hdf):
        hdf_filemap = dict([(os.path.splitext(i)[0], os.path.abspath(os.path.join(hdf, i))) for i in os.listdir(hdf) if i.lower().endswith('mzml')])
    else:
        hdf_filemap[os.path.splitext(os.path.split(hdf)[1])[0]] = os.path.abspath(hdf)

    results = GuessIterator(source, full=True, store=False)

    mass_scans = {}
    msf_scan_index = {}

    heavy_masses = results.getSILACLabels()

    sys.stderr.write('Loading Scans:\n')
    for index, i in enumerate(results):
        if index%1000 == 0:
            sys.stderr.write('.')
        if args.peptide and args.peptide.lower() != i.peptide.lower():
            continue
        specId = str(i.rawId)
        if specId in mass_scans:
            continue
        shift = 0.0
        for mod in i.mods:
            if 'Label' in mod[3]:
                shift += float(mod[2])
        mass_scans[specId] = {'id': specId, 'peptide': i.peptide, 'mod_peptide': i.modifiedPeptide, 'mass_shift': shift, 'rt': i.rt, 'charge': i.charge, 'modifications': i.getModifications()}

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
        # import pdb; pdb.set_trace();
        scan = msf_scans.get((specId, i.get('peptide')))
        fname = os.path.splitext(scan.file)[0]
        if fname not in hdf_filemap:
            if skip:
                continue
            sys.stderr.write('{0} not found in filemap. Filemap is {1}.'.format(fname, hdf_filemap))
            return 1
        i['file'] = fname
        i['ms1'] = scan.ms1_scan.title
        try:
            raw_files[fname].append((i.get('rt'), i))
        except KeyError:
            raw_files[fname] = [(i.get('rt'), i)]
        msf_scan_index[(specId, i.get('peptide'))] = {'mass': scan.mass, 'charge': scan.charge, 'ms1': scan.ms1_scan.title}

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
                if v == 'Heavy' and 'heavy_rt' in images:
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images.get('heavy_rt'), i))
                elif v == 'Light' and 'light_rt' in images:
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images.get('light_rt'), i))
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
        worker = Worker(queue=in_queue, results=result_queue, raw_name=filepath, heavy_masses=heavy_masses,
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
            scan = msf_scan_index.get((v.get('id'), v.get('peptide')))
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
