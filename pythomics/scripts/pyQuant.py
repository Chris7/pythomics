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
from matplotlib import pyplot as plt
from multiprocessing import Process, Queue
from profilestats import profile
from Queue import Empty
from scipy import integrate, optimize
from matplotlib.backends.backend_pdf import PdfPages
from cStringIO import StringIO
import urllib
import base64

from pythomics.templates import CustomParser
from pythomics.proteomics.parsers import GuessIterator

RESULT_ORDER = (('peptide', 'Peptide'), ('heavy', 'Heavy'), ('light', 'Light'), ('h_l', 'Heavy/Light'), ('modifications', 'Modifications'),
                ('charge', 'Charge'), ('ms1', 'MS1 Spectrum ID'), ('scan', 'MS2 Spectrum ID'), ('light_precursor', 'Light Precursor M/Z'),
                ('heavy_precursor', 'Heavy Precursor M/Z'),
                ('light_clusters', 'Light Peaks Identified'), ('heavy_clusters', 'Heavy Peaks Identified'), ('light_rt_dev', 'Light Retention Time Deviation'),
                ('heavy_rt_dev', 'Heavy Retention Time Deviation'), ('light_residuals', 'Light Intensity Residuals of Fit'),
                ('heavy_residuals', 'Heavy Intensity Residuals of Fit'), ('cluster_ratio_deviation', 'Deviation of Heavy Clusters From Light'))


parser = CustomParser(description=description)
parser.add_processed_ms()
parser.add_argument('--mzml', help="The mzML files for the raw data or directory containing them. Defaults to the directory of the processed file.", type=str)
parser.add_argument('--filemap', help="By default, the mzML file is assumed to be a derivation of the scan filename. " \
                                     "If this is not true, this file provides a correct mapping.")
parser.add_argument('--rt-width', help="The width of the retention time window to search for. Default: 1.0 minute", type=float, default=1.0)
parser.add_argument('--precision', help="The precision for storing m/z values. Defaults to 6 decimal places.", type=int, default=6)
parser.add_argument('--no-sparse', help="Don't store data as sparse entities. Major memory savings at the expense of speed.", action='store_false')
parser.add_argument('--no-temp', help="This will disable saving of the intermediate files for SILAC analysis.", action='store_false')
parser.add_argument('--debug', help="This will output debug information and graphs.", action='store_true')
parser.add_argument('--skip', help="If true, skip scans with missing files in the mapping.", action='store_true')
parser.add_argument('--peptide', help="The peptide to limit quantification to.", type=str)
parser.add_argument('--html', help="Output a HTML table summary.", action='store_true')
parser.add_argument('--resume', help="Will resume from the last run. Only works if not directing output to stdout.", action='store_true')
parser.add_argument('-o', '--out', nargs='?', help='Stuff', type=str)

class Worker(Process):
    def __init__(self, queue=None, results=None, precision=6, raw_name=None, heavy_masses=None, temp=False, sparse=False, debug=False, html=False):
        super(Worker, self).__init__()
        self.precision = precision
        self.sparse = sparse
        self.queue=queue
        self.results = results
        if heavy_masses is None:
            self.heavy_masses = {'R': 10.008269, 'K': 8.0142}
        else:
            self.heavy_masses = heavy_masses
        self.hydrogen = 1.00794
        self.raw_name = raw_name
        self.filename = os.path.split(self.raw_name)[1]
        self.temp = temp
        self.rt_width = 1
        self.rt_tol = 0.2 # for fitting
        self.debug = debug
        self.html = html
        self.last_clean = 0 # for cleaning up scans before our current RT window
        sys.stderr.write('Thread spawned: Sparse: {0}\n'.format(self.sparse))

    def findMicro(self, df, pos):
        """
            We want to find the zeros around our given index value to see the boundaries of our isotopic clusters.
        """
        df_empty_index = df[df==0].index
        right = df_empty_index.searchsorted(df.index[pos])
        left = right-1
        return (df.index.searchsorted(df_empty_index[left])-1,
                df.index.searchsorted(df_empty_index[right]))

    def merge_list(self, starting_list):
        final_list = []
        for i,v in enumerate(starting_list[:-1]):
            if set(v)&set(starting_list[i+1]):
                starting_list[i+1].extend(list(set(v) - set(starting_list[i+1])))
            else:
                final_list.append(v)
        final_list.append(starting_list[-1])
        return final_list


    def findEnvelope(self, df, start_mz=None, ppm=5, charge=2, debug=False):
        neutron = 1.008701
        # returns the envelope of isotopic peaks as well as micro envelopes  of each individual cluster
        spacing = neutron/float(charge)
        start_mz = start_mz/float(charge)
        tolerance = ppm/1000000.0
        ppm_dist = []

        non_empty = df[df>0].dropna()
        non_empty_ind = non_empty.index.searchsorted(start_mz)
        start = non_empty.index[non_empty_ind]
        start_index = df.index.searchsorted(start)
        right_envelope = [start_index]
        # micro means return the 'micro envelope' which is the envelope of each isotopic cluster, start with our beginning isotope
        micro_envelopes = [self.findMicro(df, start_index)]
        pos = non_empty_ind+1
        isotope_index = 1
        offset = spacing
        df_len = non_empty.shape[0]
        last_displacement = None
        while pos < df_len:
            # search for the ppm error until it rises again, we select the minima and if this minima is
            # outside our ppm error, we stop the expansion of our isotopic cluster
            current_loc = non_empty.index[pos]
            displacement = abs(abs(start-current_loc)-offset)/current_loc
            if debug:
                print start, current_loc, displacement, last_displacement, displacement > last_displacement, last_displacement < tolerance, isotope_index, offset
            if last_displacement is not None and displacement > last_displacement:
                if last_displacement < tolerance:
                    # search out and create the micro envelope for this
                    micro_index = df.index.searchsorted(last_loc)
                    micro_envelopes.append(self.findMicro(df, micro_index))
                    right_envelope.append(micro_index)
                    ppm_dist.append(last_displacement)
                    isotope_index += 1
                    offset = spacing*isotope_index
                    displacement = abs(abs(start-current_loc)-offset)/current_loc
                else:
                    break
            last_displacement = displacement
            last_loc = current_loc
            pos += 1

        #combine any overlapping micro envelopes
        final_micros = self.merge_list(micro_envelopes)
        return {'envelope': right_envelope, 'micro_envelopes': final_micros, 'ppms': ppm_dist}

    def getScan(self, ms1, dense=True):
        import numpy as np
        if ms1 not in self.data:
            scan = self.raw.getScan(ms1)
            scan_vals = np.array(scan.scans)
            res = pd.Series(scan_vals[:, 1].astype(np.uint16), index=np.round(scan_vals[:, 0], self.precision), name=scan.rt, dtype='uint16')
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
        residuals = sum((fit-Y)**2)
        return residuals

    def run(self):
        # first, load our series in
        import time
        args = parser.parse_args()
        source = self.raw_name
        # if os.path.exists(pickle):
        #     res = np.load(pickle)
        #     scan_rt_map, rt_index_map, df_map = res['scan_rt_map'].tolist(), res['rt_index_map'].tolist(), res['df_map'].tolist()
        #     self.data = [pd.Series(i).to_sparse(fill_value=0) if self.sparse else pd.Series(i) for i in res['data']]
        if True:
            self.raw = GuessIterator(source, full=False, store=False, ms_filter=1)
            self.data = {}
            self.data2 = pd.DataFrame(dtype='uint16')
            scan_rt_map = {}
            rt_index_map = {}
            df_map = {}
            precision = self.precision
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
                # data_index = len(data)
                rt_index_map[rt] = scan_title
                # df_map[scan_title] = data_index
                # data.append(pd.Series(dict([(np.round(mz, precision), int(intensity)) for mz, intensity in i.scans]), dtype='uint16', name=rt))
                # if self.sparse:
                #     data[-1] = data[-1].to_sparse(fill_value=0)
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
                precursor = scan.get('mass')+self.hydrogen*(charge-1)
                df = self.getScan(ms1)
                shift = scan_info['mass_shift']
                peptide = scan_info['peptide']
                if self.debug:
                    sys.stderr.write('on ms {0} {1} {2} {3}\n'.format(ms1, rt, precursor, scan_info))
                if shift == 0:
                    for label, label_mass in self.heavy_masses.iteritems():
                        shift += peptide.lower().count(label.lower())*label_mass
                    light_precursor = precursor
                    heavy_precursor = precursor+shift
                else:
                    light_precursor = precursor-shift
                    heavy_precursor = precursor

                quant = light_precursor != heavy_precursor

                light = self.findEnvelope(df, start_mz=light_precursor, charge=charge, ppm=5)
                if quant:
                    heavy = self.findEnvelope(df, start_mz=heavy_precursor, charge=charge, ppm=5)
                    # check for overlap of heavy/light
                    common_peaks = set(light['envelope']).intersection(set(heavy['envelope']))
                    for common_peak in common_peaks:
                        # light yields to heavy
                        common_index = light['envelope'].index(common_peak)
                        del light['envelope'][common_index]
                        del light['micro_envelopes'][common_index]
                        del light['ppms'][common_index-1]


                start_mzs = [df.index[i[0]] for i in light['micro_envelopes']]
                end_mzs = [df.index[i[-1]] for i in light['micro_envelopes']]
                light_clusters = len(start_mzs)
                light_query = ' | '.join(['(index >= {0} & index <= {1})'.format(start, end) for start, end in zip(start_mzs, end_mzs)])

                if quant:
                    heavy_start_mzs = [df.index[i[0]] for i in heavy['micro_envelopes']]
                    heavy_end_mzs = [df.index[i[-1]] for i in heavy['micro_envelopes']]
                    heavy_clusters = len(heavy_start_mzs)
                    heavy_query = ' | '.join(['(index >= {0} & index <= {1})'.format(start, end) for start, end in zip(heavy_start_mzs, heavy_end_mzs)])

                # PLOTTING FOR DEBUG

                if self.debug or self.html:
                    if self.debug:
                        pdf = PdfPages('{2}_{0}_{1}_fix.pdf'.format(peptide, ms1, self.raw_name))
                        plt.figure()

                    light_x = []
                    light_y = []
                    smz = 0
                    for i in light['micro_envelopes']:
                        if smz == 0:
                            smz = i[0]
                        val_ind = df.index[i[0]]
                        light_x.append(val_ind)
                        light_y.append(0)
                        for j in xrange(i[0], i[1]+1):
                            val = df.iloc[j]
                            val_ind = df.index[j]
                            light_x.append(val_ind)
                            light_y.append(val)

                    if quant:
                        heavy_x = []
                        heavy_y = []
                        for i in heavy['micro_envelopes']:
                            for j in xrange(i[0], i[1]+1):
                                val = df.iloc[j]
                                val_ind = df.index[j]
                                heavy_x.append(val_ind)
                                heavy_y.append(val)
                    ax = df.iloc[smz:j].to_dense().plot()
                    fig = ax.get_figure()
                    if self.debug:
                        pdf.savefig(figure=fig)
                    fig.clf()


                    plt.figure()
                    ax = plt.plot(light_x, light_y, color='b')
                    if quant:
                        ax = plt.plot(heavy_x, heavy_y, color='r')
                    fig = ax[0].get_figure()
                    if self.debug:
                        pdf.savefig(figure=fig)
                    if self.html:
                        fname = '{2}_{0}_{1}_{3}_clusters.png'.format(peptide, ms1, self.filename, scanId)
                        fig.savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                        html_images['clusters'] = os.path.join(self.html['rel'], fname)
                    fig.clf()


                scan_query = ' | '.join([light_query, heavy_query] if quant else [light_query])
                try:
                    # use the RT to condense our search
                    t = time.time()
                    # print time.time()-t, 'time to get ', len(scan_list), 'scans'
                    scan_block = pd.concat((self.getScan(i, dense=False).dropna().to_frame().query(scan_query) for i in rt_index_map[rt-self.rt_width:rt+self.rt_width]), axis=1)
                    t = time.time()
                    # print 'concat time is', time.time()-t
                    # print 'scan block is,', scan_block.shape, scan_block.memory_usage(index=True).sum()
                except:
                    sys.stderr.write('ERROR ON {0}\n{1}\n{2}\n{3}\n'.format(ms1, scan_info, scan_query, traceback.format_exc()))
                    return
                scan_block.dropna(how='all', inplace=True)

                light_data = scan_block.query(light_query).fillna(0).sum()

                if quant:
                    heavy_data = scan_block.query(heavy_query).fillna(0).sum()

                if self.debug:
                    from matplotlib import cm
                    from mpl_toolkits.mplot3d.axes3d import Axes3D
                    fig = plt.figure(figsize=(14,6))
                    res = scan_block[(scan_block.values > 100) == True].dropna(how='all').fillna(0)

                    X,Y = np.meshgrid(res.index.values, res.columns.values)
                    Z = res.values.T

                    # surface_plot with color grading and color bar
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    p = ax.plot_surface(X, Y, Z, rstride=50, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    cb = fig.colorbar(p, shrink=0.5)
                    pdf.savefig(figure=ax.get_figure())

                del scan_block

                amp = light_data.max()
                #params = np.array([rt, 0.1, light_data.quantile([0.1]).iloc[0]])
                params = np.array([rt, 0.1])

                opt = optimize.minimize(self.gfit, params, args=(light_data.index.values, light_data.values, amp, rt), method='Nelder-Mead',
                                        options={'xtol': 1e-8})

                X = light_data.index.values


                mu = opt.x[0]
                light_rt_dev = abs(mu-rt)
                opt_std = abs(opt.x[1])
                #opt_offset = opt.x[2]
                opt_offset = 0
                light_residuals = opt.fun
                fit = lambda t : (amp-opt_offset)*np.exp(-(t-mu)**2/(2*opt_std**2))
                light_int = integrate.quad(fit, mu-6*opt_std, mu+6*opt_std)[0]
                if self.debug or self.html:
                    plt.figure()
                    ax = light_data.plot(color='b', title=str(light_int))
                    pd.Series(dict(zip(X, fit(X)+opt_offset))).plot(color='r')
                    plt.plot([rt, rt], [0, (amp-opt_offset)], color='k')
                    plt.plot([mu, mu], [0, (amp-opt_offset)], color='r')
                    if self.debug:
                        pdf.savefig(figure=ax.get_figure())
                    if self.html:
                        fname = '{2}_{0}_{1}_{3}_light_rt.png'.format(peptide, ms1, self.filename, scanId)
                        ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                        html_images['light_rt'] = os.path.join(self.html['rel'], fname)
                    ax.get_figure().clf()

                if quant:
                    amp = heavy_data.max()
                    #params = np.array([rt, 0.1, heavy_data.quantile([0.1]).iloc[0]])
                    params = np.array([rt, 0.1])

                    opt = optimize.minimize(self.gfit, params, args=(heavy_data.index.values, heavy_data.values, amp, rt), method='Nelder-Mead',
                                            options={'xtol': 1e-8})

                    X = heavy_data.index.values

                    mu = opt.x[0]
                    heavy_rt_dev = abs(mu-rt)
                    opt_std = abs(opt.x[1])
                    #opt_offset = opt.x[2]
                    opt_offset = 0
                    heavy_residuals = opt.fun
                    fit = lambda t : (amp-opt_offset)*np.exp(-(t-mu)**2/(2*opt_std**2))
                    heavy_int = integrate.quad(fit, mu-6*opt_std, mu+6*opt_std)[0]
                    light_rats = [df.iloc[v]/df.iloc[light['envelope'][i]] for i,v in enumerate(light['envelope'][1:])]
                    heavy_rats = [df.iloc[v]/df.iloc[heavy['envelope'][i]] for i,v in enumerate(heavy['envelope'][1:])]
                    cluster_rat = sum([(i-v) for i, v in zip(light_rats, heavy_rats[:-1])])
                    if self.debug or self.html:
                        plt.figure()
                        ax = heavy_data.plot(color='b', title=str(heavy_int))
                        pd.Series(dict(zip(X, fit(X)+opt_offset))).plot(color='r')
                        plt.plot([rt, rt], [0, (amp-opt_offset)], color='k')
                        plt.plot([mu, mu], [0, (amp-opt_offset)], color='r')
                        if self.debug:
                            pdf.savefig(figure=ax.get_figure())
                        if self.html:
                            fname = '{2}_{0}_{1}_{3}_heavy_rt.png'.format(peptide, ms1, self.filename, scanId)
                            ax.get_figure().savefig(os.path.join(self.html['full'], fname), format='png', dpi=50)
                            html_images['heavy_rt'] = os.path.join(self.html['rel'], fname)
                        ax.get_figure().clf()
                        plt.close('all')


                if self.debug:
                    pdf.close()

                self.results.put({'peptide': peptide, 'heavy': heavy_int if quant else 'NA', 'light': light_int,
                                  'scan': scanId, 'light_clusters': light_clusters,
                                  'h_l': heavy_int/light_int if quant and light_int else 'NA',
                                  'heavy_clusters': heavy_clusters if quant else 'NA', 'ms1': ms1,
                                  'charge': charge, 'modifications': mods, 'light_precursor': light_precursor,
                                  'heavy_precursor': heavy_precursor if quant else 'NA',
                                  'light_residuals': light_residuals, 'heavy_residuals': heavy_residuals if quant else 'NA',
                                  'cluster_ratio_deviation': cluster_rat if quant else 'NA', 'light_rt_dev': light_rt_dev,
                                  'heavy_rt_dev': heavy_rt_dev if quant else 'NA', 'html_info': html_images})
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
    source = args.processed
    threads = args.p
    skip = args.skip
    out = args.out
    hdf = args.mzml
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
        hdf_filemap[hdf] = os.path.abspath(hdf)

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
        mass_scans[specId] = {'id': specId, 'peptide': i.peptide, 'mass_shift': shift, 'rt': i.rt, 'charge': i.charge, 'modifications': i.getModifications()}

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
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images['heavy_rt'], i))
                elif v == 'Light' and 'light_rt' in images:
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images['light_rt'], i))
                elif v == 'Peptide' and 'clusters' in images:
                    out.append("""<td data-toggle="popover" data-trigger="hover click" data-content='<img src="{0}">'>{1}</td>""".format(images['clusters'], i))
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
                        temp=save_files, sparse=sparse, debug=args.debug, html=html)
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
                        res = '{0}\t{1}\n'.format(rawstr, '\t'.join([str(result[i[0]]) for i in RESULT_ORDER]))
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
                res = '{0}\t{1}\n'.format(rawstr, '\t'.join([str(result[i[0]]) for i in RESULT_ORDER]))
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
