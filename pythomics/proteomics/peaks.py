from collections import OrderedDict

import pandas as pd
from scipy import integrate

from pythomics.proteomics.config import NEUTRON
from pythomics.proteomics.peptide import Peptide

def findMicro(df, pos, ppm=None):
    """
        We want to find the boundaries of our isotopic clusters. Basically we search until our gradient
        changes, this assumes it's roughly gaussian and there is little interference
    """
    df_empty_index = df[df==0].index
    right = df_empty_index.searchsorted(df.index[pos])
    left = right-1
    return (df.index.searchsorted(df_empty_index[left])-1,
            df.index.searchsorted(df_empty_index[right]))

def merge_list(starting_list):
    final_list = []
    for i,v in enumerate(starting_list[:-1]):
        if set(v)&set(starting_list[i+1]):
            starting_list[i+1].extend(list(set(v) - set(starting_list[i+1])))
        else:
            final_list.append(v)
    final_list.append(starting_list[-1])
    return final_list

def buildEnvelope(peaks_found=None, isotopes=None, gradient=False, rt_window=None):
    isotope_index = {}
    df = pd.DataFrame(columns=rt_window if rt_window is not None else [])
    for silac_isotope, isotope_data in peaks_found.iteritems():
        if isotopes and silac_isotope not in isotopes:
            continue
        switch = False
        last_val = 0
        for envelope, micro_envelope in zip(
                isotope_data.get('envelope', []),
                isotope_data.get('micro_envelopes', [])):
            if envelope:
                iso_df = envelope['df']
                rt = iso_df.name

            if micro_envelope and micro_envelope.get('info'):
                start, end = micro_envelope['info'][0], micro_envelope['info'][1]
                selector = iso_df.index[range(start,end+1)]
                temp_series = pd.Series(iso_df[selector], index=selector)
                try:
                    series_index = isotope_index[silac_isotope]
                except KeyError:
                    series_index = temp_series.idxmax()
                    isotope_index[silac_isotope] = series_index
                temp_int = integrate.simps(temp_series.values)
                if gradient:
                    if temp_int > last_val and switch is False:
                        pass
                    elif temp_int > last_val and switch is True:
                        break
                    elif temp_int < last_val and switch is True:
                        pass
                    elif temp_int < last_val and switch is False:
                        switch = True
                        pass
                last_val = temp_int
                if series_index in df.index:
                    if rt in df and not pd.isnull(df.loc[series_index, rt]):
                        df.loc[series_index, rt] += temp_int
                    else:
                        df.loc[series_index, rt] = temp_int
                else:
                    df = df.append(pd.DataFrame(temp_int, columns=[rt], index=[series_index]))
    rt_data = df.fillna(0).sum(axis=0).sort_index()
    int_val = integrate.trapz(rt_data.values, rt_data.index.values) if not rt_data.empty else 0
    return {'data': df, 'integration': int_val, 'rt': rt_data}

def findEnvelope(df, start_mz=None, max_mz=None, ppm=5, charge=2, debug=False, heavy=False, isotope_offset=0, peptide=None):
    # returns the envelope of isotopic peaks as well as micro envelopes  of each individual cluster
    spacing = NEUTRON/float(charge)
    start_mz = start_mz/float(charge) if isotope_offset == 0 else (start_mz+isotope_offset*NEUTRON)/float(charge)
    if max_mz is not None:
        max_mz = max_mz/float(charge) if isotope_offset == 0 else (max_mz+isotope_offset*NEUTRON)/float(charge)
    tolerance = ppm/1000000.0

    non_empty = df[df>0].dropna()
    non_empty_ind = non_empty.index.searchsorted(start_mz)
    start = non_empty.index[non_empty_ind]
    attempts = 0
    env_dict, micro_dict, ppm_dict = OrderedDict(),OrderedDict(),OrderedDict()
    empty_dict = {'envelope': env_dict, 'micro_envelopes': micro_dict, 'ppms': ppm_dict}
    isotope_index = 0
    while abs(start-start_mz)/start >= tolerance:
        # check if the second peak is present
        if heavy is False or attempts >= 1:
            return empty_dict
        start_mz += spacing
        non_empty_ind = non_empty.index.searchsorted(start_mz)
        start = non_empty.index[non_empty_ind]
        if max_mz is not None and start >= max_mz:
            return empty_dict
        attempts += 1
        isotope_index += 1
    # check behind us for a peak, if we find one taller than us, quit because we're in a contaminant, or we're in a non-modified
    # version of the peptide
    con_mz = start_mz-spacing
    con_non_empty_ind = non_empty.index.searchsorted(con_mz)
    con_start = non_empty.index[non_empty_ind]
    if abs(start-start_mz)/start <= tolerance:
        # we found one
        if non_empty[start] < non_empty[con_start]:
            return empty_dict
    # we reset the isotope_index back to 1 here, so on the subsequent steps we aren't looking
    isotope_index += isotope_offset
    # find the largest intensity within our tolerance
    start_df = df.iloc[df.index.searchsorted(start-start*tolerance/2):df.index.searchsorted(start+start*tolerance/2),]
    start = start_df.idxmax()
    start_index = df.index.searchsorted(start)
    # find the largest in our tolerance
    env_dict[isotope_index] = start_index

    # micro means return the 'micro envelope' which is the envelope of each isotopic cluster, start with our beginning isotope
    micro_dict[isotope_index] = findMicro(df, start_index, ppm=ppm)
    isotope_index += 1
    pos = non_empty_ind+1
    offset = isotope_index*spacing
    df_len = non_empty.shape[0]
    last_displacement = None
    valid_locations = set([])
    while pos < df_len:
        # search for the ppm error until it rises again, we select the minima and if this minima is
        # outside our ppm error, we stop the expansion of our isotopic cluster
        current_loc = non_empty.index[pos]
        if max_mz is not None and current_loc >= max_mz:
            if not valid_locations:
                break
            displacement = last_displacement+tolerance if last_displacement is not None else tolerance*2
        else:
            displacement = abs(abs(start-current_loc)-offset)/current_loc
        if debug:
            print pos, start, current_loc, displacement, last_displacement, displacement > last_displacement, last_displacement < tolerance, isotope_index, offset
        if displacement < tolerance:
            valid_locations.add(current_loc)
        if valid_locations and displacement > last_displacement:
            # pick the largest peak within our error tolerance
            largest_loc = df[valid_locations].idxmax()
            # min_loc, max_loc = df.index.searchsorted(min(valid_locations)), df.index.searchsorted(max(valid_locations))
            # search out and create the micro envelope for this
            micro_index = df.index.searchsorted(largest_loc)
            # micro_envelopes.append((min_loc, max_loc+1))#self.findMicro(df, micro_index, charge=charge))
            micro_dict[isotope_index] = findMicro(df, micro_index, ppm=ppm)
            env_dict[isotope_index] = micro_index
            ppm_dict[isotope_index] = abs(abs(start-largest_loc)-offset)/largest_loc
            isotope_index += 1
            offset = spacing*isotope_index
            displacement = abs(abs(start-current_loc)-offset)/current_loc
            valid_locations = set([])
        elif last_displacement is not None and displacement > last_displacement and not valid_locations:
            break
        # elif not valid_locations:
        #     # we're done
        #     break
        last_displacement = displacement
        pos += 1

    #combine any overlapping micro envelopes
    #final_micros = self.merge_list(micro_dict)
    return {'envelope': env_dict, 'micro_envelopes': micro_dict, 'ppms': ppm_dict}