from collections import OrderedDict
from operator import itemgetter
import random
from collections import Counter
import math
from scipy.misc import comb

import pandas as pd
from scipy import integrate

from pythomics.proteomics.config import NEUTRON, RESIDUE_COMPOSITION, ELEMENTS
from pythomics.proteomics.peptide import Peptide


from scipy import optimize
from scipy.signal import argrelmax, argrelmin
from scipy import stats as ss
import numpy as np
import sys

ETNS = {1: {'C': .0110, 'H': 0.00015, 'N': 0.0037, 'O': 0.00038, 'S': 0.0075},
        2: {'O': 0.0020, 'S': .0421},
        4: {'S': 0.00020}}

def calculate_theoretical_distribution(peptide):
    def dio_solve(n, l=None, index=0, out=None):
        if l is None:
            l = [1,2,4]
        if out is None:
            out = [0]*len(l)
        if index != len(l):
            for i in xrange(int(n/l[index])+1):
                out[index] = i
                for j in dio_solve(n, l=l, index=index+1, out=out):
                    yield j
        else:
            if n == sum([a*b for a,b in zip(l,out)]):
                yield out

    ETN_P = {}
    element_composition = {}
    aa_counts = Counter(peptide)
    for aa, aa_count in aa_counts.items():
        for element, element_count in RESIDUE_COMPOSITION[aa].items():
            try:
                element_composition[element] += aa_count*element_count
            except KeyError:
                element_composition[element] = aa_count*element_count
    # we lose a water for every peptide bond
    peptide_bonds = len(peptide)-1
    element_composition['H'] -= peptide_bonds*2
    element_composition['O'] -= peptide_bonds
    # and gain a hydrogen for our NH3
    element_composition['H'] += 1

    total_atoms = sum(element_composition.values())
    for etn, etn_members in ETNS.items():
        p = 0.0
        for isotope, abundance in etn_members.items():
            p += element_composition.get(isotope, 0)*abundance/total_atoms
        ETN_P[etn] = p
    tp = 0
    dist = []
    while tp < 0.999:
        p=0
        for solution in dio_solve(len(dist)):
            p2 = []
            for k, i in zip(solution, [1, 2, 4]):
                petn = ETN_P[i]
                p2.append((comb(total_atoms, k)*(petn**k))*((1-petn)**(total_atoms-k)))
            p+=np.cumprod(p2)[-1]
        tp += p
        dist.append(p)
    return pd.Series(dist)

def fit_theo_dist(params, ny, ty):
    right_limit, scaling = params
    index = xrange(len(ny)) if len(ny) > len(ty) else xrange(len(ty))
    exp_dist = pd.Series(0, index=index)
    theo_dist = pd.Series(0, index=index)
    exp_dist += pd.Series(ny)
    exp_dist[int(right_limit):] = 0
    theo_dist += ty
    return ((exp_dist-theo_dist*scaling)**2).sum()

def neg_binomial((n, p), ny):
    if p >= 1 or p <= 0:
        return np.inf
    if n > len(ny) or n <= 0:
        return np.inf
    x = range(len(ny))
    guess = pd.Series(ss.nbinom.pmf(x, n, p), index=x)
    real = pd.Series(ny)
    return ((guess-real)**2).sum()

def neg_binomial2((n1,p1, prop1, n2,p2,n2_offset, prop2), ny):
    if prop1 > 1 or prop1 <0 or prop2 < 0 or prop2 > 1 or (prop1+prop2)>1:
        return np.inf
    if n2_offset < 0:
        return np.inf
    if p1 >= 1 or p2 >= 1 or p1 <= 0 or p2 <= 0:
        return np.inf
    if n2 > len(ny) or n1 > len(ny) or n2 <= 0 or n1 <= 0:
        return np.inf

    # model the left peak
    x = range(len(ny))
    fit1 = pd.Series(ss.nbinom.pmf(x, n1, p1), index=x)*prop1

    # model the right peak
    fit2 = pd.Series(ss.nbinom.pmf(x, n2, p2, loc=n2_offset), index=x)*prop2
    combined = fit1+fit2
    return ((pd.Series(ny)-combined)**2).sum()

def neg_binomial3((n1,p1, prop1, n2,p2,n2_offset, prop2), ny):
    if (prop1+prop2)>1:
        return np.inf
    if n2_offset < n1:
        return np.inf
    if prop1 > 1 or prop1 <0 or prop2 < 0 or prop2 > 1 or (prop1+prop2)>1:
        return np.inf
    if n2_offset < 0:
        return np.inf
    if p1 >= 1 or p2 >= 1 or p1 <= 0 or p2 <= 0:
        return np.inf
    if n2 > len(ny) or n1 > len(ny) or n2 <= 0 or n1 <= 0:
        return np.inf

    # model the left peak
    x = range(len(ny))
    fit1 = pd.Series(ss.nbinom.pmf(x, n1, p1), index=x)

    fit1 = (fit1/fit1.max()).fillna(0)
#     data_max = ny[fit1.idxmax()]
#     fit1 = fit1*data_max*prop1
    fit_res = ny-fit1*prop1

    fit2 = pd.Series(ss.nbinom.pmf(x, n2, p2, loc=n2_offset), index=x)
    fit2 = (fit2/fit2.max()).fillna(0)
#     data_max = fit_res[fit2.idxmax()]
#     fit2 = fit2*data_max*prop2
    fit_res -= (fit2*prop2)
#     print fit1+fit2
#     print ((ny-(fit1+fit2))**2).sum()


#     fit_res = fit_res-fit2_norm*ny
#     print 'f2',n2,p2,prop2,n2_offset,fit2_norm*ny
#     print fit_res,
#     return (fit_res**2).sum()
    return abs(fit_res).sum()

def fit_data(data, charge=1.0, peptide=None):
    spacing = NEUTRON/float(charge)
    Y = data.values
    ny = np.array(Y, dtype=float)/np.sum(Y)
    x_axis = range(len(ny))
    initial_guess = np.average(x_axis, weights=ny)
    opt_kwargs = {
        'args': (ny,),
        'method': 'Nelder-Mead',
    }
    res = optimize.minimize(neg_binomial, (initial_guess, 0.4), **opt_kwargs)
    if res.fun > 0.1 or res.success is False:
        old_res = res
        opt_kwargs.update({'method': 'Powell'})
        # try a 2 state model
        # res = optimize.minimize(neg_binomial2, (initial_guess, 0.4, 0.5, initial_guess, 0.4, 3, 0.5), **opt_kwargs)
        res = optimize.minimize(neg_binomial3, (initial_guess, 0.4, 0.5, initial_guess, 0.4, int(len(ny)/2), 0.5), tol=1e-10, **opt_kwargs)
        if res.success is False:
            res = old_res
            n1,p1 = old_res.x
        else:
            n1,p1,n1prop, n2,p2,n2_offset, n2prop = res.x
        # model the left peak
        # fit1 = pd.Series(ss.nbinom.pmf(x_axis, n1, p1), index=data.index)*n1prop
        # fit2 = pd.Series(ss.nbinom.pmf(x_axis, n2, p2, loc=n2_offset), index=data.index)*n2prop
    else:
        n1, p1 = res.x
    # take the index of our fitted max from our data and extend it out until we're at ~0.01 of the distribution
    x_axis = range(len(data))
    data_x = list(data.index)
    fit1_nonorm = pd.Series(ss.nbinom.pmf(x_axis, n1, p1), index=data_x)
    fit1 = (fit1_nonorm/fit1_nonorm.max()).fillna(0)
    tries=0
    while fit1.max() < 1 or fit1_nonorm.sum() < 0.85:
        tries+=1
        if tries > 15:
            sys.stderr.write('Failure on {}\n'.format(data))
            return {'fit': pd.Series(), 'residual': np.inf}
        data_x.append(data_x[-1]+spacing)
        x_axis.append(len(x_axis))
        fit1_nonorm = pd.Series(ss.nbinom.pmf(x_axis, n1, p1), index=data_x)
        fit1 = (fit1_nonorm/fit1_nonorm.max()).fillna(0)
    data_max = data.loc[fit1.idxmax()]
    fitted = fit1*data_max
    return {'fit': fitted, 'residual': res.fun}


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

def buildEnvelope(peaks_found=None, isotopes=None, gradient=False, rt_window=None, start_rt=None):
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

    # THis is experimental
    ndfs = []
    from scipy.signal import argrelmax
    closest = sorted([(i, abs(float(v)-start_rt)) for i,v in enumerate(df.columns)], key=itemgetter(1))[0][0]
    for index, values in df.iterrows():
        y = values.fillna(0)
        maximas = y.iloc[argrelmax(y.fillna(0).values)[0]]
        # check for peaks outside our window
        neighboring_peaks = filter(lambda x: x>0.2, maximas.index-start_rt)
        if not neighboring_peaks and len(maximas) <= 1:
            ndfs.append(pd.DataFrame(y))
            continue
        srt = sorted([(i, abs(float(v[0])-start_rt)) for i,v in enumerate(y.iteritems()) if v[1]], key=itemgetter(1))[0][0]
        if abs(srt-closest) > 2:
            continue
        # find the left/right from our peak
        peak = y.iloc[srt]
        left=0
        # find the gradient
        if any([i for i in y.iloc[-2:srt+1][-2:] if i<peak]):
            # increasing to the right
            shift = 'right'
        else:
            shift = 'left'
        for left in xrange(srt-1, -1, -1):
            val = y.iloc[left]
            if val == 0 or val > peak:
                if val == 0 or shift != 'left':
                    break
            elif shift == 'left':
                shift = None
        right = len(y)
        highest_val = peak
        for right in xrange(srt+1, len(y)):
            val = y.iloc[right]
            if val == 0 or val > peak:
                if val > highest_val:
                    highest_val = val
                if val == 0 or shift != 'right':
                    if val == 0:
                        right += 1
                    break
            elif shift == 'right':
                shift = None
                peak = highest_val
        ndfs.append(pd.DataFrame(y.iloc[left:right]))

    df = pd.concat(ndfs, axis=1).T if ndfs else df
    # THis is experimental

    rt_data = df.fillna(0).sum(axis=0).sort_index()
    rt_data.index = rt_data.index.astype(float)
    int_val = integrate.trapz(rt_data.values, rt_data.index.values) if not rt_data.empty else 0
    return {'data': df, 'integration': int_val, 'rt': rt_data}

import copy
def looper(selected=None, df=None, theo=None, index=0, out=None):
    if out is None:
        out = [0]*len(selected)
    if index != len(selected):
        for i in selected[index]:
            out[index] = i
            for j in looper(selected=selected, df=df, theo=theo, index=index+1, out=out):
                yield j
    else:
        vals = pd.Series([df[i] for i in out])
        vals = vals/vals.max()
        residual = ((theo-vals)**2).sum()
        yield (residual, copy.deepcopy(out))

def findEnvelope(df, start_mz=None, max_mz=None, ppm=5, charge=2, debug=False,
                 heavy=False, isotope_offset=0, theo_dist=None):
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
    # # check behind us for a peak, if we find one taller than us, quit because we're in a contaminant, or we're in a non-modified
    # # version of the peptide
    # con_mz = start_mz-spacing
    # con_non_empty_ind = non_empty.index.searchsorted(con_mz)
    # con_start = non_empty.index[non_empty_ind]
    # if abs(start-start_mz)/start <= tolerance:
    #     # we found one
    #     if non_empty[start] < non_empty[con_start]:
    #         return empty_dict
    # we reset the isotope_index back to 1 here, so on the subsequent steps we aren't looking
    isotope_index += isotope_offset
    # find the largest intensity within our tolerance
    start_df = df.iloc[df.index.searchsorted(start-start*tolerance/2):df.index.searchsorted(start+start*tolerance/2),]
    start = start_df.idxmax()
    start_index = df.index.searchsorted(start)
    # find the largest in our tolerance
    # env_dict[isotope_index] = start_index
    valid_locations2 = OrderedDict()
    valid_locations2[isotope_index] = list(start_df.index)

    # micro means return the 'micro envelope' which is the envelope of each isotopic cluster, start with our beginning isotope
    # micro_dict[isotope_index] = findMicro(df, start_index, ppm=ppm)
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
            valid_locations2[isotope_index] = valid_locations
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
    valid_keys = sorted(set(valid_locations2.keys()).intersection(theo_dist.keys()))
    valid_vals = [valid_locations2[i] for i in valid_keys]
    valid_theor = pd.Series([theo_dist[i] for i in valid_keys])
    valid_theor = valid_theor/valid_theor.max()
    best_locations = sorted(looper(selected=valid_vals, df=df, theo=valid_theor), key=itemgetter(0))[0][1]
    # min_loc, max_loc = df.index.searchsorted(min(valid_locations)), df.index.searchsorted(max(valid_locations))
    # search out and create the micro envelope for this
    # micro_envelopes.append((min_loc, max_loc+1))#self.findMicro(df, micro_index, charge=charge))
    for index, isotope_index in enumerate(valid_keys):
        largest_loc = best_locations[index]
        micro_index = df.index.searchsorted(largest_loc)
        micro_dict[isotope_index] = findMicro(df, micro_index, ppm=ppm)
        env_dict[isotope_index] = micro_index
        ppm_dict[isotope_index] = 0#abs(abs(start-largest_loc)-offset)/largest_loc
    return {'envelope': env_dict, 'micro_envelopes': micro_dict, 'ppms': ppm_dict}


