import math
from collections import OrderedDict
from operator import itemgetter
from collections import Counter
from scipy.misc import comb
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import savgol_filter

import pandas as pd
from scipy import integrate

from pythomics.proteomics.config import NEUTRON, RESIDUE_COMPOSITION, ELEMENTS


from scipy import optimize
from scipy.signal import argrelmax, argrelmin
from scipy import stats as ss
import numpy as np
import sys

ETNS = {1: {'C': .0110, 'H': 0.00015, 'N': 0.0037, 'O': 0.00038, 'S': 0.0075},
        2: {'O': 0.0020, 'S': .0421},
        4: {'S': 0.00020}}

_epsilon = np.sqrt(np.finfo(float).eps)

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

def findMicro(df, pos, ppm=None, start_mz=None, calc_start_mz=None, isotope=0, spacing=0, quant_method='integrate'):
    """
        We want to find the boundaries of our isotopic clusters. Basically we search until our gradient
        changes, this assumes it's roughly gaussian and there is little interference
    """
    # find the edges within our tolerance
    tolerance = ppm
    offset = spacing*isotope
    df_empty_index = df[df==0].index
    right = df_empty_index.searchsorted(df.index[pos])
    left = right-1
    left, right = (df.index.searchsorted(df_empty_index[left], side='left'),
            df.index.searchsorted(df_empty_index[right]))
    right += 1
    y = df.iloc[left:right]
    peaks, peak_centers = findAllPeaks(y, min_dist=(y.index[1]-y.index[0])*2)
    if start_mz is None:
        start_mz = df.index[pos]

    # new logic is nm
    sorted_peaks = sorted([(peaks.x[i*3:(i+1)*3], get_ppm(start_mz+offset, v)) for i,v in enumerate(peaks.x[1::3])], key=lambda x: x[1])
    fit = True

    if not filter(lambda x: x[1]<tolerance, sorted_peaks):
        if calc_start_mz is not None:
            sorted_peaks2 = sorted([(peaks.x[i*3:(i+1)*3], get_ppm(calc_start_mz+offset, v)) for i,v in enumerate(peaks.x[1::3])], key=lambda x: x[1])
            if filter(lambda x: x[1]<tolerance, sorted_peaks2):
                sorted_peaks = sorted_peaks2
            else:
                fit = False
        else:
            fit = False

    peak = sorted_peaks[0][0]
    # interpolate our mean/std to a linear range
    from scipy.interpolate import interp1d
    mapper = interp1d(y.index, range(len(y)))
    try:
        mu = mapper(peak[1])
    except:
        print 'mu', sorted_peaks, peak, y.index
        return {'int': 0, 'error': np.inf}
    try:
        std = mapper(y.index[0]+np.abs(peak[2]))-mapper(y.index[0])
    except:
        print 'std', sorted_peaks, peak, y.index
        return {'int': 0, 'error': np.inf}
    peak_gauss = (peak[0]*y.max(), mu, std)
    peak[0] *= y.max()

    lr = np.linspace(peak_gauss[1]-peak_gauss[2]*4, peak_gauss[1]+peak_gauss[2]*4, 1000)
    left_peak, right_peak = peak[1]-peak[2]*2, peak[1]+peak[2]*2
    int_val = integrate.simps(gauss(lr, *peak_gauss), x=lr)# if quant_method == 'integrate' else y[(y.index > left_peak) & (y.index < right_peak)].sum()

    return {'int': int_val if fit else 0, 'bounds': (left, right), 'params': peak, 'error': sorted_peaks[0][1]}

def gauss(x, amp, mu, std):
    return amp*np.exp(-(x - mu)**2/(2*std**2))

def gauss_ndim( xdata, *args):
    amps, mus, sigmas = args[::3], args[1::3], args[2::3]
    data = np.zeros(len(xdata))
    for amp, mu, sigma in zip(amps, mus, sigmas):
        data += gauss(xdata, amp, mu, sigma)
    return data

def gauss_func( guess, *args):
    xdata, ydata = args
    data = gauss_ndim(xdata, *guess)
    # absolute deviation as our distance metric. Empirically found to give better results than
    # residual sum of squares for this data.
    return sum(np.abs(ydata-data)**2)

def gauss_jac( guess, *args):
    # take partial differentials of each parameter as the jacobian
    xdata, ydata = args
    amps, mus, sigmas = guess[::3], guess[1::3], guess[2::3]
    jac_amps = [np.exp(-(mu)/(2*sigma**2)) for mu, sigma in zip(mus, sigmas)]
    jac_mus = [amp*np.exp(-(mu)/(2*sigma**2))*((mu)/(sigma**2)) for amp, mu, sigma in zip(amps, mus, sigmas)]
    jac_sigmas = [amp*np.exp(-(mu)/(2*sigma**2))*((mu)**2/(sigma**3)) for amp, mu, sigma in zip(amps, mus, sigmas)]
    return [j for i in zip(jac_amps, jac_mus, jac_sigmas) for j in i]

def findAllPeaks2(values, min_dist=0, filter=False):
    xdata = values.index.values.astype(float)
    ydata = values.fillna(0).values.astype(float)
    # ydata_peaks[ydata<ydata.max()*.10] = 0
    mval = ydata.max()
    ydata /= ydata.max()
    ydata_peaks = np.copy(ydata)
    if filter:
        if len(ydata) >= 5:
            ydata_peaks = gaussian_filter1d(ydata_peaks, 2, mode='constant')
            ydata_peaks[ydata_peaks<0] = 0
    peaks_found = {}
    peak_count = 0
    last_fit = None
    fit_accuracy = []
    centers = []
    run_length = 0
    peaks_found = len(argrelmax(ydata)[0])+1
    while peak_count < peaks_found:
        peak_count += 1
        # remove previous peaks, add the top residual peak as a new peak center
        ydata_new = np.copy(ydata_peaks)
        if last_fit:
            ydata_new -= gauss_ndim(xdata, *res.x)
            guess = res.x.tolist()
        else:
            guess = []
            bnds = []
        peak_index = ydata_new.argmax()
        peak_max = ydata_new[peak_index]
        # ydata_new[ydata_new<peak_max*0.37] = 0
        try:
            lb = xdata[ydata_new[:peak_index]==0][-1]
        except:
            lb = xdata[0]
        try:
            rb = xdata[peak_index:][ydata_new[peak_index:]==0][0]
        except:
            rb = xdata[-1]
        # print lb, peak_index, rb
        bnds.extend([(0, 1), (lb, rb), (0.0001, (rb-lb))])
        centers.append(xdata[peak_index])
        guess.extend([peak_max, xdata[peak_index], (rb-lb)/2])
        args = (xdata, ydata)
        opts = {'maxiter': 1000}
        res = optimize.minimize(gauss_func, guess, args, method='SLSQP', options=opts, bounds=bnds)#, jac=gauss_jac)
        n = len(xdata)
        k = len(res.x)
        bic = n*np.log(res.fun/n)+k+np.log(n)
        res.bic = bic
        fit_accuracy.append((1, bic, res, copy.deepcopy(centers)))
        last_fit = True
    # we want to maximize our BIC given our definition
    best_fits = sorted(fit_accuracy, key=itemgetter(1), reverse=True)
    return best_fits[0][2:]

def bigauss(x, amp, mu, stdl, stdr):
    sigma1 = stdl/1.177
    m1 = np.sqrt(2*np.pi)*sigma1*amp
    sigma2 = stdr/1.177
    m2 = np.sqrt(2*np.pi)*sigma2*amp
    #left side
    if isinstance(x, float):
        x = np.array([x])
    lx = x[x<=mu]
    left = m1/(np.sqrt(2*np.pi)*sigma1)*np.exp(-(lx-mu)**2/(2*sigma1**2))
    rx = x[x>mu]
    right = m2/(np.sqrt(2*np.pi)*sigma2)*np.exp(-(rx-mu)**2/(2*sigma2**2))
    return np.concatenate([left, right], axis=1)

def bigauss_ndim( xdata, *args):
    amps, mus, sigmasl, sigmasr = args[::4], args[1::4], args[2::4], args[3::4]
    data = np.zeros(len(xdata))
    for amp, mu, sigma1, sigma2 in zip(amps, mus, sigmasl, sigmasr):
        data += bigauss(xdata, amp, mu, sigma1, sigma2)
    return data

def bigauss_func( guess, *args):
    xdata, ydata = args
    if any([pd.isnull(i) for i in guess]):
        return np.inf
    data = bigauss_ndim(xdata, *guess)
    # absolute deviation as our distance metric. Empirically found to give better results than
    # residual sum of squares for this data.
    return sum(np.abs(ydata-data)**2)

def fixedMeanFit(values, peak_index=None, debug=False):
    xdata = values.index.values.astype(float)
    ydata = values.fillna(0).values.astype(float)

    ydata /= ydata.max()
    rel_peak = ydata[peak_index]
    # ydata = gaussian_filter1d(ydata, 1, mode='constant')
    # fill the zeros with the average of their left/right
    # gaussian boost
    gauss_boost = gaussian_filter1d(ydata, 0.2, mode='constant')
    ydata_original = np.copy(ydata)
    for i,v in enumerate(ydata):
        if v <= gauss_boost[i]:
            boosted = gauss_boost[i]
            interp = (ydata_original[i-1]+ydata_original[i+1])/2 if i >= 1 and i < len(ydata)-1 else 0
            ydata[i] = interp if interp > boosted else boosted
    # print ydata_original.tolist()
    peak_left, peak_right = findPeak(gaussian_filter1d(savgol_filter(ydata, 5, 3), 2, mode='constant'), peak_index)
    peaks = xdata[peak_left:peak_right]
    peak_min, peak_max = xdata[0], xdata[-1]
    bnds = [(rel_peak*0.75, 1), (peak_min, peak_max), (0, xdata[peak_index]-peak_min), (0, peak_max-xdata[peak_index])]
    # reset the fitting data to our bounds
    if debug:
        print xdata.tolist(), ydata.tolist(), peak_left, peak_index, peak_right
    if peak_index == peak_right:
        peak_index -= peak_left-1
    else:
        peak_index -= peak_left
    xdata = xdata[peak_left:peak_right]
    ydata = ydata[peak_left:peak_right]
    # top_indices = ydata>np.percentile(ydata, 80)
    # print top_indices
    # ydata = ydata[top_indices]
    # xdata = xdata[top_indices]
    # if values.name > 729.36 and values.name < 731:
    if debug:
        print xdata.tolist(), ydata.tolist(), peak_left, peak_index, peak_right
    if ydata.sum() == 0:
        return None
    if debug:
        print peak_left, peak_right, xdata, ydata, values
    average = np.average(xdata, weights=ydata)
    variance = np.sqrt(np.average((xdata-average)**2, weights=ydata))
    if variance == 0:
        # we have a singular peak if variance == 0, so set the variance to half of the x/y spacing
        if peak_index >= 1:
            variance = np.abs(xdata[peak_index]-xdata[peak_index-1])
        elif peak_index < len(xdata):
            variance = np.abs(xdata[peak_index+1]-xdata[peak_index])
        else:
            # we have only 1 data point, most RT's fall into this width
            variance = 0.05
    else:
        variance = 0.05
    if variance > xdata[peak_index]-peak_min or variance > peak_max-xdata[peak_index]:
        variance = xdata[peak_index]-peak_min
    guess = [rel_peak, xdata[peak_index], variance, variance]
    # if values.name > 729.36 and values.name < 731:
    #     print guess, bnds
    args = (xdata, ydata)
    opts = {'maxiter': 1000}
    routines = ['SLSQP', 'TNC', 'L-BFGS-B']
    routine = routines.pop()
    results = [optimize.minimize(bigauss_func, guess, args, method=routine, bounds=bnds, options=opts, tol=1e-10)]#, jac=gauss_jac)]
    # if values.name > 729.36 and values.name < 731:
    if debug:
        print results[-1]
    while not results[-1].success and routines:
        # if values.name > 729.36 and values.name < 731:
        if debug:
            print results[-1]
        routine = routines.pop()
        results.append(optimize.minimize(bigauss_func, guess, args, method=routine, bounds=bnds, options=opts))#, jac=gauss_jac)
    n = len(xdata)
    if not results[-1].success:
        res = sorted(results, key=lambda x: x.fun)[0]
    else:
        res = results[-1]
    k = len(res.x)
    bic = n*np.log(res.fun/n)+k+np.log(n)
    res.bic = bic
    return res

def findAllPeaks(values, min_dist=0, filter=False, bigauss_fit=False):
    xdata = values.index.values.astype(float)
    ydata = values.fillna(0).values.astype(float)

    ydata /= ydata.max()
    ydata_peaks = np.copy(ydata)
    if filter:
        if len(ydata) >= 5:
            # ydata_peaks = values.replace([0], np.nan).interpolate(method='index').values
            ydata_peaks = gaussian_filter1d(ydata_peaks, 3, mode='constant')
            ydata_peaks[ydata_peaks<0] = 0
    peaks_found = {}
    for peak_width in xrange(1,4):
        row_peaks = argrelmax(ydata_peaks, order=peak_width)[0]
        if not row_peaks.any():
            row_peaks = [np.argmax(ydata)]
        minima = [i for i,v in enumerate(ydata_peaks) if v == 0]
        minima.extend([i for i in argrelmin(ydata_peaks, order=peak_width)[0] if i not in minima])
        minima.sort()
        peaks_found[peak_width] = {'peaks': row_peaks, 'minima': minima}
    # collapse identical orders
    final_peaks = {}
    for peak_width in xrange(1, 3):
        if peak_width == len(peaks_found):
            final_peaks[peak_width] = peaks_found[peak_width]
            continue
        smaller_peaks, smaller_minima = peaks_found[peak_width]['peaks'],peaks_found[peak_width]['minima']
        larger_peaks, larger_minima = peaks_found[peak_width+1]['peaks'],peaks_found[peak_width+1]['minima']
        if np.array_equal(smaller_peaks, larger_peaks) and np.array_equal(smaller_minima, larger_minima):
            final_peaks[peak_width+1] = peaks_found[peak_width+1]
            if peak_width in final_peaks:
                del final_peaks[peak_width]
        else:
            final_peaks[peak_width] = peaks_found[peak_width]
    fit_accuracy = []
    for peak_width, peak_info in final_peaks.items():
        row_peaks = peak_info['peaks']
        minima = peak_info['minima']
        guess = []
        bnds = []
        last_peak = None
        for peak_num, peak_index in enumerate(row_peaks):
            next_peak = len(xdata)-1 if peak_index == row_peaks[-1] else row_peaks[peak_num+1]
            peak_min, peak_max = xdata[peak_index]-0.2, xdata[peak_index]+0.2

            peak_min = xdata[0] if peak_min < xdata[0] else peak_min
            peak_max = xdata[-1] if peak_max > xdata[-1] else peak_max
            rel_peak = ydata[peak_index]/sum(ydata[row_peaks])
            std_bounds = (0.0001, peak_max-peak_min)
            bnds.extend([(rel_peak, 1), (peak_min, peak_max), std_bounds])
            if bigauss_fit:
                bnds.extend([std_bounds])
            # find the points around it to estimate the std of the peak
            left = 0
            for i,v in enumerate(minima):
                if v >= peak_index:
                    if i != 0:
                        left = minima[i-1]
                    break
                left = v
            if last_peak is not None and left < last_peak:
                left = last_peak
            last_peak = peak_index
            right = len(xdata)
            for right in minima:
                if right > peak_index or right >= next_peak:
                    if right < len(xdata) and right != next_peak:
                        right += 1
                    break
            if right > next_peak:
                right = next_peak
            if right < peak_index:
                right = next_peak
            peak_values = ydata[left:right]
            peak_indices = xdata[left:right]
            if peak_values.any():
                average = np.average(peak_indices, weights=peak_values)
                variance = np.sqrt(np.average((peak_indices-average)**2, weights=peak_values))
                if variance == 0:
                    # we have a singular peak if variance == 0, so set the variance to half of the x/y spacing
                    if peak_index >= 1:
                        variance = np.abs(xdata[peak_index]-xdata[peak_index-1])
                    elif peak_index < len(xdata):
                        variance = np.abs(xdata[peak_index+1]-xdata[peak_index])
                    else:
                        # we have only 1 data point, most RT's fall into this width
                        variance = 0.05
            else:
                variance = 0.05
                average = xdata[peak_index]
            if variance is not None:
                guess.extend([ydata[peak_index], average, variance])
                if bigauss_fit:
                    guess.extend([variance])

        if not guess:
            average = np.average(xdata, weights=ydata)
            variance = np.sqrt(np.average((xdata-average)**2, weights=ydata))
            if variance == 0:
                variance = 0.05
            guess = [max(ydata), np.argmax(ydata), variance]
            if bigauss_fit:
                guess.extend([variance])

        args = (xdata, ydata)
        opts = {'maxiter': 1000}
        fit_func = bigauss_func if bigauss_fit else gauss_func
        routines = ['SLSQP', 'TNC', 'L-BFGS-B', 'SLSQP']
        routine = routines.pop()
        res = optimize.minimize(fit_func, guess, args, method=routine, bounds=bnds, options=opts)#, jac=gauss_jac)
        while not res.success and routines:
            routine = routines.pop()
            res = optimize.minimize(fit_func, guess, args, method=routine, bounds=bnds, options=opts)#, jac=gauss_jac)
        n = len(xdata)
        k = len(res.x)
        bic = n*np.log(res.fun/n)+k+np.log(n)
        res.bic = bic
        fit_accuracy.append((peak_width, bic, res, xdata[row_peaks]))
    # we want to maximize our BIC given our definition
    best_fits = sorted(fit_accuracy, key=itemgetter(1,0), reverse=True)
    return best_fits[0][2:]

def merge_list(starting_list):
    final_list = []
    for i,v in enumerate(starting_list[:-1]):
        if set(v)&set(starting_list[i+1]):
            starting_list[i+1].extend(list(set(v) - set(starting_list[i+1])))
        else:
            final_list.append(v)
    final_list.append(starting_list[-1])
    return final_list

def findPeak(y, srt):
    # check our SNR, if it's low, lessen our window
    left_offset = 1
    right_offset = 2
    lsrt = srt-left_offset if srt-left_offset > 0 else 0
    rsrt = srt+right_offset if srt+right_offset < len(y) else len(y)
    peak = y[srt]
    left = 0
    grad = y[lsrt:rsrt]
    shift = sum(np.sign(np.gradient(grad)))
    shift = 'left' if shift < 0 else 'right'
    ishift = shift
    slope_shifts = 0
    last_slope = -1
    for left in xrange(srt-1, -1, -1):
        val = y[left]
        grad = y[left-2:left+1]
        slope = None
        if len(grad) >= 2:
            slope = sum(np.sign(np.gradient(grad)))
            if slope < 0:
                slope = 'right'
            elif slope > 0:
                slope = 'left'
        if last_slope != -1:
            if last_slope != slope and slope != None:
                slope_shifts += 1
        last_slope = slope
        if ishift == 'right' and ishift == slope:
            break
        if ishift == 'left' and slope == 'right' and slope_shifts > 1:
            break
        if val == 0 or (val > peak and slope != 'right'):
            if val == 0 or shift != 'left':
                break
        elif shift == 'left' and slope != 'right':# slope != right logic: newsl
            shift = None
    right = len(y)
    shift = ishift
    highest_val = peak
    slope_shifts = 0
    last_slope = -1
    for right in xrange(srt+1, len(y)):
        val = y[right]
        grad = y[right:right+3]
        slope = None
        if len(grad) >= 2:
            slope = sum(np.sign(np.gradient(grad)))
            if slope < 0:
                slope = 'right'
            elif slope > 0:
                slope = 'left'
        if last_slope != -1:
            if last_slope != slope and slope != None:
                slope_shifts += 1
        last_slope = slope
        if ishift == 'left' and ishift == slope:
            break
        if ishift == 'right' and slope == 'left' and slope_shifts > 1:
            break
        if val == 0 or (val > peak and slope != 'left'):
            if val > highest_val:
                highest_val = val
            if val == 0 or shift != 'right':
                if val == 0:
                    right += 1
                break
        elif shift == 'right' and slope != 'left':
            shift = None
            peak = highest_val
    return left, right

def findValleys(y, srt):
    peak = y.iloc[srt]
    for left in xrange(srt-1, -1, -1):
        val = y.iloc[left]
        if val == 0 or val > peak:
            break
    right = len(y)
    for right in xrange(srt+1, len(y)):
        val = y.iloc[right]
        if val == 0 or val > peak:
            break
    return left, right

def buildEnvelope(peaks_found=None, isotopes=None, gradient=False, rt_window=None, start_rt=None, silac_label=None):
    isotope_index = {}
    df = pd.DataFrame(columns=rt_window if rt_window is not None else [])
    for silac_isotope, isotope_data in peaks_found.iteritems():
        if isotopes and silac_isotope not in isotopes:
            continue
        temp_data = []
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
                temp_mean = np.average(temp_series.index.values, weights=temp_series.values)
                temp_data.append((temp_mean, {'int': temp_int, 'index': series_index, 'rt': rt}))

        exclude = {}
        data = [i[0] for i in temp_data]
        if len(data) > 4:
            data_median, data_std = np.median(data), np.std(data)
            cutoff = 2*data_std
            cut_fun = lambda i: abs(i-data_median)>cutoff
            exclude = set(filter(cut_fun, data))
        for i in temp_data:
            if i[0] in exclude:
                continue
            temp_int = i[1]['int']
            series_index = i[1]['index']
            rt = i[1]['rt']
            if series_index in df.index:
                if rt in df and not pd.isnull(df.loc[series_index, rt]):
                    df.loc[series_index, rt] += temp_int
                else:
                    df.loc[series_index, rt] = temp_int
            else:
                df = df.append(pd.DataFrame(temp_int, columns=[rt], index=[series_index]))

    # THis is experimental
    if not df.empty:
        delta_df = pd.rolling_apply(df, 2, lambda x: x[1]/x[0]).fillna(0)
        ndfs = []
        indices = df.index
        for index, row in delta_df.T.iterrows():
            exclude = False
            for i,v in enumerate(row):
                if v > 3:
                    exclude = True
                    break
            ndfs.append(df.loc[indices[:i] if exclude else indices,index])
        df = pd.concat(ndfs, axis=1) if ndfs else df

        ndfs = []
        from scipy.signal import argrelmax
        closest = sorted([(i, abs(float(v)-start_rt)) for i,v in enumerate(df.columns)], key=itemgetter(1))[0][0]
        isotope_coms = []
        for index, values in df.iterrows():
            y = values.fillna(0)
            maximas = y.iloc[argrelmax(y.fillna(0).values)[0]]
            if len(maximas) == 0:
                ndfs.append(pd.DataFrame(y))
                continue
            # check for peaks outside our window
            # neighboring_peaks = filter(lambda x: abs(x)>0.3, maximas.index-start_rt)
            # if not neighboring_peaks and len(maximas) <= 1:
            #     ndfs.append(pd.DataFrame(y))
            #     continue
            peaks = sorted([(index, abs(float(index)-start_rt), v) for index, v in maximas.iteritems() if v], key=itemgetter(1))
            # choose the peak closest to our retention time
            # if any([i for i in peaks if i[1]<0.30]):
            #     peaks.sort(key=itemgetter(2), reverse=True)
            #     srt = y.index.searchsorted(peaks[0][0])
            # else:
            srt = find_nearest(y.index.values, peaks[0][0])
            # find the left/right from our peak
            # see if we're increasing to the left, otherwise assume we are to the right
            left, right = findPeak(y, srt)
            peak_data = y.iloc[left:right]
            peak_com = np.average(peak_data.index, weights=peak_data.values)
            isotope_coms.append(peak_com)
            # new logic -- valleys -- don't use
            # if abs(peak_com-start_rt) > 0.25:
                # we are too far away from our starting RT, find the local peak and use it
                # left, right = findValleys(y, srt)
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
        vals = (vals/vals.max()).fillna(0)
        residual = ((theo-vals)**2).sum()
        yield (residual, copy.deepcopy(out))

def find_nearest(array, value):
    return array[find_nearest_index(array, value)]

def find_nearest_index(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx == 0:
        return 0
    elif idx == len(array):
        return -1
    elif math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
        return idx-1
    else:
        return idx

def get_ppm(theoretical, observed):
    return np.abs(float(theoretical)-float(observed))/float(theoretical)

def findEnvelope(df, measured_mz=None, theo_mz=None, max_mz=None, precursor_ppm=5, isotope_ppm=2.5, isotope_ppms=None, charge=2, debug=False,
                 isotope_offset=0, theo_dist=None, label=None, skip_isotopes=None, last_precursor=None, quant_method='integrate', reporter_mode=False):
    # returns the envelope of isotopic peaks as well as micro envelopes  of each individual cluster
    spacing = NEUTRON/float(charge)
    start_mz = measured_mz if isotope_offset == 0 else measured_mz+isotope_offset*NEUTRON/float(charge)
    if max_mz is not None:
        max_mz = max_mz-spacing*0.9 if isotope_offset == 0 else max_mz+isotope_offset*NEUTRON*0.9/float(charge)
    if isotope_ppms is None:
        isotope_ppms = {}
    tolerance = isotope_ppms.get(0, precursor_ppm)/1000000.0
    env_dict, micro_dict, ppm_dict = OrderedDict(),OrderedDict(),OrderedDict()
    empty_dict = {'envelope': env_dict, 'micro_envelopes': micro_dict, 'ppms': ppm_dict}

    non_empty = df[df>0].dropna()
    if len(non_empty) == 0:
        return empty_dict
    start = find_nearest(non_empty.index.values, start_mz)
    attempts = 0


    isotope_index = 0
    # This is purposefully verbose to be more explicit
    while get_ppm(start_mz, start) > tolerance:
        # let's try using our theoretical mass
        start = find_nearest(non_empty.index.values, theo_mz)
        if get_ppm(theo_mz, start) > tolerance:
            # let's check our last boundary
            if last_precursor is not None:
                start = find_nearest(non_empty.index.values, last_precursor)
                if get_ppm(last_precursor, start) > tolerance:
                    # repeat all of that for the next isotopic index
                    start_mz += spacing
                    theo_mz += spacing
                    last_precursor += spacing
                    isotope_index += 1
                else:
                    start_mz = last_precursor
                    break
            else:
                start_mz += spacing
                theo_mz += spacing
                isotope_index += 1
        else:
            start_mz = theo_mz
            break
        tolerance = isotope_ppms.get(isotope_index, isotope_ppm)/1000000.0
        if isotope_index == 2 or (max_mz is not None and start >= max_mz):
            return empty_dict

    # we reset the isotope_index back to 1 here, so on the subsequent steps we aren't looking
    isotope_index += isotope_offset
    # find the com
    start = findMicro(df, find_nearest_index(df.index.values, start), ppm=tolerance, start_mz=start_mz, quant_method=quant_method)
    start_error = start['error']

    if 'params' in start:
        if start['error'] > tolerance:
            start = last_precursor if last_precursor is not None else start_mz
        else:
            start = start['params'][1]
    else:
        return empty_dict

    valid_locations2 = OrderedDict()
    valid_locations2[isotope_index] = [(isotope_index, start)]

    if not reporter_mode:
        isotope_index += 1
        pos = find_nearest_index(non_empty.index.values, start)+1
        offset = isotope_index*spacing
        df_len = non_empty.shape[0]
        last_displacement = None
        valid_locations = []
        tolerance = isotope_ppms.get(isotope_index, isotope_ppm)/1000000.0

        while pos < df_len:
            # search for the ppm error until it rises again, we select the minima and if this minima is
            # outside our ppm error, we stop the expansion of our isotopic cluster
            current_loc = non_empty.index[pos]
            if max_mz is not None and current_loc >= max_mz:
                if not valid_locations:
                    break
                displacement = last_displacement+tolerance if last_displacement is not None else tolerance*2
            else:
                displacement = get_ppm(start+offset, current_loc)
            if debug:
                print pos, start, current_loc, displacement, last_displacement, displacement > last_displacement, last_displacement < tolerance, isotope_index, offset
            if displacement < tolerance:
                valid_locations.append((displacement, current_loc))
            if valid_locations and displacement > last_displacement:
                # pick the largest peak within our error tolerance
                valid_locations2[isotope_index] = valid_locations
                isotope_index += 1
                tolerance = isotope_ppms.get(isotope_index, isotope_ppm)/1000000.0
                offset = spacing*isotope_index
                displacement = get_ppm(start+offset, current_loc)
                valid_locations = []
            elif last_displacement is not None and displacement > last_displacement and not valid_locations:
                break
            last_displacement = displacement
            pos += 1

    #combine any overlapping micro envelopes
    #final_micros = self.merge_list(micro_dict)
    valid_keys = sorted(set(valid_locations2.keys()).intersection(theo_dist.keys() if theo_dist is not None else valid_locations2.keys()))
    # This attempts to use a diophantine equation to match the clusters to the theoretical distribution of isotopes
    #valid_vals = [j[1] for i in valid_keys for j in valid_locations2[i]]
    #valid_theor = pd.Series([theo_dist[i] for i in valid_keys])
    #valid_theor = valid_theor/valid_theor.max()
    #best_locations = sorted(looper(selected=valid_vals, df=df, theo=valid_theor), key=itemgetter(0))[0][1]

    best_locations = [sorted(valid_locations2[i], key=itemgetter(0))[0][1] for i in valid_keys]

    for index, isotope_index in enumerate(valid_keys):
        if skip_isotopes is not None and isotope_index in skip_isotopes:
            continue
        largest_loc = best_locations[index]
        micro_index = find_nearest_index(df.index.values, largest_loc)
        # if micro_index == 0:
        #     pass
        isotope_tolerance = isotope_ppms.get(isotope_index, isotope_ppm)/1000000.0
        precursor_tolerance = isotope_ppms.get(0, precursor_ppm)/1000000.0
        micro_bounds = findMicro(df, micro_index, ppm=precursor_tolerance if isotope_index == 0 else isotope_tolerance,
                                 calc_start_mz=start_mz, start_mz=start, isotope=isotope_index, spacing=spacing, quant_method=quant_method)
        if isotope_index == 0:
            micro_bounds['error'] = start_error

        micro_dict[isotope_index] = micro_bounds
        env_dict[isotope_index] = micro_index
        ppm_dict[isotope_index] = micro_bounds.get('error')


    # in all cases, the envelope is going to be either monotonically decreasing, or a parabola (-x^2)
    isotope_pattern = [(isotope_index, isotope_dict['int']) for isotope_index, isotope_dict in micro_dict.items()]
    # are we monotonically decreasing?
    remove = False
    if len(isotope_pattern) > 2:
        # check if the 2nd isotope is smaller than the first. This is a classical case looking like:
        #
        #  |
        #  |  |
        #  |  |  |
        #  |  |  |  |

        if isotope_pattern[1][1] < isotope_pattern[0][1]:
            # we are, check this trend holds and remove isotopes it fails for
            for i,j in zip(isotope_pattern, isotope_pattern[1:]):
                if j[1]*0.9 > i[1]:
                    # the pattern broke, remove isotopes beyond this point
                    remove = True
                if remove:
                    env_dict.pop(j[0])
                    micro_dict.pop(j[0])
                    ppm_dict.pop(j[0])

        # check if the 2nd isotope is larger than the first. This is a case looking like:
        #
        #
        #     |  |
        #     |  |
        #  |  |  |  |

        elif isotope_pattern[1][1] > isotope_pattern[0][1]:
            shift = False
            for i,j in zip(isotope_pattern, isotope_pattern[1:]):
                if shift and j[1]*0.9 > i[1]:
                    remove = True
                elif shift is False and j[1] < i[1]*0.9:
                    if shift:
                        remove = True
                    else:
                        shift = True
                if remove:
                    env_dict.pop(j[0])
                    micro_dict.pop(j[0])
                    ppm_dict.pop(j[0])

    return {'envelope': env_dict, 'micro_envelopes': micro_dict, 'ppms': ppm_dict}
