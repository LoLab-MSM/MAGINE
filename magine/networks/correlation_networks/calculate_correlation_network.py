'''


# import pathos.multiprocessing as mp
import time
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from scipy import ndimage
from scipy.stats import chi2_contingency
from scipy.stats.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests

np.set_printoptions(linewidth=300)


global_data = None
samples = None
import scipy.special as special
from scipy.stats.stats import distributions

from numba import float64, jit
import bottleneck as bn

@jit((float64, float64, float64), nopython=True, nogil=True)
def calc_prob(n, dsq, df):
    denom = n * (n ** 2 - 1) / 6.
    rho = 1. - dsq / denom
    val = (rho + 1.0) * (1.0 - rho)
    if val == 0.0:
        return rho, np.inf
    d = df / val
    return rho, d


@jit(float64(float64, float64), nopython=True, nogil=True)
def calc_val_2(df, t):
    return df / (df + t * t)


def spearmanr_optimized(x, y):
    n = x.size
    xy = x * y
    nans = np.isnan(xy)
    x[nans] = np.nan
    y[nans] = np.nan
    n -= nans.sum()
    df = n - 2
    if df < 0:
        raise ValueError("The input must have at least 3 entries!")

    # Gets the ranks and rank differences
    rankx2 = bn.nanrankdata(x)
    ranky2 = bn.nanrankdata(y)
    dsq = bn.nansum((rankx2 - ranky2) ** 2)
    rho, d = calc_prob(n, dsq, df)
    if d == np.inf or d < 0:
        prob = 0
    else:
        t = np.sqrt(d) * rho
        value_2 = calc_val_2(df, t)
        prob = special.betainc(0.5 * df, 0.5, value_2)
    return rho, prob


def calculate_spearmanr(a, b):
    # new code
    # return spearmanr_optimized(a, b)
    if bn.anynan(a):
        return spearmanr_optimized(a, b)
    elif bn.anynan(b):
        return spearmanr_optimized(a, b)
    ar = np.apply_along_axis(bn.rankdata, 0, a)
    br = np.apply_along_axis(bn.rankdata, 0, b)
    n = a.shape[0]
    rs = np.corrcoef(ar, br, rowvar=0)
    olderr = np.seterr(divide='ignore')  # rs can have elements equal to 1
    try:
        # clip the small negative values possibly caused by rounding
        # errors before taking the square root
        t = rs * np.sqrt(((n - 2) / ((rs + 1.0) * (1.0 - rs))).clip(0))
    finally:
        np.seterr(**olderr)

    prob = 2 * distributions.t.sf(np.abs(t), n - 2)

    if rs.shape == (2, 2):
        return rs[1, 0], prob[1, 0]
    else:
        return rs, prob

def calculate_spearman(pos):
    tmp = global_data[samples[pos, :], :]
    n_nans = np.isfinite(tmp[0, :] * tmp[1, :]).sum()

    if n_nans > 2:
        # cov_v, p_value = calculate_spearmanr(tmp[0, :], tmp[1, :])
        cov_v, p_value = spearmanr(tmp[0, :], tmp[1, :], nan_policy='omit')
        return cov_v, p_value
    else:
        return np.nan, np.nan


EPS = np.finfo(float).eps


def mi_1(d):
    """
    Computes (normalized) mutual information between two 1D variate from a
    joint histogram.
    Parameters
    ----------
    x : 1D array
        first variable
    y : 1D array
        second variable
    sigma: float
        sigma for Gaussian smoothing of the joint histogram
    Returns
    -------
    nmi: float
        the computed similariy measure
    """
    bins = (16, 16)
    x, y, _ = d
    sigma = 1
    normalized = False
    jh = np.histogram2d(x, y, bins=bins)[0]

    # smooth the jh with a gaussian filter of given sigma
    ndimage.gaussian_filter(jh, sigma=sigma, mode='constant',
                            output=jh)

    # compute marginal histograms
    jh = jh + EPS
    sh = np.sum(jh)
    jh = jh / sh
    s1 = np.sum(jh, axis=0).reshape((-1, jh.shape[0]))
    s2 = np.sum(jh, axis=1).reshape((jh.shape[1], -1))

    # Normalised Mutual Information of:
    # Studholme,  jhill & jhawkes (1998).
    # "A normalized entropy measure of 3-D medical image alignment".
    # in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
    if normalized:
        mi = ((np.sum(s1 * np.log(s1)) + np.sum(s2 * np.log(s2)))
              / np.sum(jh * np.log(jh))) - 1
    else:
        mi = (np.sum(jh * np.log(jh)) - np.sum(s1 * np.log(s1))
              - np.sum(s2 * np.log(s2)))

    return mi


def mi_2(d):
    x, y, bins = d
    c_xy = np.histogram2d(x, y, bins)[0]
    g, p, dof, expected = chi2_contingency(c_xy, lambda_="log-likelihood")
    mi = 0.5 * g / c_xy.sum()
    print(mi)
    return mi


def mi_3(d):
    x, y, bins = d
    c_xy = np.histogram2d(x, y, bins)[0]
    c_x = np.histogram(x, bins)[0]
    c_y = np.histogram(y, bins)[0]

    h_x = shan_entropy(c_x)
    h_y = shan_entropy(c_y)
    h_xy = shan_entropy(c_xy)

    mi = h_x + h_y - h_xy
    return mi


def shan_entropy(c):
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    return -sum(c_normalized * np.log2(c_normalized))


def correlation_sampling(data, names, save_name, create_plots=True, pool=None):
    global global_data
    global_data = data
    total_values = len(global_data)

    # Create pairwise samples
    global samples

    samples = np.array(list(combinations(range(len(data)), 2)))
    n_samples = int(total_values * (total_values - 1) / 2)
    print("Sampling all = {}".format(n_samples))
    print('Starting to calculate spearman correlations')
    samples_range = range(0, n_samples)
    start_time = time.time()

    all_sample = []
    n_bins = 3
    for i, j in list(combinations(range(len(data)), 2)):
        all_sample.append([data[i], data[j], n_bins])

    # all_sample = [[data[i], data[j], n_bins] for i, j in
    #               combinations(range(len(data)), 2)]

    if pool is None:
        print("Running in serial")
        # x = map(calculate_spearman, samples_range)
        x = map(mi_1, all_sample)
        # x = map(calc_MI, all_sample)

    else:
        print("Running in parallel")
        # x = pool.map(calculate_spearman, samples_range)
        # x = pool.map(calc_mi, all_sample)
        x = pool.map(mi_1, all_sample)
        pool.close()
        pool.join()
    time_taken = time.time() - start_time
    print('Done with {} samples in {}'.format(n_samples, time_taken))
    x = np.array(x)
    use_mi = True

    if use_mi:
        output = pd.DataFrame()
        output['species_1'] = names[samples[:, 0]]
        output['species_2'] = names[samples[:, 1]]
        print(np.shape(x))
        print(output.shape)
        output['mi'] = x
        print(output.mi.min())
        print(output.mi.max())
        plt.hist(output.mi)
        plt.show()
        quit()
    else:
        co_var = x[:, 0]
        pvals = x[:, 1]
        nan_index = np.isfinite(pvals)
        co_var = co_var[nan_index]
        pvals = pvals[nan_index]

        samples = samples[nan_index]

        # Create output file
        output = pd.DataFrame()
        output['species_1'] = names[samples[:, 0]]
        output['species_2'] = names[samples[:, 1]]

        # adjust pvalues for multiple hypothesis testing
        adj_pvalues = multipletests(pvals=pvals, method='fdr_bh')

        output['adj_pvalue'] = adj_pvalues[1]
        output['co_var'] = co_var
        output.to_csv(
            '{}_correlation_sampling_output.csv.gz'.format(save_name),
            index=False, compression='gzip')
    # print('Saved file!')

    if create_plots:
        print("Plotting output")
        fig = plt.figure(figsize=(8, 4))
        bins = np.linspace(-1, 1, 41)
        ax1 = fig.add_subplot(121)
        ax1.hist(np.array(output.co_var), bins=bins)

        output2 = output[output['adj_pvalue'] <= 0.05]
        ax1.hist(np.array(output2.co_var), bins=bins, color='red')
        ax1.set_xlim(-1.1, 1.1)
        plt.xlabel('Correlation', fontsize=20)
        start, end = ax1.get_ylim()
        ax1.yaxis.set_ticks(np.linspace(start, end, 4))

        ax2 = fig.add_subplot(122)
        bins = np.linspace(0, 1, 21)
        ax2.hist(np.array(output.adj_pvalue), bins=bins)
        ax2.hist(np.array(output2.adj_pvalue), bins=bins, color='red')
        ax2.set_xlim(0, 1)
        plt.xlabel('Adjusted P value', fontsize=20)

        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1g'))
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1g'))
        start, end = ax2.get_ylim()
        ax2.yaxis.set_ticks(np.linspace(start, end, 4))
        ax1.tick_params(axis='y', labelsize=16)
        ax1.tick_params(axis='x', labelsize=16)
        ax2.tick_params(axis='y', labelsize=16)
        ax2.tick_params(axis='x', labelsize=16)
        plt.tight_layout(h_pad=0.15)
        _same_name = '{}_correlation_and_pvalue_histogram.png'.format(
            save_name)
        print("Saving to  {}".format(_same_name))
        plt.savefig(_same_name, bbox_inches='tight')
        plt.close()
        print("Done 2")
    return output
'''
