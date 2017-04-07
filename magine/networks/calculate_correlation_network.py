try:
    import bottleneck as bn
except ImportError:
    print("Must install bottleneck to calculate correlation network")
try:
    from rpy2.robjects.packages import importr
    from rpy2.robjects.vectors import FloatVector
except ImportError:
    print("Must install R and rpy2 to calculate correlation network")
import time
from itertools import combinations
import matplotlib.pyplot
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.special as special
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numba import float64, jit
from scipy.stats.stats import distributions

matplotlib.use('Agg')
np.set_printoptions(linewidth=300)
plt = matplotlib.pyplot
stats = importr('stats')


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


# @profile
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


# @profile

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
        cov_v, p_value = calculate_spearmanr(tmp[0, :], tmp[1, :])
        # cov_v, p_value = spearmanr(tmp[0, :], tmp[1, :], nan_policy='omit')
        return cov_v, p_value
    else:
        return np.nan, np.nan


def correlation_sampling(data, names, n_samples, save_name, sample_all=False,
                         create_plots=True):
    global global_data
    global_data = data
    total_values = len(global_data)

    # Create pairwise samples
    global samples
    if sample_all:
        samples = np.array(list(combinations(xrange(n_samples), 2)))
        total_possible_samples = total_values * (total_values - 1) / 2
        n_samples = total_possible_samples
        # print("Sampling all = {}".format(total_possible_samples))
    else:
        samples = np.array(list(combinations(xrange(len(data)), 2)))
        samples = samples[np.random.choice(len(samples), size=n_samples,
                                           replace=False), :]
        # print("Sampling = {}".format(len(samples)))

    # start_time = time.time()
    # for i in range(5000):
    #     calculate_spearman(i)
    # time_taken = time.time() - start_time
    # print('Done with {} samples in {}'.format(n_samples, time_taken))
    # quit()
    pool = mp.Pool(processes=6)
    print('Starting to calculate spearman correlations')

    samples_range = range(0, n_samples)
    start_time = time.time()
    x = pool.map(calculate_spearman, samples_range, n_samples / 100)
    pool.close()
    pool.join()
    time_taken = time.time() - start_time
    # print('Done calculating correlation')
    # return time_taken
    print('Done with {} samples in {}'.format(n_samples, time_taken))
    x = np.array(x)
    co_var = x[:, 0]
    pvals = x[:, 1]
    nan_index = np.isfinite(pvals)
    co_var = co_var[nan_index]
    pvals = pvals[nan_index]

    samples = samples[nan_index]
    print(np.shape(samples))
    # Create output file
    output = pd.DataFrame()
    output['species_1'] = names[samples[:, 0]]
    output['species_2'] = names[samples[:, 1]]

    # adjust pvalues for multiple hypothesis testing
    adj_pvalues = stats.p_adjust(FloatVector(pvals), method='BH')
    output['adj_pvalue'] = adj_pvalues
    output['co_var'] = co_var
    output.to_csv('{}_correlation_sampling_output.csv.gz'.format(save_name),
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
        plt.savefig(
                '{}_correlation_and_pvalue_histogram.png'.format(save_name),
                bbox_inches='tight')
        plt.close()
        print("Done 2")


def create_heatmap(data, savename, xlabels=None):
    data = np.nan_to_num(data)
    # data = np.log2(data)
    threshold = 2.
    data[data > threshold] = threshold
    data[data < -1 * threshold] = -1 * threshold
    name = 'bwr'

    size_of_data = np.shape(data)[1]
    length_matrix = len(data)
    x_ticks = np.linspace(.5, size_of_data - .5, size_of_data)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    im = plt.imshow(data, interpolation='nearest', aspect='auto',
                    extent=(0, size_of_data, 0, length_matrix + 1),
                    cmap=plt.get_cmap(name))

    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    if xlabels is not None:
        plt.xticks(x_ticks, xlabels, fontsize=16, rotation='90')
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax1)

    plt.savefig('{}_heatmap.png'.format(savename), bbox_inches='tight')
    plt.close()

    # distanceMatrix = dist.pdist(data, 'euclidean')
    # linkageMatrix = sch.linkage(distanceMatrix, method='centroid')
    linkageMatrix = sch.linkage(data, method='centroid')
    heatmapOrder = sch.leaves_list(linkageMatrix)

    data = data[heatmapOrder, :]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    im = plt.imshow(data, interpolation='nearest', aspect='auto',
                    extent=(0, size_of_data, 0, length_matrix + 1),
                    cmap=plt.get_cmap(name))

    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    if xlabels is not None:
        plt.xticks(x_ticks, xlabels, fontsize=16, rotation='90')
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax1)

    plt.savefig('{}_heatmap_clustered.png'.format(savename),
                bbox_inches='tight')
