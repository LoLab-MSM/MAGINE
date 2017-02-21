import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch


def create_heatmaps_pvalue_xs(data, save_name, labels=None, mark_pvalues=False,
                              trim_nodes=False, n_hits_per_time=None):
    from magine.ontology.ontology_tools import check_term_list

    if isinstance(data, str):
        data = pd.read_csv(data, index_col=0)
    assert 'pvalue' in data.columns
    assert 'enrichment_score' in data.columns
    data['genes'] = data['genes'].astype(set)

    data = data[data['ref'] >= 5]
    data.loc[data['enrichment_score'] > 20, 'enrichment_score'] = 20

    n_samples = range(len(labels))
    tmp = data[data['pvalue'] < 0.05]

    def find_n_go_terms(n_terms):
        terms = set(tmp.head(n_terms).index)
        terms = check_term_list(terms)
        return terms

    if n_hits_per_time is not None:
        tmp = pd.pivot_table(tmp, index=['GO_id', ], columns='sample_index')
        list_all_go = set()
        for i in n_samples:
            tmp = tmp.sort_values(by=('enrichment_score', i), ascending=False)
            list_of_go = find_n_go_terms(n_hits_per_time)
            count = 1
            while len(list_of_go) < n_hits_per_time:
                list_of_go = find_n_go_terms(n_hits_per_time + count)
                count += 1

            if trim_nodes:
                list_all_go = check_term_list(list_all_go)
            print(list_of_go)
            list_all_go.update(list_of_go)

        if trim_nodes:
            list_all_go = check_term_list(list_all_go)
        data = data[data['GO_id'].isin(list_all_go)]

    tmp2 = pd.pivot_table(data, index=['GO_id', 'GO_name'],
                          columns='sample_index',
                          aggfunc=lambda gene: ' '.join(gene))

    tmp = pd.pivot_table(data, index=['GO_id', 'GO_name'],
                         columns='sample_index', )
    # print(tmp2.head(10))
    enrichment_list = [('enrichment_score', i) for i in n_samples]
    tmp = tmp.sort_values(by=enrichment_list, ascending=False)

    # Want to show both enrichment as a panel as pvalie
    pvalue = tmp['pvalue'].fillna(1).as_matrix()
    enrichment_value = tmp['enrichment_score'].as_matrix()
    fig = plt.figure(figsize=(14, 7))

    # Enrichment panel
    ax1 = fig.add_subplot(111)
    plt.imshow(enrichment_value, aspect='auto',
               interpolation='nearest', cmap=plt.get_cmap('YlOrRd'),
               # extent=[0, len(labels), 0, len(tmp.index)],
               )
    y_array = np.arange(0, len(tmp.index), 1)
    x_array = np.arange(0, len(labels), 1)
    if mark_pvalues:
        # text portion
        x, y = np.meshgrid(x_array, y_array)

        for x_val, y_val in zip(x.flatten(), y.flatten()):
            if pvalue[y_val, x_val] < 0.05:
                ax1.text(x_val, y_val, '*', va='center', ha='center',
                         fontdict={
                             'color':    'black',
                             'weight':   'normal',
                             'fontsize': 18
                         })

    ax1.tick_params(axis='both', direction='out')
    ax1.set_xticks(x_array)
    ax1.set_xticklabels(labels, fontsize=16)

    if n_hits_per_time is not None:
        ax1.set_yticks(y_array)
        ax1.set_yticklabels(tmp.index, fontsize=16)
    else:
        ax1.set_yticks([])
        ax1.set_yticklabels([])

    plt.colorbar()
    plt.savefig('{}_heatplot.png'.format(save_name), bbox_inches='tight',
                transparent=True
                )
    plt.savefig('{}_heatplot.pdf'.format(save_name), bbox_inches='tight',
                transparent=True)
    plt.close()


def plot_heatmap(enrichment_array, names_col, labels, global_go,
                 start=0, stop=None, savename='tmp', out_dir='.'):
    """
    General function to create a heatmap of enrichment array data

    Parameters
    ----------
    enrichment_array : array like
        numpy array of enrichment values
    names_col : array like
        array of names
    labels : list
        list of column names
    start : int
        starting index for size of array to be shown
    stop : int or None
        stopping index for size of array to be shown
    savename : str
        name of figure to be saved as

    Returns
    -------

    """

    if stop is None:
        matrix = enrichment_array[start:, :]
    else:
        matrix = enrichment_array[start:stop, :]

    size_of_data = np.shape(enrichment_array)[1]
    length_matrix = len(matrix)
    fig = plt.figure(figsize=(14, 20))
    ax1 = fig.add_subplot(111)

    im = ax1.imshow(matrix, aspect=.25, interpolation='nearest',
                    extent=(0, size_of_data, 0, length_matrix + 1),
                    origin='lower')

    names_2 = []
    for i in names_col[start:stop]:
        names_2.append(global_go[i])

    x_ticks = np.linspace(.5, size_of_data - .5, size_of_data)
    y_ticks = np.linspace(.5, length_matrix + .5, length_matrix)

    plt.yticks(y_ticks, names_2, fontsize=16)
    if labels:
        plt.xticks(x_ticks, labels, fontsize=16, rotation='90')
    else:
        print("Provide labels")
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.savefig(os.path.join(out_dir, 'Figures', '%s.pdf' % savename),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_heatmap_cluster(tmp_array, names, savename, out_dir):
    fig = plt.figure()
    axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
    linkage = sch.linkage(tmp_array, method='centroid')
    dendrogram = sch.dendrogram(linkage, orientation='right')
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
    index = dendrogram['leaves']
    tmp_array = tmp_array[index, :]
    names_sorted = names[index]
    im = axmatrix.matshow(tmp_array, aspect='auto', origin='lower')
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    plt.colorbar(im, cax=axcolor)
    fig.savefig(os.path.join(out_dir, 'Figures',
                             '%s_clustered.pdf' % savename))
    plt.close()
    return tmp_array, names_sorted
