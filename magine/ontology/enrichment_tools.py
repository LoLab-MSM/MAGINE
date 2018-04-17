import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch


def jaccard_index(first_set, second_set):
    """
    Computes the similarity between two sets.
        https://en.wikipedia.org/wiki/Jaccard_index

    Parameters
    ----------
    first_set : set
    second_set : set



    Returns
    -------
    index : float


    References
    ----------
    .. [1] `Wikipedia entry for the Jaccard index
           <https://en.wikipedia.org/wiki/Jaccard_index>`_
    """
    set1 = set(first_set)
    set2 = set(second_set)

    return float(len(set1.intersection(set2))) / len(set1.union(set2))


def filter_rows(local_df, column, options):
    """
    Filters a pandas dataframe provides a column and filter selection.

    Parameters
    ----------
    local_df : pd.DataFrame
    column : str
    options : str, list
        Can be a single entry or a list

    Returns
    -------
    pd.DataFrame
    """
    copy_df = local_df.copy()
    valid_opts = list(local_df[column].unique())
    if isinstance(options, str):
        if options not in valid_opts:
            print('{} not in {}'.format(options, valid_opts))
        else:
            copy_df = local_df[local_df[column] == options]
    elif isinstance(options, list):
        for i in options:
            if i not in valid_opts:
                print('{} not in {}'.format(i, valid_opts))
        copy_df = local_df[local_df[column].isin(options)]
    return copy_df


def filter_dataframe(df, p_value=None, combined_score=None, db=None,
                     sample_id=None, category=None, rank=None):
    """
    Filters an enrichment array.

    This is an aggregate function that allows ones to filter an entire
    dataframe with a single function call.

    Parameters
    ----------
    df : pd.DataFrame
    p_value : float
        filters all values less than or equal
    combined_score : float
        filters all values greater than or equal
    db : str, list

    sample_id : str, list
    category : str, list
    rank : int

    Returns
    -------

    """
    copy_df = df.copy()
    if p_value is not None:
        assert isinstance(p_value, float)
        copy_df = copy_df[copy_df['adj_p_value'] <= p_value]
    if combined_score is not None:
        assert isinstance(combined_score, float)
        copy_df = copy_df[copy_df['combined_score'] >= combined_score]
    if isinstance(rank, (int, float)):
        copy_df = copy_df[copy_df['rank'] <= rank]
    if db is not None:
        copy_df = filter_rows(copy_df, 'db', db)
    if sample_id is not None:
        copy_df = filter_rows(copy_df, 'sample_id', sample_id)
    if category is not None:
        copy_df = filter_rows(copy_df, 'category', category)
    copy_df.sort_values('combined_score', inplace=True, ascending=False)
    return copy_df


def term_to_genes(df, term):
    genes = df[df['term_name'] == term]['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').as_matrix()))


def filter_based_on_words(df, words):
    if isinstance(words, str):
        words = [words]
    return df[df['term_name'].str.lower().str.contains('|'.join(words))]


def all_genes_from_df(df):
    genes = df['genes']
    return set(itertools.chain.from_iterable(genes.str.split(',').as_matrix()))


def filter_similar_terms(data, threshold=0.75, verbose=False,
                         sort_by='combined_score'):
    if sort_by == 'rank':
        ascending = True
    else:
        ascending = False
    data.sort_values(sort_by, inplace=True, ascending=ascending)
    data_copy = data.copy()

    to_remove = _terms_to_remove(data_copy, threshold, verbose)

    data_copy = data_copy[~data_copy['term_name'].isin(to_remove)]
    print("Number of rows went from {} to {}".format(data.shape[0],
                                                     data_copy.shape[0]))
    return data_copy


def find_similar_terms(term, df):
    first_genes = term_to_genes(df, term)
    rest_of_df = df[~(df['term_name'] == term)]
    array = rest_of_df[['term_name', 'genes']].as_matrix()
    dist_m = [None] * len(array)
    for j in range(len(array)):
        name = array[j, 0]
        c = jaccard_index(first_genes, set(array[j, 1].split(',')))
        dist_m[j] = [name, c]

    dist_m = pd.DataFrame(dist_m, columns=['term_name', 'similarity_score'])
    dist_m.sort_values('similarity_score', inplace=True, ascending=False)
    return dist_m


def dist_matrix(df, fig_size=(8, 8)):
    array = df[['term_name', 'genes']].as_matrix()
    n_dim = len(array)
    dist_mat = np.empty((n_dim, n_dim), dtype=float)
    lookup = dict()
    for i, ac in enumerate(array):
        first_genes = set(array[i, 1].split(','))
        for j, bc in enumerate(array):
            if i > j:
                continue
            term_name = array[j, 0]
            if term_name in lookup:
                second_genes = lookup[term_name]
            else:
                second_genes = set(array[j, 1].split(','))
                lookup[term_name] = second_genes
            c = jaccard_index(first_genes, second_genes)
            dist_mat[i, j] = c
            dist_mat[j, i] = c
    names = array[:, 0]

    return cluster_distance_mat(dist_mat, names, fig_size)


def dist_matrix2(df, fig_size=(8, 8)):
    names = df['term_name'].unique()
    n_dim = len(names)
    dist_mat = np.empty((n_dim, n_dim), dtype=float)
    for i, n1 in enumerate(names):
        first_genes = term_to_genes(df, n1)
        for j, n2 in enumerate(names):
            if i > j:
                continue
            second_genes = term_to_genes(df, n2)
            c = jaccard_index(first_genes, second_genes)
            dist_mat[i, j] = c
            dist_mat[j, i] = c
    return cluster_distance_mat(dist_mat, names, fig_size)


def cluster_distance_mat(dist_mat, names, fig_size=(8, 8)):
    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=fig_size)

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.2])
    Y = sch.linkage(dist_mat, method='average')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])

    # reorder matrix
    idx1 = Z2['leaves']
    dist_mat = dist_mat[idx1, :]
    dist_mat = dist_mat[:, idx1]
    names = names[idx1]

    # create figure
    im = axmatrix.matshow(dist_mat, aspect='auto', origin='lower',
                          cmap=plt.cm.Reds, vmin=0, vmax=1)

    # add xtick labels
    axmatrix.set_xticks(range(len(names)))
    axmatrix.set_xticklabels(names, minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    plt.xticks(rotation=90, fontsize=8)

    # add ytick labels
    axmatrix.set_yticks(range(len(names)))
    axmatrix.set_yticklabels(names, minor=False)
    axmatrix.yaxis.set_label_position('left')
    axmatrix.yaxis.tick_left()
    plt.yticks(rotation=0, fontsize=8)

    # add colorbar
    axcolor = fig.add_axes([0.94, 0.1, 0.02, 0.6])
    plt.colorbar(im, cax=axcolor)

    return fig


def remove_redundant(data, threshold=0.75, verbose=False,
                     sort_by='combined_score'):
    if sort_by == 'rank':
        ascending = True
    else:
        ascending = False
    data.sort_values(sort_by, inplace=True, ascending=ascending)

    data_copy = data.copy()
    sample_ids = list(sorted(data['sample_id'].unique()))

    to_remove = set()
    for i in sample_ids:
        tmp = data_copy[data_copy['sample_id'] == i]
        redundant_terms = _terms_to_remove(tmp, threshold, verbose)
        to_remove.update(redundant_terms)

    data_copy = data_copy[~data_copy['term_name'].isin(to_remove)]
    print("Number of rows went from {} to {}".format(data.shape[0],
                                                     data_copy.shape[0]))

    return data_copy


def _terms_to_remove(data, threshold=0.75, verbose=False):

    array = data[['term_name', 'genes']].as_matrix()
    lookup = dict()
    to_remove = set()
    for j in range(len(data) - 2):
        top_term_name = array[j, 0]
        top = set(array[j, 1].split(','))
        if verbose:
            print("Finding matches for {}".format(array[j, 0]))
        if top_term_name in to_remove:
            continue
        for i in range(j + 1, len(data)):
            term_name = array[i, 0]

            if top_term_name == term_name:
                continue
            if term_name in lookup:
                new_set = lookup[term_name]
            else:
                new_set = set(array[i, 1].split(','))
                lookup[term_name] = new_set

            if len(new_set.difference(top)) == 0:
                to_remove.add(term_name)
                continue
            score = jaccard_index(top, new_set)
            if verbose:
                print("\tScore for {} is {:.3f}".format(array[i, 0], score))
            if score > threshold:
                to_remove.add(term_name)
                if verbose:
                    print("\t\tRemoving {}".format(term_name))
    return to_remove


if __name__ == '__main__':
    a = {0, 1, 2}
    b = {0, 1, 3}
    print(b.difference(a))
    jaccard_index(a, b)
