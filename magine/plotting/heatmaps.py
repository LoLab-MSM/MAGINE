import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def heatmap_from_array(data, convert_to_log=False, yticklabels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False,
                       fig_size=(6, 4)):
    """

    Parameters
    ----------


    data : pandas.DataFrame
    convert_to_log : bool
    yticklabels : list_like
    columns : str
        Name of columns of df for pivotn
    index : str
        Name of index of df for pivot
    values : str
        Name of values of df for pivot
    cluster_col : bool
        Cluster the data using searborn.clustermap
    cluster_row : bool
        Cluster the data using searborn.clustermap
    div_colors : bool
        Use divergent colors for plotting
    fig_size : tuple
        Size of figure, passed to matplotlib/seaborn
    Returns
    -------

    """
    array = _pivot_table(data, convert_to_log, columns=columns, index=index,
                         values=values)
    if div_colors:
        pal = sns.color_palette("coolwarm", 7)
        center = 0
    else:
        pal = sns.light_palette("purple", as_cmap=True)
        center = None

    if cluster_col or cluster_row:

        fig = sns.clustermap(array, yticklabels=yticklabels, figsize=fig_size,
                             col_cluster=cluster_col, row_cluster=cluster_row,
                             center=center, cmap=pal)
    else:
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(111)
        sns.heatmap(array, ax=ax, yticklabels=yticklabels, cmap=pal,
                    center=center)

    return fig


def _pivot_table(data, convert_to_log, index, columns, values):
    d_copy = data.copy()
    if convert_to_log:
        d_copy = _log2_normalize_df(d_copy, values)

    # remove human suffix from term names
    if index == 'term_name':
        d_copy[index] = d_copy.apply(_cut_word, axis=1)

    array = pd.pivot_table(d_copy, index=index,
                           columns=columns,
                           values=values)
    array.fillna(0, inplace=True)

    array.sort_values(by=list(sorted(d_copy[columns].unique())),
                      ascending=False, inplace=True)
    return array


def _cut_word(row):
    term_name = row['term_name']
    if len(term_name.split('_Homo sapiens')):
        return term_name.split('_Homo sapiens')[0]
    else:
        return term_name


def _convert_to_log(data):
    d_copy = data.copy()
    above = d_copy > 0.
    below = d_copy < 0.
    d_copy[above] = np.log2(d_copy[above])
    d_copy[below] = -1 * np.log2(-1 * d_copy[below])
    return d_copy


def _log2_normalize_df(df, fold_change):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of fold changes
    fold_change : str
        column that contains fold change values

    Returns
    -------

    """
    tmp_df = df.copy()
    greater_than = tmp_df[fold_change] > 0
    less_than = tmp_df[fold_change] < 0
    tmp_df.loc[greater_than, fold_change] = \
        np.log2(tmp_df[greater_than][fold_change])
    tmp_df.loc[less_than, fold_change] = \
        -np.log2(-tmp_df[less_than][fold_change])
    return tmp_df
