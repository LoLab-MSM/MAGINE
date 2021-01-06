import warnings
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as sch
import seaborn as sns


def heatmap_from_array(data, convert_to_log=False, y_tick_labels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False, num_colors=7,
                       figsize=(6, 4), sort_row=None, annotate_sig=False,
                       rank_index=None,
                       linewidths=0.0, cluster_by_set=False, min_sig=0):
    """

    Parameters
    ----------
    data : magine.data.base.BaseData
    convert_to_log : bool
        Convert fold_change column to log2 scale
    y_tick_labels : list_like
    columns : str
        Name of columns of df for pivot
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
    figsize : tuple
        Size of figure, passed to matplotlib/seaborn
    sort_row : str
        Sort rows by ('index', 'mean', max')
    num_colors : int
        Number of colors for color bar
    annotate_sig : bool
        Add '*' annotation to plot for significant changed terms
    linewidths : float or None
        Add white line between plots
    cluster_by_set: bool
        Cluster by gene set column. Only works for enrichment_array
    min_sig : int
        Minimum number of significant 'index' across samples. Can be used to
        remove rows that are not significant across any sample.
    rank_index: bool
        Rank rows by index. Deprecated , plus use sort_row arg instead.
    Returns
    -------
    plt.Figure

    """
    if min_sig:
        d_copy = data.require_n_sig(columns=columns, index=index,
                                    n_sig=min_sig)
    else:
        d_copy = data.copy()
    array = d_copy.pivoter(convert_to_log, columns=columns, index=index,
                           fill_value=0.0, values=values, min_sig=min_sig)
    if not len(array):
        warnings.warn("Empty array after filtering.")
        return
    # default values to be overwritten below
    col_colors = None
    annotations = None
    col_color_map = None
    col_labels = None
    fmt = None
    linkage = None
    add_col_group = False

    if rank_index is not None:
        warnings.warn("rank_index is deprecated; use sort_row='index'",
                      DeprecationWarning)
        if rank_index:
            # assuming provided true, will sort by rank.
            array.sort_index(ascending=True, inplace=True)

    if sort_row is not None:
        if isinstance(sort_row, (list, np.ndarray)):
            array = array.reindex(sort_row)
        elif sort_row not in ('index', 'max', 'mean', 'min', 'sum'):
            raise ValueError("Can sort rows by 'index' name or 'max', 'min',"
                             "'mean' of values")

    # rank by index or cluster by term column
    if sort_row == 'index':
        array.sort_index(ascending=True, inplace=True)
    elif sort_row == 'mean':
        new_index = array.mean(axis=1).sort_values(ascending=False).index
        array = array.reindex(new_index)
    elif sort_row == 'max':
        new_index = array.max(axis=1).sort_values(ascending=False).index
        array = array.reindex(new_index)
    elif sort_row == 'min':
        new_index = array.max(axis=1).sort_values(ascending=False).index
        array = array.reindex(new_index)
    elif sort_row == 'sum':
        new_index = array.sum(axis=1).sort_values(ascending=True).index
        array = array.reindex(new_index)
    elif isinstance(sort_row, list):
        array = array.reindex(sort_row)
    if cluster_by_set and "genes" in d_copy.columns:
        # clustering will be based on jaccard index of terms
        dist_mat, names = d_copy.calc_dist(level='sample')
        linkage = sch.linkage(dist_mat, method='average')
        # Add row cluster flag in case user didn't set
        cluster_row = True
        array = array.reindex(names)

    # set coloring scheme for heatmap
    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("purple", as_cmap=True)
        center = None

    # Group together by columns if provided
    if isinstance(columns, (list, tuple)) and len(columns) == 2:
        add_col_group = True
        col_labels, col_colors, col_color_map = _set_col_colors(array)

    # check annotations exist
    if annotate_sig:
        annotate_sig, annotations, fmt = _get_sig_annotations(array, d_copy,
                                                              columns,
                                                              index, min_sig)
    cluster_args = dict(method='complete', metric='correlation')
    if cluster_row or cluster_col or add_col_group:
        fig = sns.clustermap(array, cmap=pal, center=center,
                             yticklabels=y_tick_labels,
                             col_colors=col_colors,
                             col_cluster=cluster_col,
                             row_cluster=cluster_row,
                             row_linkage=linkage,
                             figsize=figsize, linewidths=linewidths,
                             annot=annotations, fmt=fmt,
                             **cluster_args
                             )

        # We need to reorder the annotations if we cluster
        if annotate_sig:
            if cluster_row:
                annotations = annotations[fig.dendrogram_row.reordered_ind]
            if cluster_col:
                annotations = annotations[:, fig.dendrogram_col.reordered_ind]
            # Only need figure for dendrogram ordering, not actual plot.
            # Can probably do this without the plotting interface, but this
            # seems to do the job for now.

            plt.close()

            # make final figure
            fig = sns.clustermap(array, cmap=pal, center=center,
                                 yticklabels=y_tick_labels,
                                 col_colors=col_colors,
                                 col_cluster=cluster_col,
                                 row_cluster=cluster_row,
                                 row_linkage=linkage,
                                 figsize=figsize, linewidths=linewidths,
                                 annot=annotations, fmt=fmt,
                                 **cluster_args
                                 )

        # add labels to column colors
        if add_col_group:
            fig = _add_column_color_groups(d_copy, fig, col_color_map,
                                           col_labels, columns)

        # add clustered columns and rows, basically allows us to extract out
        # the clusters from the figure, if we wanted to do something with them
        # ie run enrichment analysis.
        if cluster_col:
            col_cltrs = sch.fcluster(fig.dendrogram_col.linkage, t=2,
                                     criterion='maxclust')
            col_cltrs = col_cltrs[fig.dendrogram_col.reordered_ind]
            col_clusters = dict()
            for i in sorted(set(col_cltrs)):
                cols = fig.data2d.columns.values[col_cltrs == i]
                col_clusters[i] = fig.data2d[cols]
            fig.col_clusters = col_clusters
        if cluster_row:
            row_cltrs = sch.fcluster(fig.dendrogram_row.linkage, t=2,
                                     criterion='maxclust')
            row_cltrs = row_cltrs[fig.dendrogram_row.reordered_ind]
            row_clusters = dict()
            for i in sorted(set(row_cltrs)):
                row_clusters[i] = fig.data2d.loc[row_cltrs == i].index.values
            fig.row_clusters = row_clusters
        fig.ax_heatmap.set_ylabel('')
        fig.ax_heatmap.set_xlabel('')
        # plt.subplots_adjust(right=0.7, top=1.5)
    else:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        sns.heatmap(array, ax=ax, yticklabels=y_tick_labels, cmap=pal,
                    center=center, annot=annotations, fmt=fmt,
                    linewidths=linewidths)
        ax.set_ylabel('')
        ax.set_xlabel('')
    return fig


def heatmap_by_terms(data, term_labels, term_sets, colors=None, min_sig=None,
                     convert_to_log=False, y_tick_labels='auto',
                     columns='sample_id', index='identifier',
                     values='fold_change', linewidths=0,
                     cluster_row=False, cluster_col=False, div_colors=False,
                     num_colors=21, figsize=None, annotate_sig=False,
                     **kwargs):
    """

    Parameters
    ----------
    data : pd.DataFrame
    term_labels : list_like
        List of labels for grouping
    term_sets : list_like
        List of list like that create the terms
    colors : list_like
        Colors for plotting, if not provided it will be created
    min_sig : int
        Number of sign
    convert_to_log : bool
    y_tick_labels : list_like
    columns : str
        Name of columns of df for pivotn
    index : str
        Name of index of df for pivot
    values : str
        Name of values of df for pivot
    cluster_col : bool
        Cluster the data using searborn.clustermap
    cluster_row : bool
        Cluster rows
    div_colors : bool
        Use divergent colors for plotting
    figsize : tuple
        Size of figure, passed to matplotlib/seaborn
    num_colors : int
        Number of colors for color bar
    annotate_sig : bool
        Add '*' annotation to plot for significant changed terms
    linewidths : float or None
        Add white line between plots
    min_sig : int
        Minimum number of significant 'index' across samples. Can be used to
        remove rows that are not significant across any sample.

    Returns
    -------
    plt.Figure

    """

    if len(term_labels) != len(term_sets):
        raise AssertionError("Number of term_labels must "
                             "equal number of term_sets")
    # default values to be overwritten below
    annotations = None
    fmt = None
    add_col_group = False
    tmp_d = data.copy()
    if index == 'label':
        id_to_label = dict()
        for i, j in tmp_d[['identifier', 'label']].values:
            if i not in id_to_label:
                id_to_label[i] = set()
            id_to_label[i].add(j)
        new_term_sets = [
            set(chain.from_iterable([id_to_label[j] for j in i
                                     if j in id_to_label]))
            for i in term_sets
        ]
    else:
        new_term_sets = term_sets
    all_items = set(chain.from_iterable(new_term_sets))
    tmp_d = tmp_d.loc[tmp_d[index].isin(all_items)]

    # pivot datatable
    array = tmp_d.pivoter(convert_to_log, columns=columns, index=index,
                          fill_value=0.0, values=values, min_sig=min_sig)
    if not len(array):
        warnings.warn("Empty array after filtering.")
        return
    if colors is None:
        colors = sns.color_palette("Paired", n_colors=len(term_labels))
    else:
        if len(colors) != len(new_term_sets):
            raise AssertionError("Number of colors must "
                                 "equal number of term_labels")
    vals = set(array.index.values)
    final_sorted = sorted(new_term_sets[0].intersection(vals))
    added = set(final_sorted)
    # create colors for each

    row_colors = [colors[0] for _ in added]
    to_remove = set()
    for term, color, cname in zip(new_term_sets[1:], colors[1:],
                                  term_labels[1:]):
        added_any = False
        for i in sorted(term.intersection(vals)):
            if i not in added:
                added_any = True
                row_colors.append(color)
                final_sorted.append(i)
                added.add(i)
        if not added_any:
            to_remove.add(cname)
    # only keep indexes that are in the provided sets
    array = array[array.index.isin(final_sorted)]

    # resort according to color
    array = array.reindex(final_sorted)
    if isinstance(columns, list) and len(columns) == 2:
        add_col_group = True
        col_labels, col_colors, col_color_map = _set_col_colors(array)
    else:
        col_labels, col_colors, col_color_map = None, None, None
    # set colors map for heatmap
    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("red", n_colors=len(new_term_sets),
                                as_cmap=True)
        center = None

    if annotate_sig:
        annotate_sig, annotations, fmt = _get_sig_annotations(array, data,
                                                              columns,
                                                              index, min_sig)
    cluster_args = dict(method='single', metric='correlation')
    fig = sns.clustermap(array,
                         yticklabels=y_tick_labels,
                         figsize=figsize, linewidths=linewidths,
                         row_colors=row_colors, col_colors=col_colors,
                         col_cluster=cluster_col, row_cluster=cluster_row,
                         cmap=pal, center=center,
                         annot=annotations, fmt=fmt,
                         **cluster_args
                         )
    if annotate_sig:
        if cluster_col:
            annotations = annotations[:, fig.dendrogram_col.reordered_ind]
        if cluster_row:
            annotations = annotations[fig.dendrogram_row.reordered_ind]
        plt.close()
        fig = sns.clustermap(array,
                             yticklabels=y_tick_labels,
                             figsize=figsize, linewidths=linewidths,
                             row_colors=row_colors, col_colors=col_colors,
                             col_cluster=cluster_col, row_cluster=cluster_row,
                             cmap=pal, center=center,
                             annot=annotations, fmt=fmt,
                             **cluster_args
                             )
    for color, label in zip(colors, term_labels):
        if label in to_remove:
            continue
        fig.ax_row_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    fig.ax_row_dendrogram.legend(loc=0, ncol=1)
    if add_col_group:
        fig = _add_column_color_groups(tmp_d, fig, col_color_map, col_labels,
                                       columns)
    fig.ax_heatmap.set_ylabel('')
    fig.ax_heatmap.set_xlabel('')
    plt.subplots_adjust(right=0.7)
    return fig


def cluster_distance_mat(dist_mat, names, figsize=(8, 8)):
    """ Creates heatmap from distance matrix.

    Parameters
    ----------
    dist_mat : np.array
        Distance matrix array.
    names : list_like
        Names of ticks for distance matrix
    figsize : tuple
        Size of figure, passed to matplotlib

    Returns
    -------

    """
    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=figsize)

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


def _set_col_colors(array):
    col_labels = array.columns.levels[0]
    labels = list(array.columns.levels[1])
    col_color_map = sns.color_palette("Dark2", len(col_labels))
    col_colors = [col_color_map[i] for i in array.columns.codes[0]]
    array.columns = [labels[i] for i in array.columns.codes[1]]
    return col_labels, col_colors, col_color_map


def _add_column_color_groups(data, fig, colors, color_labels, columns):
    for color, label in zip(colors, color_labels):
        fig.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    plt.setp(fig.ax_col_dendrogram.yaxis.get_majorticklabels(), rotation=0,
             fontsize=16)
    fig.ax_col_dendrogram.legend(loc="center", fontsize=12, ncol=2,
                                 bbox_to_anchor=(0.5, 1., 0.5, 0.5))
    v_line_list = []
    prev = 0
    for i in color_labels:
        n_samples = len(data[data[columns[0]] == i][columns[1]].unique())
        prev += n_samples
        v_line_list.append(prev)

    fig.fig.axes[2].vlines(v_line_list, *fig.fig.axes[2].get_ylim())
    fig.fig.axes[3].vlines(v_line_list, *fig.fig.axes[3].get_ylim())
    return fig


def _get_sig_annotations(arr, dat, columns, index, min_sig):
    # Have to rank by column for this to work
    if 'significant' in dat.columns:
        tmp2 = dat.pivoter(False, columns=columns, index=index,
                           values='significant', fill_value=0, min_sig=min_sig)
        tmp2 = tmp2.reindex(arr.index)
        tmp2[tmp2 > 0] = True
        tmp2 = tmp2.replace(0, '')
        tmp2 = tmp2.replace(False, '')
        tmp2 = tmp2.replace(True, '+')
        return True, tmp2.values, ''
    else:
        print("To annotate please add a significant column to data")
        return False, None, None
