import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import seaborn as sns


def heatmap_from_array(data, convert_to_log=False, y_tick_labels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False, num_colors=7,
                       fig_size=(6, 4), rank_index=False, annotate_sig=False,
                       linewidths=0.0, cluster_by_set=False):
    """

    Parameters
    ----------


    data : magine.data.experimental_data.ExpreimentalData
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
        Cluster the data using searborn.clustermap
    div_colors : bool
        Use divergent colors for plotting
    fig_size : tuple
        Size of figure, passed to matplotlib/seaborn
    rank_index : bool
        Order by index.
    num_colors : int
        Number of colors for color bar
    annotate_sig : bool
        Add '*' annotation to plot for significant changed terms
    linewidths : float or None
        Add white line between plots
    cluster_by_set: bool
        Cluster by gene set column. Only works for enrichment_array

    Returns
    -------

    """
    array = data.pivoter(convert_to_log, columns=columns, index=index,
                         fill_value=0.0, values=values)
    annotations = None
    fmt = None

    if rank_index:
        array.sort_index(ascending=True, inplace=True)

    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("purple", as_cmap=True)
        center = None

    linkage = None

    if cluster_by_set and "genes" in data.columns:
        # generate structure of graph based on jaccard index
        dist_mat, names = data.calc_dist(level='sample')
        linkage = sch.linkage(dist_mat, method='average')

        # Add row cluster flag in case user didn't set
        cluster_row = True
        array = array.reindex(names)

    if cluster_row or cluster_col:
        fig = sns.clustermap(array,
                             yticklabels=y_tick_labels,
                             row_linkage=linkage,
                             col_linkage=None,
                             row_cluster=cluster_row,
                             col_cluster=cluster_col,
                             linewidths=linewidths,
                             figsize=fig_size,
                             center=center,
                             cmap=pal
                             )
        if annotate_sig:
            annotate_sig, annotations, fmt = get_sig_annotations(array, data,
                                                                 columns,
                                                                 index)
            if not annotate_sig:
                return fig
            if cluster_by_set or cluster_row:
                annotations = annotations[fig.dendrogram_row.reordered_ind]
            if cluster_col:
                annotations = annotations[:, fig.dendrogram_col.reordered_ind]
            plt.close()
            fig = sns.clustermap(array,
                                 yticklabels=y_tick_labels,
                                 row_linkage=linkage,
                                 col_linkage=None,
                                 row_cluster=cluster_row,
                                 col_cluster=cluster_col,
                                 linewidths=linewidths,
                                 figsize=fig_size,
                                 center=center,
                                 cmap=pal,
                                 annot=annotations,
                                 fmt=fmt
                                 )

    else:
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(111)
        if annotate_sig:
            annotate_sig, annotations, fmt = get_sig_annotations(array, data,
                                                                 columns,
                                                                 index)

        sns.heatmap(array, ax=ax, yticklabels=y_tick_labels, cmap=pal,
                    center=center, annot=annotations, fmt=fmt,
                    linewidths=linewidths)
    return fig


def heatmap_by_category(data,
                        convert_to_log=False,
                        figsize=(6, 4),
                        columns=('category', 'sample_id'),
                        index='term_name', cluster_by_set=False,
                        values='combined_score', linewidths=0,
                        num_colors=11, div_colors=False,
                        annotate_sig=False, cluster_row=False):
    array = data.pivoter(convert_to_log, columns=columns,
                         index=index, fill_value=0.0,
                         values=values)
    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("purple", as_cmap=True)
        center = None

    color_labels = array.columns.levels[0]
    col_labels = list(array.columns.levels[1])
    colors = sns.color_palette("Dark2", len(color_labels))
    row_colors = [colors[i] for i in array.columns.labels[0]]
    array.columns = [col_labels[i] for i in array.columns.labels[1]]
    annotations = None
    fmt = None
    linkage = None
    if cluster_by_set and "genes" in data.columns:
        # generate structure of graph based on jaccard index
        dist_mat, names = data.calc_dist(level='sample')
        linkage = sch.linkage(dist_mat, method='average')

        # Add row cluster flag in case user didn't set
        cluster_row = True
        array = array.reindex(names)

    if cluster_row or cluster_by_set:

        fig = sns.clustermap(array, cmap=pal, center=center,
                             row_linkage=linkage,
                             yticklabels='auto',
                             col_colors=row_colors,
                             col_cluster=False, row_cluster=cluster_row,
                             figsize=figsize, linewidths=linewidths,
                             annot=annotations, fmt=fmt)

        if annotate_sig:
            annotate_sig, annotations, fmt = get_sig_annotations(array, data,
                                                                 columns,
                                                                 index)
            # makes sure it is possible
            if annotate_sig:
                annotations = annotations[fig.dendrogram_row.reordered_ind]
                plt.close()
                fig = sns.clustermap(array, cmap=pal, center=center,
                                     row_linkage=linkage,
                                     yticklabels='auto',
                                     col_colors=row_colors,
                                     col_cluster=False, row_cluster=cluster_row,
                                     figsize=figsize, linewidths=linewidths,
                                     annot=annotations, fmt=fmt)
    else:
        fig = sns.clustermap(array, cmap=pal, center=center,
                             row_linkage=linkage,
                             yticklabels='auto', col_colors=row_colors,
                             col_cluster=False, row_cluster=cluster_row,
                             figsize=figsize, linewidths=linewidths,
                             annot=annotations, fmt=fmt)

    for color, label in zip(colors, color_labels):
        fig.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    plt.setp(fig.ax_col_dendrogram.yaxis.get_majorticklabels(), rotation=0,
             fontsize=16)
    fig.ax_col_dendrogram.legend(loc="center", ncol=3, fontsize=12)
    v_line_list = []
    prev = 0
    for i in color_labels:
        n_samples = len(data[data['category'] == i]['sample_id'].unique())
        prev += n_samples
        v_line_list.append(prev)

    fig.fig.axes[2].vlines(v_line_list, *fig.fig.axes[2].get_ylim())
    fig.fig.axes[3].vlines(v_line_list, *fig.fig.axes[3].get_ylim())
    return fig


def heatmap_by_terms(data, terms, color_labels, colors=None,
                     convert_to_log=False, y_tick_labels='auto',
                     columns='sample_id', index='term_name', cluster_col=False,
                     values='combined_score', div_colors=False, num_colors=7,
                     fig_size=(6, 4), annotate_sig=False):
    # pivot datatable
    array = data.pivoter(convert_to_log, columns=columns, index=index,
                         fill_value=0.0, values=values)

    if colors is None:
        colors = sns.color_palette("Dark2", len(terms))
    vals = set(array.index.values)
    final_sorted = sorted(terms[0].intersection(vals))
    added = set(final_sorted)
    # create colors for each
    row_colors = [colors[0] for _ in added]

    for term, color in zip(terms, colors):
        for i in sorted(term.intersection(vals)):
            if i not in added:
                row_colors.append(color)
                final_sorted.append(i)
                added.add(i)

    # only keep indexes that are in the provided sets
    array = array[array.index.isin(final_sorted)]

    # resort according to color
    array = array.reindex(final_sorted)

    labels = None
    fmt = None

    if annotate_sig:
        annotate_sig, annotations, fmt = get_sig_annotations(array, data,
                                                             columns,
                                                             index)

    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("red", as_cmap=True)
        center = None

    fig = sns.clustermap(array,
                         yticklabels=y_tick_labels,
                         row_colors=row_colors,
                         figsize=fig_size,
                         col_cluster=cluster_col,
                         row_cluster=False,
                         center=center,
                         cmap=pal,
                         annot=labels, fmt=fmt
                         )
    for color, label in zip(colors, color_labels):
        fig.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    fig.ax_col_dendrogram.legend(loc="center", ncol=6)

    return fig


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


def get_sig_annotations(arr, dat, columns, index):
    # Have to rank by column for this to work
    if 'significant' in dat.columns:
        tmp2 = dat.pivoter(False, columns=columns, index=index,
                           values='significant', fill_value=0,
                           min_sig=False)
        tmp2 = tmp2.reindex(arr.index)
        tmp2[tmp2 > 0] = True
        tmp2 = tmp2.replace(False, '')
        tmp2 = tmp2.replace(True, '*')
        return True, tmp2.values, ''
    else:
        print("To annotate please add a significant_flag column to data")
        return False, None, None
