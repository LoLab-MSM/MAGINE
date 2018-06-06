import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import seaborn as sns


def heatmap_from_array(data, convert_to_log=False, y_tick_labels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False, num_colors=7,
                       fig_size=(6, 4), rank_index=False, annotate_sig=False,
                       linewidths=0.0):
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
    Returns
    -------

    """
    array = data.pivoter(convert_to_log, columns=columns, index=index,
                         fill_value=0.0, values=values)
    labels = None
    fmt = None

    if rank_index:
        array.sort_index(ascending=True, inplace=True)

    if annotate_sig:
        # Have to rank by column for this to work
        if 'significant_flag' in data.columns:
            tmp2 = data.pivoter(False, columns=columns, index=index,
                                fill_value=False, values='significant_flag',
                                min_sig=False)

            tmp2 = tmp2.reindex(array.index)
            tmp2[tmp2 > 0] = True
            tmp2 = tmp2.replace(False, '')
            tmp2 = tmp2.replace(True, '*')

            labels = tmp2.values
            fmt = ''
        else:
            annotate_sig = False
            print("To annotate please add a significant_flag column to data")

    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
        center = 0
    else:
        pal = sns.light_palette("purple", as_cmap=True)
        center = None

    if cluster_col or cluster_row:
        fig = sns.clustermap(array, yticklabels=y_tick_labels,
                             figsize=fig_size, col_cluster=cluster_col,
                             row_cluster=cluster_row, center=center, cmap=pal)

        if cluster_row and annotate_sig:
            labels = labels[fig.dendrogram_row.reordered_ind]
        if cluster_col and annotate_sig:
            labels = labels[:, fig.dendrogram_col.reordered_ind]
        if (cluster_col or cluster_row) and annotate_sig:
            plt.close()
            fig = sns.clustermap(array, yticklabels=y_tick_labels,
                                 figsize=fig_size, col_cluster=cluster_col,
                                 row_cluster=cluster_row, center=center,
                                 cmap=pal, annot=labels, fmt=fmt)

    else:
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(111)
        sns.heatmap(array, ax=ax, yticklabels=y_tick_labels, cmap=pal,
                    center=center, annot=labels, fmt=fmt,
                    linewidths=linewidths)

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
        # Have to rank by column for this to work
        if 'significant_flag' in data.columns:
            tmp2 = data.pivoter(False, columns=columns, index=index,
                                fill_value=False, min_sig=False,
                                values='significant_flag')

            tmp2 = tmp2.reindex(array.index)
            tmp2[tmp2 > 0] = True
            tmp2 = tmp2.replace(False, '')
            tmp2 = tmp2.replace(True, '*')

            labels = tmp2.values
            fmt = ''
        else:
            print("To annotate please add a significant_flag column to data")

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
