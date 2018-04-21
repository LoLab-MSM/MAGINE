import matplotlib.pyplot as plt
import seaborn as sns

from magine.data.tools import pivot_table


def heatmap_from_array(data, convert_to_log=False, y_tick_labels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False, num_colors=7,
                       fig_size=(6, 4), rank_index=False, annotate_sig=False):
    """

    Parameters
    ----------


    data : pandas.DataFrame
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
    Returns
    -------

    """
    array = pivot_table(data, convert_to_log, columns=columns, index=index,
                        fill_value=0.0, values=values)
    labels = None
    fmt = None

    if rank_index:
        array.sort_index(ascending=True, inplace=True)

    if annotate_sig:
        # Have to rank by column for this to work
        if 'significant_flag' in data.columns:
            tmp2 = pivot_table(data, False, columns=columns,
                               index=index, fill_value=False,
                               values='significant_flag', min_sig=False)

            tmp2 = tmp2.reindex(array.index)
            tmp2[tmp2 > 0] = True
            tmp2 = tmp2.replace(False, '')
            tmp2 = tmp2.replace(True, '*')

            labels = tmp2.as_matrix()
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
                    center=center, annot=labels, fmt=fmt)

    return fig


def heatmap_by_terms(data, terms, colors, color_labels, convert_to_log=False,
                     y_tick_labels='auto', cluster_col=False,
                     columns='sample_id', index='term_name',
                     values='combined_score', div_colors=False, num_colors=7,
                     fig_size=(6, 4), annotate_sig=False):
    array = pivot_table(data, convert_to_log, columns=columns, index=index,
                        fill_value=0.0, values=values)

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

    array = array.reindex(final_sorted)

    labels = None
    fmt = None

    if annotate_sig:
        # Have to rank by column for this to work
        if 'significant_flag' in data.columns:
            tmp2 = pivot_table(data.copy(), False, columns=columns,
                               index=index, fill_value=False,
                               values='significant_flag', min_sig=False)

            tmp2 = tmp2.reindex(array.index)
            tmp2[tmp2 > 0] = True
            tmp2 = tmp2.replace(False, '')
            tmp2 = tmp2.replace(True, '*')

            labels = tmp2.as_matrix()
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
