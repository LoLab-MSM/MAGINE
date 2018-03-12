import matplotlib.pyplot as plt
import seaborn as sns

from magine.data.tools import pivot_table


def heatmap_from_array(data, convert_to_log=False, yticklabels='auto',
                       cluster_row=False, cluster_col=False,
                       columns='sample_id', index='term_name',
                       values='combined_score', div_colors=False, num_colors=7,
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
    array = pivot_table(data, convert_to_log, columns=columns, index=index,
                        values=values)
    if div_colors:
        pal = sns.color_palette("coolwarm", num_colors)
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
