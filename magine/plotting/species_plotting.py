import os
import re
import time
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathos.multiprocessing as mp
import plotly.graph_objs as plotly_graph
from plotly.offline import plot

import magine.html_templates.html_tools as ht
from magine.data.formatter import pivot_table_for_export, log2_normalize_df

gene_index = 'gene'
protein = 'protein'
metabolites = 'metabolites'
meta_index = 'compound'
species_type = 'species_type'
sample_id = 'sample_id'
fold_change = 'treated_control_fold_change'
flag = 'significant_flag'
cm = plt.get_cmap('jet')


def write_table_to_html(data, save_name='index', out_dir=None,
                        run_parallel=False, exp_data=None):
    """
    Creates a table of all the plots of all genes for each GO term.

    Uses last calculated enrichment array.

    Parameters
    ----------
    data : pandas.DataFrame
    save_name : str
        name of html output file
    out_dir : str, optional
        output path for all plots
    run_parallel : bool
        Create plots in parallel
    exp_data : magine.data.datatypes.ExperimentalData
    Returns
    -------

    """

    list_of_terms = list(data['term_name'].unique())
    fig_dict, to_remove = plot_genes_by_ont(data=data,
                                            list_of_terms=list_of_terms,
                                            save_name=save_name,
                                            out_dir=out_dir,
                                            exp_data=exp_data,
                                            run_parallel=run_parallel
                                            )

    for i in fig_dict:
        data.loc[data['term_name'] == i, 'term_name'] = fig_dict[i]

    data = data[~data['term_name'].isin(to_remove)]

    ht.write_single_table(data, 'MAGINE GO analysis', save_name)
    html_out = save_name + '_filter'
    ht.write_filter_table(data, html_out)


def plot_genes_by_ont(data, list_of_terms, save_name, out_dir=None,
                      exp_data=None, run_parallel=False, plot_type='plotly'):
    """ Creates a figure for each GO term in data

    Data should be a result of running calculate_enrichment.
    This function creates a plot of all proteins per term if a term is
    significant and the number of the reference set is larger than 5 and
    the total number of species measured is less than 100.


    Parameters
    ----------
    data : pandas.DataFrame
        previously ran enrichment analysis
    list_of_terms : list_list

    save_name : str
        name to save file
    out_dir : str
        output path for file
    exp_data : magine.ExperimentalData
        data to plot
    run_parallel : bool
        To run in parallel using pathos.multiprocessing
    plot_type : str
        plotly or matplotlib

    Returns
    -------
    out_array : dict
        dict where keys are pointers to figure locations
    """

    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if not os.path.exists(os.path.join(out_dir, 'Figures')):
            os.mkdir(os.path.join(out_dir, 'Figures'))
    data = data.copy()
    figure_locations = {}
    plots_to_create = []
    to_remove = set()
    assert plot_type == ('plotly' or 'matplotlib')
    # filter data by significance and number of references
    if len(list_of_terms) == 0:
        print("No significant GO terms!!!")
        return figure_locations, to_remove
    # here we are going to iterate through all sig GO terms and create
    # a list of plots to create. For the HTML side, we need to point to
    # a location

    # create plot of genes over time
    for n, i in enumerate(list_of_terms):
        # want to plot all species over time
        index = data['term_name'] == i

        name = data[index]['term_name'].unique()

        if len(name) > 0:
            name = name[0]

        gene_set = set()
        genes = data[index]['genes']
        for g in genes:
            if isinstance(g, list):
                each = g
            else:
                each = g.split(',')

            gene_set.update(set(each))

        if plot_type == 'matplotlib':
            # too many genes isn't helpful on plots, so skip them
            if len(gene_set) > 100:
                figure_locations[i] = '<a>{0}</a>'.format(name)
                continue
        local_save_name = 'Figures/{0}_{1}'.format(n, save_name)
        if out_dir is not None:
            local_save_name = '{0}/{1}'.format(out_dir, local_save_name)

        local_save_name = local_save_name.replace(':', '')
        out_point = '<a href="{0}.html">{1}</a>'.format(local_save_name, name)
        figure_locations[i] = out_point

        title = "{0} : {1}".format(str(i), name)
        local_df = exp_data.data[
            exp_data.data[gene_index].isin(list(gene_set))].copy()
        p_input = (local_df, list(gene_set), 'gene', local_save_name, '.',
                   title, plot_type)

        plots_to_create.append(p_input)

    print("Starting to create plots for each ontology term")
    _make_plots(plots_to_create, plot_list_of_genes, run_parallel)

    return figure_locations, to_remove


def write_table_to_html_with_figures(data, exp_data, save_name='index',
                                     out_dir=None, run_parallel=True):
    # create plots of everything
    if isinstance(data, str):
        data = pd.read_csv(data)

    fig_dict, to_remove = plot_genes_by_ont(
        data, save_name, out_dir, exp_data, run_parallel=run_parallel
    )

    for i in fig_dict:
        data.loc[data['GO_id'] == i, 'GO_name'] = fig_dict[i]

    data = data[~data['GO_id'].isin(to_remove)]

    tmp = pivot_table_for_export(data)

    html_out = save_name
    if out_dir is not None:
        html_out = os.path.join(out_dir, html_out)
    print("Saving to : {}".format(html_out))

    html_out = save_name + '_filter'
    if out_dir is not None:
        html_out = os.path.join(out_dir, html_out)

    ht.write_filter_table(tmp, html_out)


def plot_dataframe(exp_data, html_filename, out_dir='proteins',
                   plot_type='plotly', type_of_species='protein',
                   run_parallel=False):
    """
    Creates a plot of all proteins

    Parameters
    ----------
    exp_data : pandas.DataFrame
    html_filename : str
    out_dir: str, path
        Directory that will contain all proteins
    plot_type : str
        plotly or matplotlib output
    type_of_species : str
        proteins or metabolites
    run_parallel : bool
        create plots in parallel
    Returns
    -------

    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if type_of_species == 'protein':
        idx_key = 'gene'
        plot_func = plot_list_of_genes
    elif type_of_species == 'metabolites':
        idx_key = 'compound'
        plot_func = plot_list_of_metabolites
    else:
        print('type_of_species can only be "protein" or "metabolites"')
        return
    local_data = exp_data[exp_data['species_type'] == type_of_species].copy()

    assert idx_key in local_data.dtypes,\
        '{} is not in your species_type column.\n'.format(idx_key)

    species_to_plot = local_data[idx_key].unique()

    print("Plotting {} {}".format(len(species_to_plot), type_of_species))
    figure_locations = {}
    list_of_plots = []

    suffix = 'html' if plot_type == 'plotly' else 'pdf'

    for i in species_to_plot:
        save_name = re.sub('[/_.]', '', i)
        list_of_plots.append(
            (
            local_data, [i], type_of_species, save_name, out_dir, i, plot_type)
        )
        figure_locations[i] = '<a href="{0}/{1}.{2}">{1}</a>'.format(out_dir,
                                                                     save_name,
                                                                     suffix)

    _make_plots(list_of_plots, plot_func, run_parallel)

    # Place a link to the species for each key
    for i in figure_locations:
        local_data.loc[exp_data[idx_key] == i, idx_key] = figure_locations[i]

    if type_of_species == 'protein':
        cols = ['gene', 'treated_control_fold_change', 'protein',
                'p_value_group_1_and_group_2', 'time', 'data_type',
                'significant_flag',
                ]
    elif type_of_species == 'metabolites':

        cols = ['compound',
                'treated_control_fold_change',
                'p_value_group_1_and_group_2', 'significant_flag',
                'data_type', 'time',  # 'time_points',
                ]
        if 'compound_id' in local_data.columns:
            cols.insert(2, 'compound_id')

    local_data = local_data[cols]
    ht.write_filter_table(local_data, html_filename)


def _make_plots(plots_to_make, plot_func, parallel=False):
    if parallel:
        st2 = time.time()
        pool = mp.Pool()
        pool.map_async(plot_func, plots_to_make)
        pool.close()
        pool.join()
        end2 = time.time()
        print("parallel time = {}".format(end2 - st2))
        print("Done creating plots for each GO term")

    else:
        st1 = time.time()
        map(plot_func, plots_to_make)
        end1 = time.time()
        print("sequential time = {}".format(end1 - st1))


def plot_list_of_metabolites(dataframe, list_of_metab=None,
                             species_type='gene',
                             save_name='test', out_dir=None, title=None,
                             plot_type='plotly', image_format='pdf'):
    """

    Parameters
    ----------
    dataframe: pandas.DataFrame
        magine formatted dataframe
    list_of_metab: list
        List of genes to be plotter
    save_name: str
        Filename to be saved as
    out_dir: str
        Path for output to be saved
    title: str
        Title of plot, useful when list of genes corresponds to a GO term
    plot_type : str
        Use plotly to generate html output or matplotlib to generate pdf
    image_format : str
        pdf or png, only used if plot_type="matplotlib"

    Returns
    -------

    """
    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    if list_of_metab is None:
        dataframe, list_of_metab, species_type, save_name, out_dir, title, plot_type = dataframe

    if 'sample_id' not in dataframe.dtypes:
        dataframe['sample_id'] = dataframe['time']

    # gather x axis points
    x_points = sorted(dataframe[sample_id].unique())

    if isinstance(x_points[0], np.float):
        x_point_dict = {i: x_points[n] for n, i in enumerate(x_points)}
    else:
        x_point_dict = {i: n for n, i in enumerate(x_points)}
    if species_type == 'metabolites':
        idx = meta_index
    else:
        idx = gene_index

    local_df = dataframe[dataframe[idx].isin(list_of_metab)].copy()
    local_df = log2_normalize_df(local_df, fold_change=fold_change)

    n_plots = len(local_df[idx].unique())

    color_list = [cm(1. * i / n_plots) for i in range(n_plots)]
    colors = enumerate(color_list)

    if plot_type == 'matplotlib':
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_prop_cycle(plt.cycler('color', color_list))

    plotly_list = []
    names_list = []
    total_counter = 0

    group = local_df.groupby(idx)
    for name, m in group:

        index_counter = 0
        x = np.array(m[sample_id])
        if len(x) < 1:
            continue
        y = np.array(m['log2fc'])
        sig_flag = np.array(m[flag])
        index = np.argsort(x)
        x = x[index]
        y = y[index]
        s_flag = sig_flag[index]

        # x values with scaled values (only changes things if non-float
        # values are used for sample_id
        x_index = np.array([x_point_dict[ind] for ind in x])

        index_counter += 1
        total_counter += 1
        # create matplotlib plot
        if plot_type == 'matplotlib':
            label = "\n".join(wrap(name, 40))
            p = ax.plot(x_index, y, '.-', label=label)
            if len(s_flag) != 0:
                color = p[0].get_color()
                ax.plot(x_index[s_flag], y[s_flag], '^', color=color)

        # create plotly plot
        elif plot_type == 'plotly':
            c = next(colors)[1]
            plotly_list.append(_create_ploty_graph(x_index, y, name, name, c))
            if len(s_flag) != 0:
                index_counter += 1
                total_counter += 1
                plotly_list.append(_create_ploty_graph(x_index[s_flag],
                                                       y[s_flag], name, name,
                                                       c, marker='x-open-dot'))
        names_list.append([name, index_counter])

    if plot_type == 'matplotlib':
        _save_matplotlib_output(ax, save_name, out_dir, image_format,
                                x_point_dict, x_points)

    elif plot_type == 'plotly':

        _save_ploty_output(out_dir, save_name, total_counter, n_plots,
                           names_list, x_point_dict, title, x_points,
                           plotly_list)


def plot_list_of_genes(df, genes=None, save_name='test',
                       out_dir=None, title=None, plot_type='plotly',
                       image_format='pdf'):
    """

    Parameters
    ----------
    df: pandas.DataFrame
        magine formatted dataframe
    genes: list
        List of genes to be plotter
    save_name: str
        Filename to be saved as
    out_dir: str
        Path for output to be saved
    title: str
        Title of plot, useful when list of genes corresponds to a GO term
    plot_type : str
        Use plotly to generate html output or matplotlib to generate pdf
    image_format : str
        pdf or png, only used if plot_type="matplotlib"

    Returns
    -------

    """

    if genes is None:
        df, genes, species_type, save_name, out_dir, title, plot_type = df

    ldf = df.copy(deep=True)
    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    if 'sample_id' not in ldf.dtypes:
        ldf['sample_id'] = ldf['time']

    # gather x axis points
    x_points = sorted(ldf[sample_id].unique())
    if len(x_points) == 0:
        return
    if isinstance(x_points[0], np.float):
        x_point_dict = {i: x_points[n] for n, i
                        in enumerate(x_points)}
    else:
        x_point_dict = {i: n for n, i
                        in enumerate(x_points)}

    ldf = ldf[ldf[gene_index].isin(genes)].copy(deep=True)
    ldf = log2_normalize_df(ldf, fold_change=fold_change)

    n_plots = len(ldf[gene_index].unique())
    num_colors = len(ldf[protein].unique())
    color_list = [cm(1. * i / num_colors) for i in range(num_colors)]

    if plot_type == 'matplotlib':
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_prop_cycle(plt.cycler('color', color_list))

    colors = enumerate(color_list)

    plotly_list = []
    names_list = []
    total_counter = 0
    group = ldf.groupby(gene_index)
    for i, j in group:
        name = i
        group2 = j.groupby(protein)
        index_counter = 0
        for n, m in group2:

            x = np.array(m[sample_id])
            if len(x) < 1:
                continue
            y = np.array(m['log2fc'])
            sig_flag = np.array(m[flag])
            index = np.argsort(x)
            x = x[index]
            y = y[index]
            s_flag = sig_flag[index]

            # x values with scaled values (only changes things if non-float
            # values are used for sample_id
            x_index = np.array([x_point_dict[ind] for ind in x])

            index_counter += 1
            total_counter += 1

            # create matplotlib plot
            if plot_type == 'matplotlib':
                label = "\n".join(wrap(n, 40))
                p = ax.plot(x_index, y, '.-', label=label)
                if len(s_flag) != 0:
                    color = p[0].get_color()
                    ax.plot(x_index[s_flag], y[s_flag], '^', color=color)

            # create plotly plot
            elif plot_type == 'plotly':
                c = next(colors)[1]
                plotly_list.append(_create_ploty_graph(x_index, y, n, n, c))
                if len(s_flag) != 0:
                    index_counter += 1
                    total_counter += 1
                    plotly_list.append(_create_ploty_graph(x_index[s_flag],
                                                           y[s_flag], n, n,
                                                           c,
                                                           marker='x-open-dot'))
        names_list.append([name, index_counter])

    if plot_type == 'matplotlib':
        _save_matplotlib_output(ax, save_name, out_dir, image_format,
                                x_point_dict, x_points)

    elif plot_type == 'plotly':
        _save_ploty_output(out_dir, save_name, total_counter, n_plots,
                           names_list, x_point_dict, title, x_points,
                           plotly_list)


def _save_matplotlib_output(ax, save_name, out_dir, image_format, x_point_dict,
                            x_points, ):
    ax.set_xlim(min(x_point_dict.values()) - 2, max(x_point_dict.values()) + 2)
    ax.set_xticks(sorted(x_point_dict.values()))
    ax.set_xticklabels(x_points, rotation=90)
    plt.ylabel('log$_2$ Fold Change')

    plt.axhline(y=np.log2(1.5), linestyle='--')
    plt.axhline(y=-np.log2(1.5), linestyle='--')

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='best', ncol=3,
                    bbox_to_anchor=(1.01, 1.0))

    tmp_savename = "{}.{}".format(save_name, image_format)
    if out_dir is not None:
        tmp_savename = os.path.join(out_dir, tmp_savename)

    plt.savefig(tmp_savename, bbox_extra_artists=(lgd,),
                bbox_inches='tight')


def _save_ploty_output(out_dir, save_name, total_counter, n_plots, names_list,
                       x_point_dict, title, x_points, plotly_list):
    true_list = [True] * total_counter
    scroll_list = [dict(args=['visible', true_list],
                        label='All',
                        method='restyle')]

    prev = 0
    # making all false except group defined by protein name
    for i in range(n_plots):
        t_row = [False] * total_counter
        for j in range(prev, prev + names_list[i][1]):
            t_row[j] = True
        prev += names_list[i][1]
        scroll = dict(args=['visible', t_row],
                      label=names_list[i][0], method='restyle')
        scroll_list.append(scroll)

    update_menu = list([dict(x=-0.05,
                             y=1,
                             yanchor='top',
                             buttons=scroll_list, )])
    ticks = np.sort(list(x_point_dict.values()))
    min_tick = np.min(ticks)
    max_tick = np.max(ticks)
    layout = plotly_graph.Layout(
        title=title,
        showlegend=True,
        xaxis=dict(title='Sample index',
                   range=[min_tick, max_tick],
                   showticklabels=True,
                   ticktext=x_points,
                   tickmode='array',
                   tickvals=ticks,
                   ),
        yaxis=dict(title='log2fc'),
        hovermode="closest",
        updatemenus=update_menu
    )

    fig = plotly_graph.Figure(data=plotly_list, layout=layout)
    tmp_savename = "{}.html".format(save_name)
    if out_dir is not None:
        tmp_savename = os.path.join(out_dir, tmp_savename)

    x = plot(fig, filename=tmp_savename, auto_open=False,
             include_plotlyjs=False, output_type='div')
    ht.format_ploty(x, tmp_savename)


def _create_ploty_graph(x, y, label, enum, color, marker='circle'):
    """
    Creates a single scatter plot
    Parameters
    ----------
    x : list_like
    y : list_like
    label : str
    enum : int
    color : str
    marker : str

    Returns
    -------

    """
    l_color = 'rgba({},{},{},1.)'.format(color[0], color[1], color[2])
    if marker != 'circle':
        mode = 'markers'
        show = False
        size = 12
    else:
        mode = 'lines+markers'
        show = True
        size = 8
    legend = 'group_{}'.format(enum)

    g = plotly_graph.Scatter(
            x=x,
            y=y,
            hoveron='text',
            name=label,
            visible=True,
            mode=mode,
        legendgroup=legend,
            showlegend=show,
            line=dict(color=l_color),
            marker=dict(symbol=marker,
                        size=size,
                        color=l_color),
    )
    return g
