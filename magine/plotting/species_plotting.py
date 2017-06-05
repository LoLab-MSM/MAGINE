import os
import time
from ast import literal_eval
from textwrap import wrap

import matplotlib
import numpy as np
import pandas as pd
import pathos.multiprocessing as mp
import plotly
import plotly.graph_objs as plotly_graph
from plotly.offline import plot

from magine.data.formatter import pivot_tables_for_export

matplotlib.use('Agg')
import matplotlib.pyplot as plt


plotly.plotly.sign_in(username='james.ch.pino',
                      api_key='BnUcJSpmPcMKZg0yEFaL')

gene = 'gene'
protein = 'protein'
metabolites = 'metabolites'
species_type = 'species_type'
sample_id = 'sample_id'
fold_change = 'treated_control_fold_change'
flag = 'significant_flag'


def create_gene_plots_per_go(data, save_name, out_dir, exp_data,
                             run_parallel=False, plot_type='plotly'):
    """ Creates a figure for each GO term in data

    Data should be a result of running calculate_enrichment.
    This function creates a plot of all proteins per term if a term is
    significant and the number of the reference set is larger than 5 and
    the total number of species measured is less than 100.


    Parameters
    ----------
    data : pandas.DataFrame
        previously ran enrichment analysis
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

    if isinstance(data, str):
        data = pd.read_csv(data)

    assert plot_type == ('plotly' or 'matplotlib')
    # get list of all terms
    list_of_go_terms = data['GO_id'].unique()

    # filter data by significance and number of references
    list_of_sig_go = data[(data['ref'] >= 5)
                          &
                          (data['ref'] <= 2000)
                          &
                          (data['pvalue'] < 0.05)]['GO_id'].unique()

    # here we are going to iterate through all sig GO terms and create
    # a list of plots to create. For the HTML side, we need to point to
    # a location
    figure_locations = {}
    plots_to_create = []
    to_remove = set()
    # create plot of genes over time
    for n, i in enumerate(list_of_go_terms):

        # want to plot all species over time
        index = data['GO_id'] == i

        name = data[index]['GO_name'].unique()

        if len(name) > 0:
            name = name[0]

        # want to only plot significant species
        if i not in list_of_sig_go:
            to_remove.add(i)
            figure_locations[i] = '<a>{0}</a>'.format(name)
            continue

        gene_set = set()
        genes = data[index]['genes']
        for g in genes:
            if isinstance(g, list):
                each = g
            else:
                each = literal_eval(g)

            for j in each:
                gene_set.add(j)
        if plot_type == 'matplotlib':
            # too many genes isn't helpful on plots, so skip them
            if len(gene_set) > 100:
                figure_locations[i] = '<a>{0}</a>'.format(name)
                continue
        out_point = '<a href="Figures/go_{0}_{1}.html">{2}</a>'
        out_point = out_point.format(i, save_name, name).replace(':', '')

        figure_locations[i] = out_point

        local_save_name = '{0}/Figures/go_{1}_{2}'.format(out_dir, i,
                                                          save_name)
        local_save_name = local_save_name.replace(':', '')
        title = "{0} : {1}".format(str(i), name)
        local_df = exp_data.data[exp_data.data[gene].isin(list(gene_set))]
        p_input = (local_df, list(gene_set), local_save_name, '.', title,
                   plot_type)
        plots_to_create.append(p_input)

    # return figure_locations, to_remove

    print("Starting to create plots for each GO term")
    # just keeping this code just in case using pathos is a bad idea
    # ultimately, using matplotlib is slow.

    if run_parallel:
        st2 = time.time()
        pool = mp.Pool()
        pool.map_async(plot_list_of_genes2, plots_to_create)
        # pool.map(plot_list_of_genes2, plots_to_create)
        pool.close()
        pool.join()
        end2 = time.time()
        print("parallel time = {}".format(end2 - st2))
        print("Done creating plots for each GO term")

    else:
        st1 = time.time()
        for i in plots_to_create:
            plot_list_of_genes2(i)
        end1 = time.time()
        print("sequential time = {}".format(end1 - st1))

    return figure_locations, to_remove


def plot_dataframe(exp_data, html_filename, out_dir='proteins',
                   plot_type='plotly'):
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

    Returns
    -------

    """
    if os.path.exists(out_dir):
        pass
    else:
        os.mkdir(out_dir)

    genes = exp_data[exp_data['species_type'] == 'protein'].copy()
    genes_to_plot = genes['gene'].unique()
    print("Plotting {} genes".format(len(genes_to_plot)))
    figure_locations = {}
    for i in genes_to_plot:
        plot_list_of_genes2(genes, list_of_genes=[i], save_name=i,
                            out_dir=out_dir, title=i, plot_type=plot_type)
        if plot_type == 'plotly':
            out_point = '<a href="{0}/{1}.html">{1}</a>'.format(out_dir, i)
            figure_locations[i] = out_point
        else:
            out_point = '<a href="{0}/{1}.pdf">{1}</a>'.format(out_dir, i)
            figure_locations[i] = out_point

    print(figure_locations)
    print(exp_data.head(20))
    for i in figure_locations:
        exp_data.loc[exp_data['gene'] == i, 'gene'] = figure_locations[i]
    print(exp_data.head(20))
    genes_out, meta_out = pivot_tables_for_export(exp_data)
    from magine.html_templates.html_tools import write_filter_table
    write_filter_table(genes_out, html_filename, 'genes')

    # quit()
    # proteins = pd.DataFrame(figure_locations, columns=['Genes'])
    # print(proteins.head(10))
    # quit()


def plot_list_of_genes2(dataframe, list_of_genes=None, save_name='test',
                        out_dir=None, title=None, plot_type='plotly',
                        image_format='pdf'):
    """

    Parameters
    ----------
    dataframe: pandas.DataFrame
        magine formatted dataframe
    list_of_genes: list
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
    from magine.html_templates.html_tools import format_ploty
    if out_dir is None:
        pass
    elif not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if list_of_genes is None:
        dataframe, list_of_genes, save_name, out_dir, title, plot_type = dataframe

    if 'sample_id' not in dataframe.dtypes:
        dataframe['sample_id'] = dataframe['time']
    x_points = sorted(dataframe[sample_id].unique())

    if isinstance(x_points[0], np.float):
        x_point_dict = {i: x_points[n] for n, i
                        in enumerate(x_points)}
    else:
        x_point_dict = {i: 2 ** (n + 1) for n, i
                        in enumerate(x_points)}

    local_df = dataframe[dataframe[gene].isin(list_of_genes)].copy()
    crit_1 = local_df[fold_change] > 0
    crit_2 = local_df[fold_change] < 0
    local_df.loc[crit_1, 'log2fc'] = np.log2(
        local_df[crit_1][fold_change].astype(np.float64))
    local_df.loc[crit_2, 'log2fc'] = -np.log2(
        -local_df[crit_2][fold_change].astype(np.float64))
    n_genes = len(local_df[gene].unique())
    group = local_df.groupby(gene)

    cm = plt.get_cmap('jet')
    num_colors = len(local_df[protein].unique())

    color_list = [cm(1. * i / num_colors) for i in range(num_colors)]

    if plot_type == 'matplotlib':
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_prop_cycle(plt.cycler('color', color_list))

    colors = enumerate(color_list)

    plotly_list = []
    names_list = []
    total_counter = 0
    for i, j in group:
        name = i

        group2 = j.groupby(protein)
        index_counter = 0
        for n, m in group2:

            x_index = []
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
            for ind in x:
                x_index.append(x_point_dict[ind])
            x_index = np.array(x_index)

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
        plt.xlim(min(x_point_dict.values()) - 2,
                 max(x_point_dict.values()) + 2)
        ax.set_xticks(sorted(x_point_dict.values()))
        ax.set_xticklabels(x_points)
        plt.ylabel('log$_2$ Fold Change')
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.axhline(y=np.log2(1.5), linestyle='--')
        plt.axhline(y=-np.log2(1.5), linestyle='--')

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='best',
                        bbox_to_anchor=(1.01, 1.0))
        if out_dir is None:
            tmp_savename = "{}.{}".format(save_name, image_format)
        else:
            tmp_savename = os.path.join(out_dir, "{}.{}".format(save_name,
                                                                image_format))
        # print("Saving {}".format(tmp_savename))
        plt.savefig(tmp_savename, bbox_extra_artists=(lgd,),
                    bbox_inches='tight')

    elif plot_type == 'plotly':
        true_list = [True] * total_counter
        scroll_list = [dict(args=['visible', true_list],
                            label='All',
                            method='restyle')]

        prev = 0
        # making all false except group defined by protein name
        for i in range(n_genes):
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

        layout = plotly_graph.Layout(
                title=title,
                showlegend=True,
                xaxis=dict(title='Sample index',
                           range=[min(x_point_dict.values()),
                                  max(x_point_dict.values())],
                           showticklabels=True,
                           ticktext=x_points,
                           tickmode='array',
                           tickvals=np.sort(x_point_dict.values()),
                           ),
                yaxis=dict(title='log2fc'),
                hovermode="closest",
                updatemenus=update_menu
        )

        fig = plotly_graph.Figure(data=plotly_list, layout=layout)
        if out_dir is None:
            tmp_savename = "{}.html".format(save_name)
        else:
            tmp_savename = os.path.join(out_dir, "{}.html".format(save_name))

        x = plot(fig, filename=tmp_savename, auto_open=False,
                 include_plotlyjs=False, output_type='div')

        format_ploty(x, tmp_savename)
        # return tmp_savename


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
    legendgroup = 'group_{}'.format(enum)

    g = plotly_graph.Scatter(
            x=x,
            y=y,
            hoveron='text',
            name=label,
            visible=True,
            mode=mode,
            legendgroup=legendgroup,
            showlegend=show,
            line=dict(color=l_color),
            marker=dict(symbol=marker,
                        size=size,
                        color=l_color),
    )
    return g
