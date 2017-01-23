import os
import warnings
from ast import literal_eval

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from magine.html_templates.html_tools import write_single_table
from magine.networks.cytoscape_view import RenderModel
from magine.networks.go_network_generator import GoNetworkGenerator
from magine.ontology.ontology_analysis import GoAnalysis


class Analyzer:
    """
    MAGINE analyzer
    Class to perform entire MAGINE pipeline.
    Creates a network
    Performs gene enrichment analysis
    Combines the two
    Requires cytoscape session to be opened if you want the output networks
    """

    def __init__(self, experimental_data, network=None, species='hsa',
                 metric='pvalue', output_directory='tmp', save_name='tmp',
                 build_network=False):
        self.exp_data = experimental_data
        self.species = species
        self.save_name = save_name
        self.out_dir = output_directory
        self.go = GoAnalysis
        self.metric = metric
        self.network = None
        if build_network:
            if network is not None:
                warnings.warn("Warning : Passing build_network=True "
                              "and a network file! ", RuntimeWarning)
            from magine.networks.network_generator import build_network
            self.build_network = build_network
            self.generate_network(save_name)
        elif network is not None:
            self.network = network
            self.go_net_gen = None
        if self.network is not None:

            self.go_net_gen = GoNetworkGenerator(self.species, self.network,
                                                 self.out_dir)
        self.html_names = []
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        if not os.path.exists(os.path.join(self.out_dir, 'Figures')):
            os.mkdir(os.path.join(self.out_dir, 'Figures'))
        if not os.path.exists(os.path.join(self.out_dir, 'Network_files')):
            os.mkdir(os.path.join(self.out_dir, 'Network_files'))

    def generate_network(self, save_name):
        """
        Generates network with default parameters

        Parameters
        ----------
        save_name : str
            name of network to be saved

        Returns
        -------

        """
        proteins = self.exp_data.list_sig_proteins
        self.network = self.build_network(proteins, num_overlap=1,
                                          save_name=save_name,
                                          species=self.species)

    def run_go(self, data_type='proteomics', slim=False, metric=None,
               run_type='all', aspect='P', create_go_networks=False):
        """
        Runs enrichment analysis
        Parameters
        ----------
        data_type : str, {'proteomics','rnaseq'}
            data type to run enrichment analysis
        slim : bool
        metric : {'enrichment', 'pvalue', 'fraction'}
            metric to score by
        run_type : {'proteomics', 'down', 'changed', 'all'}
        aspect : {'F','B','C'}
            Can be any combination of GO aspects.
            F: Molecular function
            B: Biological process
            C: Cellular component
        create_go_networks: bool
            create a GO level network

        Returns
        -------

        """
        if metric is None:
            metric = self.metric
            assert metric is not None, 'Must provide metric'

        if data_type == 'proteomics':
            up_reg = list(self.exp_data.proteomics_up_over_time[0:2])
            print(up_reg)
            down_reg = self.exp_data.proteomics_down_over_time
            sig_changes = self.exp_data.proteomics_over_time
            labels = list(self.exp_data.protomics_time_points)
            save_name = self.save_name + '_proteomics'
        if data_type == 'rnaseq':
            up_reg = self.exp_data.rna_up_over_time
            down_reg = self.exp_data.rna_down_over_time
            sig_changes = self.exp_data.rna_over_time
            labels = list(self.exp_data.rna_time_points)
            save_name = self.save_name + '_rnaseq'
        if slim:
            save_name += '_slim'
            go = self.go(species=self.species, output_directory=self.out_dir,
                         metric=metric, slim='goslim_pir',
                         experimental_data=self.exp_data)
        else:
            go = self.go(species=self.species, output_directory=self.out_dir,
                         metric=metric, experimental_data=self.exp_data)

        if run_type == 'up' or run_type == 'all':
            save_name_up = save_name + '_up' + '_' + aspect
            self.html_names.append(save_name_up)
            # '''
            go.analyze_data(up_reg, aspect=aspect, savename=save_name_up,
                            labels=labels, analyze=True)
            go.export_to_html(labels, html_name=save_name_up)
            # go.print_ranked_over_time(create_plots=False, number=5)
            if create_go_networks:
                self.create_go_network_and_render(go, savename=save_name_up)
                # '''
        if run_type == 'down' or run_type == 'all':
            save_name_down = save_name + '_down' + '_' + aspect
            self.html_names.append(save_name_down)
            # '''
            go.analyze_data(down_reg, aspect=aspect, savename=save_name_down,
                            labels=labels, analyze=True)
            # go.export_to_html(labels, html_name=save_name_down)
            # go.print_ranked_over_time(create_plots=False, number=5)
            if create_go_networks:
                self.create_go_network_and_render(go, savename=save_name_down)
                # '''

        if run_type == 'changed' or run_type == 'all':
            save_name_both = save_name + '_changed' + '_' + aspect
            self.html_names.append(save_name_both)
            # '''
            go.analyze_data(sig_changes, aspect=aspect,
                            savename=save_name_both, labels=labels,
                            analyze=True)
            # go.export_to_html(labels, html_name=save_name_both)
            # go.print_ranked_over_time(create_plots=False, number=5)
            if create_go_networks:
                self.create_go_network_and_render(go, savename=save_name_both)
                # '''

        return go

    def run_single_go(self, data_type='proteomics', fold_change='up',
                      slim=False, metric=None, aspect='P',
                      create_go_networks=False):
        """
        Runs enrichment analysis
        Parameters
        ----------
        data_type : str, {'proteomics','rnaseq'}
            data type to run enrichment analysis
        slim : bool
        metric : {'enrichment', 'pvalue', 'fraction'}
            metric to score by
        aspect : {'F','B','C'}
            Can be any combination of GO aspects.
            F: Molecular function
            B: Biological process
            C: Cellular component
        create_go_networks: bool
            create a GO level network

        Returns
        -------

        """
        if metric is None:
            metric = self.metric
            assert metric is not None, 'Must provide metric'

        if data_type == 'proteomics':
            if fold_change == 'up':
                data = self.exp_data.proteomics_up_over_time
            elif fold_change == 'down':
                data = self.exp_data.proteomics_down_over_time
            elif fold_change == 'both':
                data = self.exp_data.proteomics_over_time
            labels = list(self.exp_data.protomics_time_points)
        elif data_type == 'rnaseq':
            if fold_change == 'up':
                data = self.exp_data.rna_up_over_time
            elif fold_change == 'down':
                data = self.exp_data.rna_down_over_time
            elif fold_change == 'both':
                data = self.exp_data.rna_over_time
            labels = list(self.exp_data.rna_time_points)
        save_name = '_'.join([self.save_name, data_type, fold_change, aspect])

        if slim:
            save_name += '_slim'
            slim_name = 'goslim_pir'
        else:
            slim_name = None

        filename = '<a href="{}.html">link</a>'.format(
            self.out_dir + '/' + save_name)

        html_dict = {'DataType': data_type,
                     'FoldChange': fold_change,
                     'fileName': filename,
                     'aspect': aspect,
                     'slim': slim}

        go = self.go(species=self.species, output_directory=self.out_dir,
                     metric=metric, slim=slim_name,
                     experimental_data=self.exp_data)

        go.analyze_data(data, aspect=aspect, savename=save_name,
                        labels=labels, analyze=True)
        go.export_to_html(labels, html_name=save_name)
        # go.print_ranked_over_time(create_plots=False, number=5)
        if create_go_networks:
            self.create_go_network_and_render(go, savename=save_name)

        return go, html_dict

    def create_go_network_and_render(self, go, savename):
        """
        Creates GO level network using provided GO rankings

        Parameters
        ----------
        go : list of array_like
            list of arrays of go rankings per sample
        savename : str
            name of file to save GO level network

        Returns
        -------

        """

        all_timepoints = []
        list_of_timepoints = []

        for n, i in enumerate(go.top_hits):
            remove_zeros(i, n)
            list_of_timepoints.append(i)
            for each in i.keys():
                all_timepoints.append(each)

        if self.go_net_gen is None:
            warnings.warn("Warning : You must provide Analyzer with a network "
                          "or pass the build_network flag"
                          "or generate a network in order to generate"
                          " a go network", RuntimeError)

        tall = self.go_net_gen.create_network_from_list(
            list_of_go_terms=all_timepoints,
            save_name='{0}_network'.format(savename),
            threshold=0)
        if len(tall.nodes()) == 0:
            return
        # Do these in reverse so we see where the signal started
        for n, each in enumerate(reversed(list_of_timepoints)):
            paint_graph(tall, each, colors[n])

        save_name = os.path.join(self.out_dir, 'Network_files',
                                 '{0}_all_colored_{1}.graphml'.format(savename,
                                                                      self.metric))
        nx.nx.write_graphml(tall, save_name)
        rm = RenderModel(tall)
        rm.visualize_by_list_of_time(create_names(len(go.top_hits)),
                                     prefix=savename,
                                     out_dir=self.out_dir)

    def create_go_network_from_list(self, go, terms, savename):

        all_timepoints = []
        list_of_timepoints = []
        for n, i in enumerate(go.top_hits):
            remove_zeros(i, n)
            list_of_timepoints.append(i)
            for each in i.keys():
                if each in terms:
                    all_timepoints.append(each)
        tall = self.go_net_gen.create_network_from_list(
            list_of_go_terms=all_timepoints,
            save_name='{0}_network'.format(savename), threshold=0)

        # Do these in reverse so we see where the signal started
        for n, each in enumerate(reversed(list_of_timepoints)):
            paint_graph(tall, each, colors[n])

        save_name = os.path.join(self.out_dir, 'Network_files',
                                 '{0}_all_colored_{1}.graphml'.format(
                                     savename, self.metric))
        nx.nx.write_graphml(tall, save_name)
        rm = RenderModel(tall)
        rm.visualize_by_list_of_time(
            create_names(len(self.exp_data.protomics_time_points)),
            prefix=savename,
            out_dir='.')

    def create_html_report(self, d_type):
        if d_type == 'proteomics':
            out = proteomics_html.format(*tuple(self.html_names))
        if d_type == 'rna':
            out = rna_html.format(*tuple(self.html_names))
        if d_type == 'all':
            print(len(self.html_names[0:9]))
            out = proteomics_html.format(*tuple(self.html_names[0:9]))
            out2 = rna_html.format(*tuple(self.html_names[9:]))
            out += out2
        output_file = os.path.join(self.out_dir,
                                   self.save_name + '_all_analysis.html')
        with open(output_file, 'w') as f:
            f.write(out)
        pass

    def run_all(self, data_type='proteomics'):
        """ performs all analysis

        :return:
        """

        if data_type == 'proteomics' or data_type == 'all':
            # self.run_go(data_type='proteomics', aspect='P')
            # self.run_go(data_type='proteomics', aspect='F')
            # self.run_go(data_type='proteomics', aspect='C')
            self.run_go(data_type='proteomics', slim=True, aspect='P')
            # self.run_go(data_type='proteomics', slim=True, aspect='F')
            # self.run_go(data_type='proteomics', slim=True, aspect='C')
            # if data_type == 'rnaseq' or data_type == 'all':
            # self.run_go(data_type='rnaseq', aspect='P')
            # self.run_go(data_type='rnaseq', aspect='F')
            # self.run_go(data_type='rnaseq', aspect='C')
            # self.run_go(data_type='rnaseq', slim=True, aspect='P')
            # self.run_go(data_type='rnaseq', slim=True, aspect='F')
            # self.run_go(data_type='rnaseq', slim=True, aspect='C')

        self.create_html_report(data_type)

    def run_proteomics_go(self):
        """
        Updated function to run all proteomics
        Returns
        -------

        """

        list_of_go_dict = []

        for i in ['up', 'down', 'both']:
            go, x = self.run_single_go(data_type='proteomics',
                                       fold_change=i, slim=True,
                                       metric='enrichment', aspect='P',
                                       create_go_networks=True)
            list_of_go_dict.append(x)
            go, x = self.run_single_go(data_type='proteomics',
                                       fold_change=i, slim=False,
                                       metric='enrichment', aspect='P',
                                       create_go_networks=False)
            list_of_go_dict.append(x)
            go, x = self.run_single_go(data_type='rnaseq',
                                       fold_change=i, slim=False,
                                       metric='enrichment', aspect='P',
                                       create_go_networks=False)
            list_of_go_dict.append(x)

        # self.create_selected_go_network(
        #     '20161026_enrichment_proteomics_up_P_all_data.csv',
        #     'proteomics_up')
        # self.create_selected_go_network(
        #     '20161026_enrichment_rnaseq_up_P_all_data.csv',
        #     'rnaseq_up')

        df = pd.DataFrame(list_of_go_dict)
        write_single_table(df, 'index.html', 'MAGINE Output')

    def create_selected_go_network(self, file_name, save_name):
        data = pd.read_csv(file_name, converters={"genes": literal_eval},
                           index_col=0)
        data['genes'] = data['genes'].astype(set)
        labels = self.exp_data.timepoints
        data = data[data['ref'] > 6]

        # data = data[data['depth'] < 10]
        data = data.sort_values('enrichment_score', ascending=False)
        tmp = data[data['pvalue'] < 0.05]

        tmp = pd.pivot_table(tmp, index=['GO_id', ], columns='sample_index')
        n_samples = range(len(labels))
        list_all_go = set()
        for i in n_samples:
            tmp = tmp.sort_values(by=('enrichment_score', i), ascending=False)
            list_of_go = set(tmp.head(8).index)
            # list_of_go = check_term_list(list_of_go)
            list_all_go.update(list_of_go)
            # g = find_disjunction_common_ancestor(list_of_terms=list_of_go)
            # export_to_dot(g, save_name='{}_graph'.format(i), view=False)

        # g = find_disjunction_common_ancestor(list_of_terms=list_all_go)
        # export_to_dot(g, save_name='all_go_graph', view=False)
        # trim_nodes = True
        trim_nodes = True
        if trim_nodes:
            from magine.ontology.ontology_tools import check_term_list
            list_all_go = check_term_list(list_all_go)
            data = data[data['GO_id'].isin(list_all_go)]

        data.loc[data['enrichment_score'] > 20, 'enrichment_score'] = 20

        tmp = pd.pivot_table(data, index=['GO_id', 'GO_name'],
                             columns='sample_index')
        tmp2 = pd.pivot_table(data, index=['GO_id'], columns='sample_index')
        enrichment_list = [('enrichment_score', i) for i in n_samples]
        tmp = tmp.sort_values(by=enrichment_list, ascending=False)

        # Want to show both enrichment as a panel as pvalie
        pvalue = tmp['pvalue'].fillna(1).as_matrix()
        enrichment_value = tmp['enrichment_score'].as_matrix()
        fig = plt.figure(figsize=(12, 6))

        # Enrichment panel
        ax1 = fig.add_subplot(111)
        plt.imshow(enrichment_value, aspect='auto',
                   interpolation='nearest', cmap=plt.get_cmap('Blues'))

        # text portion
        y_array = range(0, len(tmp.index), 1)
        x_array = range(0, len(labels), 1)

        x, y = np.meshgrid(x_array, y_array)

        for x_val, y_val in zip(x.flatten(), y.flatten()):
            if pvalue[y_val, x_val] < 0.05:
                ax1.text(x_val, y_val, '*', va='center', ha='center')

        ax1.tick_params(axis='both', direction='out')
        ax1.set_xticks(range(len(labels)))
        ax1.set_xticklabels(labels, fontsize=16)

        ax1.set_yticks(range(len(tmp.index)))
        ax1.set_yticklabels(tmp.index, fontsize=16)
        plt.colorbar()
        plt.tight_layout(w_pad=2.8, h_pad=.6)
        plt.savefig('{}_heatplot.png'.format(save_name), bbox_inches='tight')
        plt.close()

        if self.go_net_gen is None:
            warnings.warn("Warning : You must provide Analyzer with a network "
                          "or pass the build_network flag"
                          "or generate a network in order to generate"
                          " a go network", RuntimeError)
        tall = self.go_net_gen.create_network_from_list(
            list_of_go_terms=list_all_go,
            save_name='{0}_network'.format(save_name),
            threshold=0)
        if len(tall.nodes()) == 0:
            print('No nodes')
            quit()

        score_array = tmp2.copy()

        x = score_array['enrichment_score'].fillna(0)

        for i in tall.nodes():
            tmp = tall.node[i]['go']
            tmp = tmp[:2] + ':' + tmp[2:]
            values = x.loc[tmp]
            tall.node[i]['color'] = 'red'
            for n, time in enumerate(n_samples):
                t = 'time_{0:04d}'.format(n)
                tall.node[i][t] = float(values[time])

        savename = os.path.join(self.out_dir, 'Network_files',
                                '{0}_all_colored_{1}.graphml'.format(
                                    save_name, self.metric))

        nx.nx.write_graphml(tall, savename)
        size_of_data = len(n_samples)
        rm = RenderModel(tall)
        rm.visualize_by_list_of_time(create_names(size_of_data),
                                     prefix=save_name,
                                     labels=self.exp_data.timepoints,
                                     out_dir=self.out_dir)
headers = {'DataType',
           'FoldChange',
           'fileName',
           'aspect'}


def remove_zeros(top_hits, value):
    """ removes any GO scores that are zero enriched

    :param top_hits: dictionary of GO enrichments
    :param value: index to find the score
    :return:
    """
    to_remove = []
    for i in top_hits:
        if top_hits[i][value] == 0.0:
            to_remove.append(i)
    for i in to_remove:
        top_hits.pop(i)


def paint_graph(graph, time_point, color):
    """ sets color attribute according to time

    :param graph: networkx graph
    :param time_point: dictionary of enrichment scores
    :param color:
    :return:
    """
    for i in graph.nodes():
        tmp = graph.node[i]['go']
        tmp = tmp[:2] + ':' + tmp[2:]
        if tmp in time_point:
            graph.node[i]['color'] = color  # 'red'
            for n, time in enumerate(time_point[tmp]):
                t = 'time_{0:04d}'.format(n)
                graph.node[i][t] = float(time)
                # graph.node[i][t] = 25#float(time)


def create_names(n):
    names = []
    for i in range(n):
        names.append('time_{0:04d}'.format(i))
        print('time_{0:04d}'.format(i))
    return names


colors = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6",
          "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762",
          "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6",
          "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A",
          "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
          "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018",
          "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
          "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C",
          "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1",
          "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459",
          "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA",
          "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81", "#575329",
          "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1",
          "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
          "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
          "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534",
          "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72",
          "#6A3A4C", ]

proteomics_html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Analysis</title>
</head>
<body>
<h1> Proteomics GO enrichment analysis</h1>
<p>Up regulated genes only</p>
<a href="{0}.html">Up regulated biological processes</a><br>
<a href="{3}.html">Up regulated molecular function</a><br>
<a href="{6}.html">Up regulated cellular component</a><br>
<p>Down regulated only</p>
<a href="{1}.html">Down regulated biological processes</a><br>
<a href="{4}.html">Down regulated molecular function</a><br>
<a href="{7}.html">Down regulated cellular component</a><br>
<p>Both up and down regulated genes.</p>
<a href="{2}.html">Absolute change biological processes</a><br>
<a href="{5}.html">Absolute change molecular function</a><br>
<a href="{8}.html">Absolute change cellular component</a><br>
"""
prot_slim_html = """
<h1>GO enrichment analysis with slim terms</h1>
<p>GO slim enrichment analysis </p>
<p>Up regulated genes only</p>
<a href="{9}.html">Up regulated biological processes</a><br>
<a href="{12}.html">Up regulated molecular function</a><br>
<a href="{15}.html">Up regulated cellular component</a><br>
<p>Down regulated only</p>
<a href="{10}.html">Down regulated biological processes</a><br>
<a href="{13}.html">Down regulated molecular function</a><br>
<a href="{16}.html">Down regulated cellular component</a><br>
<p>Both up and down regulated genes.</p>
<a href="{11}.html">Absolute change biological processes</a><br>
<a href="{14}.html">Absolute change molecular function</a><br>
<a href="{17}.html">Absolute change cellular component</a><br>
"""

rna_html = """
<h1>RNA seq GO enrichment analysis</h1>
<p>Up regulated genes only</p>
<a href="{0}.html">Up regulated biological processes</a><br>
<a href="{3}.html">Up regulated molecular function</a><br>
<a href="{6}.html">Up regulated cellular component</a><br>
<p>Down regulated only</p>
<a href="{1}.html">Down regulated biological processes</a><br>
<a href="{4}.html">Down regulated molecular function</a><br>
<a href="{7}.html">Down regulated cellular component</a><br>
<p>Both up and down regulated genes.</p>
<a href="{2}.html">Absolute change biological processes</a><br>
<a href="{5}.html">Absolute change molecular function</a><br>
<a href="{8}.html">Absolute change cellular component</a><br>
</body>
</html>
"""

network_images = """
<p>Network images from up regulated genes.
    <br>
    <a href="prot_up_time_0_wpr.pdf">Time 30 min</a>
    <br>
    <a href="prot_up_time_1_wpr.pdf">Time 1 hour</a>
    <br>
    <a href="prot_up_time_2_wpr.pdf">Time 6 hour</a>
    <br>
    <a href="prot_up_time_3_wpr.pdf">Time 12 hour</a>
    <br>
    <a href="prot_up_time_4_wpr.pdf">Time 48 hour</a>
    <br>
</p>
<p>Network images from down regulated genes.
    <br>
    <a href="prot_down_time_0_wpr.pdf">Time 30 min</a>
    <br>
    <a href="prot_down_time_1_wpr.pdf">Time 1 hour</a>
    <br>
    <a href="prot_down_time_2_wpr.pdf">Time 6 hour</a>
    <br>
    <a href="prot_down_time_3_wpr.pdf">Time 12 hour</a>
    <br>
    <a href="prot_down_time_4_wpr.pdf">Time 48 hour</a>
    <br>
</p>
<p>Network images from significantly changes genes.
    <br>
    <a href="prot_changes_time_0_wpr.pdf">Time 30 min</a>
    <br>
    <a href="prot_changes_time_0_wpr.pdf">Time 1 hour</a>
    <br>
    <a href="prot_changes_time_0_wpr.pdf">Time 6 hour</a>
    <br>
    <a href="prot_changes_time_0_wpr.pdf">Time 12 hour</a>
    <br>
    <a href="prot_changes_time_0_wpr.pdf">Time 48 hour</a>
    <br>
</p>
"""
