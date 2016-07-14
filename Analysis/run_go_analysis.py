import os

import networkx as nx

from Analysis.cytoscape_layout import RenderModel
from Analysis.go_analysis import GoAnalysis
from Analysis.go_network_generator import GoNetworkGenerator


class Analyzer():
    def __init__(self, experimental_data, network, species='hsa', metric='pvalue', output_directory='tmp'):
        self.exp_data = experimental_data
        self.network = network
        self.species = species
        self.out_dir = output_directory
        self.go = GoAnalysis
        self.render_model = RenderModel
        self.go_net_gen = GoNetworkGenerator
        self.metric = metric
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

    def run_go(self, data_type='proteomics', slim=False):
        if slim:
            go = self.go(species=self.species, output_directory='Figures', metric=self.metric, slim='goslim_pir',
                         experimental_data=self.exp_data)
        else:
            go = self.go(species=self.species, output_directory='Figures', metric=self.metric,
                         experimental_data=self.exp_data)

        if data_type == 'proteomics':
            up_reg = self.exp_data.proteomics_up_over_time
            down_reg = self.exp_data.proteomics_down_over_time
            sig_changes = self.exp_data.proteomics_over_time
            labels = list(self.exp_data.protomics_time_points)
            save_name = 'proteomics'

        save_name_up = save_name + '_up'
        save_name_down = save_name + '_down'
        save_name_both = save_name + '_changed'
        go.analysis_data(up_reg, aspect='P', savename=save_name_up, labels=labels, analyze=True)
        go.export_to_html(labels, x=labels, html_name='%s_slim' % save_name_up)
        self.create_up_reg_go_network(go, savename=save_name_up)

        go.analysis_data(down_reg, aspect='P', savename=save_name_down, labels=labels, analyze=True)
        go.export_to_html(labels, x=labels, html_name='%s_slim' % save_name_down)
        self.create_up_reg_go_network(go, savename=save_name_down)

        go.analysis_data(sig_changes, aspect='P', savename=save_name_both, labels=labels, analyze=True)
        go.export_to_html(labels, x=labels, html_name='%s_slim' % save_name_both)
        self.create_up_reg_go_network(go, savename=save_name_both)

    def create_up_reg_go_network(self, go, savename):

        all_timepoints = []
        list_of_timepoints = []
        for n, i in enumerate(go.top_hits):
            remove_zeros(i, n)
            list_of_timepoints.append(i)
            for each in i.keys():
                all_timepoints.append(each)
        print(all_timepoints)
        tall = self.go_net_gen.create_network_from_list(list_of_go_terms=all_timepoints,
                                                        save_name='{0}_network'.format(savename),
                                                        threshold=0)

        # Do these in reverse so we see where the signal started
        # paint_graph(tall, time_point_6, 'orange')
        colors = ['pink', 'yellow', 'green', 'lightblue', 'red']
        for n, each in enumerate(reversed(list_of_timepoints)):
            paint_graph(tall, each, colors[n])

        save_name = '{0}_all_colored'.format(savename)
        nx.nx.write_graphml(tall, 'Figures/{0}_{1}.graphml'.format(save_name, self.metric))
        rm = RenderModel(tall)
        rm.visualize_by_list_of_time(create_names(5), prefix=savename, directory='Figures')

    def run_all(self):
        self.run_go(data_type='proteomics')
        self.run_go(data_type='proteomics', slim=True)


def remove_zeros(top_hits, value):
    to_remove = []
    for i in top_hits:
        if top_hits[i][value] == 0.0:
            to_remove.append(i)
    for i in to_remove:
        top_hits.pop(i)


def paint_graph(graph, time_point, color):
    for i in graph.nodes():
        tmp = graph.node[i]['go']
        tmp = tmp[:2] + ':' + tmp[2:]
        if tmp in time_point:
            graph.node[i]['color'] = color
            for n, time in enumerate(time_point[tmp]):
                t = 'time_{0}'.format(n)
                graph.node[i][t] = float(time)


def create_names(n):
    names = []
    for i in range(n):
        names.append('time_{0}'.format(i))
    return names
