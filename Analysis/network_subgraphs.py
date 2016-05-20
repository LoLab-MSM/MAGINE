import networkx as nx
import pygraphviz as pyg

from magine.network_tools import compress_edges
import os


class NetworkSubgraphs:
    """ Class to create subgraphs of larger network

    """
    def __init__(self, network, exp_data):
        self.network = network
        self.nodes = self.network.nodes()
        self.exp_data = exp_data

    def generate_shortest_paths_between_two_proteins(self, protein_1, protein_2):
        """
        Generates a graph based on all shortest paths between two proteins
        :param protein_1:
        :param protein_2:
        :return: graph
        """
        graph = pyg.AGraph(directed=True)
        if not nx.has_path(self.network, protein_1, protein_2):
            print("Path does NOT exist between %s and %s" % (protein_1, protein_2))
        else:
            for path in nx.all_shortest_paths(self.network, protein_1, protein_2):
                self.add_edges_from_path(graph, path)
        if not nx.has_path(self.network, protein_2, protein_1):
            print("Path does NOT exist between %s and %s" % (protein_2, protein_1))
        else:
            for path in nx.all_shortest_paths(self.network, protein_2, protein_1):
                self.add_edges_from_path(graph, path)
        #self.create_legend(graph)
        self.write(graph)
        graph.write("%s_and_%s.dot" % (protein_1, protein_2))
        graph.draw("%s_and_%s.pdf" % (protein_1, protein_2), prog='dot')
        return graph

    def generate_shortest_paths_between_lists_of_proteins(self, protein_list,savename):
        """
        Generates a graph based on all shortest paths between two proteins
        :param protein_list:
        :param savename:
        :return: graph
        """
        graph = pyg.AGraph(directed=True)
        for protein_1 in protein_list:
            if not protein_1 in self.network.nodes():
                continue
            for protein_2 in protein_list:
                if not protein_2 in self.network.nodes():
                    continue
                if protein_1 == protein_2:
                    continue
                else:
                    if not nx.has_path(self.network, protein_1, protein_2):
                        print("Path does NOT exist between %s and %s" % (protein_1, protein_2))
                    else:
                        for path in nx.all_shortest_paths(self.network, protein_1, protein_2):
                            self.add_edges_from_path(graph, path)
                    if not nx.has_path(self.network,  protein_2,protein_1):
                            print("Path does NOT exist between %s and %s" % (protein_2, protein_1))
                    else:
                        for path in nx.all_shortest_paths(self.network, protein_2, protein_1):
                            self.add_edges_from_path(graph, path)
        # self.create_legend(graph)
        self.write(graph)
        graph.write("%s.dot" % savename)
        graph.draw("%s.pdf" % savename, prog='dot')
        return graph

    #TODO create a function that finds all neighbors of each proteins and links them
    # THis would find more than the shortest path
    def generate_all_paths_between_two_proteins(self, protein_1, protein_2):
        """
        Generates a graph based on all shortest paths between two proteins
        :param protein_1:
        :param protein_2:
        :return: graph
        """
        graph = pyg.AGraph(directed=True)
        if not nx.has_path(self.network, protein_1, protein_2):
            print("Path does NOT exist between %s and %s" % (protein_1, protein_2))
            return graph
        if not nx.has_path(self.network, protein_2, protein_1):
            print("Path does NOT exist between %s and %s" % (protein_2, protein_1))
            return graph
        graph = pyg.AGraph(directed=True)
        for path in nx.all_shortest_paths(self.network, protein_1, protein_2):
            self.add_edges_from_path(graph, path)
        for path in nx.all_shortest_paths(self.network, protein_2, protein_1):
            self.add_edges_from_path(graph, path)
        self.create_legend(graph)
        graph.write("%s_and_%s.dot" % (protein_1, protein_2))
        graph.draw("%s_and_%s.pdf" % (protein_1, protein_2), prog='dot')
        return graph

    def generate_shortest_paths_from_protein_to_protein(self, protein_1, protein_2):
        """
        Generates a graph from all the shortests paths FROM protein 1 to protein 2
        :param protein_1: string
        :param protein_2: string
        :return:
        """
        graph = pyg.AGraph(directed=True)
        if not nx.has_path(self.network, protein_1, protein_2):
            print("Path does NOT exist between %s and %s" % (protein_1, protein_2))
            return graph
        short_distance = nx.shortest_path_length(self.network, protein_1, protein_2)
        print("Shortest path length from %s to %s = %i" % (protein_1, protein_2, short_distance))

        for path in nx.all_shortest_paths(self.network, protein_1, protein_2):
            self.add_edges_from_path(graph, path)
        #self.create_legend(graph)
        graph.write("%s_to_%s.dot" % (protein_1, protein_2))
        graph.draw("%s_to_%s.pdf" % (protein_1, protein_2), prog='dot')
        return graph

    def generate_paths_between_proteins_list(self, protein_list, save_name):
        """ Generate s graph from paths from a list of pair of proteins

        :param protein_list:
        :param save_name:
        :return:
        """
        protein_pair_graph = pyg.AGraph(directed=True)
        for pair in protein_list:
            print(pair)
            graph = self.generate_shortest_paths_from_protein_to_protein(pair[0],pair[1])
            for node in graph.nodes():
                print(node)
                protein_pair_graph.add_node(node, **self.network.node[node])
            for edge in graph.edges():
                protein_pair_graph.add_edge(edge[0],edge[1],**self.network.edge[edge[0]][edge[1]])
        protein_pair_graph.draw("%s.pdf" % save_name, prog='dot')
        return protein_pair_graph


    def generate_upstream_network_of_protein(self, protein_1, savename='test',compress=False):
        """
        Finds all linear pathways leading up to a protein.
        :param protein_1: string
        :param savename: string
        :return: graph
        """
        paths = []
        for i in self.nodes:
            if i in self.nodes:
                if nx.has_path(self.network, i, protein_1):
                    path = nx.shortest_path(self.network, i, protein_1)
                    paths.append(path)
        temporary_graph = pyg.AGraph(directed=True)
        for path in paths:
            cont = True
            for protein in list(path):
                if protein not in self.exp_data.sign_proteins_omics+self.exp_data.metabolomics:
                    cont = False
            if not cont:
                continue
            self.add_edges_from_path(temporary_graph, path)
        if compress:
            temporary_graph = compress_edges(temporary_graph)
        #self.create_legend(temporary_graph)
        temporary_graph.write('%s.dot' % savename)
        temporary_graph.draw('%s.pdf' % savename, prog='dot')
        return temporary_graph

    def generate_downstream_network_of_protein(self, network, protein_1, exp_data, savename='test', compress=False):
        """
        Finds all linear pathways leading up to a protein.
        :param protein_1: string
        :param savename: string
        :return: graph
        """
        paths = []
        nodes = network.nodes()
        for i in nodes:
            if i in nodes:
                if nx.has_path(network, protein_1,i):
                    path = nx.shortest_path(network, protein_1,i)
                    paths.append(path)
        temporary_graph = pyg.AGraph(directed=True)
        for path in paths:
            cont = True
            for protein in list(path):
                if protein not in self.exp_data.measured:
                    cont = False
            if not cont:
                continue
            self.add_edges_from_path(temporary_graph, path)
        if compress:
            temporary_graph = compress_edges(temporary_graph)
        #self.create_legend(temporary_graph)
        temporary_graph.write('%s.dot' % savename)
        temporary_graph.draw('%s.pdf' % savename, prog='dot')
        return temporary_graph

    def write(self,network):
        dict_of_types = {'activation': 'onormal',
                         'indirect effect': 'odiamondodiamond',
                         'expression': 'normal',
                         'inhibition': 'tee',
                         'binding/association': 'curve',
                         'phosphorylation': 'dot',
                         'missing interaction': 'odiamond',
                         'compound': 'dotodot',
                         'dissociation': 'diamond',
                         'ubiquitination': 'oldiamond',
                         'state change': 'teetee',
                         'dephosphorylation': 'onormal',
                         'repression': 'obox',
                         'glycosylation': 'dot'}
        activators = ['activation', 'expression','phosphorylation']
        inhibitors = ['inhibition','repression','dephosphorylation']
        physical_contact = ['binding/association','dissociation','state change']
        chemical = ['compound','glycosylation','ubiquitination']
        for i in network.edges():
            n = network.get_edge(i[0],i[1])
            edge_type = str(i.attr['interactionType'])
            cont = False
            for j in activators:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'normal'
                    cont = True
            if cont:
                continue
            for j in inhibitors:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'tee'
                    cont = True
            if cont:
                continue
            for j in physical_contact:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'diamond'
                    cont = True
            if cont:
                continue
            for j in chemical:
                if edge_type.startswith(j):
                    n.attr['arrowhead'] = 'dotodot'
                    cont = True
            if cont:
                continue
            else:
                print(n,n.attr)





    def paint_network_overtime(self, graph, list_of_lists, color_list,savename):
        """

        :param graph:
        :param list_of_lists:
        :param color_list:
        :return:
        """
        if len(list_of_lists) != len(color_list):
            print('Length of list of data must equal len of color list')
        string = 'convert '
        tmp_graph = graph.copy()
        for n, i in enumerate(list_of_lists):
            graph2 = paint_network(tmp_graph, i, color_list[n])
            self.write(graph2)
            graph2.draw('%s_%04i.png' % (savename, n), prog='dot')
            string += '%s_%04i.png ' % (savename, n)
        string += '  %s.pdf' % savename
        print(string)
        os.system(string)

    def add_edges_from_path(self, graph, path):
        """
        Adds a path to a graph by extracting the edges and edge attributes
        from the ddn
        :param graph:
        :param path:
        :return: None
        """
        previous = None
        for protein in list(path):
            if previous is None:
                previous = protein
                continue
            else:
                graph.add_edge(previous, protein, **self.network.edge[previous][protein])
                previous = protein


    def create_legend(self, graph):
        """
        adds a legend to a graph
        :param graph:
        :return:
        """
        dict_of_types = {'activation': 'onormal',
                         'indirect effect': 'odiamondodiamond',
                         'expression': 'normal',
                         'inhibition': 'tee',
                         'binding/association': 'curve',
                         'phosphorylation': 'dot',
                         'missing interaction': 'odiamond',
                         'compound': 'dotodot',
                         'dissociation': 'diamond',
                         'ubiquitination': 'oldiamond',
                         'state change': 'teetee',
                         'dephosphorylation': 'onormal',
                         'repression': 'obox',
                         'glycosylation': 'dot'}
        subgraph = []
        len_dic = len(dict_of_types)
        for n, i in enumerate(dict_of_types):
            subgraph.append(n)
            subgraph.append(n + len_dic)
            graph.add_node(n, label="")
            graph.add_node(n + len_dic, label="")
            graph.add_edge(n, n + len_dic, dir='both', arrowhead=dict_of_types[i], arrowtail="none", label=i)
        graph.add_subgraph(subgraph, name='cluster_legend', rank="LR")


def paint_network(graph, list_to_paint, color):
    """
    Paints a graph given a list of nodes and a color for that list
    :param graph: pygraphvix.AGraph
    :param list_to_paint: list
    :param color: string
    :return:
    """
    tmp_g = graph.copy()
    nodes1 = tmp_g.nodes()
    for i in list_to_paint:
        if i in nodes1:
            n = tmp_g.get_node(i)
            n.attr['measured'] = 'True'
            n.attr['color'] = 'black'
            n.attr['fillcolor'] = color
            n.attr['style'] = 'filled'
    return tmp_g