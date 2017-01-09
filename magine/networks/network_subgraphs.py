import itertools
import os
import warnings

import networkx as nx

from magine.network_tools import compress_edges, export_to_dot


class NetworkSubgraphs:
    def __init__(self, network, exp_data=None):
        """
        Generates network subgraphs

        Parameters
        ----------
        network : networkx.DiGraph
        exp_data : magine.data_handler.ExperimentalData
        """
        self.network = network
        self.nodes = set(self.network.nodes())
        self.exp_data = exp_data

    def _add_edges_from_path(self, graph, path):
        """
        Adds a path to a graph

        Parameters
        ----------
        graph : networkx.DiGraph
            graph to add paths
        path : list_like
            list of species that create a path

        """
        previous = None
        for protein in list(path):
            if previous is None:
                previous = protein
                continue
            else:
                graph.add_node(previous, **self.network.node[previous])
                graph.add_node(protein, **self.network.node[protein])
                graph.add_edge(previous, protein,
                               **self.network.edge[previous][protein])
                previous = protein

    def shortest_paths_between_two_proteins(self, node_1, node_2, draw=False):
        """
        Generates a graph based on all shortest paths between two species


        Parameters
        ----------
        node_1 : str
            name of first species
        node_2 : str
            name of second species
        draw : bool
            create an image of returned network


        Returns
        -------
        graph : networkx.DiGraph


        Example
        -------
        >>> from networkx import DiGraph
        >>> from magine.networks.network_subgraphs import NetworkSubgraphs
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
        >>> path_a_d = net_sub.shortest_paths_between_two_proteins('a','d')
        >>> path_a_d.edges()
        [('a', 'd')]
        >>> path_a_e = net_sub.shortest_paths_between_two_proteins('a','e')

        """
        # graph = pyg.AGraph(directed=True)
        graph = nx.DiGraph()
        direction_1, direction_2 = True, True
        if not nx.has_path(self.network, node_1, node_2):
            direction_1 = False
        else:
            for path in nx.all_shortest_paths(self.network, node_1, node_2):
                self._add_edges_from_path(graph, path)
        if not nx.has_path(self.network, node_2, node_1):
            direction_2 = False
        else:
            for path in nx.all_shortest_paths(self.network, node_2, node_1):
                self._add_edges_from_path(graph, path)
        if not direction_1 and not direction_2:
            warnings.warn("Warning : No paths between {} and {}. Returning "
                          "None".format(node_1, node_2), RuntimeWarning)
            return None

        nx.write_gml(graph, "%s_and_%s.gml" % (node_1, node_2))
        if draw:
            export_to_dot(graph, format='png', engine='dot', view=False)
            graph.draw("%s_and_%s.pdf" % (node_1, node_2), prog='dot')
        return graph

    def shortest_paths_between_lists_of_proteins(self, species_list,
                                                 save_name='tmp',
                                                 single_path=False,
                                                 draw=False):
        """

        Parameters
        ----------
        species_list : list_like
            list of species
        save_name : str
            name to save
        single_path : bool
            use single shortest path if True, else use all shortest paths
        draw : bool
            create a dot generated figure

        Returns
        -------
        graph : networkx.DiGraph
            graph containing paths between species list provided


        Example
        -------
        >>> from networkx import DiGraph
        >>> from magine.networks.network_subgraphs import NetworkSubgraphs
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
        >>> path_a_d = net_sub.shortest_paths_between_lists_of_proteins(['a','c','d'])
        Path does NOT exist between c and a
        Path does NOT exist between d and a
        Path does NOT exist between d and c
        >>> path_a_d.edges()
        [('a', 'b'), ('a', 'd'), ('c', 'd'), ('b', 'c')]
        """

        graph = nx.DiGraph()
        nodes = self.nodes
        tmp_protein_list = set(species_list)
        for i in species_list:
            if i not in nodes:
                print('{} not in network'.format(i))
                tmp_protein_list.remove(i)
        node_pairs = itertools.combinations(tmp_protein_list, 2)
        for node_1, node_2 in node_pairs:
            if not nx.has_path(self.network, node_1, node_2):
                print(
                    "Path does NOT exist between %s and %s" % (node_1, node_2))
            else:
                if single_path:
                    self._add_edges_from_path(graph, nx.shortest_path(
                        self.network, node_1, node_2))
                else:
                    for path in nx.all_shortest_paths(self.network, node_1,
                                                      node_2):
                        self._add_edges_from_path(graph, path)
            if not nx.has_path(self.network, node_2, node_1):
                print("Path does NOT exist between %s and %s" % (node_2,
                                                                 node_1))
            else:
                if single_path:
                    self._add_edges_from_path(graph, nx.shortest_path(
                        self.network, node_2, node_1))
                else:
                    for path in nx.all_shortest_paths(self.network, node_2,
                                                      node_1):
                        self._add_edges_from_path(graph, path)
        if draw:
            export_to_dot(graph, save_name=save_name)
        nx.write_gml(graph, "{}.gml".format(save_name))
        return graph

    def upstream_network_of_specie(self, species_1, include_list=None,
                                   save_name='test', compress=False,
                                   draw=False, ):
        """
        Generate network of all upstream species of provides species


        Parameters
        ----------
        species_1 : str
            species name
        save_name : str
            name to save gml file
        compress : bool
            compress the graph by making species in linear path a single node
        draw : bool
            create figure of graph
        include_list : list_like
            list of species that must be in path in order to consider a path
        Returns
        -------
        graph : networkx.DiGraph


        Example
        -------
        >>> from networkx import DiGraph
        >>> from magine.networks.network_subgraphs import NetworkSubgraphs
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
        >>> upstream_d = net_sub.upstream_network_of_specie('d')
        >>> upstream_d.edges()
        [('a', 'd'), ('c', 'd'), ('b', 'c'), ('e', 'd')]
        >>> upstream_c = net_sub.upstream_network_of_specie('c')
        >>> upstream_c.edges()
        [('a', 'b'), ('b', 'c')]

        """
        if include_list is None:
            include_list = self.nodes

        graph = nx.DiGraph()
        for i in self.nodes:
            if i in include_list:
                if nx.has_path(self.network, i, species_1):
                    path = nx.shortest_path(self.network, i, species_1)
                    self._add_edges_from_path(graph, path)

        if compress:
            graph = compress_edges(graph)

        if draw:
            export_to_dot(graph, save_name)
        nx.write_gml(graph, "{}.gml".format(save_name))

        return graph

    def downstream_network_of_specie(self, species_1, include_list=None,
                                     save_name='test', compress=False,
                                     draw=False, ):
        """
        Generate network of all downstream species of provides species


        Parameters
        ----------
        species_1 : str
            species name
        save_name : str
            name to save gml file
        compress : bool
            compress the graph by making species in linear path a single node
        draw : bool
            create figure of graph
        include_list : list_like
            list of species that must be in path in order to consider a path
        Returns
        -------
        graph : networkx.DiGraph


        Example
        -------
        >>> from networkx import DiGraph
        >>> from magine.networks.network_subgraphs import NetworkSubgraphs
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
        >>> upstream_d = net_sub.upstream_network_of_specie('d')
        >>> upstream_d.edges()
        [('a', 'd'), ('c', 'd'), ('b', 'c'), ('e', 'd')]
        >>> upstream_c = net_sub.upstream_network_of_specie('c')
        >>> upstream_c.edges()
        [('a', 'b'), ('b', 'c')]

        """
        if include_list is None:
            include_list = self.nodes

        graph = nx.DiGraph()
        for i in self.nodes:
            if i in include_list:
                if nx.has_path(self.network, species_1, i):
                    path = nx.shortest_path(self.network, species_1, i)
                    self._add_edges_from_path(graph, path)

        if compress:
            graph = compress_edges(graph)

        if draw:
            export_to_dot(graph, save_name)
        nx.write_gml(graph, "{}.gml".format(save_name))

        return graph

    # deprecated function
    def paint_network_overtime(self, graph, list_of_lists, color_list,
                               savename):
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
            write(graph2)
            graph2.draw('%s_%04i.png' % (savename, n), prog='dot')
            string += '%s_%04i.png ' % (savename, n)
        string += '  %s.pdf' % savename
        print(string)
        os.system(string)



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


def write(network):
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
    activators = ['activation', 'expression', 'phosphorylation']
    inhibitors = ['inhibition', 'repression', 'dephosphorylation']
    physical_contact = ['binding/association', 'dissociation', 'state change']
    chemical = ['compound', 'glycosylation', 'ubiquitination']
    for i in network.edges():
        n = network.get_edge(i[0], i[1])
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
            print(n, n.attr)


def create_legend(graph):
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
        graph.add_edge(n, n + len_dic, dir='both', arrowhead=dict_of_types[i],
                       arrowtail="none", label=i)
    graph.add_subgraph(subgraph, name='cluster_legend', rank="LR")
