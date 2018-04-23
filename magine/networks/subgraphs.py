import itertools
from functools import partial

import networkx as nx

import magine.networks.dev_tools as nt
import magine.networks.exporters
import magine.networks.utils


class Subgraph(object):
    def __init__(self, network, exp_data=None, pool=None):
        """
        Generates network subgraphs

        Parameters
        ----------
        network : networkx.DiGraph
        exp_data : magine.data.datatypes.ExperimentalData
        
        """
        self.network = network
        self.nodes = set(self.network.nodes())
        self.exp_data = exp_data
        self._paths = None

        if pool is None:
            self.map = map
        else:
            self.map = pool.map

    def find_paths(self, source, target):
        try:
            return self._paths[source][target]
        except:
            return None

    def shortest_paths_between_two_proteins(self, node_1, node_2,
                                            bidirectional=False,
                                            single_path=False,
                                            draw=False, image_format='png'):
        """
        Generates a graph based on all shortest paths between two species


        Parameters
        ----------
        node_1 : str
            name of first species
        node_2 : str
            name of second species
        bidirectional : bool
            If you want to search bidirectionally
        single_path : bool
            If you only want a single shortest path
        draw : bool
            create an image of returned network
        image_format : str, optional
            If draw=True you can pass an image format. (pdf, png, svg). 
            default=png

        Returns
        -------
        graph : networkx.DiGraph


        Examples
        --------
        >>> from networkx import DiGraph
        >>> from magine.networks.subgraphs import Subgraph
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = Subgraph(g)
        >>> path_a_d = net_sub.shortest_paths_between_two_proteins('a','d')
        >>> path_a_d.edges()
        [('a', 'd')]
        >>> path_a_e = net_sub.shortest_paths_between_two_proteins('a','e')

        """

        for i in (node_1, node_2):
            if i not in self.nodes:
                print("{} is not in network!")
                return

        _paths = list()
        _paths.append(_nx_find_path(self.network, node_1, node_2,
                                    single_path=single_path))
        if bidirectional:
            _paths.append(_nx_find_path(self.network, node_2, node_1,
                                        single_path=single_path))
        graph = self._list_paths_to_graph(_paths)
        if draw:
            save_name = "%s_and_%s" % (node_1, node_2)
            self._save_or_draw(graph, save_name, draw, image_format)

        return graph

    def shortest_paths_between_lists(self, species_list, save_name=None,
                                     single_path=False, draw=False, pool=None,
                                     image_format='png', max_length=None,
                                     include_only=None):
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
        image_format : str
            dot acceptable output formats, (pdf, png, etc)
        pool : multiprocessing.Pool
            If it it provided, it uses its map function to run this function.
        max_length : int
            Max length for path between any 2 species
        include_only : list_like
            List of species that must be present
        Returns
        -------
        graph : networkx.DiGraph
            graph containing paths between species list provided


        Examples
        --------
        >>> from networkx import DiGraph
        >>> from magine.networks.subgraphs import Subgraph
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), ('e', 'd')])
        >>> net_sub = Subgraph(g)
        >>> path_a_d = net_sub.shortest_paths_between_lists(['a','c','d'])
        >>> path_a_d.edges()
        [('a', 'b'), ('a', 'd'), ('c', 'd'), ('b', 'c')]
        """

        tmp_species_list = list(self._check_node(species_list))

        if pool is not None:
            paths = pool.map(
                partial(_find_nx_path,
                        network=self.network,
                        single_path=single_path),
                itertools.combinations(tmp_species_list, 2)
            )
            pool.close()
        else:
            paths = map(partial(_find_nx_path,
                                network=self.network,
                                single_path=single_path),
                        itertools.combinations(tmp_species_list, 2))

        if max_length is not None:
            new_paths = []
            for p in paths:
                local_p = []
                for i in p:
                    if len(i) <= max_length:
                        local_p.append(i)
                new_paths.append(local_p)
            paths = new_paths[:]

        graph = self._list_paths_to_graph(paths)
        if include_only is not None:
            graph = self._include_only(graph, include_only)
        if save_name is not None:
            self._save_or_draw(graph, save_name, draw, image_format)
        return graph

    def neighbors(self, node, up=True, down=True, max_dist=1,
                  include_only=None):
        if max_dist > 3:
            print("Max distance is 3. Big networks")
            max_dist = 3
        sg = nx.DiGraph()

        if node not in self.nodes:
            return sg

        def _add_node(n):
            sg.add_node(n, **self.network.node[n])

        def _add_edge(i, j):
            sg.add_edge(i, j, **self.network[i][j])

        _add_node(node)

        def _get_upstream(new_node):
            upstream = self.network.predecessors(new_node)
            for i in upstream:
                if i != new_node:
                    _add_node(i)
                    _add_edge(i, new_node)
            return upstream

        def _get_downstream(new_node):
            downstream = self.network.successors(new_node)
            for i in downstream:
                if i != new_node:
                    _add_node(i)
                    _add_edge(new_node, i)
            return downstream

        up_layer = [node]
        down_layer = [node]
        for _ in range(0, max_dist):
            if up:
                up_layer = {i for n in up_layer for i in _get_upstream(n)}
            if down:
                down_layer = {i for n in down_layer for i in
                              _get_downstream(n)}
        if include_only is not None:
            sg = self._include_only(sg, include_only)

        return sg

    def neighbors_of_list(self, list_start, up_stream=True, down_stream=True,
                          max_dist=1, include_only=None):
        new_g = nx.DiGraph()
        for start in list_start:
            sg = self.neighbors(start, up_stream, down_stream, max_dist)
            new_g = nx.compose(new_g, sg)

        if include_only is not None:
            new_g = self._include_only(new_g, include_only)

        return new_g

    def upstream_network_of_specie(self, species_1, include_list=None,
                                   save_name=None, compress=False,
                                   draw=False):
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


        Examples
        --------
        >>> from networkx import DiGraph
        >>> from magine.networks.subgraphs import Subgraph
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), \
        ('e', 'd')])
        >>> net_sub = Subgraph(g)
        >>> upstream_d = net_sub.upstream_network_of_specie('d')
        >>> upstream_d.edges()
        [('a', 'd'), ('c', 'd'), ('b', 'c'), ('e', 'd')]
        >>> upstream_c = net_sub.upstream_network_of_specie('c')
        >>> upstream_c.edges()
        [('a', 'b'), ('b', 'c')]

        """
        if include_list is None:
            include_list = set(self.nodes.copy())
        else:
            include_list = self._check_node(include_list)

        graph = self._list_paths_to_graph(
            [_nx_find_path(self.network, i, species_1, single_path=True)
             for i in self.nodes if i in include_list])

        if compress:
            graph = nt.compress_edges(graph)
        if save_name is not None:
            self._save_or_draw(graph, save_name, draw)

        return graph

    def downstream_network_of_specie(self, species_1, include_list=None,
                                     save_name=None, compress=False,
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


        Examples
        --------
        >>> from networkx import DiGraph
        >>> from magine.networks.subgraphs import Subgraph
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'),\
         ('e', 'd')])
        >>> net_sub = Subgraph(g)
        >>> downstream_d = net_sub.downstream_network_of_specie('d')
        >>> downstream_d.edges()
        []
        >>> downstream_c = net_sub.downstream_network_of_specie('c')
        >>> downstream_c.edges()
        [('c', 'd')]

        """
        if include_list is None:
            include_list = set(self.nodes.copy())
        else:
            include_list = self._check_node(include_list)

        graph = self._list_paths_to_graph(
            [_nx_find_path(self.network, species_1, i, single_path=True)
             for i in self.nodes if i in include_list])

        if compress:
            graph = nt.compress_edges(graph)

        if save_name is not None:
            self._save_or_draw(graph, save_name, draw)

        return graph

    def expand_neighbors(self, network, nodes=None, up_stream=False,
                         down_stream=False, include_list=None,
                         expand_all=False):
        if include_list is None:
            include_list = set(self.nodes.copy())
        else:
            include_list = self._check_node(include_list)
        if expand_all:
            nodes = set(network.nodes())
        elif isinstance(nodes, str):
            nodes = [nodes]
        elif not isinstance(nodes, (list, set)):
            print("Must provide node, list of nodes, or expand_all=True")
            return network

        new_net = network.copy()

        for i in nodes:
            if up_stream:
                upstream = self.network.predecessors(i)
                up_nodes = [g for g in upstream if g in include_list]
                for n in up_nodes:
                    new_net.add_edge(n, i, **self.network[n][i])

            if down_stream:
                downstream = self.network.successors(i)
                up_nodes = [g for g in downstream if g in include_list]
                for n in up_nodes:
                    new_net.add_edge(i, n, **self.network[i][n])
        return new_net

    def measured_networks_over_time(self, graph, colors, prefix):
        """ Adds color to a network over time
        
        Parameters
        ----------
        graph : nx.DiGraph
        colors : list
            List of colors for time points
        prefix : str
            Prefix for image files
            

        Returns
        -------

        """
        measured_list = []
        for i, j in sorted(self.exp_data.sig_species_over_time.items()):
            measured_list.append(j)
        magine.networks.utils.paint_network_overtime(graph, measured_list,
                                                     colors, prefix,
                                                     self.exp_data.proteomics_sample_ids)

    def measured_networks_over_time_up_down(self, graph, prefix,
                                            color_up='tomato',
                                            color_down='lightblue'):
        """

        Parameters
        ----------
        graph : nx.DiGraph

        prefix : str
            Prefix for image files
        color_up : str
        
        color_down : str


        Returns
        -------

        """
        labels = self.exp_data.proteomics_sample_ids
        up_measured_list = []
        down_measured_list = []
        for i, j in sorted(self.exp_data.sig_species_over_time.items()):
            up_measured_list.append(j)

        for i, j in sorted(self.exp_data.sig_species_down_over_time.items()):
            down_measured_list.append(j)
        magine.networks.utils.paint_network_overtime_up_down(graph,
                                                             list_up=up_measured_list,
                                                             list_down=down_measured_list,
                                                             save_name=prefix,
                                                             color_down=color_down,
                                                             color_up=color_up,
                                                             labels=labels)

    def _check_node(self, node_list):
        """
        Checks to see if list of nodes is in the graph, removes them if not.
        Parameters
        ----------
        node_list : list_like

        Returns
        -------
        set
        """
        node_list = set(node_list)
        missing_nodes = set()
        for i in node_list:
            if i not in self.nodes:
                missing_nodes.add(i)
        if len(missing_nodes) != 0:
            print("Warning : {} do not exist in graph\n"
                  "Removing from list".format(len(missing_nodes)))
            node_list.difference_update(missing_nodes)
        return sorted(node_list)

    def _list_paths_to_graph(self, paths):
        graph = nx.DiGraph()
        current_nodes = set(graph.nodes())
        _edges = set()

        def _add_node(node):
            if node not in current_nodes:
                graph.add_node(node, **self.network.node[node])
                current_nodes.add(node)

        def _add_edge(i, j):
            if (i, j) not in _edges:
                graph.add_edge(i, j, **self.network.edge[i][j])
                _edges.add((i, j))

        # flatten paths into list of lists
        paths = [n for x in paths for n in x if n]

        for path in paths:
            previous = path[0]
            _add_node(previous)
            for protein in path[1:]:
                _add_node(protein)
                _add_edge(previous, protein)
                previous = protein
        if len(graph.nodes()) == 0:
            return None
        return graph

    @staticmethod
    def _include_only(network, include_list):
        assert isinstance(include_list, (list, set))
        sg = network.copy()
        all_nodes = set(sg.nodes())
        not_found = all_nodes.difference(set(include_list))
        sg.remove_nodes_from(not_found)
        if len(sg.nodes()) == 0:
            print("Warning: no nodes were found in include_only list! "
                  "Network doesn't contain any nodes!")
        # nt.delete_disconnected_network(sg)
        return sg

    @staticmethod
    def _save_or_draw(graph, save_name, draw, img_format='png'):
        nx.write_gml(graph, "{}.gml".format(save_name))
        if draw:
            graph = magine.networks.exporters.format_to_directions(graph)
            magine.networks.exporters.export_to_dot(graph, save_name=save_name,
                                                    image_format=img_format)


def _find_nx_path(node, network, single_path):
    node1, node2 = node

    l1 = _nx_find_path(network, node1, node2, single_path)
    l2 = _nx_find_path(network, node2, node1, single_path)
    return l1 + l2


def _nx_find_path(network, node1, node2, single_path=False):
    try:
        if single_path:
            return [nx.shortest_path(network, node1, node2)]
        else:
            return [p for p in nx.all_shortest_paths(network, node1, node2)]

    except nx.NetworkXNoPath:
        return []


if __name__ == '__main__':
    net = nx.DiGraph()
    net.add_node('X', label='XX', db='test1')
    net.add_node('Y', label='YY', db='test2')
    net.add_edge('X', 'B', intType='here')
    net.add_edge('X', 'C', intType='here')
    net.add_edge('B', 'Y', intType='here')
    net.add_edge('C', 'Y', intType='here')
    net.add_edge('D', 'Y', intType='here')
    net.add_edge('R', 'D', intType='here')
    net.add_edge('X', 'R', intType='here')

    net = nx.read_gpickle('background_network.p.gz')

    net = net.subgraph(
        net.neighbors('BAX') + ['BAX'] + net.neighbors('TP53') + ['TP53'])
    # print(len(net.edges()))
    # quit()
    # print(x.neighbors('X').nodes())
    # print(x.downstream_network_of_specie('D').nodes())

    # print(x.path_between_2('X', ['Y', 'B', 'C']))
    # quit()

    # print(x.shortest_paths_between_two_proteins('X', 'Y', bidirectional=True, single_path=False).edges())
    x = Subgraph(net)
    import sys

    print(sys.getsizeof(x._paths))
    print(x._paths['BAX'])
    quit()
    print(x.find_paths('X', 'Y'))
    print(x.find_paths('Y', 'X'))

    g = x.shortest_paths_between_lists(['X', 'Y', 'B', 'D'], single_path=True,
                                       draw=False, save_name='test')

    print(g.nodes())
    # g = x.shortest_paths_between_lists(['X', 'Y'], single_path=False,  draw=True, save_name='test')
    # print(g.nodes(data=True))

