import itertools

import networkx as nx

import magine.network_tools as nt


class NetworkSubgraphs(object):
    def __init__(self, network, exp_data=None):
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
        graph = nx.DiGraph()
        direction_1, direction_2 = True, True
        for i in (node_1, node_2):
            if i not in self.nodes:
                print("{} is not in network!")
                return

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
            raise RuntimeWarning("No paths between {} and {}. Returning "
                                 "None".format(node_1, node_2))
            return None

        if draw is not None:
            save_name = "%s_and_%s.pdf" % (node_1, node_2)
            self._save_or_draw(graph, save_name, draw)

        return graph

    def shortest_paths_between_lists(self, species_list, save_name=None,
                                     single_path=False, draw=False,
                                     image_format='png'):
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
        >>> path_a_d = net_sub.shortest_paths_between_lists(['a','c','d'])
        Path does NOT exist between c and a
        Path does NOT exist between d and a
        Path does NOT exist between d and c
        >>> path_a_d.edges()
        [('a', 'b'), ('a', 'd'), ('c', 'd'), ('b', 'c')]
        """

        graph = nx.DiGraph()

        tmp_protein_list = self._check_node(species_list)
        node_pairs = itertools.combinations(tmp_protein_list, 2)
        warning_message = "Path does NOT exist between {} and {}"
        for node_1, node_2 in node_pairs:
            if not nx.has_path(self.network, node_1, node_2):
                print(warning_message.format(node_1, node_2))
            else:
                if single_path:
                    self._add_edges_from_path(
                            graph,
                            nx.shortest_path(self.network, node_1, node_2)
                    )
                else:
                    for path in nx.all_shortest_paths(self.network, node_1,
                                                      node_2):
                        self._add_edges_from_path(graph, path)
            if not nx.has_path(self.network, node_2, node_1):
                print(warning_message.format(node_2, node_1))
            else:
                if single_path:
                    self._add_edges_from_path(
                            graph,
                            nx.shortest_path(self.network, node_2, node_1)
                    )
                else:
                    for path in nx.all_shortest_paths(self.network, node_2,
                                                      node_1):
                        self._add_edges_from_path(graph, path)

        if save_name is not None:
            self._save_or_draw(graph, save_name, draw, image_format)
        return graph

    def upstream_network_of_specie(self, species_1, include_list=None,
                                   save_name=None, compress=False,
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
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'), \
        ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
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

        graph = nx.DiGraph()
        for i in self.nodes:
            if i in include_list:
                if nx.has_path(self.network, i, species_1):
                    path = nx.shortest_path(self.network, i, species_1)
                    self._add_edges_from_path(graph, path)

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


        Example
        -------
        >>> from networkx import DiGraph
        >>> from magine.networks.network_subgraphs import NetworkSubgraphs
        >>> g = DiGraph()
        >>> g.add_edges_from([('a','b'),('b','c'), ('c', 'd'), ('a', 'd'),\
         ('e', 'd')])
        >>> net_sub = NetworkSubgraphs(g)
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

        graph = nx.DiGraph()
        for i in self.nodes:
            if i in include_list:
                if nx.has_path(self.network, species_1, i):
                    self._add_edges_from_path(
                            graph,
                            nx.shortest_path(self.network, species_1, i)
                    )

        if compress:
            graph = nt.compress_edges(graph)

        if save_name is not None:
            self._save_or_draw(graph, save_name, draw)

        return graph

    def measured_networks_over_time(self, graph, colors, prefix):
        """ Adds color to a network over time
        
        Parameters
        ----------
        graph : pygraphviz.AGraph
        colors : list
            List of colors for time points
        prefix : str
            Prefix for image files
            

        Returns
        -------

        """
        nt.paint_network_overtime(graph, self.exp_data.sig_species_over_time,
                                  colors, prefix,
                                  self.exp_data.protomics_time_points)

    def measured_networks_over_time_up_down(self, graph, prefix,
                                            color_up='tomato',
                                            color_down='lightblue'):
        """

        Parameters
        ----------
        graph : pygraphviz.AGraph

        prefix : str
            Prefix for image files
        color_up : str
        
        color_down : str


        Returns
        -------

        """
        up_species = self.exp_data.sig_species_up_over_time
        down_species = self.exp_data.sig_species_down_over_time
        labels = self.exp_data.protomics_time_points
        nt.paint_network_overtime_up_down(graph, list_up=up_species,
                                          list_down=down_species,
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
                  "Removing from list".format(missing_nodes))
            print(len(node_list))
            node_list.difference_update(missing_nodes)
            print(len(node_list))
            return node_list
        else:
            return node_list
    
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

    @staticmethod
    def _save_or_draw(graph, save_name, draw, img_format):
        nx.write_gml(graph, "{}.gml".format(save_name))
        if draw:
            graph = nt._format_to_directions(graph)
            nt.export_to_dot(graph, save_name=save_name,
                             image_format=img_format)
