import itertools
import os

import networkx as nx
from orangecontrib.bio import go

from magine.network_tools import export_to_dot


class GoNetworkGenerator:
    """ Generates GO networks from molecular networks.

    Nodes are GO terms.
    Edges between nodes are determined based on molecular network graphs.

    """

    def __init__(self, organism='hsa', network=None, directory='tmp'):
        self.ontology = go.Ontology()
        self.annotations = go.Annotations(organism, ontology=self.ontology)
        if network is None:
            print("Please provide networkx network to map GO onto")
        self.network = network
        if network is not None:
            self.nodes = set(self.network.nodes())
            self.edges = set(self.network.edges())
        self.go_network = None
        self.molecular_network = None
        self.out_dir = directory
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
            os.mkdir(os.path.join(self.out_dir, 'Network_files'))

    def count_neighbors(self, term_a, term_b):
        """
        Calculate the number of direct edges between two GO terms species


        Parameters
        ----------
        term_a : list_like
            list of species
        term_b : list_like
            list of species

        Returns
        -------
        int, int, list_like
            number of edges from A to B
            number of edges from B to A
            genes responsible for edges
        """
        term_a = set(term_a)
        term_b = set(term_b)
        genes_in_both_go = set()
        # calculate edges between A and B
        a_to_b, gene_in = self._determine_edges(term_1=term_a, term_2=term_b)
        genes_in_both_go.update(gene_in)

        # calculate edges between B and A
        b_to_a, gene_in = self._determine_edges(term_1=term_b, term_2=term_a)
        genes_in_both_go.update(gene_in)

        return a_to_b, b_to_a, genes_in_both_go

    def _determine_edges(self, term_1, term_2):
        """
        calculate the number of neighbors that connect between two terms

        Parameters
        ----------
        term_1
        term_2

        Returns
        -------

        """
        counter = 0
        genes_in_go = set()
        term_1_good = {i for i in term_1 if i in self.nodes}
        term_2_good = {i for i in term_2 if i in self.nodes}

        for i in term_1_good:
            for j in term_2_good:
                if i == j:
                    continue
                if (i, j) in self.edges:
                    genes_in_go.add(i)
                    genes_in_go.add(j)
                    counter += 1
        return counter, genes_in_go

    # @profile
    def create_network_from_list(self, network=None, list_of_go_terms=None,
                                 save_name=None, draw=False, threshold=0):
        """
        Creates a GO level network from list of GO terms

        Parameters
        ----------
        network : nx.DiGraph
            molecular level network
        list_of_go_terms : list_list
            list of GO terms
        save_name : str
            name to save network
        draw : bool
            create a go_graph of network
        threshold : int
            integer threshold of number of neighbors between two GO terms
            default = 0

        Returns
        -------
        networkx.DiGraph

        """
        # if provided a new network, nodes and edges will be used from here
        if network is not None:
            self.network = network
            self.nodes = set(self.network.nodes())
            self.edges = set(self.network.edges())

        go_graph = nx.DiGraph()
        molecular_network = nx.DiGraph()
        gene_annotations_dict = dict()
        go_names = dict()
        all_genes = set()
        list_of_go_terms = set(list_of_go_terms)
        for i in list_of_go_terms:
            gene_annotations_dict[i] = set(self.annotations.get_all_genes(i))
            # go_names[i] = "\n".join(wrap(self.ontology[i].name, 20))
            go_names[i] = self.ontology[i].name

        for i in itertools.combinations(list_of_go_terms, 2):
            term1 = i[0]
            term2 = i[1]
            term_1 = set(gene_annotations_dict[term1])
            term_2 = set(gene_annotations_dict[term2])
            label_1 = go_names[term1]
            label_2 = go_names[term2]
            go_1 = str(term1).replace(':', '')
            go_2 = str(term2).replace(':', '')

            a_to_b, b_to_a, genes_in_edges = self.count_neighbors(term_1,
                                                                  term_2)

            # add to graph if at least one edge found between terms
            if a_to_b or b_to_a:
                x = self.network.subgraph(genes_in_edges)
                molecular_network.add_edges_from(x.edges(data=True))
                for gene in genes_in_edges:
                    if gene in all_genes:
                        if gene in term_1:
                            names = molecular_network.node[gene]['go'].split(
                                ',')
                            if go_1 not in names:
                                molecular_network.node[gene][
                                    'go'] += ',' + go_1
                                molecular_network.node[gene][
                                    'goName'] += ',' + label_1

                        if gene in term_2:
                            names = molecular_network.node[gene]['go'].split(
                                ',')
                            if go_2 not in names:
                                molecular_network.node[gene][
                                    'go'] += ',' + go_2
                                molecular_network.node[gene][
                                    'goName'] += ',' + label_2
                    else:
                        if gene in term_1:
                            molecular_network.node[gene]['go'] = go_1
                            molecular_network.node[gene]['goName'] = label_1
                        if gene in term_2:
                            molecular_network.node[gene]['go'] = go_2
                            molecular_network.node[gene]['goName'] = label_2
                for g in genes_in_edges:
                    all_genes.add(g)
            else:
                print('No edges between {} and {}'.format(label_1, label_2))

            # if a_to_b > threshold or b_to_a > threshold:
            #     go_graph.add_node(label_1, go=go_1, label=label_1)
            #     go_graph.add_node(label_2, go=go_2, label=label_2)
            #     if a_to_b > b_to_a:
            #         go_graph.add_edge(label_2, label_1, label=str(a_to_b + b_to_a), weight=a_to_b + b_to_a,
            #                        weightAtoB=b_to_a, weightBtoA=a_to_b)
            #     else:
            #         go_graph.add_edge(label_1, label_2, label=str(a_to_b + b_to_a), weight=a_to_b + b_to_a,
            #                        weightAtoB=a_to_b, weightBtoA=b_to_a)
            # else:
            #     pass

            if a_to_b > threshold:
                go_graph.add_node(label_1, go=go_1, label=label_1, )
                go_graph.add_node(label_2, go=go_2, label=label_2, )
                go_graph.add_edge(label_1, label_2, label=str(a_to_b),
                                  weight=a_to_b)

            if b_to_a > threshold:
                go_graph.add_node(label_1, go=go_1, label=label_1, )
                go_graph.add_node(label_2, go=go_2, label=label_2, )
                go_graph.add_edge(label_2, label_1, label=str(b_to_a),
                                  weight=b_to_a)

        self.molecular_network = molecular_network

        nx.write_gml(molecular_network,
                     os.path.join(self.out_dir, 'Network_files',
                                  '{0}_subgraph.gml'.format(save_name)))

        nx.write_graphml(go_graph,
                         os.path.join(self.out_dir, 'Network_files',
                                      '{0}.graphml'.format(save_name)))

        if draw:
            s_name = os.path.join(self.out_dir, 'Network_files',
                                  '{0}'.format(save_name))
            export_to_dot(go_graph, s_name)
            export_to_dot(molecular_network,
                          '{}_molecular_network'.format(save_name))
        return go_graph


if __name__ == '__main__':
#    ddn = nx.read_gml(
#        'C:\Users\James Pino\PycharmProjects\Network_projects\Cisplatin_project\Network_files\ddn3.gml')
    # test_list = ['GO:0008219', 'GO:0006281', 'GO:0008283', 'GO:0043066', 'GO:0043065', 'GO:1902175', 'GO:0006805',
    #              'GO:0006766']
    test_list = ['GO:1902175', 'GO:0006805', 'GO:0006766', 'GO:0015893',
                 'GO:0006936']

    gnc = GoNetworkGenerator('hsa', ddn)
    # import cProfile
    # cProfile.run("gnc.create_network_from_list(ddn, test_list, 'xeno', draw=False)", sort=1)
    gnc.create_network_from_list(ddn, test_list, 'xeno', draw=False)
    # create_go_network,save_name='death_dnaRepair_proliferation')

    # create_go_network('GO:0043066', 'GO:0043065', 'GO:1902175', 'apoptosis')
