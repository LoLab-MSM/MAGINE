import itertools
import os

import networkx as nx
import numpy as np
import pygraphviz as pyg
from orangecontrib.bio import go


class GoNetworkGenerator:
    """ Generates GO networks from molecular networks.

    Nodes are GO terms, edges between nodes are determined based on molecular network graphs.

    """

    def __init__(self, organism='hsa', network=None, directory='tmp'):
        self.ontology = go.Ontology()
        self.annotations = go.Annotations(organism, ontology=self.ontology)
        if network is None:
            print("Please provide networkx network to map GO onto")
        self.network = network
        if network is not None:
            self.nodes = set(self.network.nodes())
        self.go_network = None
        self.molecular_network = None
        self.out_dir = directory
        if os.path.exists(self.out_dir):
            pass
        else:
            os.mkdir(self.out_dir)

    def calculate_terms_between_a_b(self, term_a, term_b):
        """ calculates the number of edges between species in two lists

        :param term_a: proteins list A
        :param term_b: protein list B
        :return: number of edges from A to B, number of edges from B to A, genes responsible for edges
        """
        term_a = set(term_a)
        term_b = set(term_b)
        genes_in_both_go = set()

        path_a_to_b, gene_in = self.determine_edges(term_1=term_a, term_2=term_b)
        for i in gene_in:
            genes_in_both_go.add(i)

        path_b_to_a, gene_in = self.determine_edges(term_1=term_b, term_2=term_a)
        for i in gene_in:
            genes_in_both_go.add(i)

        return path_a_to_b, path_b_to_a, genes_in_both_go

    def determine_edges(self, term_1, term_2):
        """ calculates the number of edges between two terms

        :param term_1: list of proteins
        :param term_2: list of proteins
        :return: count of edges from 1 to 2, genes with edges
        """
        counter = 0
        genes_in_go = set()
        for i in term_1:
            if i in self.nodes:
                neigh = set(self.network.neighbors(i))
                for j in neigh:
                    if j in term_2:
                        if i == j:
                            continue
                        if self.network.has_edge(i, j):
                            genes_in_go.add(i)
                            genes_in_go.add(j)
                            counter += 1
        return counter, genes_in_go

    # @profile
    def create_network_from_list(self, network=None, list_of_go_terms=None, save_name=None, draw=False, threshold=0):
        """ creates a network from lists of GO terms

        :param network: networkx graph of molecular terms, species must be same species as ontology in init
        :param list_of_go_terms:
        :param save_name:
        :param draw: boolean, generates graphviz images of network, can be time consuming
        :param threshold: int, requirment for minimum number of edges between two GO terms to be considered valid
        :return:
        """
        if network is not None:
            self.network = network
            self.nodes = set(self.network.nodes())
        list_of_go_terms = np.array(np.unique(list_of_go_terms))
        graph = nx.DiGraph()
        molecular_network_subgraph = nx.DiGraph()
        gene_annotations_dict = dict()
        go_names = dict()
        all_genes = set()
        for i in list_of_go_terms:
            gene_annotations_dict[i] = self.annotations.get_all_genes(i)
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

            a_to_b, b_to_a, genes_in_edges = self.calculate_terms_between_a_b(term_1, term_2)
            if a_to_b or b_to_a:
                x = self.network.subgraph(genes_in_edges)
                # molecular_network_subgraph.add_edges_from(x.edges())
                for i, j, s in x.edges(data=True):
                    molecular_network_subgraph.add_edge(i, j, s)
                for gene in genes_in_edges:
                    if gene in all_genes:
                        if gene in term_1:
                            molecular_network_subgraph.node[gene]['go'] += ',' + go_1
                        if gene in term_2:
                            molecular_network_subgraph.node[gene]['go'] += ',' + go_2
                    else:
                        if gene in term_1:
                            molecular_network_subgraph.node[gene]['go'] = go_1
                        if gene in term_2:
                            molecular_network_subgraph.node[gene]['go'] = go_2
                all_genes.add(g for g in genes_in_edges)
            if a_to_b > threshold or b_to_a > threshold :
                graph.add_node(label_1, go=go_1, label=label_1)
                graph.add_node(label_2, go=go_2, label=label_2)
                graph.add_edge(label_1, label_2, label=str(a_to_b+b_to_a), weight=a_to_b+b_to_a)
            else:
                print(label_1,label_2)

            # if a_to_b > threshold:
            #     graph.add_node(label_1, go=go_1, label=label_1)
            #     graph.add_node(label_2, go=go_2, label=label_2)
            #     graph.add_edge(label_1, label_2, label=str(a_to_b), weight=a_to_b)
            #
            # if b_to_a > threshold:
            #     graph.add_node(label_1, go=go_1, label=label_1)
            #     graph.add_node(label_2, go=go_2, label=label_2)
            #     graph.add_edge(label_2, label_1, label=str(b_to_a), weight=b_to_a)

        self.molecular_network = molecular_network_subgraph
        nx.write_gml(molecular_network_subgraph,
                     os.path.join(self.out_dir, '{0}_subgraph.gml'.format(save_name)))

        nx.nx.write_dot(graph,
                        os.path.join(self.out_dir, '{0}.dot'.format(save_name)))
        nx.write_graphml(graph,
                         os.path.join(self.out_dir, '{0}.graphml'.format(save_name)))

        if draw:
            g = pyg.AGraph()
            g.read(os.path.join(self.out_dir, '{0}.dot'.format(save_name)))
            g.draw(os.path.join(self.out_dir, '{0}.pdf'.format(save_name)), prog='dot')
        return graph


if __name__ == '__main__':
    ddn = nx.read_gml('/home/pinojc/git/Network_projects/Cisplatin_project/Network_files/ddn3.gml')
    test_list = ['GO:0008219', 'GO:0006281', 'GO:0008283', 'GO:0043066', 'GO:0043065', 'GO:1902175', 'GO:0006805',
                 'GO:0006766']
    # test_list = ['GO:1902175', 'GO:0006805', 'GO:0006766', 'GO:0015893', 'GO:0006936']
    gnc = GoNetworkGenerator('hsa', ddn)

    gnc.create_network_from_list(ddn, test_list, 'xeno', draw=False)
    # create_go_network,save_name='death_dnaRepair_proliferation')

    # create_go_network('GO:0043066', 'GO:0043065', 'GO:1902175', 'apoptosis')
