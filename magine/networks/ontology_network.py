import itertools
import os

import networkx as nx
import pandas as pd

from magine.networks.exporters import export_to_dot
from magine.networks.visualization.igraph_tools import create_igraph_figure


class OntologyNetworkGenerator(object):
    """ Generates ontology network from molecular networks.

    Nodes are the ontology terms
    Edges between nodes are determined based on molecular network graphs.

    """

    def __init__(self, molecular_network=None):

        self.mol_network = molecular_network
        self._nodes = None
        self._edges = None
        self.molecular_network = None

    @property
    def edges(self):
        if self._edges is None:
            self._edges = set(self.mol_network.edges())
        return self._edges

    @property
    def nodes(self):
        if self._nodes is None:
            self._nodes = set(self.mol_network.nodes())
        return self._nodes

    def _count_neighbors(self, term_a, term_b):
        """
        Calculate the number of direct edges between species of two terms


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

        for i, j in itertools.product(term_1_good, term_2_good):
            if i != j:
                if (i, j) in self.edges:
                    genes_in_go.add(i)
                    genes_in_go.add(j)
                    counter += 1

        return counter, genes_in_go

    def create_network_from_list(self, list_of_ontology_terms,
                                 ont_to_species_dict, ont_to_label_dict,
                                 save_name=None,
                                 draw=False, threshold=0, out_dir=None,
                                 merge_edges=False):
        """
        Creates a GO level network from list of GO terms

        Parameters
        ----------

        list_of_ontology_terms : list_list
            list of GO terms
        ont_to_species_dict : dict
            dictionary, where keys are ontology terms and values are species
            for that term
        ont_to_label_dict: dict
            dictionary where keys are ontology terms and values are labels
        save_name : str
            name to save network
        draw : bool
            create a go_graph of network
        threshold : int
            integer threshold of number of neighbors between two GO terms
            default = 0
        out_dir : str
            output directory
        merge_edges : bool
            merge the edges between GO nodes

        Returns
        -------
        networkx.DiGraph

        """
        if out_dir is not None:
            out_path = os.path.join(out_dir, 'Network_files')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                os.mkdir(out_path)
        else:
            out_path = '.'

        # make sure a network exists
        if self.mol_network is None:
            print("Must provide a network! Returning None")
            return None
        go_graph = nx.DiGraph()
        mol_net = nx.DiGraph()
        all_genes = set()
        list_of_go_terms = set(list_of_ontology_terms)
        sp_to_term = dict()
        sp_to_label = dict()

        def _add(gene_name, term_name, all_term, label):
            if gene_name in all_term:
                if gene_name in sp_to_term:
                    sp_to_term[gene_name].add(term_name)
                    sp_to_label[gene_name].add(label)
                else:
                    sp_to_term[gene_name] = {term_name}
                    sp_to_label[gene_name] = {label}

        for i in itertools.combinations(list_of_go_terms, 2):
            term1 = i[0]
            term2 = i[1]
            term_1 = set(ont_to_species_dict[term1])
            term_2 = set(ont_to_species_dict[term2])
            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]

            a_to_b, b_to_a, genes_in_edges = self._count_neighbors(term_1,
                                                                   term_2)

            # add to graph if at least one edge found between terms
            if a_to_b or b_to_a:
                x = self.mol_network.subgraph(genes_in_edges)
                mol_net.add_edges_from(x.edges(data=True))
                for gene in genes_in_edges:
                    _add(gene, term1, term_1, label_1)
                    _add(gene, term2, term_2, label_2)
                all_genes.update(genes_in_edges)
            else:
                print('No edges between {} and {}'.format(label_1, label_2))
            go_graph.add_node(label_1, term=term1, label=label_1)
            go_graph.add_node(label_2, term=term2, label=label_2)

            if merge_edges:
                if a_to_b > threshold or b_to_a > threshold:
                    if a_to_b > b_to_a:
                        go_graph.add_edge(label_2, label_1,
                                          label=str(a_to_b + b_to_a),
                                          weight=a_to_b + b_to_a, dir='both',
                                          weightAtoB=b_to_a, weightBtoA=a_to_b)
                    else:
                        go_graph.add_edge(label_1, label_2,
                                          label=str(a_to_b + b_to_a),
                                          weight=a_to_b + b_to_a, dir='both',
                                          weightAtoB=a_to_b, weightBtoA=b_to_a)
            else:
                if a_to_b > threshold:
                    go_graph.add_edge(label_1, label_2, label=str(a_to_b),
                                      weight=a_to_b)

                if b_to_a > threshold:
                    go_graph.add_edge(label_2, label_1, label=str(b_to_a),
                                      weight=b_to_a)
        for i in sp_to_term:
            labels = sp_to_label[i]
            terms = sp_to_term[i]
            assert len(labels) == len(terms), \
                'len(labels) should equal len(terms)'
            mol_net.node[i]['termName'] = ','.join(sorted(labels))
            mol_net.node[i]['terms'] = ','.join(sorted(terms))

        self.molecular_network = mol_net
        if save_name is not None:
            if out_dir is not None:
                save_name = os.path.join(out_path, save_name)
            out_name = '{}_subgraph.gml'.format(save_name)
            nx.write_gml(mol_net, out_name)
            nx.write_gml(go_graph, '{}.gml'.format(save_name))

            if draw:
                export_to_dot(go_graph, save_name)
                create_igraph_figure(mol_net, save_name + '_subgraph_igraph')

        return go_graph, mol_net


def visualize_go_network(go_network, data, save_name,
                         format_only=False, out_dir=None, merge=False):
    """ Renders GO network with py2cytoscape

    Parameters
    ----------
    go_network : networkx.DiGraph
        GO network created with GoNetworkGenerator
    data : pd.Dataframe
        enrichment data from GOAnalysis.calculate_enrichment
    save_name: str
        prefix to save images of GO network
    format_only: boolean
        option to return formatted network and labels for rendering
    out_dir: str
        If you want to files to be output into directory
    merge: bool
        merge images into gif

    Returns
    -------

    """

    from magine.networks.visualization.cytoscape_view import RenderModel
    from magine.plotting.heatmaps import heatmap_from_array

    if len(go_network.nodes()) == 0:
        print('No nodes')
        quit()

    labels = data['sample_index'].unique()

    score_array = pd.pivot_table(data, index=['term_name'],
                                 columns='sample_id')

    heatmap_from_array(data, cluster_row=False, convert_to_log=True,
                       index='term_name', values='combined_score',
                       columns='sample_id', div_colors=True)

    x = score_array['combined_score'].fillna(0)

    for i in go_network.nodes():
        values = x.loc[i]
        go_network.node[i]['color'] = 'red'
        go_network.node[i]['label'] = i
        for n, time in enumerate(labels):
            go_network.node[i][time] = float(values[time])

    if format_only:
        return go_network, data

    savename = os.path.join('{}_all_colored.graphml'.format(save_name))

    nx.write_graphml(go_network, savename)
    # size_of_data = len(labels)
    rm = RenderModel(go_network, layout='force-directed')
    rm.visualize_by_list_of_time(labels,
                                 prefix=save_name,
                                 out_dir=out_dir,
                                 )
    if merge:
        os.system(_s.format(save_name))


_s = "convert -delay 100 -dispose previous -loop 0 ont_network_*_formatted.png {}.gif"


def _create_names(n):
    names = []
    for i in range(n):
        names.append('time_{0:04d}'.format(i))
        print('time_{0:04d}'.format(i))
    return names
