import itertools
import os
from itertools import combinations, product

import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection

import magine.networks.utils as nt
from magine.networks.exporters import export_to_dot
from magine.networks.visualization.igraph_tools import render_igraph


class OntologyNetworkGenerator(object):
    """ Generates ontology network from molecular networks.

    Nodes are the ontology terms
    Edges between nodes are determined based on molecular network graphs.

    """

    def __init__(self, molecular_network=None):
        """

        Parameters
        ----------
        molecular_network : nx.DiGraph
        """

        self.network = molecular_network
        self._nodes = None
        self._edges = None
        self.molecular_network = None

    @property
    def edges(self):
        if self._edges is None:
            self.network.remove_edges_from(nx.selfloop_edges(self.network))
            self._edges = set(self.network.edges())
        return self._edges

    @property
    def nodes(self):
        if self._nodes is None:
            self._nodes = set(self.network.nodes())
        return self._nodes

    def create_network_from_list(self, list_of_ontology_terms,
                                 ont_to_species_dict, ont_to_label_dict,
                                 save_name=None, draw=False, out_dir=None,
                                 use_threshold=True, use_fdr=False):
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
        use_threshold : bool
            Use binomial test to check for edge significance
            Uses bh correction for multiple hypothesis testing
        out_dir : str
            output directory
        use_fdr : bool
            Use FDR correction for edge

        Returns
        -------
        go_graph : networkx.DiGraph
            Ontology level network
        mol_net : networkx.DiGraph
            Sub-network made my molecular level species

        """
        if out_dir is not None:
            out_path = os.path.join(out_dir, 'Network_files')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                os.mkdir(out_path)
        else:
            out_path = '.'

        # make sure a network exists
        if self.network is None:
            print("Must provide a network! Returning None")
            return None

        list_of_go_terms = set(list_of_ontology_terms)

        n_edges = float(len(self.edges))
        p_values = []
        # create dictionaries that add the label and terms as node attributes
        gene_to_term, gene_to_label = dict(), dict()

        def add_to_dict(genes, term):

            for gene in genes:
                if gene in gene_to_term:
                    gene_to_term[gene].add(term)
                    gene_to_label[gene].add(ont_to_label_dict[term])
                else:
                    gene_to_term[gene] = {term}
                    gene_to_label[gene] = {ont_to_label_dict[term]}

        def get_edges(set1, set2):
            possible_edges = set(product(set1, set2.difference(set1)))

            n_possible = len(possible_edges)
            edge_hits = possible_edges.intersection(self.edges)
            n_hits = len(edge_hits)
            odds_to_find = float(n_possible) / n_edges
            p_val = stats.binom_test(n_hits, n_possible, odds_to_find)
            return [edge_hits, n_hits, p_val]

        for term1, term2 in combinations(list_of_go_terms, 2):
            term_1 = set(ont_to_species_dict[term1]).intersection(self.nodes)
            term_2 = set(ont_to_species_dict[term2]).intersection(self.nodes)

            add_to_dict(term_1, term1)
            add_to_dict(term_2, term2)

            p_values.append([term1, term2] + get_edges(term_1, term_2))
            p_values.append([term2, term1] + get_edges(term_2, term_1))

        cols = ['term1', 'term2', 'edges', 'n_edges', 'p_values', ]

        df = pd.DataFrame(p_values, columns=cols)

        # FDR correction
        _, df['adj_p_values'] = fdrcorrection(df['p_values'])
        df = df.loc[df['n_edges'] > 0]

        if use_threshold:
            if use_fdr:
                df = df.loc[df['adj_p_values'] <= .05]
            else:
                df = df.loc[df['p_values'] <= .05]
        cols += ['adj_p_values']
        # create empty networks
        go_graph = nx.DiGraph()
        mol_net = nx.DiGraph()
        for term1, term2, edge, n_edges, p_value, adj_pval in df[cols].values:
            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]
            go_graph.add_node(label_1, term=term1, label=label_1)
            go_graph.add_node(label_2, term=term2, label=label_2)
            go_graph.add_edge(label_1, label_2, label=str(n_edges),
                              weight=n_edges, pvalue=p_value, adjPval=adj_pval)

            nodes = list(itertools.chain(*edge))

            nodes += ont_to_species_dict[term1]
            nodes += ont_to_species_dict[term2]

            x = self.network.subgraph(nodes).copy()

            mol_net.add_edges_from(list(x.edges(data=True)))
        exp_net = self.network.subgraph(mol_net.nodes)

        for i in mol_net.nodes:
            labels = gene_to_label[i]
            terms = gene_to_term[i]
            assert len(labels) == len(terms), \
                'len(labels) should equal len(terms)'
            mol_net.node[i]['termName'] = ','.join(sorted(labels))
            mol_net.node[i]['terms'] = ','.join(sorted(terms))
            exp_net.node[i]['termName'] = ','.join(sorted(labels))
            exp_net.node[i]['terms'] = ','.join(sorted(terms))

        self.molecular_network = mol_net
        if save_name is not None:
            if out_dir is not None:
                save_name = os.path.join(out_path, save_name)
            out_name = '{}_subgraph.gml'.format(save_name)
            nx.write_gml(mol_net, out_name)
            nx.write_gml(go_graph, '{}.gml'.format(save_name))

            if draw:
                export_to_dot(go_graph, save_name)
                render_igraph(mol_net, save_name + '_subgraph_igraph')
        return go_graph, exp_net


def create_subnetwork(df, network, terms=None, save_name=None, draw_png=False,
                      remove_isolated=False, create_only=True, merge=False,
                      out_dir=None, use_threshold=False, use_fdr=False):
    """

    Parameters
    ----------

    df : pd.DataFrame
    network : nx.DiGraph
    terms : list, optional
        List of terms to use in ont network. Default is all terms
    save_name : str
    draw_png : bool
    remove_isolated : bool
        Remove nodes that are not connected in the final graphs
    create_only : bool
        Create ontology network only, don't visualize with cytoscape
    merge : bool
        Create gif of cytoscape session
    out_dir : str
        Output directory of images
    use_threshold : bool
        Use binomial test to calculate edge relevance between terms
    use_fdr : bool
        Use BH correction on binomial test


    Returns
    -------

    """
    if terms is not None:
        df_copy = df[df['term_name'].isin(terms)].copy()
        terms = set(terms)
    else:
        df_copy = df.copy()
        terms = set(df['term_name'].values)

    # normalize enriched scores
    df_copy['combined_score'] = np.abs(df_copy['combined_score'])
    df_copy['combined_score'] = np.log2(df_copy['combined_score'])
    df_copy.loc[df['combined_score'] > 150, 'combined_score'] = 150

    labels = df_copy['sample_id'].unique()
    # create dictionary of values
    label_dict, term_dict = dict(), dict()
    for i in terms:
        genes = set(df_copy.term_to_genes(i))
        term_dict[i] = genes
        label_dict[i] = i

    ong = OntologyNetworkGenerator(molecular_network=network)

    print("Creating ontology network")
    term_g, molecular_g = ong.create_network_from_list(
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png,
        use_threshold=use_threshold, use_fdr=use_fdr
    )

    if remove_isolated:
        nt.remove_isolated_nodes(term_g)
        nt.remove_isolated_nodes(molecular_g)

    score_array = pd.pivot_table(df_copy, index=['term_name'],
                                 columns='sample_id')
    x = score_array['combined_score'].fillna(0)

    for i in term_g.nodes:
        values = x.loc[i]
        term_g.node[i]['color'] = 'red'
        term_g.node[i]['label'] = i
        for n, time in enumerate(labels):
            term_g.node[i]['sample{}'.format(time)] = float(values[time])

    if not create_only:
        from magine.networks.visualization.cytoscape_view import RenderModel
        rm = RenderModel(term_g, layout='force-directed')
        rm.visualize_by_list_of_time(labels,
                                     prefix=save_name,
                                     out_dir=out_dir,
                                     )
        _s = "convert -delay 100 -dispose previous -loop 0 " \
             "ont_network_*_formatted.png {}.gif"
        if merge:
            os.system(_s.format(save_name))
    return term_g, molecular_g
