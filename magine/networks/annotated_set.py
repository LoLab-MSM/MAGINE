import itertools
import os
import random
from itertools import combinations, product

import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection

import magine.networks.utils as nt
from magine.networks.visualization import draw_graphviz, draw_igraph


class AnnotatedSetNetworkGenerator(object):
    """ Generates ontology network from molecular networks.

    Nodes are the ontology terms
    Edges between nodes are determined based on molecular network graphs.

    """

    def __init__(self, network):
        """

        Parameters
        ----------
        network : nx.DiGraph
        """
        # make sure a network exists
        if not(isinstance(network, nx.DiGraph)):
            raise Exception("Must provide a nx.DiGraph")
        self.network = network.copy()
        # remove all self loops from graph
        self.network.remove_edges_from(nx.selfloop_edges(self.network))
        self.nodes = set(self.network.nodes)
        self.edges = set(self.network.edges)
        self.molecular_network = None

    def create_network_from_list(self, list_of_ontology_terms,
                                 ont_to_species_dict, ont_to_label_dict,
                                 save_name=None, draw=False, out_dir=None,
                                 use_threshold=True, use_fdr=False,
                                 min_edges=0):
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
        min_edges : int
            Minimum number of edges between terms

        Returns
        -------
        asn : networkx.DiGraph
            Annotated set network
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

        list_of_go_terms = set(list_of_ontology_terms)

        dd = nx.out_degree_centrality(self.network)
        x = 0
        for _, i in dd.items():
            x += i
        avg_conn = x / len(dd)
        # n_total_edges = len(self.edges)
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
            # odds_to_find = float(n_possible) / n_total_edges
            # p_val = stats.binom_test(n_hits, n_possible, odds_to_find)
            # print(stats.binom_test(n_hits, n_possible, odds_to_find),
            #       stats.binom_test(n_hits, n_possible, avg_conn))
            p_val = stats.binom_test(n_hits, n_possible, avg_conn)
            return [edge_hits, n_hits, p_val]

        for term1, term2 in combinations(list_of_go_terms, 2):
            term_1 = set(ont_to_species_dict[term1]).intersection(self.nodes)
            term_2 = set(ont_to_species_dict[term2]).intersection(self.nodes)

            add_to_dict(term_1, term1)
            add_to_dict(term_2, term2)
            # print(get_edges(term_1, term_2)[1:],
            # gather_stats(term_1, term_2, self.nodes, self.edges)[1:])

            # p_values.append(
            #     [term1, term2] + gather_stats(term_1, term_2, self.nodes,
            #                                   self.edges))
            # p_values.append(
            #     [term2, term1] + gather_stats(term_2, term_1, self.nodes,
            #                                   self.edges))

            p_values.append([term1, term2] + get_edges(term_1, term_2))
            p_values.append([term2, term1] + get_edges(term_2, term_1))

        cols = ['term1', 'term2', 'edges', 'n_edges', 'p_value', ]

        df = pd.DataFrame(p_values, columns=cols)
        df = df.loc[df['n_edges'] > min_edges].copy()
        # FDR correction
        _, df['adj_p_value'] = fdrcorrection(df['p_value'])

        if use_threshold:
            if use_fdr:
                df = df.loc[df['adj_p_value'] <= .05]
            else:
                df = df.loc[df['p_value'] <= .05]
        cols += ['adj_p_value']
        # create empty networks
        asn = nx.DiGraph()
        mol_net = nx.DiGraph()
        for term1, term2, edge, n_edges, p_value, adj_pval in df[cols].values:

            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]
            asn.add_node(label_1, term=term1, label=label_1)
            asn.add_node(label_2, term=term2, label=label_2)
            asn.add_edge(label_1, label_2, label=str(n_edges), weight=n_edges,
                         pvalue=p_value, adjPval=adj_pval)

            nodes = list(itertools.chain(*edge))

            nodes += ont_to_species_dict[term1]
            nodes += ont_to_species_dict[term2]

            x = self.network.subgraph(nodes).copy()

            mol_net.add_edges_from(list(x.edges(data=True)))

        # addes missing edges between that maybe absent due to threshold, want
        # to keep in network
        for i, j, d in self.network.subgraph(mol_net.nodes).edges(data=True):
            mol_net.add_edge(i, j, **d)

        for i in mol_net.nodes:
            labels = gene_to_label[i]
            terms = gene_to_term[i]
            if len(labels) != len(terms):
                raise AssertionError('len(labels) should equal len(terms)')
            mol_net.node[i]['termName'] = ','.join(sorted(labels))
            mol_net.node[i]['terms'] = ','.join(sorted(terms))

        self.molecular_network = mol_net
        if save_name is not None:
            if out_dir is not None:
                save_name = os.path.join(out_path, save_name)
            nx.write_gml(mol_net, '{}_subgraph.gml'.format(save_name))
            nx.write_gml(asn, '{}.gml'.format(save_name))

            if draw:
                draw_graphviz(asn, save_name=save_name)
                draw_igraph(mol_net, save_name + '_subgraph_igraph')
        return asn, mol_net


def gather_stats(set1, set2, backgroud_nodes, edges):
    n_set1 = len(set1)
    n_set2 = len(set2.difference(set1))
    total_nodes = n_set1 + n_set2

    n_tests = 100
    odds = np.zeros(n_tests)
    nodes_list = list(backgroud_nodes)
    n_possible = range(len(backgroud_nodes))
    for i in range(n_tests):
        s = [nodes_list[j] for j in random.sample(n_possible, total_nodes)]
        odds[i] = len(set(product(s[:n_set1], s[n_set1:])).intersection(edges))

    # calculates edges between two sets
    possible_edges = set(product(set1, set2.difference(set1)))
    edge_hits = possible_edges.intersection(edges)
    n_hits = len(edge_hits)

    # sort to find probability of finding random sets connectivity
    odds.sort()
    first_index = np.searchsorted(odds, n_hits, side='right')
    probability = float(n_tests - first_index) / n_tests
    p_val = stats.binom_test(n_hits, len(possible_edges), probability)

    return [edge_hits, len(edge_hits), p_val]


def create_subnetwork(df, network, terms=None, save_name=None, draw_png=False,
                      remove_isolated=False, use_cytoscape=False, merge=False,
                      out_dir=None, use_threshold=False, use_fdr=False,
                      min_edges=0):
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
    use_cytoscape : bool
        Create ontology network only, don't visualize with cytoscape
    merge : bool
        Create gif of cytoscape session
    out_dir : str
        Output directory of images
    use_threshold : bool
        Use binomial test to calculate edge relevance between terms
    use_fdr : bool
        Use BH correction on binomial test
    min_edges : int
        Minimum number of edges between two terms

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
    # df_copy['combined_score'] = np.log2(df_copy['combined_score'])
    # df_copy.loc[df['combined_score'] > 150, 'combined_score'] = 150

    labels = df_copy['sample_id'].unique()
    # create dictionary of values
    label_dict, term_dict = dict(), dict()
    for i in terms:
        genes = set(df_copy.term_to_genes(i))
        term_dict[i] = genes
        label_dict[i] = i

    ong = AnnotatedSetNetworkGenerator(network=network)

    print("Creating ontology network")
    term_g, molecular_g = ong.create_network_from_list(
        terms, term_dict, label_dict, save_name=save_name, draw=draw_png,
        use_threshold=use_threshold, use_fdr=use_fdr, min_edges=min_edges
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
            fmt_label = 'sample{}'.format(time.replace('_', ''))
            term_g.node[i][fmt_label] = float(values[time])
    labels = [i.replace('_', '') for i in labels]
    if save_name:
        nx.write_graphml(term_g, "{}_ags_network.graphml".format(save_name))
    if use_cytoscape:
        from magine.networks.visualization.cytoscape import RenderModel
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
