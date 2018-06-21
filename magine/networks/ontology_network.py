import itertools
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection

import magine.networks.utils as nt
from magine.networks.exporters import export_to_dot
from magine.networks.visualization.igraph_tools import render_igraph
from magine.plotting.heatmaps import heatmap_from_array


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
                                 use_threshold=True):
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

        def get_edges(set1, set2, background):
            non_overlap = set2.difference(set1)
            possible_edges = set(itertools.product(set1, non_overlap))
            edge_hits = possible_edges.intersection(background)
            n_hits = len(edge_hits)
            return edge_hits, n_hits

        for term1, term2 in itertools.combinations(list_of_go_terms, 2):
            term_1 = set(ont_to_species_dict[term1]).intersection(self.nodes)
            term_2 = set(ont_to_species_dict[term2]).intersection(self.nodes)
            add_to_dict(term_1, term1)
            add_to_dict(term_2, term2)

            # calculate how many possible edge combinations there are
            total = 0
            for i in term_1:
                for j in term_2:
                    # need to not count an edge if one of the species
                    # is in both
                    # g1 = (a, b)
                    # g3 = (b, d)
                    # we count a-d, and b-d, not a-b.
                    if i != j:
                        total += 1

            edges_1_to_2, a_to_b = get_edges(term_1, term_2, self.edges)
            edges_2_to_1, b_to_a = get_edges(term_2, term_1, self.edges)

            odds_to_find = float(total) / n_edges
            p_values.append([term1, term2, edges_1_to_2, a_to_b,
                             stats.binom_test(a_to_b, total, odds_to_find)])

            p_values.append([term2, term1, edges_2_to_1, b_to_a,
                             stats.binom_test(b_to_a, total, odds_to_find)])

        cols = ['term1', 'term2', 'edges', 'n_edges', 'p_values']

        df = pd.DataFrame(p_values, columns=cols)
        # FDR correction
        _, df['adj_p_values'] = fdrcorrection(df['p_values'])

        if use_threshold:
            df = df.loc[df['adj_p_values'] < .05]

        # create empty networks
        go_graph = nx.DiGraph()
        mol_net = nx.DiGraph()
        for term1, term2, edge, n_edges, p_value in df[cols].values:
            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]
            go_graph.add_node(label_1, term=term1, label=label_1)
            go_graph.add_node(label_2, term=term2, label=label_2)
            go_graph.add_edge(label_1, label_2, label=str(n_edges),
                              weight=n_edges, pvalue=p_value)

            nodes = list(itertools.chain(*edge))

            nodes += ont_to_species_dict[term1]
            nodes += ont_to_species_dict[term2]

            x = self.network.subgraph(nodes).copy()

            mol_net.add_edges_from(list(x.edges(data=True)))

        for i in mol_net.nodes:
            labels = gene_to_label[i]
            terms = gene_to_term[i]
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
                render_igraph(mol_net, save_name + '_subgraph_igraph')
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

    labels = data['sample_id'].unique()

    score_array = pd.pivot_table(data, index=['term_name'],
                                 columns='sample_id')

    heatmap_from_array(data, cluster_row=False, convert_to_log=True,
                       index='term_name', values='combined_score',
                       columns='sample_id', div_colors=False)

    plt.savefig('{}_heatplot.png'.format(save_name), bbox_inches='tight')

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
    _s = "convert -delay 100 -dispose previous -loop 0 " \
         "ont_network_*_formatted.png {}.gif"
    if merge:
        os.system(_s.format(save_name))


def create_subnetwork(df, network, terms=None, save_name=None, draw_png=False,
                      remove_isolated=False, create_only=True, merge=False,
                      out_dir=None):
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
        use_threshold=True
    )

    if remove_isolated:
        nt.remove_isolated_nodes(term_g)
        nt.remove_isolated_nodes(molecular_g)

    fig = heatmap_from_array(df_copy, cluster_row=False, convert_to_log=True,
                             index='term_name', values='combined_score',
                             columns='sample_id', div_colors=False)
    fig.savefig('{}.png'.format(save_name), bbox_inches='tight')

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

# def fdrcorrection(p_vals):
#     """ Benjamini/Hochberg false discovery rate correction
#
#     Parameters
#     ----------
#     p_vals : array_like
#         set of p-values of the individual tests.
#
#     Returns
#     -------
#     adj_p_values : np.array
#
#     """
#     p_vals = np.asarray(p_vals)
#     n_samples = len(p_vals)
#     ecdf_factor = np.arange(1, n_samples + 1) / float(n_samples)
#     corrected_p_vals = p_vals / ecdf_factor
#     pvals_corrected = np.minimum.accumulate(corrected_p_vals[::-1])[::-1]
#     pvals_corrected[pvals_corrected > 1] = 1
#     return pvals_corrected
