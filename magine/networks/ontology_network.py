import itertools
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats

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
        if self.network is None:
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

        n_edges = float(len(self.edges))
        p_values = []
        for term1, term2 in itertools.combinations(list_of_go_terms, 2):
            term_1 = set(ont_to_species_dict[term1]).intersection(self.nodes)
            term_2 = set(ont_to_species_dict[term2]).intersection(self.nodes)
            total = 0
            for i in term_1:
                for j in term_2:
                    if i != j:
                        total += 1
            # total_poss = len(list(itertools.product(term_1, term_2)))
            edges_1_to_2 = set(itertools.product(
                term_1, term_2.difference(term_1))
            ).intersection(self.edges)
            edges_2_to_1 = set(itertools.product(
                term_2, term_1.difference(term_2))
            ).intersection(self.edges)

            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]

            a_to_b = len(edges_1_to_2)
            b_to_a = len(edges_2_to_1)

            p_values.append([term1, term2, edges_1_to_2, a_to_b,
                             stats.binom_test(a_to_b, total, total / n_edges)])
            p_values.append([term2, term1, edges_2_to_1, b_to_a,
                             stats.binom_test(b_to_a, total, total / n_edges)])

            # add to graph if at least one edge found between terms
            if a_to_b or b_to_a:
                x = self.network.edge_subgraph(edges_1_to_2).copy()
                mol_net.add_edges_from(list(x.edges(data=True)))
                for gene in x.nodes:
                    _add(gene, term1, term_1, label_1)
                    _add(gene, term2, term_2, label_2)
                all_genes.update(x.nodes)
            # else:
            #     print('No edges between {} and {}'.format(label_1, label_2))
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
                render_igraph(mol_net, save_name + '_subgraph_igraph')
        cols = ['term1', 'term2', 'edges', 'n_edges', 'p_values']
        df = pd.DataFrame(p_values, columns=cols)
        df['sig'], df['adj_p_values'] = fdrcorrection(df['p_values'])
        new_df = df[df['adj_p_values'] < .05]
        pd.set_option('display.max_colwidth', -1)
        print(new_df)
        go_graph = nx.DiGraph()

        for term1, term2, edge, n_edges, p_value in new_df[cols].values:
            label_1 = ont_to_label_dict[term1]
            label_2 = ont_to_label_dict[term2]
            go_graph.add_node(label_1, term=term1, label=label_1)
            go_graph.add_node(label_2, term=term2, label=label_2)
            if n_edges > threshold:
                go_graph.add_edge(label_1, label_2, label=str(n_edges),
                                  weight=n_edges, pvalue=p_value)
        print(go_graph.nodes)
        print(go_graph.edges)
        print(df[['term1', 'term2']].values)
        print(new_df[['term1', 'term2']].values)

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
                      threshold=0, remove_isolated=False, create_only=True,
                      merge=False, out_dir=None):
    """

    Parameters
    ----------

    df : pd.DataFrame
    network : nx.DiGraph
    terms : list, optional
        List of terms to use in ont network. Default is all terms
    save_name : str
    draw_png : bool
    threshold : float, int
        Threshold for number of edges between two terms to consider in graph
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
        threshold=threshold
    )

    if remove_isolated:
        nt.remove_isolated_nodes(term_g)
        nt.remove_isolated_nodes(molecular_g)

    heatmap_from_array(df_copy, cluster_row=False, convert_to_log=True,
                       index='term_name', values='combined_score',
                       columns='sample_id', div_colors=False)

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


def fdrcorrection(p_vals, alpha=0.05):
    """ Benjamini/Hochberg false discovery rate correction

    Parameters
    ----------
    p_vals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate

    Returns
    -------
    rejected : np.array
        True if a hypothesis is rejected, False if not
    adj_p_values : np.array

    """
    p_vals = np.asarray(p_vals)
    original_sort = np.argsort(p_vals)
    p_vals_sorted = np.take(p_vals, original_sort)
    n_samples = len(p_vals)
    ecdf_factor = np.arange(1, n_samples + 1) / float(n_samples)

    reject = p_vals_sorted <= ecdf_factor * alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    corrected_p_vals = p_vals_sorted / ecdf_factor
    pvals_corrected = np.minimum.accumulate(corrected_p_vals[::-1])[::-1]

    pvals_corrected[pvals_corrected > 1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[original_sort] = pvals_corrected

    reject_ = np.empty_like(reject)
    reject_[original_sort] = reject
    return reject_, pvals_corrected_
