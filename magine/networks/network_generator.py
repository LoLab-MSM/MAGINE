import logging
import os
import networkx as nx
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import magine.networks.utils as nt
import magine.networks.databases as db
from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

try:
    import cPickle as pickle
except ImportError:
    import pickle
from magine.logging import get_logger

# logger = get_logger("magine.networks.network_generator", log_level=logging.INFO)
logger = get_logger(__name__, log_level=logging.INFO)


def build_network(seed_species, species='hsa', save_name=None,
                  all_measured_list=None, trim_source_sink=False,
                  use_reactome=True, use_hmdb=False,
                  use_biogrid=True, use_signor=True, verbose=False):
    """
    Construct a network from a list of gene names.

    Parameters
    ----------

    seed_species : list
        list of genes to construct network
    save_name : str, optional
        output name to save network. Will save one before and after ID
        conversion
    species : str
        species of proteins ('hsa': human, 'mmu':murine)
    all_measured_list : list
        list of all species that should be considered in network
    use_reactome : bool
        Add ReactomeFunctionalInteraction reaction to network
    use_biogrid : bool
        Add BioGrid reaction to network
    use_hmdb : bool
        Add HMDB reaction to network
        all_measured_list
    use_signor : bool
        Add SIGNOR reaction to network
    trim_source_sink : bool, optional
        Remove source and sink nodes if they are not measured in network
    verbose : bool

    Returns
    -------
    networkx.DiGraph
    """
    cm = ChemicalMapper()
    path_to_graph, node_to_path = db.load_kegg_mappings(species)

    seed_species = set(x.upper() for x in seed_species)
    updated_accession = set()
    old_accession = set()
    # This maps HMDB numbers to the new numbering format
    for i in seed_species:
        if i.startswith('HMDB'):
            if i in cm.hmdb_accession_to_main:
                old_accession.add(i)
                updated_accession.add(cm.hmdb_accession_to_main[i][0])

    seed_species.difference_update(old_accession)
    seed_species.update(updated_accession)

    seeds_in_kegg = seed_species.intersection(node_to_path)

    pathway_list = set()
    for seed in seeds_in_kegg:
        pathway_list.update(node_to_path[seed])

    graph_list = []
    for each in pathway_list:
        tmp = path_to_graph[each]
        if len(tmp.edges) == 0:
            continue
        graph_list.append(tmp)
    end_network = nt.compose_all(graph_list)

    if all_measured_list is None:
        all_measured_set = set(i.upper() for i in end_network.nodes)
    else:
        all_measured_set = set(str(x).upper() for x in all_measured_list)

    all_measured_set.update(seed_species)
    hmdb_ids = set(i for i in all_measured_set if i.startswith('HMDB'))
    updated_accession = set()
    old_accession = set()
    for i in hmdb_ids:
        if i in cm.hmdb_accession_to_main:
            all_measured_set.remove(i)
            all_measured_set.add(cm.hmdb_accession_to_main[i][0])
    networks_to_expand = []
    logger.info("Gathering networks")
    if use_hmdb:
        networks_to_expand.append(db.load_hmdb_network())

    if use_reactome:
        networks_to_expand.append(db.load_reactome_fi())

    if use_biogrid:
        networks_to_expand.append(db.load_biogrid_network())

    if use_signor:
        networks_to_expand.append(db.load_signor())
    logger.info("Merging networks")
    if len(networks_to_expand) != 0:
        entire_expansion_network = nt.compose_all(networks_to_expand)
        end_network = expand_by_db(end_network, entire_expansion_network,
                                   all_measured_set)

    logger.info("Trimming network")
    # makes all similar edge names the same
    nt.standardize_edge_types(end_network)
    # removes everything not connected to the largest graph
    end_network = nt.delete_disconnected_network(end_network)

    if trim_source_sink:
        end_network = nt.trim_sink_source_nodes(end_network, all_measured_list,
                                                remove_self_edge=True)
    if save_name is not None:
        nx.write_gml(end_network, '{}.gml'.format(save_name))
        nx.write_gpickle(end_network, '{}.p'.format(save_name))

    final_nodes = set(end_network.nodes)
    n_hits = len(seed_species.intersection(final_nodes))
    logger.info('Network has {} nodes and {} edges'.format(
        len(final_nodes), len(end_network.edges))
    )

    logger.info("Found {} of {} seed species in network".format(
        n_hits, len(seed_species))
    )
    if all_measured_list is not None:
        n_measured_hits = len(set(all_measured_list).intersection(final_nodes))
        logger.info("Found {} of {} background species in network"
                    "".format(n_measured_hits, len(all_measured_list)))

    return end_network


def expand_by_db(starting_network, expansion_source, measured_list,
                 verbose=False):
    """ add reference network to main network

    Parameters
    ----------
    starting_network : nx.DiGraph
    expansion_source : nx.DiGraph
    measured_list : list_like
    verbose : bool

    Returns
    -------
    new_graph : nx.DiGraph
    """

    new_graph = nx.DiGraph()
    measured_set = set(measured_list)
    current_nodes = set(starting_network.nodes)

    # gathers possible nodes to expand

    # must be in new nodes AND measured set
    nodes_to_check = set(expansion_source.nodes).intersection(measured_set)

    # add all current nodes, might be edges between them that are missed
    nodes_to_check.update(current_nodes)

    added_nodes = set()

    for i, j, k in expansion_source.edges(data=True):
        if i in nodes_to_check and j in nodes_to_check:
            added_nodes.add(i)
            added_nodes.add(j)
            new_graph.add_edge(i, j, **k)

    for node in added_nodes:
        new_graph.add_node(node, **expansion_source.node[node])

    new_graph = nt.compose(starting_network, new_graph)
    logger.info("\t\t\tbefore\tafter")
    logger.info("\tNodes\t{}\t{}".format(len(starting_network.nodes),
                                         len(new_graph.nodes)))

    logger.info("\tEdges\t{}\t{}".format(len(starting_network.edges),
                                         len(new_graph.edges)))

    return new_graph


def create_background_network(save_name='background_network',
                              fresh_download=False, verbose=True,
                              create_overlap=False):
    """

    Parameters
    ----------
    save_name : str
        Name of the network
    fresh_download : bool
        Download a fresh copy of the databases
    verbose: bool
        Print information about the databases
    create_overlap : bool
        Creates a figure comparing the databses
    Returns
    -------
    nx.DiGraph
    """
    logger.info("Generating background network from all databases")
    kegg_network = db.load_kegg(fresh_download=fresh_download)
    hmdb_network = db.load_hmdb_network(fresh_download=fresh_download)
    biogrid_network = db.load_biogrid_network(fresh_download=fresh_download)
    signor_network = db.load_signor(fresh_download=fresh_download)
    reactome_network = db.load_reactome_fi()
    network_list = [hmdb_network, kegg_network, biogrid_network,
                    reactome_network, signor_network]
    names = ['hmdb', 'kegg', 'biogrid', 'reactome', 'signor']

    def find_overlap(n1, n2):
        nodes1 = set(n1.nodes())
        nodes2 = set(n2.nodes())
        e1 = set(n1.edges())
        e2 = set(n2.edges())
        edge_overlap = len(e2.intersection(e1))
        node_overlap = len(nodes1.intersection(nodes2))
        return node_overlap, edge_overlap

    if create_overlap:
        db_maps = {i: j for i, j in zip(names, network_list)}
        n_dbs = len(names)

        pal = sns.light_palette("purple", as_cmap=True)
        node_mat = np.zeros((n_dbs, n_dbs), dtype=np.int)
        edge_mat = np.zeros((n_dbs, n_dbs), dtype=np.int)

        for i in range(n_dbs):
            row = db_maps[names[i]]
            for j in range(i + 1, n_dbs):
                col = db_maps[names[j]]
                n_overlap, e_overlap = find_overlap(row, col)
                node_mat[i, j] = n_overlap
                node_mat[j, i] = node_mat[i, j]
                edge_mat[i, j] = e_overlap
                edge_mat[j, i] = edge_mat[i, j]

        fig = plt.figure(figsize=(10, 4))
        ax = fig.add_subplot(121)
        plt.title("Number of node overlaps")
        sns.heatmap(node_mat, fmt='d', annot=True, linewidths=0.02, cmap=pal,
                    yticklabels=names, xticklabels=names)
        plt.yticks(rotation=0)
        ax = fig.add_subplot(122)
        plt.title("Number of edge overlaps")
        sns.heatmap(edge_mat, fmt='d', annot=True, linewidths=0.02, cmap=pal,
                    yticklabels=names, xticklabels=names)
        plt.tight_layout()
        plt.yticks(rotation=0)
        plt.subplots_adjust(wspace=.3)
        plt.savefig('compare_network_dbs.png', dpi=300, bbox_inches='tight')
        plt.close()
    logger.info("Combining networks")
    full_network = nt.compose_all(network_list)
    logger.info("Finished combining networks")
    nt.standardize_edge_types(full_network)
    full_network = nt.delete_disconnected_network(full_network)
    n_nodes = len(full_network.nodes)
    n_edges = len(full_network.edges)
    logger.info("Background network {} nodes and {} edges".format(n_nodes,
                                                                  n_edges))

    nx.write_gpickle(full_network, '{}.p.gz'.format(save_name))
    nx.write_gml(full_network, '{}.gml'.format(save_name))
    return full_network


if __name__ == '__main__':
    create_background_network(fresh_download=False, verbose=True,
                              create_overlap=True)
