import os
import networkx as nx
import magine.networks.utils as nt
from magine.networks.databases import load_hmdb_network, load_biogrid_network, \
    load_signor, load_reactome_fi, load_kegg_mappings
from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

try:
    import cPickle as pickle
except ImportError:
    import pickle

cm = ChemicalMapper()


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

    path_to_graph, node_to_path = load_kegg_mappings(species, verbose=False)

    seed_species = set(x.upper() for x in seed_species)
    updated_accession = set()
    old_accession = set()
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

    if use_hmdb:
        networks_to_expand.append(load_hmdb_network(verbose))

    if use_reactome:
        networks_to_expand.append(load_reactome_fi(verbose))

    if use_biogrid:
        networks_to_expand.append(load_biogrid_network(verbose))

    if use_signor:
        networks_to_expand.append(load_signor(verbose))

    if len(networks_to_expand) != 0:
        entire_expansion_network = nt.compose_all(networks_to_expand)
        end_network = expand_by_db(end_network, entire_expansion_network,
                                   all_measured_set)

    print("Trimming network")
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

    print('Network has {} nodes and {} edges'.format(len(final_nodes),
                                                     len(end_network.edges)))

    print("Found {} of {} seed species in network"
          "".format(n_hits, len(seed_species)))
    if all_measured_list is not None:
        n_measured_hits = len(set(all_measured_list).intersection(final_nodes))
        print("Found {} of {} background species in network"
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
    if verbose:
        print("\t\t\tbefore\tafter")
        print("\tNodes\t{}\t{}".format(len(starting_network.nodes),
                                       len(new_graph.nodes)))

        print("\tEdges\t{}\t{}".format(len(starting_network.edges),
                                       len(new_graph.edges)))

    return new_graph


def create_background_network(save_name='background_network'):
    """

    Parameters
    ----------
    save_name : str

    Returns
    -------

    """

    reactome_network = load_reactome_fi()
    kegg_network = load_all_of_kegg()
    hmdb_network = load_hmdb_network()
    biogrid_network = load_biogrid_network()
    signor_network = load_signor()

    def find_overlap(n1, n2):
        nodes1 = set(n1.nodes())
        nodes2 = set(n2.nodes())
        e1 = set(n1.edges())
        e2 = set(n2.edges())
        print("\tnode overlap = {}".format(len(nodes1.intersection(nodes2))))
        print("\tnode difference = {} | {}".format(
            len(nodes1.intersection(nodes2)),
            len(nodes2.intersection(nodes1)))
        )
        print("\tedge overlap = {}".format(len(e2.intersection(e1))))
        print("\tedge difference = {} | {}".format(len(e2.difference(e1)),
                                                   len(e2.difference(e1))))

    network_list = [hmdb_network, kegg_network, biogrid_network,
                    reactome_network, signor_network]
    names = ['hmdb', 'kegg', 'biogrid', 'reactome', 'signor']

    for i, n in zip(network_list, names):
        for j, m in zip(network_list, names):
            if n != m:
                print('{} : {}'.format(n, m))
                find_overlap(i, j)

    full_network = nt.compose_all(
        [hmdb_network, kegg_network, biogrid_network, reactome_network,
         signor_network]
    )

    nt.delete_disconnected_network(full_network)
    nt.standardize_edge_types(full_network)

    # find_overlap(reactome_network, full_network)
    n_nodes = len(full_network.nodes)
    n_edges = len(full_network.edges)
    print("Background network {} nodes and {} edges".format(n_nodes, n_edges))

    nx.write_gpickle(full_network, '{}.p.gz'.format(save_name))


if __name__ == '__main__':
    # create_background_network()
    # load_hmdb_network(create_new=True)
    load_hmdb_network(False)
