# -*- coding: utf-8 -*-
import os

import networkx as nx

import magine.mappings.maps as mapper
import magine.networks.utils as nt
from magine.data.storage import network_data_dir
from magine.mappings.gene_mapper import GeneMapper
from magine.networks.databases import load_reactome_fi
from magine.networks.databases.biogrid_interactions import \
    create_biogrid_network
from magine.networks.databases.kegg_kgml import load_kegg
from magine.networks.databases.signor import load_signor

try:
    import cPickle as pickle
except ImportError:
    import pickle


def build_network(gene_list, species='hsa', save_name=None,
                  all_measured_list=None, metabolite_list=None,
                  use_reactome=True, use_hmdb=False, use_biogrid=True,
                  use_signor=True, trim_source_sink=False):
    """
    Construct a network from a list of gene names.

    Parameters
    ----------

    gene_list : list
        list of genes to construct network
    save_name : str
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
    metabolite_list : list
        List of metabolites with HMDB ids
    trim_source_sink : bool, optional
        Remove source and sink nodes if they are not measured in network
    Returns
    -------
    networkx.DiGraph
    """

    gm = GeneMapper()

    path_to_graph, node_to_path = load_kegg(species, verbose=False)

    gene_list = set(x.upper() for x in gene_list)

    to_remove = gene_list.difference(gm.gene_name_to_kegg)
    genes_in_kegg = gene_list.intersection(gm.gene_name_to_kegg)
    pathway_list = set()

    for gene in genes_in_kegg:
        tmp_list = gm.gene_name_to_kegg[gene]
        if len(tmp_list) == 0:
            to_remove.add(gene)
        for n in tmp_list:
            if n in node_to_path:
                for j in node_to_path[n]:
                    pathway_list.add(j)

    if len(to_remove) != 0:
        print("{} species not found in KEGG".format(len(to_remove)))

    graph_list = []
    for each in pathway_list:
        tmp = path_to_graph[each]
        if len(tmp.edges) == 0:
            continue
        graph_list.append(tmp)

    end_network = nt.compose_all(graph_list)
    end_network = mapper.convert_all(end_network, species=species)

    if all_measured_list is None:
        all_measured_list = gene_list
    else:
        all_measured_list = set(str(x).upper() for x in all_measured_list)

    if use_hmdb:
        if metabolite_list is None:
            print("warning: automatic integration currently in progress.\n")
            print("Please provide a list of metabolites")
            metabolite_list = list(set(i for i in end_network.nodes()
                                       if i.startswith('HMDB')))
        end_network = expand_by_hmdb(end_network, metabolite_list)

    if use_reactome:
        end_network = expand_by_db(end_network, all_measured_list, 'reactome')

    if use_biogrid:
        end_network = expand_by_db(end_network, all_measured_list, 'biogrid')

    if use_signor:
        end_network = expand_by_db(end_network, all_measured_list, 'signor')

    print("Trimming network")
    # removes everything not connected to the largest graph
    nt.delete_disconnected_network(end_network)

    # makes all similar edge names the same
    nt.standardize_edge_types(end_network)

    if trim_source_sink:
        end_network = nt.trim_sink_source_nodes(end_network, all_measured_list,
                                                remove_self_edge=True)
    print('Network has {} nodes and {} edges'
          ''.format(len(end_network.nodes()),
                    len(end_network.edges()))
          )

    nx.write_gml(end_network, '{}.gml'.format(save_name))
    nx.write_gpickle(end_network, '{}.p'.format(save_name))

    return end_network


def expand_by_db(network, measured_list, db='reactome'):
    """
    add reference network to main network
    Parameters
    ----------
    network : nx.DiGraph
    measured_list : list
    db : str


    Returns
    -------

    """

    new_graph = nx.DiGraph()
    measured_set = set(measured_list)
    if db == 'reactome':
        network_to_add = load_reactome_fi()
        db_name = 'ReactomeFI'

    elif db == 'biogrid':
        network_to_add = create_biogrid_network()
        db_name = 'BioGrid'
    elif db == 'signor':
        network_to_add = load_signor()
        db_name = 'signor'
    else:
        print("Must provide 'reactome', 'biogrid', or 'signor'")
        return network

    edges = network_to_add.edges(data=True)

    # new list in reactome_nodes
    nodes_to_check = set(network_to_add.nodes()).intersection(measured_set)
    nodes_to_check.update(set(network.nodes()))

    added_nodes = set()

    for i, j, k in edges:
        k['databaseSource'] = db_name
        if i in nodes_to_check and j in measured_list:
            added_nodes.add(i)
            added_nodes.add(j)
            new_graph.add_edge(i, j, **k)
        elif j in nodes_to_check and i in measured_list:
            added_nodes.add(i)
            added_nodes.add(j)
            new_graph.add_edge(i, j, **k)

    for node in added_nodes:
        attr = network_to_add.node[node]
        new_graph.add_node(node, **attr)

    new_graph = nt.compose(network, new_graph)

    print("{} database".format(db_name))
    print("\t\tbefore\tafter")
    print(
        "\tNodes\t{}\t{}".format(len(network.nodes()), len(new_graph.nodes())))
    print(
        "\tEdges\t{}\t{}".format(len(network.edges()), len(new_graph.edges())))

    return new_graph


def expand_by_hmdb(graph, metabolite_list, verbose=False):
    """
    Expands a network using HMDB metabolites-protein information

    Parameters
    ----------
    graph : nx.DiGraph
    metabolite_list : list
        List of HMDB ids

    Returns
    -------
    nx.DiGraph

    Examples
    --------
    >>> g = nx.DiGraph()
    >>> g.add_edge('PNLIP', 'LIPC')
    >>> new_g = expand_by_hmdb(graph=g, \
                           metabolite_list=['HMDB42489', 'HMDB59874'], \
                           )
    Metabolites data= 2
    Metabolites not in starting network = 2
    Metabolites added to network = 1
    Metabolites not able to add = 1
    Before # of nodes = 2, edge = 1
    After # of nodes = 3, edge = 3
    Added 2 nodes and 1 edges
    Number of proteins not in network = 33
    Number of add proteins = 0
    <BLANKLINE>
    missing metabolites-protein info = 0
    >>> len(new_g.nodes())
    3
    >>> len(new_g.edges())
    3
    """
    from magine.mappings.chemical_mapper import ChemicalMapper

    cm = ChemicalMapper()

    tmp_graph = nx.DiGraph()
    start_nodes = set(graph.nodes())
    start_edges = set(graph.edges())
    metabolite_set = set(i for i in metabolite_list if i.startswith('HMDB'))
    metabolite_set = metabolite_set.intersection(cm.hmdb_accession_to_protein)
    missing_metabolites = metabolite_set.difference(start_nodes)
    count_in_network, count_not_in_network = 0, 0
    missing_edge = 0
    protein_hits = 0
    added_proteins = set()
    missing_protein_info = 0
    missing_proteins = set()

    def _add_node(node, node_type):
        if node != u'':
            attrs = dict()
            attrs['databaseSource'] = 'HMDB'
            attrs['speciesType'] = node_type
            if node_type == 'compound':
                if node in cm.hmdb_accession_to_chemical_name:
                    name = cm.hmdb_accession_to_chemical_name[node]
                    attrs['chemName'] = name

            tmp_graph.add_node(node, **attrs)

    def _add_edge(source, target):
        if source != u'' and target != u'':
            tmp_graph.add_edge(source, target,
                               interactionType='chemical',
                               databaseSource='HMDB'
                               )

    for i in metabolite_set:
        tmp_list = cm.hmdb_accession_to_protein[i]
        if tmp_list is None or len(tmp_list) == 0:
            # No information to add so skip
            continue

        # checks if metabolite is in the network
        # if it is, it can add an associated gene
        if i in start_nodes:
            count_in_network += 1
            for each in tmp_list:
                if each is None:
                    continue
                elif each not in start_nodes:
                    added_proteins.add(each)
                _add_node(each, 'gene')
                missing_edge += 1
                _add_edge(each, i)
        # if the metabolite is NOT in the network,
        # it checks to see if it has any protein relationships
        # if it does and they are new to the graph, it adds it
        else:
            count_not_in_network += 1
            for each in tmp_list:
                if each is None:
                    continue
                elif each in start_nodes:
                    _add_node(i, 'compound')
                    _add_edge(i, each)
                    missing_edge += 1
                    protein_hits += 1
                else:
                    missing_proteins.add(each)

    final_graph = nt.compose(graph, tmp_graph)
    end_nodes = set(final_graph.nodes())
    end_edges = set(final_graph.edges())
    new_nodes = end_nodes.intersection(start_nodes)
    new_edges = end_edges.intersection(start_edges)
    metabolites_added = missing_metabolites.intersection(end_nodes)
    still_missing = missing_metabolites.difference(end_nodes)
    if verbose:
        out = list()
        out.append('Metabolites data= {}'.format(len(metabolite_set)))
        out.append('Metabolites not in starting network = {}'
                   ''.format(len(missing_metabolites)))
        out.append('Metabolites added to network = {}'
                   ''.format(len(metabolites_added)))
        out.append(
            'Metabolites not able to add = {0}'.format(len(still_missing)))

        out.append('Before # nodes = {} edges = {}'.format(len(start_nodes),
                                                           len(graph.edges())))

        out.append('After # of nodes = {},'
                   ' edge = {}'.format(len(end_nodes),
                                       len(final_graph.edges())))

        out.append('Added {} nodes and {} edges'.format(len(new_nodes),
                                                        len(new_edges)))

        out.append('Number of proteins not in network'
                   ' = {}'.format(len(missing_proteins)))
        out.append('Number of add proteins = {0}'.format(len(added_proteins)))
        out.append('missing metabolites-protein info = {0}'.format(
            missing_protein_info))
        print('\n'.join(out))
    return final_graph


def create_hmdb_network():
    out_name = os.path.join(network_data_dir, 'hmdb_graph.p.gz')
    if os.path.exists(out_name):
        tmp_graph = nx.read_gpickle(out_name)
        print("HMDB network {} nodes and {} edges"
              "".format(len(tmp_graph.nodes()), len(tmp_graph.edges()))
              )
        return tmp_graph
    from magine.mappings.chemical_mapper import ChemicalMapper

    cm = ChemicalMapper()

    tmp_graph = nx.DiGraph()

    def _add_node(node, node_type):
        if node != u'':
            attrs = dict()
            attrs['databaseSource'] = 'HMDB'
            attrs['speciesType'] = node_type
            if node_type == 'compound':
                if node in cm.hmdb_accession_to_chemical_name:
                    name = cm.hmdb_accession_to_chemical_name[node]
                    attrs['chemName'] = name

            tmp_graph.add_node(node, **attrs)

    def _add_edge(source, target):
        if source != u'' and target != u'':
            tmp_graph.add_edge(source, target,
                               interactionType='chemical',
                               databaseSource='HMDB'
                               )

    for compound, genes in cm.hmdb_main_accession_to_protein.items():
        _add_node(compound, 'compound')
        for ge in genes:
            _add_node(ge, 'gene')
            _add_edge(compound, ge)
    print("HMDB network {} nodes and {} edges"
          "".format(len(tmp_graph.nodes()), len(tmp_graph.edges()))
          )
    nx.write_gpickle(tmp_graph, out_name)
    return tmp_graph


def create_background_network(save_name='background_network'):
    """

    Parameters
    ----------
    save_name : str

    Returns
    -------

    """

    from magine.networks.databases.biogrid_interactions import \
        create_biogrid_network
    from magine.networks.databases.kegg_kgml import create_all_of_kegg
    from magine.networks.databases.signor import load_signor

    reactome_network = load_reactome_fi()
    kegg_network = create_all_of_kegg()
    hmdb_network = create_hmdb_network()
    biogrid_network = create_biogrid_network()
    signor_network = load_signor()

    def find_overlap(n1, n2):
        e1 = set(n1.edges())
        e2 = set(n2.edges())
        print(len(e1), len(e2), len(e2.difference(e1)))
        e1 = set(n1.nodes())
        e2 = set(n2.nodes())
        print(len(e1), len(e2), len(e2.difference(e1)))
        for i in sorted(e1.difference(e2)):
            print(i)

    full_network = nt.compose_all(
        [hmdb_network, kegg_network, biogrid_network, reactome_network,
         signor_network]
    )

    nt.delete_disconnected_network(full_network)
    nt.standardize_edge_types(full_network)

    # find_overlap(reactome_network, full_network)
    n_nodes = len(full_network.nodes())
    n_edges = len(full_network.edges())
    print("Background network {} nodes and {} edges".format(n_nodes, n_edges))

    nx.write_gpickle(full_network, '{}.p.gz'.format(save_name))


if __name__ == '__main__':
    create_background_network()
