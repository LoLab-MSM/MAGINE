# -*- coding: utf-8 -*-
"""
File that generates networks
"""
import os

import networkx as nx
import numpy as np

import magine.mappings.maps as mapper
import magine.networks.network_tools as nt
from magine.data.storage import network_data_dir
from magine.mappings.gene_mapper import GeneMapper
from magine.networks.databases import download_all_of_kegg, load_reactome_fi

try:
    import cPickle as pickle
except:
    import pickle


def build_network(gene_list, save_name='tmp', species='hsa',
                  all_measured_list=None, use_reactome=True, use_hmdb=False,
                  metabolite_list=None):
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
        use _reactome functional interaction to expand network
    use_hmdb : bool
        use hmdb to expand network incorporating metabolites provided in
        all_measured_list
    metabolite_list : list
        List of metabolites with HMDB ids
    Returns
    -------
    networkx.DiGraph
    """
    end_network = nx.DiGraph()
    gm = GeneMapper()
    _kegg_raw_out_path = os.path.join(network_data_dir, 'KEGG')
    _kegg_node_to_pathway = os.path.join(network_data_dir,
                                         'kegg_node_to_pathway.p')

    if not os.path.exists(_kegg_raw_out_path) \
            or not os.path.exists(_kegg_node_to_pathway):
        download_all_of_kegg(species=species)

    node_to_path = pickle.load(open(_kegg_node_to_pathway, 'rb'))

    to_remove = set()
    pathway_list = set()
    gene_list = set(x.upper() for x in gene_list)

    for i in gene_list:
        gene = i.upper()
        if i not in gm.gene_name_to_kegg:
            to_remove.add(gene)
        else:
            tmp_list = gm.gene_name_to_kegg[gene]
            if len(tmp_list) == 0:
                to_remove.add(gene)
            for n in tmp_list:
                if n in node_to_path:
                    for j in node_to_path[n]:
                        pathway_list.add(str(j.replace(':', '')[4:]))
    if len(to_remove) != 0:
        print("{} not found in KEGG".format(to_remove))

    for each in pathway_list:
        tmp = nx.read_gpickle(
                os.path.join(_kegg_raw_out_path, "{}.p".format(each)))
        for n in tmp.nodes():
            if isinstance(n, float):
                print("Found the float... {}".format(each))
                tmp.remove_node(np.nan)
        if len(tmp.edges()) == 0:
            continue
        end_network = nt.compose(end_network, tmp)

    drug_dict = {}
    for i in end_network.nodes():
        if i.startswith('dr'):
            split_name = i.split(' ')
            if len(split_name) > 1:
                if split_name[1].startswith('cpd:'):
                    drug_dict[i] = split_name[1]
                    end_network.node[i]['drug'] = split_name[0]
        elif i == 'nan':
            end_network.remove_node(i)
        elif isinstance(i, float):
            end_network.remove_node(i)
    end_network = nx.relabel_nodes(end_network, drug_dict)

    end_network = mapper.convert_all(end_network, species=species)

    if use_reactome:
        if all_measured_list is None:
            end_network = expand_by_reactome(end_network, gene_list)
        else:
            all_measured_list = set(str(x).upper() for x in all_measured_list)
            end_network = expand_by_reactome(end_network, all_measured_list)

    if use_hmdb:
        print("warning: automatic integration currently in progress.\n")
        if metabolite_list is None:
            print("Please provide a list of metabolites")
            metabolite_list = list(set(i for i in end_network.nodes()
                                       if i.startswith('HMDB')))
            end_network = expand_by_hmdb(end_network, metabolite_list)
        else:
            end_network = expand_by_hmdb(end_network, metabolite_list)
    print("Trimming network")
    # removes everything not connected to the largest graph
    nt.delete_disconnected_network(end_network)
    print('Network has {} nodes and {} edges'
          ''.format(len(end_network.nodes()),
                    len(end_network.edges())
                    )
          )

    # nx.write_gml(end_network, '{}.gml'.format(save_name))
    nx.write_gpickle(end_network, '{}.p'.format(save_name))

    return end_network


def expand_by_hmdb(graph, metabolite_list):
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
    Loading class data
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
    2
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
    out = list()
    out.append('Metabolites not in starting network = {}'
               ''.format(len(missing_metabolites)))
    out.append('Metabolites added to network = {}'
               ''.format(len(metabolites_added)))
    out.append('Metabolites not able to add = {0}'.format(len(still_missing)))

    out.append('Before # of nodes = {},'
               ' edge = {}'.format(len(start_nodes), len(graph.edges())))

    out.append('After # of nodes = {},'
               ' edge = {}'.format(len(end_nodes), len(final_graph.edges())))

    out.append('Added {} nodes and {} edges'.format(len(new_nodes),
                                                    len(new_edges)))

    out.append('Number of proteins not in network'
               ' = {}'.format(len(missing_proteins)))
    out.append('Number of add proteins = {0}'.format(len(added_proteins)))
    out.append(','.join(added_proteins))
    out.append('missing metabolites-protein info = {0}'.format(
            missing_protein_info))
    print('\n'.join(out))
    return tmp_graph


def expand_by_reactome(network, measured_list):
    """
    add _reactome functional interaction network to network
    Parameters
    ----------
    network : nx.DiGraph

    measured_list : list


    Returns
    -------

    """
    print("Checking for Reactome edges to add to network.")
    new_graph = nx.DiGraph()
    measured_set = set(measured_list)

    fi_network = load_reactome_fi()
    reactome_edges = fi_network.edges(data=True)

    # new list in reactome_nodes
    nodes_to_check = set(fi_network.nodes()).intersection(measured_set)
    nodes_to_check.update(set(network.nodes()))

    def _add_node(node):
        new_graph.add_node(node,
                           databaseSource='ReactomeFI',
                           speciesType='gene')

    for i, j, k in reactome_edges:
        if i in nodes_to_check and j in measured_list:
            _add_node(i)
            new_graph.add_edge(i, j, **k)

    new_graph = nt.compose(new_graph, network, "ReactomeExpansion")

    print("Nodes before Reactome expansion = {}, after = {}".format(
            len(network.nodes()), len(new_graph.nodes()))
    )

    print("Edges before Reactome expansion = {}, after = {}".format(
            len(network.edges()), len(new_graph.edges()))
    )

    return new_graph


def create_hmdb_network():
    out_name = os.path.join(network_data_dir, 'hmdb_graph.p')
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

    for compound, genes in cm.hmdb_accession_to_protein.items():
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

    reactome_network = load_reactome_fi()
    hmdb_network = create_hmdb_network()
    kegg_network = create_all_of_kegg()
    bgn = create_biogrid_network()
    full_network = nt.compose_all(
            [reactome_network, hmdb_network, kegg_network, bgn],
            'background')

    print("Background network {} nodes and {} edges"
          "".format(len(full_network.nodes()), len(full_network.edges()))
          )

    nx.write_gpickle(full_network, '{}.p'.format(save_name))


if __name__ == '__main__':
    create_background_network()
