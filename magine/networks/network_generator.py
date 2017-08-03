# -*- coding: utf-8 -*-
"""
File that generates networks
"""
import os
from sys import modules

import networkx as nx
import numpy as np


import magine.mappings.maps as mapper
from magine.data.storage import network_data_dir
from magine.mappings.gene_mapper import GeneMapper
from magine.networks.databases import download_all_of_kegg, load_reactome_fi

from magine.networks.network_tools import delete_disconnected_network

try:
    import cPickle as pickle
except:
    import pickle



def build_network(gene_list, num_overlap=1, save_name='tmp', species='hsa',
                  overwrite=True, all_measured_list=None, use_reactome=True,
                  use_hmdb=False, metabolite_list=None):
    """
    Construct a network from a list of gene names.

    Parameters
    ----------
    gene_list : list
        list of genes to construct network
    num_overlap : int
        number of KEGG pathways that a species must be in to add pathway to
        network
    save_name : str
        output name to save network. Will save one before and after ID
        conversion
    species : str
        species of proteins ('hsa': human, 'mmu':murine)
    overwrite : bool
        remove old KEGG xml files and download new ones
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
    # gm.load()
    _kegg_raw_out_path = os.path.join(network_data_dir, 'KEGG')
    _kegg_node_to_pathway = os.path.join(network_data_dir, 'kegg_node_to_pathway.p')

    if not os.path.exists(_kegg_raw_out_path):
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
        tmp = nx.read_gml(os.path.join(_kegg_raw_out_path, "{}.gml".format(each)))
        for n in tmp.nodes():
            if isinstance(n, float):
                print("Found the float... {}".format(each))
                tmp.remove_node(np.nan)
        if len(tmp.edges()) == 0:
            continue
        end_network = nx.compose(end_network, tmp)

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

    end_network = mapper.convert_all(end_network, species=species,
                                     use_hmdb=use_hmdb)

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
        else:
            end_network = expand_by_hmdb(end_network, metabolite_list,
                                         all_measured_list)
    print("Trimming network")
    delete_disconnected_network(end_network)
    print('Number of nodes {}'.format(len(end_network.nodes())))
    print('Number of edges {}'.format(len(end_network.edges())))
    nx.write_gml(end_network, '{}.gml'.format(save_name))
    nx.write_gpickle(end_network, '{}.p'.format(save_name))

    return end_network


def expand_by_hmdb(graph, metabolite_list, all_measured):
    """
    Expands a network using HMDB metabolites-protein information

    Parameters
    ----------
    graph : nx.DiGraph
    metabolite_list : list
        List of HMDB ids
    all_measured : list
        list of all species measured, including proteins
    Returns
    -------
    nx.DiGraph

    """
    from magine.mappings.chemical_mapper import ChemicalMapper

    try:
        cm = modules['cm']
    except:
        cm = ChemicalMapper()
    tmp_graph = graph.copy()
    start_nodes = set(tmp_graph.nodes())
    start_edges = tmp_graph.edges()
    missing_metabolites = []
    for i in metabolite_list:
        if i not in start_nodes:
            if i.startswith('HMDB'):
                missing_metabolites.append(i)
            # else:
            #     print('Not an HMDB : {}'.format(i))

    count_in_network, count_not_in_network = 0, 0
    missing_edge = 0
    protein_hits = 0
    added_proteins = set()
    missing_protein_info = 0
    missing_proteins = set()
    tmp_nodes = set(tmp_graph.nodes())
    metabolite_set = set(metabolite_list)
    for i in metabolite_set:
        if i == np.nan:
            continue
        # checks if metabolite is in the network
        # if it is, it can add an associated gene
        if i in tmp_nodes:
            count_in_network += 1

            if i in cm.hmdb_accession_to_protein:
                if not isinstance(i, str):
                    print("Found and integer in hmdb")
                    print(type(i), i)
                tmp_list = cm.hmdb_accession_to_protein[i]
                if tmp_list is None:
                    missing_protein_info += 1
                    continue

                for each in tmp_list:
                    for each in tmp_list:
                        if not isinstance(each, str):
                            print("Found and integer in hmdb")
                            print(type(each), each)
                    if each is None:
                        continue
                    elif each not in tmp_nodes:
                        added_proteins.add(each)
                    missing_edge += 1
                    tmp_graph.add_node(i, speciesType='compound',
                                       databaseSource='HMDB')
                    tmp_nodes.add(i)
                    tmp_graph.add_edge(each, i, interactionType='chemical',
                                       databaseSource='HMDB')
        # if the metabolite is NOT in the network,
        # it checks to see if it has any protein relationships
        # if it does and they are new to the graph, it adds it
        else:
            count_not_in_network += 1
            if i in cm.hmdb_accession_to_protein:
                if not isinstance(i, str):
                    print("Found and integer in hmdb")
                    print(type(i), i)
                tmp_list = cm.hmdb_accession_to_protein[i]
                if tmp_list is None or len(tmp_list) == 0:
                    missing_protein_info += 1
                    continue

                for each in tmp_list:
                    if not isinstance(each, str):
                        print("Found and integer in hmdb")
                        print(type(each), each)
                    if each is None:
                        pass
                    elif each in start_nodes:
                        missing_edge += 1
                        tmp_graph.add_node(each,
                                           speciesType='gene',
                                           databaseSource='HMDB')

                        tmp_graph.add_node(i,
                                           speciesType='compound',
                                           databaseSource='HMDB')

                        tmp_graph.add_edge(i, each, interactionType='chemical',
                                           databaseSource='HMDB')
                        protein_hits += 1
                    else:
                        missing_proteins.add(each)
    end_nodes = set(tmp_graph.nodes())
    end_edges = tmp_graph.edges()
    metabolites_added = set()
    still_missing = set()
    for i in missing_metabolites:
        if i in end_nodes:
            metabolites_added.add(i)
        else:
            still_missing.add(i)

    count = 0
    all_measured = set(all_measured)
    for each in added_proteins:
        if each in all_measured:
            count += 1

    # print(missing_proteins)
    print('Metabolites not in starting network = {0}'.format(len(missing_metabolites)))
    print('Metabolites added to network = {0}'.format(len(metabolites_added)))
    print('Metabolites still not in network = {0}'.format(len(still_missing)))
    print('\n')
    print('Before number of nodes = {0}'.format(len(start_nodes)))
    print('Before number of edges = {0}'.format(len(start_edges)))
    print('\n')
    print('After number of nodes = {0}'.format(len(end_nodes)))
    print('After number of edges = {0}'.format(len(end_edges)))
    print('\n')
    print('Added nodes = {0}'.format(len(end_nodes) - len(start_nodes)))
    print('Added edges = {0}'.format(len(end_edges) - len(start_edges)))
    print('\n')
    print('Number of proteins not in KEGG = {0}'.format(len(missing_proteins)))
    print('Number of add proteins = {0}'.format(len(added_proteins)))
    print(added_proteins)
    print('Proteins addeddow that were measured = {0}'.format(count))
    print('\n')
    print('missing metabolites-protein info = {0}'.format(
        missing_protein_info))

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
    combined_network = network.copy()
    existing_nodes = set(network.nodes())
    existing_edges = set(network.edges())
    measured_set = set(measured_list)

    fi_network = load_reactome_fi()
    reactome_nodes = set(fi_network.nodes())
    reactome_edges = fi_network.edges(data=True)

    added_nodes = 0
    added_edges = 0

    nodes_to_check = set()
    for i in measured_set:
        if i not in existing_nodes:
            if i in reactome_nodes:
                nodes_to_check.add(i)
                added_nodes += 1

    for i, j, k in reactome_edges:
        if not isinstance(i, str):
            print("Found nan in reactome")
            continue
        if not isinstance(j, str):
            print("Found nan in reactome")
            continue

        if (i, j) in existing_edges:
            continue
        if i in measured_set and j in measured_set:
            combined_network.add_edge(i, j, k)
            added_edges += 1
        elif i in nodes_to_check and j in measured_list:
            combined_network.add_node(i, databaseSource='ReactomeFI',
                                      speciesType='gene')
            combined_network.add_edge(i, j, attr_dict=k)
            added_edges += 1
        elif j in nodes_to_check and i in measured_list:
            combined_network.add_node(j, databaseSource='ReactomeFI',
                                      speciesType='gene')
            combined_network.add_edge(i, j, attr_dict=k)
            added_edges += 1
    print("Found {} nodes in REACTOME not in KEGG".format(added_nodes))
    print("Added {} edges from REACTOME".format(added_edges))

    print("Nodes before Reactome expansion "
          "= {}, after = {}".format(len(existing_nodes),
                                    len(combined_network.nodes())))

    return combined_network


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
    import magine.networks.network_tools as nt

    reactome_network = load_reactome_fi()
    hmdb_network = create_hmdb_network()
    kegg_network = create_all_of_kegg()
    bgn = create_biogrid_network()
    full_network = nt.compose_all(
            [reactome_network, hmdb_network, kegg_network, bgn],
            'background')

    print("Background network "
          "{} nodes and {} edges".format(
            len(full_network.nodes()),
            len(full_network.edges()))
          )

    nx.write_gpickle(full_network, '{}.p'.format(save_name))
    # nx.write_gml(full_network, '{}.gml'.format(save_name))


def create_hmdb_network():

    out_name = os.path.join(network_data_dir, 'hmdb_graph.p')
    if os.path.exists(out_name):
        return nx.read_gpickle(out_name)
    from magine.mappings.chemical_mapper import ChemicalMapper

    try:
        cm = modules['cm']
    except:
        cm = ChemicalMapper()

    tmp_graph = nx.DiGraph()
    nodes = set()

    def _add_node(node, node_type):
        if node != u'':
            if node not in nodes:
                tmp_graph.add_node(node, speciesType=node_type,
                                   databaseSource='HMDB')
                nodes.add(node)

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

    nx.write_gpickle(tmp_graph, out_name)
    out_name = os.path.join(network_data_dir, 'hmdb_graph.gml')
    nx.write_gml(tmp_graph, out_name)
    return tmp_graph


if __name__ == '__main__':
    create_background_network()
    # import time
    # st = time.time()
    # g = nx.read_gpickle('background_network.p')
    # et = time.time()
    # print(et-st)
    # st = time.time()
    # g = nx.read_gml('background_network.gml')
    # et = time.time()
    # print(et - st)

