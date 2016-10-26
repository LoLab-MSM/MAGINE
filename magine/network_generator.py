# -*- coding: utf-8 -*-
"""
File that generates networks
"""
import os
import warnings
from sys import modules

import networkx as nx
from bioservices import KEGG

import Mappings.maps as mapper
from kgml_to_networkx_parser import kgml_to_graph

try:
    kegg = modules['kegg']
except KeyError:
    kegg = KEGG()


def find_kegg_pathways(protein_names, num_overlap=1, download=True,
                       species='hsa', overwrite=True):
    """
    find_kegg_pathways(protein_names,num_overlap,download)
    Find pathways associated with a list of genes.
    :param overwrite:
    :param download: optional, downloads the kgml for each pathway
    :param num_overlap: optional, default 1
    :param protein_names: List
    :param species:
    :return results; list of pathways which contain proteins from list
    """

    pathway_list = []
    list_not_found = []
    for i in protein_names:
        tmp = kegg.get_pathway_by_gene(i, species)
        list2 = []
        if tmp is None:
            print("%s not in KEGG?" % i)
            list_not_found.append(i)
            continue
        for each in tmp:
            list2.append(str(each))
        pathway_list.append(list2)
    overlap = {}
    for i in range(len(pathway_list)):
        for each in pathway_list[i]:
            if each in overlap:
                overlap[each] += 1
            else:
                overlap[each] = 1
    results = []
    for name, count in overlap.items():
        if count > num_overlap:
            results.append(name)
    if download:
        for i in results:
            download_kegg_pathway(i, overwrite=overwrite)
    if len(list_not_found) != 0:
        print("Not found in KEGG", list_not_found)
    return results


def delete_disconnected_network(full_graph, verbose=False):
    """

    :param full_graph:
    :param verbose:
    :return:
    """
    tmp_g = full_graph.to_undirected()
    sorted_graphs = sorted(nx.connected_component_subgraphs(tmp_g), key=len,
                           reverse=True)
    for i in range(1, len(sorted_graphs)):
        for node in sorted_graphs[i].nodes():
            if verbose:
                print("Removing disconnected node %s", node)
            full_graph.remove_node(node)


def download_kegg_pathway(pathway, overwrite=True):
    """
    Downloads a pathway from kegg
    :param overwrite:
    :param pathway: string; kegg pathway (without .xml)
    :return: None
    """
    download = True
    if os.path.exists("KEGG"):
        pass
    else:
        os.mkdir("KEGG")
    if os.path.exists(os.path.join("KEGG", '%s.xml' % pathway)):
        if overwrite:
            warnings.warn('Removing pathway to download new one, '
                          'can change behavior if you set "overwrite"=False')
        else:
            print("%s already downloaded" % pathway)
            download = False
    if download:
        file_name = kegg.get(pathway, "kgml")
        if file_name == 404:
            print("%s ended with 404 error" % pathway)
        else:
            with open(os.path.join("KEGG", '%s.xml' % pathway), 'w') as f:
                f.write(file_name)


def build_network(gene_list, num_overlap=1, save_name='tmp', species='hsa',
                  overwrite=True):
    """ Construct a network from a list of gene names.
    build_network_kegg_ids(gene_list, num_overlap ,save_name )
    :param gene_list: List of genes to construct network
    :param overwrite:
    :param species: organism type
    :param save_name: Output filename, optional, default tmp
    :param num_overlap: number of hits by protein to be added to
        the network. optional, Default 1
    """

    end_network = nx.DiGraph()
    list_of_all = find_kegg_pathways(protein_names=gene_list,
                                     num_overlap=num_overlap, species=species,
                                     overwrite=overwrite)
    list_of_graphs = []
    networks_added = []
    for each in list_of_all:
        networks_added.append(each)
        tmp, pathways_to_add, compounds = kgml_to_graph("%s.xml" % each,
                                                        species=species)
        list_of_graphs.append(tmp)
    for tmp in list_of_graphs:
        for nd in tmp.nodes():
            end_network.add_node(nd, **tmp.node[nd])
        for edge in tmp.edges():
            one, two = edge[0], edge[1]
            end_network.add_edge(one, two, **tmp.edge[one][two])
    drug_dict = {}
    for i in end_network.nodes():
        if i.startswith('dr'):
            split_name = i.split(' ')
            if len(split_name) > 1:
                if split_name[1].startswith('cpd:'):
                    drug_dict[i] = split_name[1]
                    end_network.node[i]['drug'] = split_name[0]
    end_network = nx.relabel_nodes(end_network, drug_dict)
    end_network = mapper.convert_all(end_network, species=species)
    delete_disconnected_network(end_network)
    print('Number of edges', len(end_network.edges()))
    print('Number of nodes', len(end_network.nodes()))
    nx.write_gml(end_network, '%s.gml' % save_name)
    nx.nx.nx_agraph.write_dot(end_network, '%s.dot' % save_name)

    return end_network


import cPickle as pickle

directory = '/home/pinojc/PycharmProjects/Magine/Mappings'
human_kegg_to_uniprot = pickle.load(
        open(os.path.join(directory, 'human_kegg_mapper.p'), 'rb'))


def download_all_of_kegg(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    kegg.organism = 'hsa'
    list_of_kegg_pathways = kegg.pathwayIds
    print("Number of pathways  = {}".format(len(list_of_kegg_pathways)))
    for i in list_of_kegg_pathways:
        pathway = i[5:]
        file_name = kegg.get(pathway, "kgml")
        if file_name == 404:
            print("%s ended with 404 error" % pathway)
        else:
            with open(os.path.join(output_dir, '%s.xml' % pathway), 'w') as f:
                f.write(file_name)

        graph, pathway_name, _ = kgml_to_graph('%s.xml' % pathway,
                                               output_dir,
                                               species='hsa')
        print('{} has {} nodes and {} edges'.format(pathway_name,
                                                    len(graph.nodes()),
                                                    len(graph.edges())))
        output_name = os.path.join(output_dir, '{}.gml'.format(pathway))
        nx.write_gml(graph, output_name)


def create_all_of_kegg(output_dir):
    kegg.organism = 'hsa'
    from Mappings.human_mapper import GeneMapper
    gm = GeneMapper()
    list_of_kegg_pathways = kegg.pathwayIds
    all_of_kegg = nx.DiGraph()
    species_not_in_dict = set()
    for i in list_of_kegg_pathways:
        pathway = i[5:]
        output_name = os.path.join(output_dir, '{}.gml'.format(pathway))
        graph = nx.read_gml(output_name)
        all_of_kegg = nx.compose(all_of_kegg, graph)
        for j in graph.nodes():
            if not j.startswith('hsa'):
                continue
            if j in gm.kegg_to_genename:

                continue
            else:
                species_not_in_dict.add(j)

    print(len(species_not_in_dict))
    print('{} has {} nodes and {} edges'.format("all of kegg",
                                                len(all_of_kegg.nodes()),
                                                len(all_of_kegg.edges())))
    nx.write_gml(all_of_kegg, 'all_of_kegg.gml')
