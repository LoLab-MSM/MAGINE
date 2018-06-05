# -*- coding: utf-8 -*-
"""
Mapping between data identifiers
"""

import networkx as nx
from bioservices import HGNC, KEGG, UniChem, UniProt

from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle

kegg = KEGG()
uniprot = UniProt()
hugo = HGNC()
chem = UniChem()
gm = GeneMapper()
cm = ChemicalMapper()


def convert_all(network, species='hsa', verbose=False):
    """
    Maps gene names to HGNC and kegg compound to HMDB

    Parameters
    ----------
    network : networkx.DiGraph
        network to convert mappings
    species : str
        species of network (hsa, mmu)
    verbose : str

    Returns
    -------

    """

    change_dict = dict()
    renamed_network = network.copy()
    if verbose:
        print('Started converting kegg compounds to HMDB')
    dict1 = cm.convert_kegg_nodes(renamed_network)
    change_dict.update(dict1)
    if verbose:
        print('Started converting kegg genes to HGNC')
    dict2, found_all = gm.convert_kegg_nodes(renamed_network, species=species)
    change_dict.update(dict2)

    if not found_all:
        if verbose:
            print('Started to check for miRNAs')
        dict3 = hugo_mapper(renamed_network, species=species)
        change_dict.update(dict3)
    change_dict = _check_dict_for_int(change_dict)

    renamed_network = nx.relabel_nodes(renamed_network, change_dict)
    renamed_network = drug_nodes(renamed_network)
    return renamed_network


def hugo_mapper(network, species='hsa'):
    """
    Converts all KEGG names to HGNC 
    
    Parameters
    ----------
    network : networkx.DiGraph
    species : str

    Returns
    -------
    dict
    """
    prefix = species+':'
    nodes = set(network.nodes)
    hugo_dict = {}
    not_found = set()
    for i in nodes:
        if str(i).startswith(prefix):
            tmp_name = str(i).replace(prefix, '')
            mapping = hugo.search(tmp_name)
            if 'response' in mapping:
                response = mapping['response']
                if 'numFound' in response:
                    if response['numFound'] == 0:
                        not_found.add(i)
                        continue
                    elif response['numFound'] == 1:
                        docs = response['docs'][0]
                        hugo_dict[i] = docs['symbol']
                        continue
                    else:
                        if 'symbol' in response['docs'][0]:
                            hugo_dict[i] = response['docs'][0]['symbol']
            else:
                not_found.add(i)
    if not_found != 0:
        print("{} mappings not found after HGNC mapping".format(len(not_found)))
        print("{} ".format(not_found))
    return hugo_dict


def drug_nodes(network):
    drug_dict = {}
    for i in network.nodes():
        if i.startswith('dr'):
            split_name = i.split(' ')
            if len(split_name) > 1:
                if split_name[1].startswith('cpd:'):
                    drug_dict[i] = split_name[1]
                    network.node[i]['drug'] = split_name[0]
        elif i == 'nan':
            network.remove_node(i)
        elif isinstance(i, float):
            network.remove_node(i)
    end_network = nx.relabel_nodes(network, drug_dict)
    return end_network


def _check_dict_for_int(dic):
    new_dic = dict()
    for key, value in dic.items():
        if isinstance(value, float):
            continue
        else:
            new_dic[key] = value
    return new_dic
