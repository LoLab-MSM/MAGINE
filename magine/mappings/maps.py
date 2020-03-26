# -*- coding: utf-8 -*-
"""
Mapping between data identifiers
"""

import networkx as nx

from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle

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
    verbose : bool

    Returns
    -------

    """

    change_dict = dict()
    renamed_network = network.copy()
    # fix kegg drugs to their compounds (if possible)

    new_dict = drug_nodes(renamed_network)
    renamed_network = nx.relabel_nodes(renamed_network, new_dict)

    if verbose:
        print('Started converting kegg compounds to HMDB')
    hmdb, kegg_short, chem_names = cm.convert_kegg_nodes(renamed_network)
    nx.set_node_attributes(renamed_network, kegg_short, 'keggName')
    nx.set_node_attributes(renamed_network, chem_names, 'chemName')
    print(chem_names)
    nx.set_node_attributes(renamed_network, hmdb, 'hmdbNames')
    change_dict.update(hmdb)

    if verbose:
        print('Started converting kegg genes to HGNC')
    dict2, kegg_short = gm.convert_kegg_nodes(renamed_network, species=species)
    nx.set_node_attributes(renamed_network, kegg_short, 'keggName')
    change_dict.update(dict2)
    renamed_network = nx.relabel_nodes(renamed_network, change_dict)

    return renamed_network


def drug_nodes(network):
    drug_dict = {}
    hits = set(i for i in network.nodes if i.startswith('dr'))
    for i in hits:
        split_name = i.split(' ')
        if len(split_name) > 1:
            if split_name[1].startswith('cpd:'):
                drug_dict[i] = split_name[1]
                network.node[i]['drug'] = split_name[0]
    return drug_dict


def _check_dict_for_int(dic):
    new_dic = dict()
    for key, value in dic.items():
        if isinstance(value, float):
            continue
        else:
            new_dic[key] = value
    return new_dic
