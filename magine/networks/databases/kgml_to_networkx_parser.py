# -*- coding: utf-8 -*-
import os
import xml.etree.cElementTree as ET
from sys import modules

import networkx as nx
from bioservices import KEGG

try:
    kegg = modules['kegg']
except:
    kegg = KEGG()
    kegg.TIMEOUT = 100


def kgml_to_graph(xmlfile, output_dir='KEGG', species='hsa'):
    """ Converts a kgml to networkx DiGraph

    Parameters
    ----------
    xmlfile : str
    output_dir : str
    species : str

    Returns
    -------

    """
    try:
        tree = ET.parse(os.path.join(output_dir, xmlfile))
    except TypeError:
        print("This file (%s)is messed up! Investigate" % str(xmlfile))
        return nx.DiGraph(), [], []
    pathway_local = nx.DiGraph()
    compounds_local = set()
    genes = set()
    connecting_maps = []
    name_label_dict = {}
    pathways_to_add_local = []
    organism = tree.getroot().get('org')
    pathway_name = tree.getroot().get('title')
    if organism != species:
        raise NotImplementedError("Didn't implement EC pathways yet")
    for entry in tree.getiterator('entry'):
        node_type = entry.get('type')
        name = entry.get('name')
        node_id = entry.get('id')
        graphics = entry.find('graphics')
        node_title = graphics.get('name')
        if node_title is None:
            connecting_maps.append(node_id)
            continue
        if node_type == 'map':
            connecting_maps.append(node_id)
            pathways_to_add_local.append(name[5:])
            continue
        elif node_type == 'gene':
            names = name.replace(' ', ",")
            name_label_dict[node_id] = names
            tmp = name_label_dict[node_id]
            for i in tmp.split(','):
                pathway_local.add_node(i, speciesType=node_type,
                                       databaseSource='KEGG', )
                genes.add(i)
        elif node_type == 'compound':
            names = name.replace(' ', ",")
            name_label_dict[node_id] = names
            tmp = name_label_dict[node_id]
            for i in tmp.split(','):
                pathway_local.add_node(i, speciesType=node_type,
                                       databaseSource='KEGG', )
                compounds_local.add(i)
        else:
            connecting_maps.append(node_id)

    # Add relations
    for rel in tree.getiterator('relation'):
        int_type = rel.get('type')
        e1 = rel.get('entry1')
        e2 = rel.get('entry2')
        if e1 in connecting_maps:
            continue
        if e2 in connecting_maps:
            continue
        try:
            type_of_interaction = ''
            for interaction in rel.getiterator('subtype'):
                type_of_interaction += interaction.get('name') + '_'
            type_of_interaction = type_of_interaction.rstrip('_')
        except TypeError:
            continue
        one, two = name_label_dict[e1], name_label_dict[e2]
        for i in one.split(','):
            for j in two.split(','):
                pathway_local.add_edge(i, j,
                                       databaseSource='KEGG',
                                       interactionType=type_of_interaction,
                                       intType=int_type)

    # Add reactions
    for reaction in tree.getiterator('reaction'):
        substrates = []
        products = []
        id_local = reaction.get('id')
        reaction_type = reaction.get('type')
        if id_local in connecting_maps:
            continue

        for sub in reaction.getiterator('substrate'):
            substrates.append(str(sub.get('name')))
        for prod in reaction.getiterator('product'):
            products.append(str(prod.get('name')))

        enzyme = name_label_dict[id_local]
        for i in enzyme.split(','):
            for sub in substrates:
                pathway_local.add_edge(sub, i, reactionType=reaction_type,
                                       interactionType='compound',
                                       databaseSource='KEGG',
                                       intType='Reaction')
            for prod in products:
                pathway_local.add_edge(i, prod, reactionType=reaction_type,
                                       interactionType='compound',
                                       databaseSource='KEGG',
                                       intType='Reaction')

    return pathway_local, pathway_name, set()
