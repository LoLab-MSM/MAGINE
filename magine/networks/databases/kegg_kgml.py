import os
import warnings
import xml.etree.cElementTree as ET
from sys import modules

import networkx as nx
from bioservices import KEGG

import magine.mappings.maps as mapper

try:
    kegg = modules['kegg']
except:
    kegg = KEGG()
    kegg.TIMEOUT = 100

try:
    import cPickle as pickle
except:
    import pickle

_out_path = os.path.join(os.path.dirname(__file__), 'KEGG')


def download_all_of_kegg(species='hsa'):
    """
    Downloads every KEGG pathway to provided directory


    Parameters
    ----------
    output_dir : str
        path to download pathways

    """
    print("Downloading KEGG")
    if not os.path.exists(_out_path):
        os.mkdir(_out_path)
    kegg.organism = species
    list_of_kegg_pathways = kegg.pathwayIds
    print("Number of pathways  = {}".format(len(list_of_kegg_pathways)))
    # keys are KEGG pathway ids, values are nodes
    pathway_dict = dict()
    for i in list_of_kegg_pathways:
        pathway = i[5:]
        file_name = kegg.get(pathway, "kgml")
        if file_name == 404:
            print("%s ended with 404 error" % pathway)
        else:
            with open(os.path.join(_out_path, '%s.xml' % pathway), 'w') as f:
                f.write(file_name)

        graph, pathway_name = kgml_to_graph('%s.xml' % pathway,
                                            _out_path, species='hsa')
        pathway_dict[i] = graph.nodes()
        print('{} has {} nodes and {} edges'.format(pathway_name,
                                                    len(graph.nodes()),
                                                    len(graph.edges())))
        output_name = os.path.join(_out_path, '{}.gml'.format(pathway))
        nx.write_gml(graph, output_name)

    # create a dictionary mapping species to pathways
    node_to_path = dict()
    for i in pathway_dict:
        for node in pathway_dict[i]:
            if node in node_to_path:
                node_to_path[node].add(i)
            else:
                node_to_path[node] = set()
                node_to_path[node].add(i)
    dirname = os.path.join(os.path.dirname(__file__), 'node_to_pathway.p')
    pickle.dump(node_to_path, open(dirname, 'wb'), )


def download_kegg_pathway(pathway, overwrite=True):
    """

    Parameters
    ----------
    pathway : str
        kegg pathway to download
    overwrite : bool
        overwrite file if is exists
    """
    download = True
    if os.path.exists("KEGG"):
        pass
    else:
        os.mkdir("KEGG")
    fname = os.path.join("KEGG", '{}.xml'.format(pathway))
    if os.path.exists(fname):
        if overwrite:
            warnings.warn('Removing KEGG pathway to download new one, '
                          'can change behavior if you set "overwrite"=False')
        else:
            print("%s already downloaded" % pathway)
            download = False
    if download:
        file_string = kegg.get(pathway, "kgml")
        if file_string == 404:
            print("{} ended with 404 error".format(pathway))
        else:
            with open(fname, 'w') as f:
                f.write(file_string)


def find_kegg_pathways(protein_names, num_overlap=1, download=True,
                       species='hsa', overwrite=True):
    """

    Parameters
    ----------
    protein_names : list
        list of gene names to map to KEGG pathways
    num_overlap : int
        number of species that must exist in multiple pathways in order
        to add to larger pathway
    download : bool
        Download the pathways
    species : str
        which species pathways to consider
    overwrite : bool
        overwrite previous kgml files if they exist

    Returns
    -------

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


def create_all_of_kegg(species='hsa'):
    dirname = os.path.join(os.path.dirname(__file__), 'node_to_pathway.p')
    tmp_dir = os.path.join(os.path.dirname(__file__), 'KEGG')
    if not os.path.exists(dirname):
        download_all_of_kegg(species=species)

    node_to_path = pickle.load(open(dirname, 'rb'))
    pathway_list = set()
    for n in node_to_path:
        if n in node_to_path:
            for j in node_to_path[n]:
                pathway_list.add(str(j.replace(':', '')[4:]))
    all_of_kegg = nx.DiGraph()
    for each in pathway_list:
        tmp = nx.read_gml(os.path.join(tmp_dir, "{}.gml".format(each)))
        all_of_kegg = nx.compose(all_of_kegg, tmp)

    drug_dict = {}
    for i in all_of_kegg.nodes():
        if i.startswith('dr'):
            split_name = i.split(' ')
            if len(split_name) > 1:
                if split_name[1].startswith('cpd:'):
                    drug_dict[i] = split_name[1]
                    all_of_kegg.node[i]['drug'] = split_name[0]
    all_of_kegg = nx.relabel_nodes(all_of_kegg, drug_dict)
    all_of_kegg = mapper.convert_all(all_of_kegg, species=species,
                                     use_hmdb=True)

    print('{} has {} nodes and {} edges'.format("all of kegg",
                                                len(all_of_kegg.nodes()),
                                                len(all_of_kegg.edges())))
    nx.write_gml(all_of_kegg, 'all_of_kegg.gml')
    return all_of_kegg


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

    return pathway_local, pathway_name
