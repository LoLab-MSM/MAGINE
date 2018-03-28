import gzip
import os
import xml.etree.cElementTree as element_tree

import networkx as nx
from bioservices import KEGG

import magine.networks.network_tools as nt
from magine.data.storage import network_data_dir

try:
    import cPickle as pickle
except:
    import pickle

kegg = KEGG()
kegg.TIMEOUT = 100

_kegg_raw_out_path = os.path.join(network_data_dir, 'KEGG')


def pathway_id_to_network(pathway_id, species='hsa'):
    """
    Creates a network from a KEGG pathway id

    Parameters
    ----------
    pathway_id : str
        KEGG id
    species : str

    Returns
    -------
    nx.DiGraph

    """
    kegg.organism = species
    list_of_kegg_pathways = [i[5:] for i in kegg.pathwayIds]

    assert pathway_id in list_of_kegg_pathways
    pathway = kegg.get(pathway_id, "kgml")
    graph, pathway_name = kgml_to_nx(pathway, species=species)
    return graph


def kgml_to_nx(xmlfile, species='hsa'):
    """ Converts a kgml to nx.DiGraph

    Parameters
    ----------
    xmlfile : str
    species : str

    Returns
    -------

    """
    try:
        tree = element_tree.ElementTree(element_tree.fromstring(xmlfile))
    except TypeError:
        print("This file (%s)is messed up! Investigate" % str(xmlfile))
        return nx.DiGraph(), []
    pathway_local = nx.DiGraph()
    connecting_maps = []
    name_label_dict = {}
    organism = tree.getroot().get('org')
    pathway_name = tree.getroot().get('title')

    def _add_node(node, type_species):
        pathway_local.add_node(node, speciesType=type_species,
                               databaseSource='KEGG')

    def _add_edge(source, target, type_int):
        pathway_local.add_edge(source, target,
                               databaseSource='KEGG',
                               interactionType=type_int,
                               )

    if organism != species:
        print(organism, species)
        raise NotImplementedError("Error with species not matching from KEGG")

    # get all nodes of pathway
    for entry in tree.iter('entry'):

        # get node id
        node_id = entry.get('id')
        node_type = entry.get('type')
        name = entry.get('name')
        if node_type in ('gene', 'compound'):
            names = name.split()
            name_label_dict[node_id] = names
            for i in names:
                _add_node(i, node_type)
        else:
            # add everything that is not a gene or compound to exclude
            # orthologs are from other species
            # groups is just a trick for KGML to layout species together
            # maps point to other dbs
            if node_type not in ('map', 'ortholog', 'group'):
                print("Not a gene or compound!!!")
                print(node_type, node_id)
            connecting_maps.append(node_id)

    # Add all relations of pathway
    for rel in tree.iter('relation'):
        e1 = rel.get('entry1')
        e2 = rel.get('entry2')
        if e1 in connecting_maps or e2 in connecting_maps:
            continue
        try:
            int_type = set()
            for interaction in rel.getiterator('subtype'):
                int_type.add(interaction.get('name'))
            int_type = '|'.join(int_type)
        except TypeError:
            continue

        one, two = name_label_dict[e1], name_label_dict[e2]
        for i in one:
            for j in two:
                _add_edge(i, j, int_type)

    # Add all reactions of pathway
    for reaction in tree.iter('reaction'):

        id_local = reaction.get('id')
        if id_local in connecting_maps:
            print(id_local)
            continue

        reaction_type = "chemical|reaction|{}".format(reaction.get('type'))
        subs = [str(i.get('name')) for i in reaction.getiterator('substrate')]
        prods = [str(i.get('name')) for i in reaction.getiterator('product')]
        enzyme = name_label_dict[id_local]

        for i in enzyme:
            for sub in subs:
                _add_edge(sub, i, type_int=reaction_type)

            for prod in prods:
                _add_edge(i, prod, type_int=reaction_type)

    return pathway_local, pathway_name


def download_all_of_kegg(species='hsa', verbose=False):
    """
    Downloads every KEGG pathway to provided directory

    """

    if not os.path.exists(_kegg_raw_out_path):
        os.mkdir(_kegg_raw_out_path)

    kegg.organism = species
    list_of_kegg_pathways = [i[5:] for i in kegg.pathwayIds]
    if verbose:
        print("Number of pathways  = {}".format(len(list_of_kegg_pathways)))

    # keys are KEGG pathways ids, values are kegg pathways
    kegg_dict = dict()

    for pathway_id in list_of_kegg_pathways:
        pathway = kegg.get(pathway_id, "kgml")

        if pathway == 404:
            print("%s ended with 404 error" % pathway_id)
            continue

        graph, pathway_name = kgml_to_graph(pathway, species=species)

        kegg_dict[pathway_id] = graph
        if verbose:
            print('{} has {} nodes and {} edges'.format(pathway_name,
                                                        len(graph.nodes()),
                                                        len(graph.edges())))

    # create a dictionary mapping species to pathways
    node_to_path = dict()
    for i in kegg_dict:
        for node in kegg_dict[i].nodes():
            if node in node_to_path:
                node_to_path[node].add(i)
            else:
                node_to_path[node] = set()
                node_to_path[node].add(i)

    # save kegg pathway its to networks
    n = '{}_kegg_path_ids_to_networks.p.gz'.format(species)
    save_path_id_to_graph = os.path.join(network_data_dir, n)
    save_gzip_pickle(save_path_id_to_graph, kegg_dict)

    # save nodes to pathways
    n = '{}_kegg_node_to_pathway.p.gz'.format(species)
    save_node_to_path = os.path.join(network_data_dir, n)
    save_gzip_pickle(save_node_to_path, node_to_path)


def create_all_of_kegg(species='hsa', verbose=False):
    """
    Combines all KEGG pathways into a single network

    Parameters
    ----------
    species : species
        Default 'hsa'
    verbose : bool

    Returns
    -------

    """
    p_name = os.path.join(_kegg_raw_out_path,
                          '{}_all_of_kegg.p.gz'.format(species))
    # load in network if it exists
    if os.path.exists(p_name):
        if verbose:
            print('Reading in KEGG network')
        all_of_kegg = nx.read_gpickle(p_name)
    else:
        # create the network
        path_to_graph, _ = load_kegg(species, verbose)

        # merge into single networks
        all_of_kegg = nt.compose_all(path_to_graph.values())

        # relabel notes
        if species == 'hsa':
            from magine.mappings.maps import convert_all
            all_of_kegg = convert_all(all_of_kegg, species=species)

        # save network
        nx.write_gpickle(all_of_kegg, p_name)

    if verbose:
        n_n = len(all_of_kegg.nodes())
        n_e = len(all_of_kegg.edges())
        print("KEGG network {} nodes and {} edges".format(n_n, n_e))
    return all_of_kegg


def load_kegg(species, verbose=False):
    n = '{}_kegg_path_ids_to_networks.p.gz'.format(species)
    save_path_id_to_graph = os.path.join(network_data_dir, n)
    if not os.path.exists(save_path_id_to_graph):
        download_all_of_kegg(species=species, verbose=verbose)

    n = '{}_kegg_node_to_pathway.p.gz'.format(species)
    save_node_to_path = os.path.join(network_data_dir, n)

    # load all networks
    return load_gz_p(save_path_id_to_graph), load_gz_p(save_node_to_path)


def save_gzip_pickle(file_name, obj):
    with gzip.open(file_name, 'wb') as f:
        f.write(pickle.dumps(obj, protocol=-1))


def load_gz_p(file_name):
    with gzip.open(file_name, 'rb') as f:
        data = f.read()
    try:
        return pickle.loads(data, encoding='utf-8')
    except:
        return pickle.loads(data)


if __name__ == '__main__':
    download_all_of_kegg('hsa')
    create_all_of_kegg('hsa')
