# -*- coding: utf-8 -*-
"""
File that generates networks
"""
import os
import warnings
from sys import modules

import networkx as nx
from bioservices import KEGG

import magine.mappings.maps as mapper
from magine.mappings.gene_mapper import GeneMapper
from magine.network_database_expansions.reactome \
    import ReactomeFunctionalInteraction
from magine.networks.kgml_to_networkx_parser import kgml_to_graph

try:
    import cPickle as pickle
except:
    import pickle

try:
    kegg = modules['kegg']
except KeyError:
    kegg = KEGG()


dirname = os.path.join(os.path.dirname(__file__), 'node_to_pathway.p')


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


def delete_disconnected_network(full_graph, verbose=False):
    """
    Delete disconnected parts of a provided network
    
    Parameters
    ----------
    full_graph : networkx.DiGraph
    verbose : bool
        Prints all nodes that are removed
    
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
    tmp_dir = os.path.join(os.path.dirname(__file__), 'TMP_KEGG')
    if not os.path.exists(dirname):
        _download_all_of_kegg(tmp_dir)

    node_to_path = pickle.load(open(dirname, 'rb'))

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
    # pathway_list2 = find_kegg_pathways(to_remove, download=False,
    #                                    species=species)
    # pathway_list.update(pathway_list2)
    # list_of_all = find_kegg_pathways(protein_names=gene_set,
    #                                  num_overlap=num_overlap, species=species,
    #                                  overwrite=overwrite)
    list_of_graphs = []
    networks_added = []

    for each in pathway_list:
        networks_added.append(each)
        # tmp, pathways_to_add, compounds = kgml_to_graph("{}.xml".format(each),
        #                                                 species=species)
        tmp = nx.read_gml(os.path.join(tmp_dir, "{}.gml".format(each)))
        end_network = nx.compose(end_network, tmp)
        list_of_graphs.append(tmp)
    # for tmp in list_of_graphs:
    #     for nd in tmp.nodes():
    #         end_network.add_node(nd, **tmp.node[nd])
    #     for edge in tmp.edges():
    #         one, two = edge[0], edge[1]
    #         end_network.add_edge(one, two, **tmp.edge[one][two])
    drug_dict = {}
    for i in end_network.nodes():
        if i.startswith('dr'):
            split_name = i.split(' ')
            if len(split_name) > 1:
                if split_name[1].startswith('cpd:'):
                    drug_dict[i] = split_name[1]
                    end_network.node[i]['drug'] = split_name[0]
    end_network = nx.relabel_nodes(end_network, drug_dict)
    end_network = mapper.convert_all(end_network, species=species,
                                     use_hmdb=use_hmdb)

    if all_measured_list is not None:
        all_measured_list = set(x.upper() for x in all_measured_list)
        end_network = add_reactome(end_network, all_measured_list)
    elif use_reactome:
        end_network = add_reactome(end_network, gene_list)

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
    # nx.nx.nx_agraph.write_dot(end_network, '{}.dot'.format(save_name))

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
            else:
                print('Not an HMDB : {}'.format(i))

    count_in_network, count_not_in_network = 0, 0
    missing_edge = 0
    protein_hits = 0
    added_proteins = set()
    missing_protein_info = 0
    missing_proteins = set()
    tmp_nodes = set(tmp_graph.nodes())
    metabolite_set = set(metabolite_list)
    for i in metabolite_set:
        # checks if metabolite is in the network
        # if it is, it can add an associated gene
        if i in tmp_nodes:
            count_in_network += 1

            if i in cm.hmdb_accession_to_protein:

                tmp_list = cm.hmdb_accession_to_protein[i]
                if tmp_list is None:
                    missing_protein_info += 1
                    continue

                for each in tmp_list:
                    if each is None:
                        continue
                    elif each not in tmp_nodes:
                        added_proteins.add(each)
                    missing_edge += 1
                    tmp_graph.add_node(i, speciesType='gene',
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
                tmp_list = cm.hmdb_accession_to_protein[i]
                if tmp_list is None or len(tmp_list) == 0:
                    missing_protein_info += 1
                    continue

                for each in tmp_list:
                    if each is None:
                        pass
                    elif each in start_nodes:
                        missing_edge += 1
                        tmp_graph.add_node(each,
                                           speciesType='gene',
                                           databaseSource='HMDB')

                        tmp_graph.add_node(i,
                                           speciesType='metabolite',
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
    print('Proteins added that were measured = {0}'.format(count))
    print('\n')
    print('missing metabolites-protein info = {0}'.format(
        missing_protein_info))

    return tmp_graph


def add_reactome(network, measured_list):
    """
    add _reactome functional interaction network to network
    Parameters
    ----------
    network : nx.DiGraph

    measured_list : list


    Returns
    -------

    """
    combined_network = network.copy()
    existing_nodes = set(network.nodes())
    existing_edges = set(network.edges())
    measured_set = set(measured_list)

    fi_network = load_reactome_fi()
    reactome_nodes = set(fi_network.nodes())
    reactome_edges = fi_network.edges(data=True)

    added_nodes = 0
    added_edges = 0
    measured_species_in_network = set()
    nodes_to_check = set()
    for i in measured_set:
        if i not in existing_nodes:
            if i in reactome_nodes:
                nodes_to_check.add(i)
                added_nodes += 1

    # for i in reactome_nodes:
    #     if i in measured_set:
    #         nodes_to_check.add(i)

    # for i in fi_network.nodes(data=True):
    #     if i[0] in existing_nodes:
    #         continue
    #     elif i in measured_list:
    #         measured_species_in_network.add(i)
    #         nodes_to_check.add(i)



    # for node_1, node_2 in itertools.combinations(nodes_to_check, 2):
    #     if nx.has_path(fi_network, node_1, node_2):
    #         for path in nx.all_shortest_paths(fi_network, node_1, node_2):
    #             _add_edges_from_path(fi_network, combined_network, path)
    #     if nx.has_path(fi_network, node_2, node_1):
    #         for path in nx.all_shortest_paths(fi_network, node_2, node_1):
    #             _add_edges_from_path(fi_network, combined_network, path)

    for i, j, k in reactome_edges:
        if (i, j) in existing_edges:
            continue
        if i in measured_set and j in measured_set:
            combined_network.add_edge(i, j, k)
            added_edges += 1
        elif i in nodes_to_check or j in measured_list:
            combined_network.add_edge(i, j, k)
            added_edges += 1

    print("Found {} nodes in REACTOME not in KEGG".format(added_nodes))
    print("Added {} edges from REACTOME".format(added_edges))

    print("Nodes before Reactome expansion "
          "= {}, after = {}".format(len(existing_nodes),
                                    len(combined_network.nodes())))

    return combined_network


def _add_edges_from_path(network, graph, path):
    """
    Adds a path to a graph

    Parameters
    ----------
    graph : networkx.DiGraph
        graph to add paths
    path : list_like
        list of species that create a path

    """
    previous = None
    for protein in list(path):
        if previous is None:
            previous = protein
            continue
        else:
            graph.add_node(previous, **network.node[previous])
            graph.add_node(protein, **network.node[protein])
            graph.add_edge(previous, protein,
                           **network.edge[previous][protein])
            previous = protein


def _download_all_of_kegg(output_dir, species='hsa'):
    """
    Downloads every KEGG pathway to provided directory
    
    
    Parameters
    ----------
    output_dir : str
        path to download pathways

    """
    print("Downloading KEGG")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
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
            with open(os.path.join(output_dir, '%s.xml' % pathway), 'w') as f:
                f.write(file_name)

        graph, pathway_name, _ = kgml_to_graph('%s.xml' % pathway,
                                               output_dir,
                                               species='hsa')
        pathway_dict[i] = graph.nodes()
        print('{} has {} nodes and {} edges'.format(pathway_name,
                                                    len(graph.nodes()),
                                                    len(graph.edges())))
        output_name = os.path.join(output_dir, '{}.gml'.format(pathway))
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
    pickle.dump(node_to_path, open(dirname, 'wb'), )


def _create_all_of_kegg(output_dir, species='hsa'):
    kegg.organism = species
    from magine.mappings.gene_mapper import GeneMapper
    from magine.mappings.chemical_mapper import ChemicalMapper
    gm = GeneMapper()
    cm = ChemicalMapper()
    list_of_kegg_pathways = kegg.pathwayIds
    all_of_kegg = nx.DiGraph()
    species_not_in_dict = set()
    name_dict = {}
    for i in list_of_kegg_pathways:
        pathway = i[5:]
        output_name = os.path.join(output_dir, '{}.gml'.format(pathway))
        graph = nx.read_gml(output_name)
        all_of_kegg = nx.compose(all_of_kegg, graph)
        for j in graph.nodes():
            if j.startswith(species):
                if j in gm.kegg_to_gene_name:
                    name_dict[j] = str(gm.kegg_to_gene_name[j][0])
                    continue
                else:
                    species_not_in_dict.add(j)
            elif j.startswith('cpd'):
                if j in cm.kegg_to_hmdb_accession:
                    name_dict[j] = str(cm.kegg_to_hmdb_accession[j][0])
                    continue
                else:
                    species_not_in_dict.add(j)
            elif j.startswith('dr'):
                split_name = j.split(' ')
                if len(split_name) > 1:
                    if split_name[1].startswith('cpd:'):
                        name_dict[j] = split_name[1]
                        graph.node[j]['drug'] = str(split_name[0])
                    else:
                        species_not_in_dict.add(j)
                else:
                    species_not_in_dict.add(j)
            else:
                species_not_in_dict.add(j)

    print(len(species_not_in_dict))
    print('{} has {} nodes and {} edges'.format("all of kegg",
                                                len(all_of_kegg.nodes()),
                                                len(all_of_kegg.edges())))
    all_of_kegg = nx.relabel_nodes(all_of_kegg, name_dict)
    nx.write_gml(all_of_kegg, 'all_of_kegg.gml')



def load_reactome_fi():
    """
    Load reactome functional interaction network
    Returns
    -------

    """
    path = os.path.join(os.path.dirname(__file__), '..',
                        'network_database_expansions'
                        'reactome',
                        'reactome_fi.gml')
    if not os.path.exists(path):
        print("Downloading Reactome Functional interaction network!")
        rfi = ReactomeFunctionalInteraction()
        rfi.parse_network()
    if not os.path.exists(path):
        print('Failed!')
        quit()
    return nx.read_gml(path)

if __name__ == '__main__':
    _download_all_of_kegg('DELETE')
    # _create_all_of_kegg('DELETE')
