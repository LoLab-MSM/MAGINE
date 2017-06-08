from sys import modules

from magine.mappings.chemical_mapper import ChemicalMapper

try:
    cm = modules['cm']
except:
    cm = ChemicalMapper()


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
