from Mappings.chemical_mapper import ChemicalMapper


def expand_by_hmdb(graph, metabolite_list, all_measured):
    """ Expands a network using HMDB metabolites-protein information


    :param graph: Networkx network
    :param metabolite_list:  List of HMDB ids
    :param all_measured: list of all species measured
    :return:
    """
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
    added_proteins = []
    missing_protein_info = 0
    missing_proteins = set()
    tmp_nodes = set(tmp_graph.nodes())
    for i in metabolite_list:
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
                        added_proteins.append(each)
                    missing_edge += 1
                    tmp_graph.add_node(i, speciesType='gene',
                                       databaseSource='HMDB')
                    tmp_nodes.add(i)
                    tmp_graph.add_edge(each, i)

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
                        tmp_graph.add_node(each, speciesType='gene',
                                           databaseSource='HMDB')
                        tmp_graph.add_node(i, speciesType='metabolite',
                                           databaseSource='HMDB')
                        tmp_graph.add_edge(i, each)
                        protein_hits += 1
                    else:
                        missing_proteins.add(each)
    end_nodes = set(tmp_graph.nodes())
    end_edges = tmp_graph.edges()
    metabolites_added = []
    still_missing = []
    for i in missing_metabolites:
        if i in end_nodes:
            metabolites_added.append(i)
        else:
            still_missing.append(i)

    count = 0
    all_measured = set(all_measured)
    for each in added_proteins:
        if each in all_measured:
            count += 1

    # print(missing_proteins)
    print('Metabolites not in starting network = {0}'.format(
        len(missing_metabolites)))
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
    print(
            'missing metabolites-protein info = {0}'.format(
                missing_protein_info))

    return tmp_graph
