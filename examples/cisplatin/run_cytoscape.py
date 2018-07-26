import multiprocessing

import networkx as nx

from exp_data import exp_data
from magine.enrichment import load_enrichment_csv
from magine.networks.ontology_network import create_subnetwork
from magine.networks.utils import trim_sink_source_nodes
from magine.networks.visualization.igraph_tools import paint_network_overtime, \
    render_igraph

if __name__ == '__main__':
    e_array = load_enrichment_csv('Data/cisplatin_enrichment.csv.gz',
                                  index_col=0)

    reactome_only = e_array.filter_multi(
        p_value=0.05,  # only sig pvalues
        combined_score=1.0,  # score threshold of positive values
        db='Reactome_2016',  # Only reactome db
        category='proteomics_up',  # from this category
        rank=100
    )

    reactome_only['term_name'] = reactome_only['term_name'].str.split(
        '_').str.get(0)

    # not_useful = ['metabolism', 'gene expression', 'developmental biology',
    #               'influenza infection', 'hiv infection',
    #               'late phase of hiv life cycle',
    #               'metabolism of proteins', 'disease',
    #               'infectious disease', 'immune system',
    #               'metabolism of amino acids and derivatives',
    #               'major pathway of rrna processing in the nucleolus',
    #               'influenza life cycle',
    #               'processing of capped intron-containing pre-mrna',
    #               'mrna splicing - major pathway',
    #               'mrna splicing - minor pathway',
    #               'innate immune system', 'cell-cell communication',
    #               'diseases of signal transduction', 'mrna splicing'
    #               ]
    #
    # reactome_only = reactome_only[~reactome_only['term_name'].isin(not_useful)]
    #
    # reactome_only.filter_by_minimum_sig_columns(columns='sample_id',
    #                                             min_terms=2,
    #                                             inplace=True)
    #
    # reactome_only.remove_redundant(level='sample', threshold=.5, inplace=True)
    # reactome_only.remove_redundant(level='dataframe', threshold=.5,
    #                                inplace=True)
    hits = [
        #     'cell cycle',
        'dna repair',
        'apoptosis',
        'transcriptional regulation by tp53',
        'g2/m checkpoints',
        'm phase'
    ]

    explore = ['apoptosis', 'cell cycle']
    # ASAP1,EXOC4,GBF1,RAB11FIP3
    subset_2_hits = reactome_only.loc[reactome_only['term_name'].isin(hits)]
    # print(subset_2_hits.shape)
    # print(subset_2_hits)
    # quit()
    network = nx.read_gpickle('Networks/cisplatin_network_w_attributes.p')

    # sg = Subgraph(network=network)
    #
    # new_g = sg.paths_between_two_lists(
    #     subset_2_hits.term_to_genes('cell cycle'),
    #     subset_2_hits.term_to_genes('apoptosis'),
    #     max_length=4  # allows two nodes to be inbetween
    # )
    #
    # print(len(new_g.nodes))
    # print(len(new_g.edges))
    #
    # nx.write_gml(new_g, 'expaned_subgraph2.gml')
    #
    # to_remove = set()
    # for i, j, d in network.edges(data=True):
    #     if d['databaseSource'] == 'ReactomeFI':
    #         to_remove.add((i, j))
    #
    # network.remove_edges_from(to_remove)

    ont_network, mol_net = create_subnetwork(subset_2_hits, network,
                                             save_name='ont_network_testing',
                                             merge=True,
                                             create_only=False,
                                             use_threshold=False
                                             )

    mol_net_trim = trim_sink_source_nodes(mol_net, [])
    quit()
    pool = multiprocessing.Pool(4)

    new_g = sg.paths_between_list(mol_net_trim.nodes, max_length=2, pool=pool)
    print(len(mol_net.nodes))
    print(len(mol_net.edges))

    print(len(new_g.nodes))
    print(len(new_g.edges))

    new_nodes = set(new_g.nodes)
    for i in mol_net_trim.nodes:
        if i in new_nodes:
            new_g.node[i]['terms'] = mol_net_trim.node[i]['terms']
            new_g.node[i]['termName'] = mol_net_trim.node[i]['termName']

    nx.write_gml(new_g, 'expaned_subgraph.gml')

    # exp_data.plot_species(mol_net.nodes, )
    quit()
    for i in mol_net.nodes():
        mol_net.node[i]['color'] = 'white'
    quit()

    # render_igraph(mol_net, save_name='mol_network', layout='graphopt',
    #               bbox=[1000, 1000], cluster=False)

    render_igraph(mol_net, save_name='ont_network_igraph', layout='graphopt',
                  bbox=[1000, 1000], cluster=False)

    quit()

    # net_sub = Subgraph(network)
    #
    new_g = net_sub.expand_neighbors(mol_net, nodes=mol_net.nodes(),
                                     up_stream=True, down_stream=True,
                                     include_list=exp_data.list_metabolites)
    paint_network_overtime(mol_net,
                           exp_data=exp_data,
                           save_name='metabolites',
                           color_list='red',
                           layout='auto',
                           fig_height=1000, fig_width=1000,
                           cluster=False)
