import networkx as nx

from magine.networks.network_generator import build_network
from magine.networks.network_subgraphs import NetworkSubgraphs

list_of_proteins = ['CASP3', 'CASP6',
                    'FAS', 'FADD', 'CASP8',  # DISC
                    'CFLAR',  # FLIP
                    'BFAR',  # BAR
                    'BAD',
                    # pro-apoptotic
                    'BID',
                    'PMAIP1',  # NOXA
                    'MCL1', 'BCL2', 'BCL2L1',  # anti-apoptotic
                    'BAX', 'BAK1',  # effector proteins
                    'DIABLO', 'CYCS',
                    'PARP1', 'APAF1', 'XIAP',
                    ]

network = build_network(list_of_proteins, save_name='earm_network',
                        use_reactome=True,
                        overwrite=False)

to_remove = set()
for i, j, data in network.edges(data=True):

    if 'predicted' in data:
        if data['predicted'] == 'True':
            to_remove.add((i, j))
    elif 'indirect effect' in data['interactionType']:
        to_remove.add((i, j))
    elif '?' == data['interactionType']:
        to_remove.add((i, j))
for i in to_remove:
    # print(i)
    network.remove_edge(i[0], i[1])

# The network should contain roughly ~2700 species and ~12K nodes
# We do this to expand all relationships that could possible be found

# Now we are going to concentrate back into the species we started with.
# Now we are going to have all connections between the species


x = NetworkSubgraphs(network)
sub = x.shortest_paths_between_lists(list_of_proteins,
                                     single_path=True,
                                     save_name='trimmed_network',
                                     draw=True)
print("Network has {} nodes".format(sub.number_of_nodes()))
print("Network has {} edgess".format(sub.number_of_edges()))
nx.write_gml(sub, 'trimmed.gml')
