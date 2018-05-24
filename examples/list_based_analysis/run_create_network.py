from magine.networks.network_generator import build_network
from magine.networks.subgraphs import Subgraph

# add an example passing a file with genes and a save name


list_of_proteins = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR',
                    'BAD', 'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX',
                    'BAK1', 'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']

network = build_network(list_of_proteins,
                        save_name='example_earm_network',
                        use_reactome=True)

print("Network has {} nodes".format(network.number_of_nodes()))
print("Network has {} edgess".format(network.number_of_edges()))

#  There are many edge types described by KEGG and Reactome
#  We have the option to remove the ones that have the least confidence
#  These include "predicted" and "?" (Reactome edges)
#  and "indirect effect" (KEGG edges)

to_remove = set()
for i, j, data in network.edges(data=True):

    if 'predicted' in data:
        if data['predicted'] == 'True':
            to_remove.add((i, j))
    elif 'indirect effect' in data['interactionType']:
        to_remove.add((i, j))
    elif 'missing interaction' in data['interactionType']:
        to_remove.add((i, j))
    elif '?' == data['interactionType']:
        to_remove.add((i, j))
for i in to_remove:
    network.remove_edge(i[0], i[1])

print("Removing edges with less detail")
print("Network has {} nodes".format(network.number_of_nodes()))
print("Network has {} edgess".format(network.number_of_edges()))

# The network should contain roughly ~2700 species and ~12K nodes
# We do this to expand all relationships that could possible be found


# Now we are going to concentrate back into the species we started with.
# Now we are going to have all connections between the species


x = Subgraph(network)
sub = x.paths_between_list(list_of_proteins,
                           single_path=False,
                           save_name='trimmed_network',
                           draw=True,
                           image_format='pdf')
print("Network has {} nodes".format(sub.number_of_nodes()))
print("Network has {} edgess".format(sub.number_of_edges()))
