import networkx as nx

from magine.networks.go_network_generator import GoNetworkGenerator

network = nx.read_gml('cisplatin_example.gml')

go_hits = ['GO:0097193', 'GO:0072331']
gng = GoNetworkGenerator(network=network)
t_all = gng.create_network_from_list(list_of_go_terms=go_hits,
                                     save_name='{0}_network'.format(
                                         'specific_go_network'),
                                     threshold=0, draw=False)
from magine.networks.network_tools import export_to_dot

export_to_dot(t_all, 'dna-apoptosis')
