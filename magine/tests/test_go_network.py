import os

import networkx as nx

from magine.networks.go_network_generator import GoNetworkGenerator

dir_name = os.path.dirname(__file__)


def test_go_network_gen():
    local_path = os.path.join(dir_name, 'Network_files', 'sample_network.gml')
    g = nx.read_gml(local_path)
    gng = GoNetworkGenerator(network=g)
    list_go = ['GO:1902175', 'GO:0006805', 'GO:0006766', 'GO:0015893',
               'GO:0006936']
    gng.create_network_from_list(list_of_go_terms=list_go, save_name='test',
                                 out_dir='out_dir_gonetwork', draw=True)
    gng.create_network_from_list(list_of_go_terms=list_go,
                                 save_name='test_merged_edge',
                                 merge_edges=True,
                                 out_dir='out_dir_gonetwork', draw=True)


def test_network_required():
    list_go = ['GO:1902175', 'GO:0006805', 'GO:0006766', 'GO:0015893',
               'GO:0006936']
    gng = GoNetworkGenerator()
    assert gng.create_network_from_list(list_of_go_terms=list_go) is None


if __name__ == '__main__':
    test_go_network_gen()
    test_network_required()
