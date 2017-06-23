import os
import matplotlib
matplotlib.use("Agg")
import networkx as nx
from nose.tools import raises
from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.tests.sample_experimental_data import exp_data

_path = os.path.join(os.path.dirname(__file__), 'Network_files',
                     'sample_network.gml')


class TestNetworkSubgraphs(object):
    def setUp(self):

        self.network = nx.read_gml(_path)
        self.net_sub = NetworkSubgraphs(self.network, exp_data)

    def tearDown(self):
        self.network = None
        self.net_sub = None

    def test_downstream_nodes(self):
        """Test finding downstream nodes."""
        down_nodes = {'MTMR14', 'GABARAPL1', 'HMDB03850', 'GABARAPL2', 'EEA1',
                      'ATG16L1', 'ATG16L2', 'WIPI1', 'GABARAP', 'WIPI2',
                      'ATG12', 'ATG5', 'ZFYVE1'}
        for i in self.net_sub.downstream_network_of_specie('MTMR14').nodes():
            assert i in down_nodes

    def test_upstream_nodes(self):
        """Test finding upstream nodes."""
        up_nodes = {'MTMR4', 'MTMR14', 'MTMR3', 'ZFYVE1', 'HMDB03850'}
        for i in self.net_sub.upstream_network_of_specie('ZFYVE1').nodes():
            assert i in up_nodes
            
    @raises(RuntimeWarning)
    def test_path_between_two_does_not_exist(self):
        start = 'HSPA9'
        end = 'ZFYVE1'
        g = self.net_sub.shortest_paths_between_two_proteins(start, end)
        assert g is None

    def test_path_between_two(self):
        start = 'AKT1'
        end = 'BAX'
        nodes = {'BCL2L1', 'PTK2', 'HMDB04249', 'BAX', 'MAPK14', 'AKT1',
                 'MAPK11', 'MAPK12', 'MAPK13', 'RAF1', 'MAP2K4',
                 'PIK3CA', 'PIK3CB', 'MAP3K1', 'PIK3CD', 'PIK3R3',
                 'PIK3R2', 'PIK3R1', 'TP53', 'CASP3', 'PTEN', 'PAK1',
                 'PAK2', 'BCL2'}
        g = self.net_sub.shortest_paths_between_two_proteins(start, end)
        assert set(g.nodes()) == nodes

    def test_paint_over_time(self):
        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.shortest_paths_between_lists(list_2, draw=False,
                                                      single_path=True,
                                                      )
        colors = ['red'] * len(exp_data.timepoints)
        self.net_sub.measured_networks_over_time(g, colors, prefix='colored')

    def test_paint_over_time_up_down(self):
        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.shortest_paths_between_lists(list_2, draw=False,
                                                      single_path=True,
                                                      )

        self.net_sub.measured_networks_over_time_up_down(g, 'colored_updown')

    def test_paths_from_list(self):
        """Test finding paths from a list."""
        up_nodes = {'MTMR4', 'ZFYVE1', 'MTMR14', 'MTMR3', 'WIPI1'}
        edges = {('MTMR14', 'HMDB03850'),
                 ('HMDB03850', 'ZFYVE1'),
                 ('HMDB03850', 'WIPI1'),
                 ('MTMR4', 'HMDB03850'),
                 ('MTMR3', 'HMDB03850')}

        g = self.net_sub.shortest_paths_between_lists(up_nodes, draw=False,
                                                      save_name='ns_test')

        assert set(g.edges()) == edges

        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.shortest_paths_between_lists(list_2, draw=True,
                                                      single_path=False,
                                                      save_name='smaller_list')
        nodes = {'TP53', 'CASP3', 'MAP3K1', 'BAX', 'MAPK10', 'MAPK8', 'MAPK9',
                 'BCL2'}
        print(g.nodes())
        assert set(g.nodes()) == nodes


if __name__ == '__main__':
    t = TestNetworkSubgraphs()
    t.setUp()
    t.test_paths_from_list()

