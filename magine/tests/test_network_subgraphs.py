import os

import networkx as nx
from nose.tools import raises

from magine.networks.network_subgraphs import NetworkSubgraphs

_path = os.path.join(os.path.dirname(__file__), 'Network_files',
                     'sample_network.gml')

from magine.tests.sample_experimental_data import exp_data


class TestSolver(object):
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

    def test_paths_from_list(self):
        """Test finding paths from a list."""
        up_nodes = {'MTMR4', 'ZFYVE1', 'MTMR14', 'MTMR3', 'WIPI1'}
        edges = {('MTMR14', u'HMDB03850'),
                 (u'HMDB03850', 'ZFYVE1'),
                 (u'HMDB03850', 'WIPI1'),
                 ('MTMR4', u'HMDB03850'),
                 ('MTMR3', u'HMDB03850')}
        g = self.net_sub.shortest_paths_between_lists(up_nodes, draw=True,
                                                      save_name='ns_test')
        assert set(g.edges()) == edges

        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.shortest_paths_between_lists(list_2, draw=True,
                                                      single_path=True,
                                                      save_name='smaller_list')
        nodes = {'TP53', 'CASP3', u'MAP3K1', u'BAX', u'MAPK10', u'BCL2'}
        assert set(g.nodes()) == nodes

    @raises(RuntimeWarning)
    def test_path_between_two_does_not_exist(self):
        start = 'HSPA9'
        end = 'ZFYVE1'
        g = self.net_sub.shortest_paths_between_two_proteins(start, end)
        assert g is None

    def test_path_between_two(self):
        start = 'AKT1'
        end = 'BAX'
        nodes = {u'BCL2L1', u'PTK2', u'HMDB04249', 'BAX', u'MAPK14', 'AKT1',
                 u'MAPK11', u'MAPK12', u'MAPK13', u'RAF1', u'MAP2K4',
                 u'PIK3CA', u'PIK3CB', u'MAP3K1', u'PIK3CD', u'PIK3R3',
                 u'PIK3R2', u'PIK3R1', u'TP53', u'CASP3', u'PTEN', u'PAK1',
                 u'PAK2', u'BCL2'}
        g = self.net_sub.shortest_paths_between_two_proteins(start, end,
                                                             draw=True,
                                                             )
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
