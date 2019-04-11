import os

import networkx as nx
from nose.tools import ok_

from magine.networks.subgraphs import Subgraph
from magine.tests.sample_experimental_data import exp_data

_path = os.path.join(os.path.dirname(__file__), 'Network_files',
                     'sample_network.gml')


class TestSubgraphs(object):
    def __init__(self):
        self.network = nx.read_gml(_path)
        self.net_sub = Subgraph(self.network, exp_data)

    def test_downstream_nodes(self):
        """Test finding downstream nodes."""
        down_nodes = {'MTMR14', 'GABARAPL1', 'HMDB03850', 'GABARAPL2', 'EEA1',
                      'ATG16L1', 'ATG16L2', 'WIPI1', 'GABARAP', 'WIPI2',
                      'ATG12', 'ATG5', 'ZFYVE1'}
        for i in self.net_sub.downstream_of_node('MTMR14').nodes:
            ok_(i in down_nodes)

    def test_upstream_nodes(self):
        """Test finding upstream nodes."""
        up_nodes = {'MTMR4', 'MTMR14', 'MTMR3', 'ZFYVE1', 'HMDB03850'}
        for i in self.net_sub.upstream_of_node('ZFYVE1').nodes:
            ok_(i in up_nodes)

    def test_path_between_two_does_not_exist(self):
        start = 'HSPA9'
        end = 'ZFYVE1'
        g = self.net_sub.paths_between_pair(start, end, bidirectional=True)
        ok_(len(g.nodes) == 0)

    def test_path_between_two(self):
        start = 'AKT1'
        end = 'BAX'
        nodes = {'BCL2L1', 'PTK2', 'HMDB04249', 'BAX', 'MAPK14', 'AKT1',
                 'MAPK11', 'MAPK12', 'MAPK13', 'RAF1', 'MAP2K4',
                 'PIK3CA', 'PIK3CB', 'MAP3K1', 'PIK3CD', 'PIK3R3',
                 'PIK3R2', 'PIK3R1', 'TP53', 'CASP3', 'PTEN', 'PAK1',
                 'PAK2', 'BCL2'}
        g = self.net_sub.paths_between_pair(start, end,
                                            bidirectional=True)
        ok_(set(g.nodes()) == nodes)

    def test_paint_over_time(self):
        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.paths_between_list(list_2, draw=False,
                                            single_path=True,
                                            )
        colors = ['red'] * len(exp_data.sample_ids)
        # self.net_sub.measured_networks_over_time(g, colors, prefix='colored')

    def test_paint_over_time_up_down(self):
        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.paths_between_list(list_2, draw=False,
                                            single_path=True,
                                            )

        # self.net_sub.measured_networks_over_time_up_down(g, 'colored_updown')

    def test_paths_from_list(self):
        """Test finding paths from a list."""
        up_nodes = {'MTMR4', 'ZFYVE1', 'MTMR14', 'MTMR3', 'WIPI1'}
        edges = {('MTMR14', 'HMDB03850'),
                 ('HMDB03850', 'ZFYVE1'),
                 ('HMDB03850', 'WIPI1'),
                 ('MTMR4', 'HMDB03850'),
                 ('MTMR3', 'HMDB03850')}

        g = self.net_sub.paths_between_list(up_nodes, draw=False,
                                            save_name='ns_test')

        ok_(set(g.edges()) == edges)

        list_2 = {'CASP3', 'BAX', 'TP53'}
        g = self.net_sub.paths_between_list(list_2, draw=False,
                                            single_path=False,
                                            save_name='smaller_list')
        nodes = {'TP53', 'CASP3', 'CDKN1A', 'MAP3K1', 'BAX', 'MAPK10', 'MAPK8',
                 'MAPK9', 'BCL2'}
        ok_(set(g.nodes) == nodes)

    def test_neighbor(self):
        gene = 'BAX'
        g = self.net_sub.neighbors(gene, True, False, 1)
        layer_1_up = {'BAX', 'BCL2', 'BCL2L1', 'BCL2L11', 'BID', 'MAPK10',
                      'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MAPK8',
                      'MAPK9', 'PRNP', 'SIRT1', 'TP53'}
        ok_(set(g.nodes) == layer_1_up)

        g = self.net_sub.neighbors(gene, True, False, 2)
        layer_2_up = {'BCL2L1', 'HSPA2', 'PMAIP1', 'HSPA6', 'MIR32', 'MIR30D',
                      'PPP2R3B', 'PPP2R3C', 'HSPA8', 'RAPGEF4', 'AKT1', 'AKT2',
                      'AKT3', 'ARAF', 'FCGR1A', 'RAPGEF3', 'MIR30A', 'MIR15A',
                      'EP300', 'HMDB03125', 'MAP3K1', 'PRNP', 'MAP3K7',
                      'PPP2R1B', 'MAP3K5', 'HSPA1A', 'MIR30C2', 'MIR30C1',
                      'RELA', 'KRAS', 'CTSH', 'DDIT3', 'CTSL', 'CTSO',
                      'PPP2R2D', 'CTSB', 'CTSC', 'CTSD', 'CTSF', 'PPP2R2B',
                      'MIR125B1', 'CTSZ', 'HMDB00902', 'MIR30B', 'MDM2',
                      'CTSS',
                      'MDM4', 'CTSV', 'CTSW', 'PPP2R3A', 'PPP2CA', 'PPP2CB',
                      'PPP2R5B', 'BMPR2', 'LAT', 'ATR', 'HRAS', 'EIF2AK4',
                      'TRAF3', 'PINK1', 'IL1R1', 'TRAF5', 'MIR223',
                      'MIR143', 'JUN', 'MAPK8IP1', 'CHEK2', 'CHEK1',
                      'SIRT1', 'NR4A1', 'FOS', 'PRKCD', 'RHOA', 'PTPN7',
                      'PTPN5', 'TP53', 'GZMB', 'TNF', 'PPP2R5D',
                      'PPP2R5E', 'CDC42', 'RALB', 'RALA', 'PPP2R5C',
                      'MAP2K1', 'PPP2R1A', 'MIR34A', 'ATM', 'RAP1B',
                      'RAP1A', 'PRKDC', 'MAPK14', 'TGFB3', 'MAPK10',
                      'MAPK11', 'MAPK12', 'MAPK13', 'RAF1', 'CREBBP',
                      'HMDB04947', 'MECOM', 'MAPK3', 'MAPK1', 'DUSP9',
                      'DUSP8', 'DUSP5', 'NFKB1', 'MAPK8', 'MAPK9',
                      'DUSP1', 'DUSP3', 'DUSP2', 'CREB1', 'PPP2R2C',
                      'DUSP4', 'USP7', 'NFKB2', 'MAPK8IP2', 'FOXO3',
                      'DUSP6', 'HSPA1B', 'EIF2AK1', 'EIF2AK3', 'EIF2AK2',
                      'BMPR1B', 'TP73', 'STAT5B', 'STAT5A', 'RIPK2',
                      'ERN1', 'MIR125A', 'MAP2K6', 'TGFB1', 'TGFB2',
                      'PTK2', 'ACVR1', 'MIR30E', 'FCGR3A', 'BAX', 'MCL1',
                      'FCGR3B', 'MAP2K3', 'BMPR1A', 'BAD', 'CTSK', 'MAP2K7',
                      'MAP2K4', 'BBC3', 'HSPA1L', 'RCHY1', 'CREB3L2',
                      'BCL2L11', 'CREB3L3', 'BID', 'SOD1', 'NRAS',
                      'CREB3L1', 'CREB5', 'TRAF1', 'TRAF2', 'NOS3',
                      'RAC2', 'RAC3', 'TRAF6', 'RAC1', 'PRKCA', 'CREB3',
                      'PRKCB', 'PRKCE', 'FCGR2C', 'FCGR2A', 'PRKCZ',
                      'HMDB03747', 'IL1RAP', 'DUSP16', 'DUSP10', 'PTPRR',
                      'MIR125B2', 'CASP8', 'PPP2R5A', 'PPP2R2A', 'MAP3K11',
                      'DUSP7', 'BRAF', 'CREB3L4', 'ATF4', 'ATF6B', 'ATF2',
                      'BCL2'}

        new_nodes = set(g.nodes)

        ok_(set(g.nodes) == layer_2_up)

        g = self.net_sub.neighbors(gene, False, True, 1)
        layer_1_down = {'CASP3', 'CYCS', 'BAX', 'CAPN2', 'CAPN1', 'BCL2'}
        ok_(set(g.nodes) == layer_1_down)

        g = self.net_sub.neighbors(gene, False, True, 2)
        layer_2_down = {'ACTG1', 'SPTAN1', 'CYCS', 'SPTA1', 'CASP12',
                        'CASP9', 'NLRP1', 'BIK', 'MAP3K1', 'TLN2', 'DCC',
                        'TLN1', 'CDK5R1', 'CAPN2', 'CAPN1', 'APAF1',
                        'PARP4', 'DFFA', 'PARP2', 'PARP1', 'STK3', 'STK4',
                        'BAX', 'PARP3', 'CASP7', 'TP53', 'CASP3', 'BAD',
                        'BAK1', 'ACTB', 'PAK1', 'PAK2', 'BCL2'}

        ok_(set(g.nodes) == layer_2_down)

    def test_expand(self):
        g = nx.DiGraph()
        g.add_edge('BCL2L1', 'BAX')
        g.add_edge('MAPK14', 'BAX')

        ng = self.net_sub.expand_neighbors(network=g, nodes=['BAX'],
                                           downstream=True)

        trues = {
            'BCL2L1', 'CASP3', 'CAPN2', 'BAX', 'MAPK14', 'CYCS',
            'CAPN1', 'BCL2'
        }

        ok_(trues == set(ng.nodes))

        includes = ['CASP3', 'CAPN2']

        ng = self.net_sub.expand_neighbors(g, 'BAX', downstream=True,
                                           include_only=includes)

        trues = {
            'BCL2L1', 'CASP3', 'CAPN2', 'BAX', 'MAPK14',
        }
        ok_(trues == set(ng.nodes))


if __name__ == '__main__':
    t = TestSubgraphs()
    t.test_neighbor()
