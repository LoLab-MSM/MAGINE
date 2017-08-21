import os

import networkx as nx

from magine.networks.network_subgraphs import NetworkSubgraphs
from magine.tests.sample_experimental_data import exp_data

_path = os.path.join(os.path.dirname(__file__), 'Network_files',
                     'sample_network.gml')


# def test():
#     return

class TestSubgraphs(object):
    def __init__(self):
        self.network = nx.read_gml(_path)
        self.net_sub = NetworkSubgraphs(self.network, exp_data)

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

    def test_path_between_two_does_not_exist(self):
        start = 'HSPA9'
        end = 'ZFYVE1'
        g = self.net_sub.shortest_paths_between_two_proteins(start, end,
                                                             bidirectional=True)
        assert g is None

    def test_path_between_two(self):
        start = 'AKT1'
        end = 'BAX'
        nodes = {'BCL2L1', 'PTK2', 'HMDB04249', 'BAX', 'MAPK14', 'AKT1',
                 'MAPK11', 'MAPK12', 'MAPK13', 'RAF1', 'MAP2K4',
                 'PIK3CA', 'PIK3CB', 'MAP3K1', 'PIK3CD', 'PIK3R3',
                 'PIK3R2', 'PIK3R1', 'TP53', 'CASP3', 'PTEN', 'PAK1',
                 'PAK2', 'BCL2'}
        g = self.net_sub.shortest_paths_between_two_proteins(start, end,
                                                             bidirectional=True)
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
        nodes = {'TP53', 'CASP3', 'CDKN1A', 'MAP3K1', 'BAX', 'MAPK10', 'MAPK8',
                 'MAPK9', 'BCL2'}
        assert set(g.nodes()) == nodes

    def test_neighbor(self):
        gene = u'BAX'
        g = self.net_sub.neighbors(gene, True, False, 1)
        layer_1_up = [u'BCL2L1', u'SIRT1', u'BCL2L11', u'BID', u'PRNP', u'BAX',
                      u'MAPK14', u'MAPK10', u'MAPK11', u'MAPK12', u'MAPK13',
                      u'TP53', u'MAPK8', u'MAPK9', u'BCL2']

        assert g.nodes() == layer_1_up

        g = self.net_sub.neighbors(gene, True, False, 2)
        layer_2_up = [u'BCL2L1', u'HSPA2', u'PMAIP1', u'HSPA6', u'MIR32',
                      u'MIR30D', u'PPP2R3B', u'PPP2R3C', u'HSPA8', u'RAPGEF4',
                      u'AKT1', u'AKT2', u'AKT3', u'ARAF', u'FCGR1A', u'RAPGEF3',
                      u'MIR30A', u'MIR15A', u'EP300', u'HMDB03125', u'MAP3K1',
                      u'PRNP', u'MAP3K7', u'PPP2R1B', u'MAP3K5', u'HSPA1A',
                      u'MIR30C2', u'MIR30C1', u'RELA', u'KRAS', u'CTSH',
                      u'DDIT3', u'CTSL', u'CTSO', u'PPP2R2D', u'CTSB', u'CTSC',
                      u'CTSD', u'CTSF', u'PPP2R2B', u'MIR125B1', u'CTSZ',
                      u'HMDB00902', u'MIR30B', u'MDM2', u'CTSS', u'MDM4',
                      u'CTSV', u'CTSW', u'PPP2R3A', u'PPP2CA', u'PPP2CB',
                      u'PPP2R5B', u'BMPR2', u'LAT', u'ATR', u'HRAS', u'EIF2AK4',
                      u'TRAF3', u'PINK1', u'IL1R1', u'TRAF5', u'MIR223',
                      u'MIR143', u'JUN', u'MAPK8IP1', u'CHEK2', u'CHEK1',
                      u'SIRT1', u'NR4A1', u'FOS', u'PRKCD', u'RHOA', u'PTPN7',
                      u'PTPN5', u'TP53', u'GZMB', u'TNF', u'PPP2R5D',
                      u'PPP2R5E', u'CDC42', u'RALB', u'RALA', u'PPP2R5C',
                      u'MAP2K1', u'PPP2R1A', u'MIR34A', u'ATM', u'RAP1B',
                      u'RAP1A', u'PRKDC', u'MAPK14', u'TGFB3', u'MAPK10',
                      u'MAPK11', u'MAPK12', u'MAPK13', u'RAF1', u'CREBBP',
                      u'HMDB04947', u'MECOM', u'MAPK3', u'MAPK1', u'DUSP9',
                      u'DUSP8', u'DUSP5', u'NFKB1', u'MAPK8', u'MAPK9',
                      u'DUSP1', u'DUSP3', u'DUSP2', u'CREB1', u'PPP2R2C',
                      u'DUSP4', u'USP7', u'NFKB2', u'MAPK8IP2', u'FOXO3',
                      u'DUSP6', u'HSPA1B', u'EIF2AK1', u'EIF2AK3', u'EIF2AK2',
                      u'BMPR1B', u'TP73', u'STAT5B', u'STAT5A', u'RIPK2',
                      u'ERN1', u'MIR125A', u'MAP2K6', u'TGFB1', u'TGFB2',
                      u'PTK2', u'ACVR1', u'MIR30E', u'FCGR3A', u'BAX',
                      u'FCGR3B',
                      u'MAP2K3', u'BMPR1A', u'BAD', u'CTSK', u'MAP2K7', u'MCL1',
                      u'MAP2K4', u'BBC3', u'HSPA1L', u'RCHY1', u'CREB3L2',
                      u'BCL2L11', u'CREB3L3', u'BID', u'SOD1', u'NRAS',
                      u'CREB3L1', u'CREB5', u'TRAF1', u'TRAF2', u'NOS3',
                      u'RAC2', u'RAC3', u'TRAF6', u'RAC1', u'PRKCA', u'CREB3',
                      u'PRKCB', u'PRKCE', u'FCGR2C', u'FCGR2A', u'PRKCZ',
                      u'HMDB03747', u'IL1RAP', u'DUSP16', u'DUSP10', u'PTPRR',
                      u'MIR125B2', u'CASP8', u'PPP2R5A', u'PPP2R2A', u'MAP3K11',
                      u'DUSP7', u'BRAF', u'CREB3L4', u'ATF4', u'ATF6B', u'ATF2',
                      u'BCL2']
        assert g.nodes() == layer_2_up

        g = self.net_sub.neighbors(gene, False, True, 1)
        layer_1_down = [u'CASP3', u'CYCS', u'BAX', u'CAPN2', u'CAPN1', u'BCL2']
        assert g.nodes() == layer_1_down

        g = self.net_sub.neighbors(gene, False, True, 2)
        layer_2_down = [u'ACTG1', u'SPTAN1', u'CYCS', u'SPTA1', u'CASP12',
                        u'CASP9', u'NLRP1', u'BIK', u'MAP3K1', u'TLN2', u'DCC',
                        u'TLN1', u'CDK5R1', u'CAPN2', u'CAPN1', u'APAF1',
                        u'PARP4', u'DFFA', u'PARP2', u'PARP1', u'STK3', u'STK4',
                        u'BAX', u'PARP3', u'CASP7', u'TP53', u'CASP3', u'BAD',
                        u'BAK1', u'ACTB', u'PAK1', u'PAK2', u'BCL2']

        assert g.nodes() == layer_2_down


if __name__ == '__main__':
    t = TestSubgraphs()
    t.test_neighbor()
