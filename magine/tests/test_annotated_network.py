import os

import networkx as nx

from magine.enrichment.enrichr import Enrichr
from magine.networks.annotated_set import create_subnetwork


def test_ont_grouping():
    e = Enrichr()

    list_1 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']

    list_2 = ['DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']

    df = e.run_samples([list_1, list_2], 'GO_Biological_Process_2017')
    df = df.sig.copy()
    _path = os.path.join(os.path.dirname(__file__), 'Network_files',
                         'sample_network.gml')
    network = nx.read_gml(_path)
    create_subnetwork(df, network, use_threshold=True, save_name='test')
    create_subnetwork(df, network, use_threshold=True, use_fdr=True)
    create_subnetwork(df, network, use_threshold=True, use_fdr=True,
                      out_dir='del')
