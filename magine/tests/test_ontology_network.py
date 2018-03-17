import os

import networkx as nx

from magine.networks.ontology_network import OntologyNetworkGenerator
from magine.ontology.enrichr import Enrichr


def test_ont_grouping():
    e = Enrichr()

    list_2 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
    df = e.run(list_2, 'GO_Biological_Process_2017')
    df['genes'] = df['genes'].str.split(',')
    df['termID'] = df['term_id'].str.replace(':', '')
    term_dict = dict(zip(df['termID'], df['genes']))
    label_dict = dict(zip(df['termID'], df['term_name']))
    term_list = list(df.head(3)['termID'])

    _path = os.path.join(os.path.dirname(__file__), 'Network_files',
                         'sample_network.gml')
    network = nx.read_gml(_path)
    ont = OntologyNetworkGenerator(network)

    ont.create_network_from_list(term_list, term_dict, label_dict,
                                 save_name='enrichr_network', draw=False)
