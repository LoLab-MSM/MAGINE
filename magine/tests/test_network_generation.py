import networkx as nx

import magine.networks.network_generator as ng


def test_expand_by_hmdb():
    g = nx.DiGraph()
    g.add_edge('PNLIP', 'LIPC')
    new_g = ng.expand_by_hmdb(graph=g,
                              metabolite_list=['HMDB42489', 'HMDB59874']
                              )
    assert len(new_g.nodes) == 3
    assert len(new_g.edges) == 3


def test_build_network():
    ng.build_network(['BAX', 'TP53', 'JAK1', 'BAD'],
                     save_name='sample_network', species='hsa',
                     all_measured_list=['CASP3', 'EGFR'],
                     use_hmdb=True, use_reactome=True)


if __name__ == '__main__':
    test_build_network()
