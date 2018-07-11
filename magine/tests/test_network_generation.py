from magine.networks.network_generator import build_network


def test_build_network():
    build_network(['BAX', 'TP53', 'JAK1', 'BAD', 'HMDB42489', 'HMDB59874'],
                  save_name='sample_network', species='hsa',
                  all_measured_list=['CASP3', 'EGFR'],
                  use_hmdb=True, use_reactome=True)


if __name__ == '__main__':
    test_build_network()
