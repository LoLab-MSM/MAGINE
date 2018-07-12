from magine.networks.network_generator import build_network


def test_build_network():
    build_network(['BAX', 'TP53', 'JAK1', 'BAD', 'HMDB0042489', 'HMDB0059874'],
                  save_name=None, species='hsa',
                  all_measured_list=['CASP3', 'EGFR'])


if __name__ == '__main__':
    test_build_network()
