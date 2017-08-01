from magine.networks.network_generator import build_network


def test_build_network():
    graph = build_network(['BAX', 'TP53', 'JAK1', 'BAD'], num_overlap=1,
                          save_name='sample_network', species='hsa',
                          overwrite=False, all_measured_list=['CASP3', 'EGFR'],
                          use_hmdb=False, use_reactome=True
                          )
    for i in graph.nodes():
        # if isinstance(i, float):
        #     print(i)
        if len(i.split(':')) > 1:
            print(i)


# _download_all_of_kegg('test')
if __name__ == '__main__':
    test_build_network()
