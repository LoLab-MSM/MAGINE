from magine.networks.network_generator import build_network


def test_build_network():
    graph = build_network(['Bax', 'tp53', 'JAK1', 'bad'], num_overlap=1,
                          save_name='sample_network', species='hsa',
                          overwrite=False, all_measured_list=['casp3', 'egfr'],
                          use_hmdb=False, use_reactome=True
                          )
    for i in graph.nodes():
        if len(i.split(':')) > 1:
            print(i)


# _download_all_of_kegg('test')
if __name__ == '__main__':
    test_build_network()
