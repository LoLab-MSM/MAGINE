from magine.networks.network_generator import build_network, create_all_of_kegg


def test_build_network():
    graph = build_network(['Bax', 'tp53', 'JAK1', 'bad'], num_overlap=1, save_name='tmp', species='hsa',
                          overwrite=False)
    for i in graph.nodes():
        if len(i.split(':')) > 1:
            print(i)


# download_all_of_kegg('test')
create_all_of_kegg('test')
