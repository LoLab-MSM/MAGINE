import os

import networkx as nx

from magine.networks.cytoscape_view import RenderModel


def test_by_time():
    ddn = nx.nx.read_graphml(
        'Network_files/test_all_colored_enrichment.graphml')
    rm = RenderModel(ddn, style='Directed', )  # layout='force-directed')
    rm.visualize_by_list_of_time(
        ['time_0000', 'time_0001', 'time_0002', 'time_0003'],
        labels=[0, 2, 3, 4])
    files_that_should_exist = [
        'tmp/Figures/go_network_tmp_time_0_formatted.png',
        'tmp/Figures/go_network_tmp_time_1_formatted.png',
        'tmp/Figures/go_network_tmp_time_2_formatted.png',
        'tmp/Figures/go_network_tmp_time_3_formatted.png',
    ]
    for i in files_that_should_exist:
        assert os.path.exists(i)


test_by_time()
