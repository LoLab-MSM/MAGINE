"""

import os

import networkx as nx

from magine.networks.visualization.cytoscape_view import RenderModel


def test_by_time():
    return
    # This test must be ran with a session of cytoscape open
    # For now this will be passed in the travis builds.
    dir_name = os.path.dirname(__file__)
    full_path = os.path.join(dir_name,
                             'Network_files',
                             'test_all_colored_enrichment.graphml')
    ddn = nx.nx.read_graphml(full_path)
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


if __name__ == '__main__':
    test_by_time()
"""
