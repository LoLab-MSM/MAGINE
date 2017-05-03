import os

import pandas as pd

from magine.plotting.species_plotting import plot_list_of_genes2


def test_plotly():
    """
    tests plotly offline graph creation from a pandas.dataframe
    Returns
    -------

    """
    # d_name = os.path.join(os.path.dirname(__file__), 'Data',
    #                       'large_example.csv')
    d_name = os.path.join(os.path.dirname(__file__), 'Data', 'example_apoptosis.csv')

    d = pd.read_csv(d_name)

    # list_species = ['CAV1', 'GPRC5A', 'PLEC']
    list_species = ['PARP1', 'TP53']
    plot_list_of_genes2(d, list_of_genes=list_species,
                        out_dir='TEST_DELETE', save_name='test_me2',
                        plot_type='plotly',
                        )


def test_matplotlib():
    """
    test matplotlib graph creation from a pandas.dataframe
    Returns
    -------

    """
    # d_name = os.path.join(os.path.dirname(__file__), 'Data',
    #                       'large_example.csv')
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')

    d = pd.read_csv(d_name)

    # list_species = ['CAV1', 'GPRC5A', 'PLEC']
    list_species = ['PARP1', 'TP53']
    plot_list_of_genes2(d, list_of_genes=list_species,
                        out_dir='DELETE', save_name='test_me',
                        plot_type='matplotlib')


if __name__ == '__main__':
    test_plotly()
    test_matplotlib()
