import sys
for i in sys.path:
    print(i)
from magine.plotting.species_plotting import plot_list_of_genes2
import pandas as pd
import os

"""
def test_plotly():
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'large_example.csv')
    d_name = os.path.join(os.path.dirname(__file__), 'Data', 'example_apoptosis.csv')

    d = pd.read_csv(d_name)

    # list_species = ['CAV1', 'GPRC5A', 'PLEC']
    list_species = ['PARP1', 'TP53']
    import time
    st = time.time()
    plot_list_of_genes2(d, list_of_genes=list_species,
                        out_dir='TEST_DELETE', save_name='test_me2',
                        plot_type='plotly')
    print(time.time()-st)


def test_matplotlib():
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'large_example.csv')
    d_name = os.path.join(os.path.dirname(__file__), 'Data', 'example_apoptosis.csv')

    d = pd.read_csv(d_name)

    # list_species = ['CAV1', 'GPRC5A', 'PLEC']
    list_species = ['PARP1', 'TP53']
    plot_list_of_genes2(d, list_of_genes=list_species,
                        out_dir='DELETE', save_name='test_me',
                        plot_type='matplotlib')


if __name__ == '__main__':
    test_plotly()
"""
