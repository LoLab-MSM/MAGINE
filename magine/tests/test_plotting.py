import os
import pandas as pd

from magine.plotting.venn_diagram_maker import create_venn3, create_venn2
from magine.plotting.volcano_plots import volcano_plot
from magine.plotting.species_plotting import plot_dataframe, \
    plot_list_of_genes2


def test_volcano():
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')
    data = pd.read_csv(d_name)
    volcano_plot(data=data, save_name='volcano_test')
    volcano_plot(data=data, save_name='volcano_test',
                 bh_criteria=True, y_range=[0, 5], x_range=[-10, 10])


def test_plotly():
    """
    tests plotly offline graph creation from a pandas.dataframe
    Returns
    -------

    """
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')

    d = pd.read_csv(d_name)

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
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')

    d = pd.read_csv(d_name)
    list_species = ['PARP1', 'TP53']
    plot_list_of_genes2(d, list_of_genes=list_species,
                        out_dir='DELETE', save_name='test_me',
                        plot_type='matplotlib')


def test_plot_df():
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')
    d = pd.read_csv(d_name)
    d['compound_id'] = d['compound']
    d['time_points'] = d['time']
    plot_dataframe(d, 'test_df_out', 'test_df_out', 'plotly')


x = ['A', 'B', 'C', 'D']
y = ['C', 'D', 'E', 'F']
z = ['D', 'E', 'F', 'N', 'Z']


def test_venn_2():
    create_venn2(x, y, 'X', 'Y', 'test_1')


def test_venn_3():
    create_venn3(x, y, z, 'X', 'Y', 'z', 'test_1')


if __name__ == '__main__':
    test_plotly()
    test_matplotlib()
    test_plot_df()
