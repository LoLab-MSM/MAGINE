import os
import pandas as pd
import tempfile
from magine.plotting.venn_diagram_maker import create_venn3, create_venn2
from magine.plotting.volcano_plots import volcano_plot
from magine.plotting.species_plotting import plot_dataframe, \
    plot_list_of_genes2


class TestSpeciesPlotting(object):
    def setUp(self):
        d_name = os.path.join(os.path.dirname(__file__), 'Data',
                              'example_apoptosis.csv')
        self.data = pd.read_csv(d_name)
        self.data['compound_id'] = self.data['compound']
        self.data['time_points'] = self.data['time']
        self.out_dir = tempfile.mkdtemp()

    def test_plotly(self):
        """
        tests plotly offline graph creation from a pandas.dataframe
        Returns
        -------

        """
        list_species = ['PARP1', 'TP53']
        plot_list_of_genes2(self.data, list_of_genes=list_species,
                            out_dir=self.out_dir, save_name='test_me2',
                            plot_type='plotly',
                            )

    def test_matplotlib(self):
        """
        test matplotlib graph creation from a pandas.dataframe
        Returns
        -------

        """
        list_species = ['PARP1', 'TP53']
        plot_list_of_genes2(self.data, list_of_genes=list_species,
                            out_dir=self.out_dir, save_name='test_me',
                            plot_type='matplotlib')

    def test_plot_df(self):
        plot_dataframe(exp_data=self.data,
                       species_type='proteins',
                       html_filename='test',
                       out_dir=self.out_dir,
                       plot_type='plotly'
                       )

    def test_volcano(self):
        volcano_plot(data=self.data, save_name='volcano_test')
        volcano_plot(data=self.data, save_name='volcano_test',
                     bh_criteria=True, y_range=[0, 5], x_range=[-10, 10])


class TestVennDiagram(object):
    def setUp(self):
        self.x = ['A', 'B', 'C', 'D']
        self.y = ['C', 'D', 'E', 'F']
        self.z = ['D', 'E', 'F', 'N', 'Z', 'A']
        self.out_dir = tempfile.mkdtemp()

    def test_venn_2(self):
        create_venn2(self.x, self.y, 'X', 'Y', 'test_1')

    def test_venn_3(self):
        create_venn3(self.x, self.y, self.z, 'X', 'Y', 'z', 'test_1')

