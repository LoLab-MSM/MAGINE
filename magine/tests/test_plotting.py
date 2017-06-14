import os
import pandas as pd
import tempfile
from magine.plotting.venn_diagram_maker import create_venn3, create_venn2
from magine.plotting.volcano_plots import volcano_plot
import magine.plotting.species_plotting as plotter


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
        plotter.plot_list_of_genes(self.data, list_of_genes=list_species,
                                   out_dir=self.out_dir,
                                   save_name='OUT_test_plotly',
                                   plot_type='plotly',
                                   )

    def test_matplotlib(self):
        """
        test matplotlib graph creation from a pandas.dataframe
        Returns
        -------

        """
        list_species = ['PARP1', 'TP53']
        plotter.plot_list_of_genes(self.data, list_of_genes=list_species,
                                   out_dir=self.out_dir,
                                   save_name='OUT_test_matplotlib',
                                   plot_type='matplotlib')

    def test_plot_df_genes(self):
        plotter.plot_dataframe(exp_data=self.data,
                               html_filename='OUT_test_plot_df_genes',
                               out_dir=self.out_dir,
                               plot_type='plotly'
                               )

    def test_plot_df_metabolites(self):
        plotter.plot_dataframe(exp_data=self.data,
                               type_of_species='metabolites',
                               html_filename='OUT_test_plot_df_metabolites',
                               out_dir=self.out_dir,
                               plot_type='plotly'
                               )
        plotter.plot_dataframe(exp_data=self.data,
                               type_of_species='metabolites',
                               html_filename='OUT_test_plot_df_metabolites2',
                               out_dir=self.out_dir,
                               plot_type='matplotlib'
                               )

    def test_plot_list_of_metabolites(self):
        ex_list = ['HMDB1', 'HMDB2']
        plotter.plot_list_of_metabolites(dataframe=self.data,
                                         list_of_metab=ex_list,
                                         out_dir=self.out_dir,
                                         save_name='metabolites'
                                         )

    def test_volcano(self):
        volcano_plot(data=self.data, save_name='volcano_test',
                     out_dir=self.out_dir)
        volcano_plot(data=self.data, save_name='volcano_test',
                     out_dir=self.out_dir, bh_criteria=True, y_range=[0, 5],
                     x_range=[-10, 10])


class TestVennDiagram(object):
    def setUp(self):
        self.x = ['A', 'B', 'C', 'D']
        self.y = ['C', 'D', 'E', 'F']
        self.z = ['D', 'E', 'F', 'N', 'Z', 'A']
        self.out_dir = tempfile.mkdtemp()

    def test_venn_2(self):
        create_venn2(self.x, self.y, 'X', 'Y', 'test_1', out_dir=self.out_dir)

    def test_venn_3(self):
        create_venn3(self.x, self.y, self.z, 'X', 'Y', 'z', 'test_1',
                     out_dir=self.out_dir)
