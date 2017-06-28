import os
import shutil
import tempfile

import pandas as pd

from magine.data.datatypes import ExperimentalData


class TestExpData(object):
    def setUp(self):
        self._dir = os.path.join(os.path.dirname(__file__), 'Data')
        self.exp_data = ExperimentalData('example_apoptosis.csv', self._dir)
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.exp_data = None
        shutil.rmtree(self.out_dir)

    def test_load_from_df(self):
        df = pd.read_csv(os.path.join(self._dir, 'test_data.csv'))
        df['compound_id'] = df['compound']
        exp_data = ExperimentalData(df)

    def test_table(self):
        df = pd.read_csv(os.path.join(self._dir, 'test_data.csv'))
        df['compound_id'] = df['compound']
        exp_data = ExperimentalData(df)
        exp_data.create_table_of_data(unique=True, save_name='unique_table')
        exp_data.create_table_of_data(unique=False, sig=True,
                                      save_name='sig_unique_table')

    def test_plot_list(self):
        x = list(self.exp_data.list_proteins)
        self.exp_data.plot_list_of_genes(x, 'del_test', out_dir=self.out_dir)

    def test_html_output(self):
        self.exp_data.plot_all_proteins(html_file_name='del',
                                        out_dir=self.out_dir)

    def test_volcano(self):
        self.exp_data.volcano_analysis(out_dir=self.out_dir)

    def test_time_series_volcano(self):
        self.exp_data.time_series_volcano('label_free', 'test_label_free',
                                          out_dir=self.out_dir, bh_critera=True)

    def test_plot_all_metabolites(self):
        self.exp_data.plot_all_metabolites('metab', out_dir=self.out_dir,
                                           plot_type='matplotlib')

    def test_list_metabolites(self):
        l = ['HMDB2', 'HMDB1']
        self.exp_data.plot_list_of_metabolites(l, save_name='metab',
                                               out_dir=self.out_dir,
                                               plot_type='matplotlib'
                                               )

    def test_histogram(self):
        self.exp_data.create_histogram_measurements('label_free',
                                                    save_name='lf_hist',
                                                    y_range=[0, 100],
                                                    out_dir=self.out_dir)
