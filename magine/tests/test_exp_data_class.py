import os
import shutil
import tempfile

import pandas as pd

from magine.data.experimental_data import ExperimentalData, load_data_csv


class TestExpData(object):
    def setUp(self):
        self._dir = os.path.join(os.path.dirname(__file__), 'Data')
        self.exp_data = load_data_csv(os.path.join(self._dir,
                                                   'example_apoptosis.csv'))
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.exp_data = None
        shutil.rmtree(self.out_dir)

    def test_load_from_df(self):
        df = pd.read_csv(os.path.join(self._dir, 'test_data.csv'))
        df['compound_id'] = df['compound']
        exp_data = ExperimentalData(df)

    def test_protein(self):

        assert self.exp_data.proteins.list == {'AHR', 'PARP4', 'ADORA1', 'BAX',
                                               'TP53', 'PARP1', 'ADRA1A',
                                               'ADORA2A', 'CASP3', 'AGTR2'}
        assert self.exp_data.proteins.sig.list == {'ADRA1A', 'BAX', 'TP53',
                                                   'AGTR2', 'PARP4', 'PARP1',
                                                   'CASP3'}

    def test_gene(self):

        assert self.exp_data.genes.list == {'ADORA1', 'PARP4', 'AKT1', 'CASP3',
                                            'ADRA1A', 'AIF1', 'PARP1', 'AGTR2',
                                            'BAX', 'AKT2', 'ADORA2A', 'TP53',
                                            'AHR'}

        assert self.exp_data.genes.sig.list == {'TP53', 'AIF1', 'AKT1', 'BAX',
                                                'CASP3', 'PARP4', 'PARP1',
                                                'AGTR2', 'ADRA1A'}

    def test_rna(self):
        assert self.exp_data.rna.list == {'AIF1', 'AKT1', 'AKT2'}
        assert self.exp_data.rna.sig.list == {'AIF1', 'AKT1'}

    def test_compounds(self):
        assert self.exp_data.compounds.sig.list == {'HMDB2', 'HMDB1'}
        assert self.exp_data.compounds.list == {'HMDB2', 'HMDB1'}

    def test_species(self):
        print(self.exp_data.species.list)

    def test_plot_list(self):
        x = self.exp_data.rna.sig.list
        self.exp_data.plot_list_of_genes(x, 'del_test', self.out_dir)

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
