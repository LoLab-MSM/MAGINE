import os
import shutil
import tempfile

import matplotlib.pyplot as plt
import pandas as pd

import magine.plotting.wordcloud_tools as wt
from magine.enrichment import load_enrichment_csv
from magine.plotting.venn_diagram_maker import create_venn3, create_venn2


class TestVennDiagram(object):
    def setUp(self):
        self.x = ['A', 'B', 'C', 'D']
        self.y = ['C', 'D', 'E', 'F']
        self.z = ['D', 'E', 'F', 'N', 'Z', 'A']
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.exp_data = None
        shutil.rmtree(self.out_dir)

    def test_venn_2(self):
        create_venn2(self.x, self.y, 'X', 'Y',
                     os.path.join(self.out_dir, 'test_1'))
        plt.close()

    def test_venn_3(self):
        create_venn3(self.x, self.y, self.z, 'X', 'Y', 'z',
                     os.path.join(self.out_dir, 'test_1'))
        plt.close()


def test_filter():
    df = load_enrichment_csv(os.path.join(os.path.dirname(__file__),
                                          'Data', 'enrichr_test_enrichr.csv'))
    df2 = df.copy()
    df['db'] = 'kegg'
    df2['db'] = 'reactome'
    df2['sample_id'] = 2

    df = pd.concat([df, df2])

    x = wt.create_wordcloud(df)
    out_dir = tempfile.mkdtemp()

    x.plot(save_name=os.path.join(out_dir, 'test_wc'))
    plt.close()
    wt.word_cloud_from_array(df, sample_ids=[1, 2])
    plt.close()
