import os
import shutil
import tempfile

import matplotlib.pyplot as plt

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
