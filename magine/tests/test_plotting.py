import os

import pandas as pd

from magine.plotting.volcano_plots import volcano_plot


def test_volcano():
    d_name = os.path.join(os.path.dirname(__file__), 'Data',
                          'example_apoptosis.csv')
    data = pd.read_csv(d_name)
    volcano_plot(data=data, save_name='volcano_test', out_dir='del')
    volcano_plot(data=data, save_name='volcano_test', out_dir='del',
                 bh_criteria=True, y_range=[0, 5], x_range=[-10, 10])
