import os

import pandas as pd

import magine.html_templates.html_tools as html_tools
from magine.tests.sample_experimental_data import exp_data


def test_filter():
    dirname = os.path.join(os.path.dirname(__file__),
                           'Data', 'test_proteomics_up_all_data.csv')

    df = pd.read_csv(dirname)
    df = df[df['aspect'] == 'biological_process']
    df = pd.pivot_table(df,
                        index=['GO_id', 'GO_name', 'depth', 'ref', 'slim',
                               'aspect'],
                        columns='sample_index')

    html_tools.write_filter_table(df, 'html_writer', 'test')

    html_tools.write_single_table(df, 'html_writer2', 'test2')

"""
def test_create_plots_per_go():
    dirname = os.path.join(os.path.dirname(__file__),
                           'Data', 'proteomics_up_enrichment_array.csv')

    html_tools.write_table_to_html_with_figures(dirname, exp_data,
                                                out_dir='DEL',
                                                run_parallel=False)
"""

if __name__ == '__main__':
    test_filter()
    # test_create_plots_per_go()
