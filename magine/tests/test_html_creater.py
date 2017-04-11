from magine.html_templates.html_tools import write_filter_table, \
    write_single_table
import pandas as pd
import os


def test_filter():
    dirname = os.path.join(os.path.dirname(__file__),
                           'Data', 'test_proteomics_up_all_data.csv')

    df = pd.read_csv(dirname)
    df = df[df['aspect'] == 'biological_process']
    df = pd.pivot_table(df,
                        index=['GO_id', 'GO_name', 'depth', 'ref', 'slim',
                               'aspect'],
                        columns='sample_index')

    write_filter_table(df, 'html_writer', 'test')

    write_single_table(df, 'html_writer2', 'test2')
