import sys

import pandas as pd

import magine.enrichment.tools as et
from magine.data.datatypes import ExperimentalData
from magine.data.formatter import process_raptr_zip
from magine.enrichment.enrichr import run_enrichment_for_project


def return_table_from_model(save_name, project_name, category, dbs):
    cols = ['term_name', 'combined_score', 'adj_p_value', 'rank', 'genes',
            'n_genes', 'sample_id', 'db', 'category']

    project_name = project_name.split(',')
    category = category.split(',')
    dbs = dbs.split(',')
    if len(project_name) > 1:
        cols.insert(0, 'project_name')

    df = pd.read_csv(save_name)
    if isinstance(project_name, str):
        df = df[df['project_name'] == project_name]
    elif isinstance(project_name, list):
        df = df[df['project_name'].isin(project_name)]

    df = df.filter(category__in=category)
    df = df.filter(db__in=dbs)
    df = pd.DataFrame(list(df.values()))[cols]

    df = et.filter_dataframe(df,
                             p_value=0.2,
                             combined_score=0,
                             db=dbs,
                             category=category)

    tmp_table = df.copy()
    tmp_table['genes'] = tmp_table['genes'].str.split(',').str.join(', ')

    tmp_table = tmp_table[cols]
    data = tmp_table.to_dict('split')
    template_vars = {"data": data}
    return template_vars


def run(zip_file, save_name):
    df = process_raptr_zip(zip_file)
    exp_data = ExperimentalData(df)
    run_enrichment_for_project(exp_data, save_name)



_help = """
please provide a savename, project_name

All are expected to be comma delimited if more than one choice
"""



if __name__ == '__main__':
    print(sys.argv)
    # assert len(sys.argv) == 2, _help
    run(sys.argv[1], sys.argv[2])
    # return_table_from_model(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
