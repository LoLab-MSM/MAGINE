import sys

from magine.data.datatypes import ExperimentalData
from magine.data.formatter import process_raptr_zip
from magine.enrichment.enrichr import run_enrichment_for_project


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
