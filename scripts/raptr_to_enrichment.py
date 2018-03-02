import sys

from magine.data.datatypes import ExperimentalData
from magine.data.formatter import process_raptr_zip
from magine.ontology.enrichr import run_enrichment_for_project


def run(zip_file, save_name):
    df = process_raptr_zip(zip_file)
    exp_data = ExperimentalData(df)
    run_enrichment_for_project(exp_data, save_name)


if __name__ == '__main__':
    print(sys.argv)
    assert len(sys.argv) == 3, 'please provide a zipfile and a savename'

    run(sys.argv[1], sys.argv[2])
