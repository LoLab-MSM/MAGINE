"""
import tempfile

from magine.ontology.ontology_analysis import GoAnalysis
from magine.tests.sample_experimental_data import exp_data


def test_html():
    out_dir = tempfile.mkdtemp()
    go = GoAnalysis(species='hsa', output_directory=out_dir,
                    verbose=False,
                    experimental_data=exp_data,
                    # reference=exp_data.list_species
                    )

    enrich_array = go.calculate_enrichment(exp_data.proteomics_up_over_time,
                            exp_data.timepoints)
    enrich_array.to_csv('proteomics_up_enrichment_array.csv')

    go.write_table_to_html('DEL', out_dir)


if __name__ == '__main__':
    test_html()
"""
