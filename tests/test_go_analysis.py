import sys

sys.path.append(".")

from magine.ontology.ontology_analysis import GoAnalysis
from sample_experimental_data import exp_data


def test_html():
    go = GoAnalysis(species='hsa', output_directory='html_output',
                    verbose=False,
                    experimental_data=exp_data)

    df = go.create_enrichment_array(exp_data.proteomics_up_over_time,
                                    exp_data.timepoints,
                                    )
    quit()
    go.write_table_to_html(save_name='index')

    go.write_table_to_html(df, save_name='index2')


if __name__ == '__main__':
    test_html()
