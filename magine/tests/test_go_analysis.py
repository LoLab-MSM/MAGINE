import sys

sys.path.append(".")

from magine.ontology.ontology_analysis import GoAnalysis
from magine.tests.sample_experimental_data import exp_data
from magine.ontology.enrichment_calculation import MagineGO, \
    create_dicts_through_orange


def test_create_dicts():
    create_dicts_through_orange(species='hsa', rev="5.2795", rev_ass='1.353')
    m_go = MagineGO()


def test_html():
    go = GoAnalysis(species='hsa', output_directory='html_output',
                    verbose=False,
                    experimental_data=exp_data,
                    # reference=exp_data.list_species
                    )

    df = go.calculate_enrichment(exp_data.proteomics_up_over_time,
                                 exp_data.timepoints,
                                 )

    print(df.sort_values(by='enrichment_score', ascending=False).head(10))

    # go.write_table_to_html(save_name='index')

    # go.write_table_to_html(df, save_name='index2')


if __name__ == '__main__':
    test_html()
