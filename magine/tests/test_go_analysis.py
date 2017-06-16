from magine.ontology.ontology_analysis import GoAnalysis
from magine.tests.sample_experimental_data import exp_data

"""
def test_create_dicts():
    create_dicts_through_orange(species='hsa', rev="5.2795", rev_ass='1.353')
    m_go = MagineGO()
"""

def test_html():
    go = GoAnalysis(species='hsa', output_directory='html_output',
                    verbose=False,
                    experimental_data=exp_data,
                    # reference=exp_data.list_species
                    )

    enrich_array = go.calculate_enrichment(exp_data.proteomics_up_over_time,
                            exp_data.timepoints)
    enrich_array.to_csv('proteomics_up_enrichment_array.csv')

    go.write_table_to_html('DEL', 'DEL')


if __name__ == '__main__':
    test_html()
