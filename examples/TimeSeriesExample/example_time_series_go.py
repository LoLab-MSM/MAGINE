import pandas as pd

from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer
from magine.ontology.ontology_tools import filter_ontology_df

if __name__ == '__main__':
    data = pd.read_csv('Data/norris_et_al_2017_cisplatin_data.csv.gz',
                       low_memory=False)
    exp_data = ExperimentalData(data)

    magine = Analyzer(experimental_data=exp_data,
                      save_name='cisplatin_go_analysis')

    # Runs go analysis
    # magine.run_proteomics_go()
    # magine.run_go_and_create_html('cisplatin_main')
    df = pd.read_csv('tmp/cisplatin_go_analysis_proteomics_up_all_data.csv')
    df.sort_values(by='pvalue', ascending=True, inplace=True)
    cols = ['GO_id', 'GO_name', 'enrichment_score', 'pvalue']
    print(df[cols].head(10))

    sig_results = filter_ontology_df(df,
                                     pvalue=0.05,
                                     go_aspects='biological_process',
                                     trim_nodes=True,
                                     n_hits_per_time=10,
                                     variable_of_interest='pvalue'
                                     )

    sig_results.sort_values(by='pvalue', ascending=True,
                            inplace=True)
    print(sig_results[cols].head(10))
