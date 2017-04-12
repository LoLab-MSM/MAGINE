import pandas as pd
from magine import ExperimentalData
from magine import Analyzer

if __name__ == '__main__':
    data = pd.read_csv('Data/norris_et_al_2017_cisplatin_data.csv.gz',
                       low_memory=False)
    exp_data = ExperimentalData(data)

    magine = Analyzer(experimental_data=exp_data,
                      save_name='cisplatin_go_analysis')

    # magine.generate_network(save_name='cisplatin_example')

    # Runs go analysis
    magine.run_go_and_create_html('cisplatin_main')
    quit()
    magine.create_selected_go_network(
            'cisplatin_go_analysis_proteomics_up_all_data.csv',
            save_name='cisplatin_magine_network',
            visualize=True
    )

    graph = nx.read_gml('cisplatin_example.gml')
    magine = Analyzer(exp_data, save_name='cisplatin_go_analysis',
                      network=graph)

    df = filter_ontology_df('cisplatin_go_analysis_proteomics_up_all_data.csv',
                            n_hits_per_time=5,
                            go_aspects='biological_process',
                            )
    df2 = filter_ontology_df(
        'cisplatin_go_analysis_proteomics_up_all_data.csv',
        n_hits_per_time=5,
        trim_nodes=True,
        go_aspects='biological_process',
        )
    create_heatmaps_pvalue_xs(df,
                              save_name='cisplatin_prot_up_not_filtered',
                              metric='enrichment_score',
                              mark_pvalues=True)
    create_heatmaps_pvalue_xs(df2,
                              save_name='cisplatin_prot_up_filtered',
                              metric='enrichment_score',
                              mark_pvalues=True)
