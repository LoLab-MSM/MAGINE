from  magine.plotting.heatmaps import create_heatmaps_pvalue_xs
from magine.ontology.ontology_tools import filter_ontology_df
import os

d_name = os.path.join(os.path.dirname(__file__), 'Data',
                      'test_df_all_data.csv')
df = filter_ontology_df(d_name,
                        trim_nodes=True,
                        n_hits_per_time=5,
                        additional_ids_to_include=['GO:0016358'])
create_heatmaps_pvalue_xs(df,
                          save_name='test_heatplot',
                          mark_pvalues=True,
                          )
