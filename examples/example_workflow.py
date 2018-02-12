from magine.data.datatypes import ExperimentalData
from magine.html_templates.html_tools import workflow_template
from magine.networks.network_generator import build_network
from magine.ontology.enrichr import Enrichr
from magine.plotting.species_plotting import plot_dataframe

if __name__ == '__main__':
    data = ExperimentalData('Data/norris_et_al_2017_cisplatin_data.csv.gz')

    # want data organized
    # data
    # proteins
    # metabolomics

    plot_dataframe(data.metabolites, 'metabolites', 'metabolites',
                   type_of_species='metabolites')

    plot_dataframe(data.proteins, 'genes', 'protein',
                   type_of_species='protein')

    e = Enrichr(exp_data=data)
    html_1 = e.run_key_dbs(data.proteomics_over_time,
                           data.proteomics_time_points,
                           save_name='proteomics_changed', create_html=True)

    html_2 = e.run_key_dbs(data.rna_over_time, data.rna_time_points,
                           save_name='rnaseq_changed', create_html=True)

    gene_html = 'genes.html'
    metabolite_html = 'metabolites.html'
    save_name = 'index_example_workflow'
    template_vars = {"title": 'Example',
                     "gene_table": gene_html,
                     "meta_table": metabolite_html,
                     "gene_up": 'Proteomics sign. changed',
                     "gene_up_table": html_1,
                     "rna_up": 'RNAseq sign. changed',
                     "rna_up_table": html_2,

                     }

    html_out = workflow_template.render(template_vars)

    with open('{}.html'.format(save_name), 'w') as f:
        f.write(html_out)

    quit()
    build_network(data.list_sig_proteins,
                  save_name='example_network',
                  all_measured_list=data.list_species,
                  metabolite_list=data.list_metabolites,
                  use_reactome=True, use_hmdb=True
                  )
