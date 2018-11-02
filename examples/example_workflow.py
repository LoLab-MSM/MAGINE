from magine.data.experimental_data import ExperimentalData
from magine.enrichment.enrichr import Enrichr
from magine.html_templates.html_tools import workflow_template
from magine.networks.network_generator import build_network


if __name__ == '__main__':
    data = ExperimentalData('Data/norris_et_al_2017_cisplatin_data.csv.gz')

    # want data organized
    # data
    # proteins
    # metabolomics

    data.compounds.plot_all('metabolites', 'metabolites')

    data.proteins.plot_all('genes', 'protein')

    e = Enrichr()
    html_1 = e.run_samples(data.proteins.sig.by_sample,
                           data.proteins.sig.sample_ids,
                           save_name='proteomics_changed', create_html=True)

    html_2 = e.run_samples(data.rna.sig.by_sample, data.rna.sig.sample_ids,
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

    build_network(data.species.sig.id_list,
                  save_name='example_network',
                  all_measured_list=data.species.id_list,
                  use_reactome=True, use_hmdb=True
                  )
