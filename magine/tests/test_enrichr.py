import tempfile

from nose.tools import ok_

from magine.data.experimental_data import ExperimentalData
from magine.enrichment.enrichr import Enrichr, clean_drug_dbs, clean_tf_names, \
    get_background_list, get_libraries, run_enrichment_for_project
from magine.tests.sample_experimental_data import exp_data

e = Enrichr()


def test_single_run():
    e.print_valid_libs()
    list_2 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
    df = e.run(list_2, 'GO_Biological_Process_2017')
    terms = df['term_name']
    ok_(len(terms) == 185)


def test_libaries():
    get_libraries()


def test_project():
    slimmed = exp_data.species.copy()
    slimmed = slimmed.loc[slimmed.source.isin(['label_free', 'silac'])]
    slimmed = slimmed.loc[slimmed.sample_id.isin(['Time_1', 'Time_2', ])]
    slimmed = ExperimentalData(slimmed)
    run_enrichment_for_project(slimmed, 'test',
                               databases=['KEGG_2016'])
    slimmed = exp_data.species.copy()
    slimmed = slimmed.loc[slimmed.source.isin(['label_free'])]
    slimmed = slimmed.loc[slimmed.sample_id.isin(['Time_1', 'Time_2', ])]
    slimmed = ExperimentalData(slimmed)
    run_enrichment_for_project(slimmed, 'test',
                               databases=['KEGG_2016'])


def test_get_gene_set_lib():
    get_background_list('Phosphatase_Substrates_from_DEPOD')


def test_clean_drug_dbs():
    list_2 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
    df = e.run(list_2, ['Drug_Perturbations_from_GEO_2014',
                        'LINCS_L1000_Chem_Pert_up'
                        ])
    df = clean_drug_dbs(df)
    ok_(len(df['term_name']) == 6243)

    ok_(len(df.sig['term_name']) == 249)


def test_multi_sample():
    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'], save_name='enrichr_test')
    df2.to_csv('enrichment.csv')

    ok_(df2.shape == (110, 11))


def test_multi_sample_plotting():
    up = exp_data.genes.sig.up_by_sample
    out_dir = tempfile.mkdtemp()
    e.run_samples(up, ['1', '2', '3'],
                  gene_set_lib=['Human_Phenotype_Ontology',
                            'MGI_Mammalian_Phenotype_2017'],
                  save_name='enrichr_test',
                  exp_data=exp_data,
                  create_html=True,
                  pivot=True,
                  out_dir=out_dir)


def test_set_of_dbs():
    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'],
                        gene_set_lib=['KEGG_2016', 'NCI-Nature_2016'],
                        save_name='t')
    ok_(df2.shape == (120, 11))


def test_tf_names():
    df = e.run(['BAX', 'BCL2', 'MCL1'], ['ARCHS4_TFs_Coexp', 'ChEA_2016'])
    tfs = clean_tf_names(df)
    for i in tfs['term_name']:
        ok_('_' not in i)


if __name__ == '__main__':
    test_single_run()
