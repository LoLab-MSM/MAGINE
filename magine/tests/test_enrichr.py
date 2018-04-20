from magine.enrichment.enrichr import Enrichr
from magine.tests.sample_experimental_data import exp_data

e = Enrichr()


def test_single_run():

    list_2 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
    df = e.run(list_2, 'GO_Biological_Process_2017')
    terms = df['term_name']
    assert len(terms) == 185


def test_multi_sample():

    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists,
                        ['1', '2', '3']
                        , save_name='enrichr_test')
    assert df2.shape == (85, 18)


def test_multi_sample_plotting():
    up = exp_data.proteomics_up_by_sample_id
    e.run_samples(up, ['1', '2', '3'],
                  save_name='enrichr_test',
                  exp_data=exp_data,
                  create_html=True,
                  out_dir='html_output2')


def test_set_of_dbs():
    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_sample_set_of_dbs(lists, ['1', '2', '3'],
                                  databases=['KEGG_2016', 'NCI-Nature_2016'],
                                  save_name='t')
    print(df2.shape)
    assert df2.shape == (74, 21)
