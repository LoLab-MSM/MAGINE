from magine.enrichment.enrichr import Enrichr, clean_tf_names
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
    df2 = e.run_samples(lists, ['1', '2', '3'], save_name='enrichr_test')
    assert df2.shape == (65, 10)


def test_multi_sample_plotting():
    up = exp_data.genes.sig.up_by_sample
    e.run_samples(up, ['1', '2', '3'],
                  save_name='enrichr_test',
                  exp_data=exp_data,
                  create_html=True,
                  out_dir='html_output2')


def test_set_of_dbs():
    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'],
                        database=['KEGG_2016', 'NCI-Nature_2016'],
                        save_name='t')
    assert df2.shape == (41, 10)


def test_tf_names():
    df = e.run(['BAX', 'BCL2', 'MCL1'], ['ARCHS4_TFs_Coexp', 'ChEA_2016'])
    tfs = clean_tf_names(df)
    for i in tfs['term_name']:
        assert '_' not in i


if __name__ == '__main__':
    test_single_run()
