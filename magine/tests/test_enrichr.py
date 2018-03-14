from magine.ontology.enrichr import Enrichr
from magine.tests.sample_experimental_data import exp_data


def test_enrichr_go():
    e = Enrichr()
    list_2 = ['CASP3', 'CASP6', 'FAS', 'FADD', 'CASP8', 'CFLAR', 'BFAR', 'BAD',
              'BID', 'PMAIP1', 'MCL1', 'BCL2', 'BCL2L1', 'BAX', 'BAK1',
              'DIABLO', 'CYCS', 'PARP1', 'APAF1', 'XIAP']
    df = e.run(list_2, 'GO_Biological_Process_2017')
    terms = df['term_name']
    assert len(terms) == 185


def test_multi_sample():
    e = Enrichr()
    lists = [['BAX', 'BCL2', 'CASP3'],
             ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists,
                        ['1', '2', '3']
                        , save_name='enrichr_test')
    assert df2.shape == (85, 21)


def test_multi_sample_plotting():
    e = Enrichr()
    up = exp_data.proteomics_up_by_sample_id
    df2 = e.run_samples(up,
                        ['1', '2', '3'],
                        save_name='enrichr_test',
                        exp_data=exp_data,
                        create_html=True,
                        out_dir='html_output2')

