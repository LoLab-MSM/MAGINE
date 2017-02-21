from magine.ontology.ontology_analysis import GoAnalysis
from sample_experimental_data import exp_data



def test_go_heatmap():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp')
    go.analyze_data([['BAX', 'BAK1', 'BOK'],
                     ['BCL2L1', 'BCL2', 'MCL1'], ],
                    # ['BAD', 'BID', 'BBC3'],
                      #['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2','MCL1', 'BAD', 'BID', 'BBC3']],
                    labels=['1', '2'], savename='test_11', )
    go.export_to_html(labels=[0, 1, ], html_name='test_out1')


def test_print_hit_by_index():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp', verbose=True)
    go.analyze_data([
        ['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD', 'BID', 'BBC3']],
            labels=['1', '2'], savename='test_11')
    go.retrieve_top_ranked(0, number=10)


def test_print_top_hits():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp2', verbose=False,)
    # go.analyze_data([['ATR', 'CHEK1', 'CHEK2', 'ATM', 'TP53', 'BAX', 'CYCS', 'CASP3', 'CASP8','FAS']],
    #                 aspect='P', labels=['1'], savename='test_11', analyze=False)

    go.analyze_data([['BAX', 'BAK1', 'BOK',
                       'BCL2L1', 'BCL2', 'MCL1',
                       'BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD',
                       'BID', 'BBC3']],
                    labels=['1', '2'], savename='test_11',
                    )
    go.export_to_html(labels=['1'], html_name='test')
    go.print_ranked_over_time(savename='yolo', labels=None, number=9,
                              create_plots=True)


def test_html():
    go = GoAnalysis(species='hsa', output_directory='html_output',
                    verbose=False,
                    experimental_data=exp_data)

    df = go.create_enrichment_array(exp_data.proteomics_up_over_time[:2],
                                    exp_data.timepoints[:2],
                                    save_name='test_df')
    go.write_table_to_html(save_name='index')

    go.write_table_to_html(df, save_name='index2')


if __name__ == '__main__':
    test_html()
#test_print_top_hits()
    # test_print_top_hits()
