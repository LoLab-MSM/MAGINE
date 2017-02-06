from magine.ontology.ontology_analysis import GoAnalysis




def test_go_heatmap():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp', slim='goslim_pir')
    go.analyze_data([['BAX', 'BAK1', 'BOK'],
                      ['BCL2L1', 'BCL2', 'MCL1'],],
                      #['BAD', 'BID', 'BBC3'],
                      #['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2','MCL1', 'BAD', 'BID', 'BBC3']],
                     aspect='P',labels=['1', '2'], savename='test_11', analyze=True)
    go.export_to_html(labels=[0, 1, ], html_name='test_out1')


def test_print_hit_by_index():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp', verbose=True)
    go.analysize_data([
        ['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD', 'BID', 'BBC3']],
        aspect='P', labels=['1', '2'], savename='test_11', analyze=False)
    go.retrieve_top_ranked(0, number=10)


def test_print_top_hits():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp2', verbose=False, save_png=True)
    # go.analyze_data([['ATR', 'CHEK1', 'CHEK2', 'ATM', 'TP53', 'BAX', 'CYCS', 'CASP3', 'CASP8','FAS']],
    #                 aspect='P', labels=['1'], savename='test_11', analyze=False)

    go.analyze_data([['BAX', 'BAK1', 'BOK',
                       'BCL2L1', 'BCL2', 'MCL1',
                       'BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD',
                       'BID', 'BBC3']],
                     #
                     aspect=None, labels=['1', '2'], savename='test_11',
                     analyze=True)
    go.export_to_html(labels=['1'], html_name='test')
    #go.print_ranked_over_time(savename='yolo', labels=None, number=9, create_plots=True)


#test_print_top_hits()
test_print_top_hits()
