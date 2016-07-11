from Analysis.go_analysis import GoAnalysis




def test_go_heatmap():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp2', slim='goslim_pir')
    go.analysis_data([['BAX', 'BAK1', 'BOK'],
                      ['BCL2L1', 'BCL2', 'MCL1'],],
                      #['BAD', 'BID', 'BBC3'],
                      #['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2','MCL1', 'BAD', 'BID', 'BBC3']],
                     aspect='P',labels=['1', '2'], savename='test_11', analyze=True)
    go.export_to_html(labels=[0, 1,], x=[0,1], html_name='test_out1')


def test_print_hit_by_index():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp2', verbose=True)
    go.analysis_data([
        ['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD', 'BID', 'BBC3']],
        aspect='P', labels=['1', '2'], savename='test_11', analyze=False)
    go.retrieve_top_ranked(0, number=10)


def test_print_top_hits():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go = GoAnalysis(species='hsa', output_directory='tmp2', verbose=False)
    go.analysis_data([['BAX', 'BAK1', 'BOK'],
                      ['BCL2L1', 'BCL2', 'MCL1'],
                      ['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2', 'MCL1', 'BAD', 'BID', 'BBC3']],
                     aspect='P', labels=['1', '2'], savename='test_11', analyze=False)
    go.print_ranked_over_time(savename='yolo', labels=None, number=20, create_plots=True)


test_print_top_hits()
