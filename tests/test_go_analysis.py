from Analysis.go_analysis import GoAnalysis

go = GoAnalysis(species='hsa',output_directory='tmp2',slim='goslim_pir')


def test_go_heatmap():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go.analysis_data([['BAX', 'BAK1', 'BOK'],
                      ['BCL2L1', 'BCL2', 'MCL1'],],
                      #['BAD', 'BID', 'BBC3'],
                      #['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2','MCL1', 'BAD', 'BID', 'BBC3']],
                     aspect='P',labels=['1', '2'], savename='test_11', analyze=True)
    go.export_to_html(labels=[0, 1,], x=[0,1], html_name='test_out1')

test_go_heatmap()
