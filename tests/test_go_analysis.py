from Analysis.go_analysis import GoAnalysis

go = GoAnalysis(species='hsa')


def test_go_heatmap():
    """
    tests the ability to create a heatmap from list of genes
    :return None
    """
    go.analysis_data([['BAX', 'BAK1', 'BOK'],
                      ['BCL2L1', 'BCL2', 'MCL1'],
                      ['BAD', 'BID', 'BBC3'],
                      ['BAX', 'BAK1', 'BOK', 'BCL2L1', 'BCL2','MCL1', 'BAD', 'BID', 'BBC3']],
                     aspect='P',labels=['1', '2'], savename='test_1', analyze=True)


test_go_heatmap()
