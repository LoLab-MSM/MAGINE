from magine.ontology.enrichr import Enrichr


def test_enrichr_go():
    e = Enrichr()
    g_list = ['PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2',
              'P2RX7', 'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1',
              'ANAPC16', 'TMCC1', 'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2',
              'PLBD2', 'LARP7', 'TECPR2', 'ZNF302', 'CUX1', 'MOB2', 'CYTH2',
              'SEC22C', 'EIF4E3', 'ROBO2', 'ADAMTS9-AS2', 'CXXC1', 'LINC01314',
              'ATF7', 'ATP5F1']

    df = e.run(g_list, 'GO_Biological_Process_2017')
    terms = df['term_name']
    assert len(terms) == 98
