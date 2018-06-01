import os

import pandas as pd

import magine.plotting.wordcloud_tools as wt
from magine.enrichment import load_enrichment_csv


def test_filter():
    df = load_enrichment_csv(os.path.join(os.path.dirname(__file__),
                                          'Data', 'enrichr_test_enrichr.csv'))
    df2 = df.copy()
    df['db'] = 'kegg'
    df2['db'] = 'reactome'
    df2['sample_id'] = 2

    df = pd.concat([df, df2])

    x = wt.create_wordcloud(df)

    x.plot('test_wc')

    wt.word_cloud_from_array(df,
                             sample_ids=[1, 2]
                             )


if __name__ == '__main__':
    test_filter()
