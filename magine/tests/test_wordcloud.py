import os

import pandas as pd

from magine.plotting.wordcloud_tools import create_wordcloud


def test_filter():
    dirname = os.path.join(os.path.dirname(__file__),
                           'Data', 'enrichr_test_enrichr.csv')

    df = pd.read_csv(dirname)
    x, d = create_wordcloud(df)

    x.plot('test_wc')


if __name__ == '__main__':
    test_filter()
