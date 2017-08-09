import os

import pandas as pd
from nose.tools import raises

from magine.data.formatter import load_csv, load_directory

pd.set_option('display.width', 10000)
pd.set_option('max_colwidth', 100)

full_path = os.path.join(os.path.dirname(__file__), 'Data', 'hilic')


def test_load_directory():
    df = load_directory(full_path, dtype='metabolites')
    print(df.head(10))
    pivot_df = pd.pivot_table(df, columns='time',
                              index=['compound', 'compound_id', 'name', ])
    print(pivot_df.shape)
    assert pivot_df.shape == (2, 28)


@raises(OSError)
def test_directory_doesnt_exist():
    load_directory('DNE')


def test_load_csv():
    df = load_csv(full_path, 'hilic_1.csv')
    print(df.shape)
    assert df.shape == (3, 50)


if __name__ == '__main__':
    test_directory_doesnt_exist()
