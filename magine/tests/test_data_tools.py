import pandas as pd
from nose.tools import raises

from magine.data import Data


class ConcentrationData(Data):
    def __init__(self, *args, **kwargs):
        super(ConcentrationData, self).__init__(*args, **kwargs)
        self._index = 'protein'

    @property
    def _constructor(self):
        return ConcentrationData


def load_concentration(file):
    dataframe = pd.read_csv(file)
    return ConcentrationData(dataframe)


def test_pivot():
    index = 'protein'
    values = 'treated_control_fold_change'
    columns = 'time_points'
    x = [
        {values: 1, index: 'x', columns: '1'},
        {values: -2, index: 'i', columns: '2'},
        {values: -4, index: 'b', columns: '1'},
        {values: -4, index: 'b', columns: '2'},
    ]

    d = ConcentrationData(x)
    df = d.pivoter(True, values=values, columns=columns)

    assert df.shape == (3, 2)


@raises(AssertionError)
def test_min_raises():
    index = 'protein'
    values = 'treated_control_fold_change'
    columns = 'time_points'
    x = [
        {values: 1, index: 'x', columns: '1'},
        {values: -2, index: 'i', columns: '2'},
        {values: -4, index: 'b', columns: '1'},
        {values: -4, index: 'b', columns: '2'},
    ]

    d = ConcentrationData(x)
    d.filter_by_minimum_sig_columns(index=index, columns=columns)


def test_filter_min():
    index = 'protein'
    values = 'treated_control_fold_change'
    columns = 'time_points'
    flag = 'significant_flag'
    x = [
        {values: 1, index: 'x', columns: '1', flag: True, 'x': 2},
        {values: -2, index: 'i', columns: '2', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '1', flag: True, 'x': 1},
        {values: -4, index: 'b', columns: '2', flag: True, 'x': 1},
    ]

    d = ConcentrationData(x)
    df = d.filter_by_minimum_sig_columns(columns=columns, min_terms=1)
    assert df.shape[0] == 4
    df = d.filter_by_minimum_sig_columns(columns=columns, min_terms=2)
    assert df.shape[0] == 2
    df = d.filter_by_minimum_sig_columns(columns=columns, min_terms=3)
    assert df.shape[0] == 0

    x = [
        {values: 1, index: 'x', columns: '1', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '1', flag: True, 'x': 1},
        {values: -2, index: 'i', columns: '2', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '2', flag: True, 'x': 1},
        {values: -4, index: 'b', columns: '2', flag: False, 'x': 2},
    ]

    d = ConcentrationData(x)
    df = d.filter_by_minimum_sig_columns(index=[index, 'x'], columns=columns,
                                         min_terms=2)

    assert df.shape == (3, 5)


def test_log_normal():
    x = [['a', 2],
         ['b', -2],
         ['c', -16],
         ['d', 16]]

    df = ConcentrationData(x, columns=['name', 'value'])
    new_df = df.log2_normalize_df('value')

    assert new_df[new_df['name'] == 'a']['value'].values == [1.]
    assert new_df[new_df['name'] == 'c']['value'].values == [-4.]
