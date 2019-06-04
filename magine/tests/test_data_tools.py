import pandas as pd
from nose.tools import raises, ok_

from magine.data.base import BaseData


class ConcentrationBaseData(BaseData):
    def __init__(self, *args, **kwargs):
        super(ConcentrationBaseData, self).__init__(*args, **kwargs)
        self._index = 'protein'

    @property
    def _constructor(self):
        return ConcentrationBaseData


def load_concentration(file):
    dataframe = pd.read_csv(file)
    return ConcentrationBaseData(dataframe)


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

    d = ConcentrationBaseData(x)
    df = d.pivoter(True, values=values, columns=columns)

    ok_(df.shape == (3, 2))


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

    d = ConcentrationBaseData(x)
    d.require_n_sig(index=index, columns=columns)


def test_filter_min():
    index = 'protein'
    values = 'treated_control_fold_change'
    columns = 'time_points'
    flag = 'significant'
    x = [
        {values: 1, index: 'x', columns: '1', flag: True, 'x': 2},
        {values: -2, index: 'i', columns: '2', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '1', flag: True, 'x': 1},
        {values: -4, index: 'b', columns: '2', flag: True, 'x': 1},
    ]

    d = ConcentrationBaseData(x)
    df = d.require_n_sig(columns=columns, n_sig=1)
    ok_(df.shape[0] == 4)
    df = d.require_n_sig(columns=columns, n_sig=2)
    ok_(df.shape[0] == 2)
    df = d.require_n_sig(columns=columns, n_sig=3)
    ok_(df.shape[0] == 0)

    x = [
        {values: 1, index: 'x', columns: '1', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '1', flag: True, 'x': 1},
        {values: -2, index: 'i', columns: '2', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '2', flag: True, 'x': 1},
        {values: -4, index: 'b', columns: '2', flag: False, 'x': 2},
    ]

    d = ConcentrationBaseData(x)

    df = d.require_n_sig(index=index, columns=columns,
                         n_sig=2)

    ok_(df.shape == (3, 5))


def test_log_normal():
    x = [['a', 2],
         ['b', -2],
         ['c', -16],
         ['d', 16]]

    df = ConcentrationBaseData(x, columns=['name', 'value'])
    new_df = df.log2_normalize_df('value')

    ok_(new_df[new_df['name'] == 'a']['value'].values == [1.])
    ok_(new_df[new_df['name'] == 'c']['value'].values == [-4.])
