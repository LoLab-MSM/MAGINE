import pandas as pd
from nose.tools import raises

import magine.data.tools as tools


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

    d = pd.DataFrame(x)
    df = tools.pivot_table(d, True, index=index, values=values,
                           columns=columns)
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

    d = pd.DataFrame(x)
    tools.filter_by_minimum_sig_columns(d, index=index, columns=columns)


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

    d = pd.DataFrame(x)
    df = tools.filter_by_minimum_sig_columns(d, index=index, columns=columns,
                                             min_terms=1)
    assert df.shape[0] == 4

    df = tools.filter_by_minimum_sig_columns(d, index=index, columns=columns,
                                             min_terms=2)

    assert df.shape[0] == 2

    df = tools.filter_by_minimum_sig_columns(d, index=index, columns=columns,
                                             min_terms=3)

    assert df.shape[0] == 0

    x = [
        {values: 1, index: 'x', columns: '1', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '1', flag: True, 'x': 1},
        {values: -2, index: 'i', columns: '2', flag: True, 'x': 2},
        {values: -4, index: 'b', columns: '2', flag: True, 'x': 1},
        {values: -4, index: 'b', columns: '2', flag: False, 'x': 2},
    ]

    d = pd.DataFrame(x)
    df = tools.filter_by_minimum_sig_columns(d, index=[index, 'x'],
                                             columns=columns, min_terms=2)

    assert df.shape == (3, 5)


def test_log_normal():
    x = [['a', 2],
         ['b', -2],
         ['c', -16],
         ['d', 16]]
    df = pd.DataFrame(x, columns=['name', 'value'])
    new_df = tools.log2_normalize_df(df, 'value')

    assert new_df[new_df['name'] == 'a']['value'].values == [1.]
    assert new_df[new_df['name'] == 'c']['value'].values == [-4.]
