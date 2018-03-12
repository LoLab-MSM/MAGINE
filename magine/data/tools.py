import numpy as np
import pandas as pd


def pivot_table(data, convert_to_log, index, columns, values):
    d_copy = data.copy()
    if convert_to_log:
        d_copy = log2_normalize_df(d_copy, values)

    array = pd.pivot_table(d_copy, index=index,
                           columns=columns,
                           values=values)
    array.fillna(0, inplace=True)

    array.sort_values(by=list(sorted(d_copy[columns].unique())),
                      ascending=False, inplace=True)
    return array


def log2_normalize_df(df, column):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of fold changes
    column : str
        column that contains fold change values

    Returns
    -------

    """
    tmp_df = df.copy()
    greater_than = tmp_df[column] > 0
    less_than = tmp_df[column] < 0
    tmp_df.loc[greater_than, column] = np.log2(tmp_df[greater_than][column])
    tmp_df.loc[less_than, column] = -np.log2(-tmp_df[less_than][column])
    return tmp_df


if __name__ == '__main__':
    index = 'protein'
    values = 'treated_control_fold_change'
    columns = 'time_points'
    x = [
        {'treated_control_fold_change': 1, 'protein': 'x', 'time_points': '1'},
        {'treated_control_fold_change': -2, 'protein': 'i',
         'time_points': '2'},
        {'treated_control_fold_change': -4, 'protein': 'b',
         'time_points': '1'},
        {'treated_control_fold_change': -4, 'protein': 'b',
         'time_points': '2'},
    ]
    d = pd.DataFrame(x)
    print(d['treated_control_fold_change'])
    print(d.head(10))
    df = pivot_table(d, True, index=index, values=values, columns=columns)
