import numpy as np
import pandas as pd


def pivot_table(data, convert_to_log, index, columns, values, fill_value=None,
                min_sig=0):
    """

    Parameters
    ----------
    data : pandas.DataFrame
    convert_to_log : bool
        Convert values column to log2
    index : str
        Index for pivot table
    columns : str
        Columns to pivot
    values : str
        Values of pivot table
    fill_value : float, optional
        Fill pivot table nans with
    min_sig : int
        Required number of significant terms to keep in a row, default 0
    Returns
    -------

    """
    d_copy = data.copy()
    if convert_to_log:
        d_copy = log2_normalize_df(d_copy, values)

    if min_sig:
        assert isinstance(min_sig, int)
        if 'significant_flag' not in d_copy.columns:
            print('In order to filter based on minimum sig figs, '
                  'please add a column of signficant terms with '
                  'a tag of "significant"')

        d_copy = filter_by_minimum_sig_columns(d_copy, index=index,
                                               columns=columns,
                                               min_terms=min_sig)

    array = pd.pivot_table(d_copy, index=index, fill_value=fill_value,
                           columns=columns, values=values)
    if isinstance(columns, list):
        x = sorted(tuple(map(tuple, d_copy[columns].values)))
        array.sort_values(by=x, ascending=False, inplace=True)
    elif isinstance(columns, str):
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


# TODO create test and example
def filter_by_minimum_sig_columns(data_frame, index, columns,
                                  min_terms=3):
    tmp_array = data_frame.copy()

    sig = pd.pivot_table(tmp_array,
                         index=index,
                         fill_value=0,
                         values='significant_flag',
                         columns=columns
                         )

    cols_to_check = list(tmp_array[columns].unique())
    sig = sig[cols_to_check]
    sig = sig[cols_to_check].T.sum()
    sig = sig[sig > min_terms]
    if isinstance(index, list):
        keepers = {i[0] for i in sig.index.values}
        return tmp_array[tmp_array[index[0]].isin(keepers)].copy()
    elif isinstance(index, str):
        keepers = {i for i in sig.index.values}
        return tmp_array[tmp_array[index].isin(keepers)].copy()
    else:
        print("Index is not a str or a list. What is it?")


def filter_by_minimum_sig_columns_old(data_frame, index, values, columns,
                                      min_terms=3):
    tmp_array = data_frame.copy()

    sig = pd.pivot_table(tmp_array[tmp_array['significant_flag']],
                         index=index,
                         fill_value=0.0,
                         values=values,
                         columns=columns
                         )

    cols_to_check = list(tmp_array[columns].unique())
    sig = (sig[cols_to_check].T == 0.).sum()
    sig = sig[sig < min_terms]
    if isinstance(index, list):
        keepers = {i[0] for i in sig.index.values}
        return tmp_array[tmp_array[index[0]].isin(keepers)].copy()
    elif isinstance(index, str):
        keepers = {i for i in sig.index.values}
        return tmp_array[tmp_array[index].isin(keepers)].copy()
    else:
        print("Index is not a str or a list. What is it?")


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
