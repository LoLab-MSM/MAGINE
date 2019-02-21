import numpy as np


def log2_normalize_df(df, column, new_col=None):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of fold changes
    column : str
        column that contains fold change values
    new_col : str, None
        Name of new column to place fold_changes values. if not provided it
        will overwrite column values with new values
    Returns
    -------

    """
    tmp_df = df.copy()
    greater_than = tmp_df[column] > 0
    less_than = tmp_df[column] < 0
    if new_col is not None:
        if not isinstance(new_col, str):
            raise AssertionError()
        tmp_df.loc[greater_than, new_col] = np.log2(
            tmp_df[greater_than][column])
        tmp_df.loc[less_than, new_col] = -np.log2(-tmp_df[less_than][column])
    else:
        tmp_df.loc[greater_than, column] = np.log2(
            tmp_df[greater_than][column])
        tmp_df.loc[less_than, column] = -np.log2(-tmp_df[less_than][column])
    return tmp_df
