import numpy as np


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