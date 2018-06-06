import numpy as np
import pandas as pd


class Data(pd.DataFrame):
    _metadata = ['added_property']
    _index = False

    def __init__(self, *args, **kwargs):
        super(Data, self).__init__(*args, **kwargs)
        # self._index = None

    @property
    def _constructor(self):
        return Data

    def pivoter(self, convert_to_log, columns, values, index=None,
                fill_value=None, min_sig=0):
        """

        Parameters
        ----------
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
        d_copy = self.copy()
        if index is None:
            index = self._index

        if convert_to_log:
            d_copy = log2_normalize_df(d_copy, values)

        if min_sig:
            assert isinstance(min_sig, int)
            if 'significant_flag' not in d_copy.columns:
                print('In order to filter based on minimum sig figs, '
                      'please add a column of signficant terms with '
                      'a tag of "significant"')

            d_copy = self.filter_by_minimum_sig_columns(index=index,
                                                        columns=columns,
                                                        min_terms=min_sig)

        array = pd.pivot_table(d_copy, index=index,
                               fill_value=fill_value, columns=columns,
                               values=values)

        if isinstance(columns, list):
            x = sorted(tuple(map(tuple, d_copy[columns].values)))
            array.sort_values(by=x, ascending=False, inplace=True)
        elif isinstance(columns, str):
            array.sort_values(by=list(sorted(d_copy[columns].unique())),
                              ascending=False, inplace=True)
        return array

    def filter_by_minimum_sig_columns(self, columns, index=None, min_terms=3,
                                      inplace=False):
        if index is None:
            index = self._index
        # create safe copy of array
        new_data = self.copy()

        # get list of columns
        cols_to_check = list(new_data[columns].unique())
        if 'significant' in new_data.columns:
            flag = 'significant'
        elif 'significant_flag' in new_data.columns:
            flag = 'significant_flag'
        else:
            flag = None
        assert flag in new_data.columns, 'Requires significant_flag column'
        # pivot
        sig = pd.pivot_table(new_data,
                             index=index,
                             fill_value=0,
                             values=flag,
                             columns=columns
                             )[cols_to_check]

        # convert everything thats not 0 to 1
        sig[sig > 0] = 1
        sig = sig[sig.T.sum() >= min_terms]
        if isinstance(index, list):
            keepers = {i[0] for i in sig.index.values}
            new_data = new_data[new_data[index[0]].isin(keepers)]
        elif isinstance(index, str):
            keepers = {i for i in sig.index.values}
            new_data = new_data.loc[new_data[index].isin(keepers)]
        else:
            print("Index is not a str or a list. What is it?")

        if inplace:
            self._update_inplace(new_data)
        else:
            return new_data


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
