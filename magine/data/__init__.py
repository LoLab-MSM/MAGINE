import pandas as pd

from magine.data.tools import log2_normalize_df


class Data(pd.DataFrame):
    _metadata = ['added_property']
    _index = False

    def __init__(self, *args, **kwargs):
        super(Data, self).__init__(*args, **kwargs)
        # self._index = None

    @property
    def _constructor(self):
        return Data

    def pivoter(self, convert_to_log, columns, values, fill_value=None,
                min_sig=0):
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
        if convert_to_log:
            d_copy = log2_normalize_df(d_copy, values)

        if min_sig:
            assert isinstance(min_sig, int)
            if 'significant_flag' not in d_copy.columns:
                print('In order to filter based on minimum sig figs, '
                      'please add a column of signficant terms with '
                      'a tag of "significant"')

            d_copy = self.filter_by_minimum_sig_columns(columns=columns,
                                                        min_terms=min_sig)

        array = pd.pivot_table(d_copy, index=self._index,
                               fill_value=fill_value, columns=columns,
                               values=values)

        if isinstance(columns, list):
            x = sorted(tuple(map(tuple, d_copy[columns].values)))
            array.sort_values(by=x, ascending=False, inplace=True)
        elif isinstance(columns, str):
            array.sort_values(by=list(sorted(d_copy[columns].unique())),
                              ascending=False, inplace=True)
        return array

    def filter_by_minimum_sig_columns(self, columns, min_terms=3):
        # create safe copy of array
        tmp_array = self.copy()

        # get list of columns
        cols_to_check = list(tmp_array[columns].unique())

        assert 'significant_flag' in tmp_array, 'Requires significant_flag column'
        # pivot
        sig = pd.pivot_table(tmp_array,
                             index=self._index,
                             fill_value=0,
                             values='significant_flag',
                             columns=columns
                             )[cols_to_check]

        # convert everything thats not 0 to 1
        sig[sig > 0] = 1
        sig = sig[sig.T.sum() >= min_terms]
        if isinstance(self._index, list):
            keepers = {i[0] for i in sig.index.values}
            return tmp_array[tmp_array[self._index[0]].isin(keepers)].copy()
        elif isinstance(self._index, str):
            keepers = {i for i in sig.index.values}
            return tmp_array[tmp_array[self._index].isin(keepers)].copy()
        else:
            print("Index is not a str or a list. What is it?")


class EnrichmentData(Data):

    def __init__(self, *args, **kwargs):
        super(EnrichmentData, self).__init__(*args, **kwargs)
        self._index = 'term_name'

    @property
    def _constructor(self):
        return EnrichmentData


class ConcentrationData(Data):

    def __init__(self, *args, **kwargs):
        super(ConcentrationData, self).__init__(*args, **kwargs)
        self._index = 'protein'

    @property
    def _constructor(self):
        return ConcentrationData


def load_enrichment_csv(file):
    dataframe = pd.read_csv(file)
    return EnrichmentData(dataframe)


def load_concentration(file):
    dataframe = pd.read_csv(file)
    return ConcentrationData(dataframe)
