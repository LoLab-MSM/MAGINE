import gzip
import itertools
import os

import pandas as pd

from magine.data.storage import id_mapping_dir

try:
    import cPickle as pickle
except:  # python3 doesnt have cPickle
    import pickle


class ChemicalMapper(object):
    """
    converts chemical species ids
    Database was creating using HMDB


    Attributes
    ----------
    hmdb_accession_to_chemical_name : dict
    chemical_name_to_hmdb_accession : dict
    hmdb_accession_to_kegg : dict
    kegg_to_hmdb_accession : dict
    hmdb_accession_to_protein : dict
    synonyms_to_hmdb : dict

    """

    valid_columns = ['kegg_id', 'name', 'accession', 'chebi_id', 'inchikey',
                     'chemspider_id', 'biocyc_id', 'synonyms', 'iupac_name',
                     'pubchem_compound_id', 'protein_associations',
                     'ontology', 'drugbank_id', 'chemical_formula',
                     'smiles', 'metlin_id', 'average_molecular_weight',
                     'secondary_accessions'
                     ]

    def __init__(self):

        self.database = None
        self.hmdb_accession_to_chemical_name = dict()
        self.chemical_name_to_hmdb_accession = dict()
        self.hmdb_accession_to_kegg = dict()
        self.kegg_to_hmdb_accession = dict()
        self.synonyms_to_hmdb = dict()
        self.drugbank_to_hmdb = dict()
        self.hmdb_accession_to_protein = dict()
        self.hmdb_main_accession_to_protein = dict()

        self._instance_filename = os.path.join(id_mapping_dir,
                                               'hmdb_instance.p.gz')

        try:
            self.reload()
            print('Loading class data')
        except:
            print('Initializing Chemical mapping')
            self.load()

    def load(self):

        _filename = os.path.join(id_mapping_dir, 'hmdb_dataframe.csv.gz')

        if not os.path.exists(_filename):
            from magine.mappings.databases.download_libraries import HMDB
            HMDB()
        hmdb_database = pd.read_csv(
            _filename, low_memory=False, encoding='utf-8',
            converters={
                "protein_associations": lambda x: x.split("|"),
                "cellular_locations": lambda x: x.split("|"),
                "biofunction": lambda x: x.split("|"),
                "synonyms": lambda x: x.split("|"),
                # "secondary_accessions": lambda x: x.split("|"),
            }
        )
        self.database = hmdb_database.where((pd.notnull(hmdb_database)), None)

        self.database['main_accession'] = self.database['accession']
        sub_db = self.database[
            self.database['secondary_accessions'].str.contains('|', na=False)]
        new_df = tidy_split(sub_db, 'secondary_accessions', '|')
        new_df['accession'] = new_df['secondary_accessions']
        self.database = pd.concat([self.database, new_df])

        self.hmdb_accession_to_chemical_name = self._to_dict("accession",
                                                             "name")
        self.chemical_name_to_hmdb_accession = self._to_dict("name",
                                                             "main_accession")
        self.hmdb_accession_to_kegg = self._to_dict("accession", "kegg_id")

        self.kegg_to_hmdb_accession = self._to_dict("kegg_id",
                                                    "main_accession")
        self.hmdb_accession_to_protein = self._from_list_dict("accession",
                                                       "protein_associations")

        self.hmdb_main_accession_to_protein = \
            self._from_list_dict("main_accession", "protein_associations")

        self.drugbank_to_hmdb = self._to_dict('drugbank_id', 'main_accession')
        self.synonyms_to_hmdb = None

        self.save()

    def _to_dict(self, key, value):
        """ creates a dictionary with a list of values for each key

        Parameters
        ----------
        key : str
        value : str

        Returns
        -------
        dict

        """
        return {k: sorted(set(list(v))) for k, v in self.database.groupby(key)[value]}

    def _from_list_dict(self, key, value):
        return {k: sorted(set(list(itertools.chain.from_iterable(v.tolist()))))
                for k, v in self.database.groupby(key)[value]}

    def convert_to_dict_from_list(self, key, value):
        """ creates a dictionary from hmdb with a list of values for each key

        Parameters
        ----------
        key : str
        value : str

        Returns
        -------
        dict

        """

        tmp_dict = {}
        for k, v in self.database.groupby(key)[value]:
            v = v.tolist()[0]
            if isinstance(v, list):
                if len(v) > 0:
                    if v[0] is not None:
                        tmp_dict[k] = v

        return tmp_dict

    def check_synonym_dict(self, term, format_name):
        """ checks hmdb database for synonyms and returns formatted name

        Parameters
        ----------
        term : str
        format_name : str

        Returns
        -------
        dict


        Examples
        --------
        >>> cm = ChemicalMapper()
        >>> cm.check_synonym_dict(term='dodecene', format_name='main_accession')
        ['HMDB0000933', 'HMDB0059874']

        """
        synonyms = self.database.copy()
        synonyms['synonyms'] = synonyms['synonyms'].apply(','.join)
        synonyms['synonyms'] = synonyms['synonyms'].str.lower()

        hits = synonyms[synonyms['synonyms'].str.contains(term.lower())]
        matches = sorted(set(hits[format_name].values))
        return matches

    def print_info(self):
        """ print information about the dataframe

        Returns
        -------

        """
        print('Number of HMDB accessions = {0}'.format(
            len(self.database['accession'].unique())))
        print('Number of unique KEGG ids = {0}'.format(
            len(self.hmdb_accession_to_kegg.keys())))
        print('Number of HMDB to KEGG mappings = {0}'.format(
            len(self.kegg_to_hmdb_accession.values())))

    def save(self):
        """ save class instance

        Returns
        -------

        """
        print('Saving class data')
        with gzip.open(self._instance_filename, 'wb') as f:
            f.write(pickle.dumps(self.__dict__, protocol=0))

    def reload(self):
        """ loads class instance

        Returns
        -------

        """
        try:
            with gzip.open(self._instance_filename, 'rb') as f:
                data = f.read()
            self.__dict__ = pickle.loads(data)
        except:
            with gzip.open(self._instance_filename, 'rb') as f:
                data = f.read()
            self.__dict__ = pickle.loads(data, encoding='latin1')


def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


if __name__ == "__main__":
    cm = ChemicalMapper()
    # print(list(cm.hmdb_accession_to_protein)[:10])
    # print(cm.hmdb_accession_to_protein['HMDB00005'])
    print(cm.check_synonym_dict(term='dodecene', format_name='main_accession'))
