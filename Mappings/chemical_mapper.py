import os

import pandas as pd

directory = os.path.dirname(__file__)


class ChemicalMapper:
    """
    converts chemical species
    Datbase was creating by pulling down everything from HMDB
    """
    categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id', 'biocyc_id', 'synonyms',
                  'pubchem_compound_id', 'iupac_name', 'protein_associations']

    def __init__(self):
        hmdb_database = pd.read_pickle(os.path.join(directory, 'hmdb_dataframe.p'))
        self.database = hmdb_database
        self.hmdb_accession_to_chemical_name = self.convert_to_dict("accession", "name")
        self.chemical_name_to_hmdb_accession = self.convert_to_dict("name", "accession")
        self.hmdb_accession_to_kegg = self.convert_to_dict("accession", "kegg_id")
        self.kegg_to_hmdb_accession = self.convert_to_dict("kegg_id", "accession")
        self.hmdb_accession_to_protein = self.convert_to_dict("accession", "protein_associations")
        self.synonyms_to_hmdb = None

    def convert_to_dict(self, key, value):
        """ creates a dictionary from hmdb with a list of values for each key

        :param key:
        :param value:
        :return:
        """
        return {k: list(v) for k, v in self.database.groupby(key)[value]}

    def check_synonym_dict(self, term, format_name):
        """ checks hmdb database for synonyms and returns formatted name

        :param term:
        :param format_name:
        :return:
        """
        for index, row in self.database.iterrows():
            each = row.synonyms
            if type(each) == list:
                if term in each:
                    return self.database.iloc[index][format_name]
        return None

    def print_info(self):
        """ print information about the dataframe

        :return:
        """

        print('Number of HMDB accessions = {0}'.format(len(self.database['accession'].unique())))
        print('Number of unique KEGG ids = {0}'.format(len(self.hmdb_accession_to_kegg.keys())))
        print('Number of HMDB to KEGG mappings = {0}'.format(len(self.kegg_to_hmdb_accession.values())))