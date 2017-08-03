try:
    import cPickle as pickle
except:  # python3 doesnt have cPickle
    import pickle
import os
from ast import literal_eval
from magine.data.storage import id_mapping_dir
import pandas as pd


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
        self.hmdb_accession_to_chemical_name = {}
        self.chemical_name_to_hmdb_accession = {}
        self.hmdb_accession_to_kegg = {}
        self.hmdb_accession_to_secondary_acc = {}
        self.kegg_to_hmdb_accession = {}
        self.hmdb_acc_to_secondary_acc = {}
        self.synonyms_to_hmdb = {}
        self.drugbank_to_hmdb = {}
        _filename = os.path.join(id_mapping_dir, 'hmdb_dataframe.csv')

        if not os.path.exists(_filename):
            from magine.mappings.databases.download_libraries import HMDB
            HMDB()
        hmdb_database = pd.read_csv(
                _filename, low_memory=False, encoding='utf-8',
                converters={
                    "protein_associations": lambda x: x.split("|"),
                    "cellular_locations":   lambda x: x.split("|"),
                    "biofunction":          lambda x: x.split("|"),
                    "synonyms":             lambda x: x.split("|"),
                    "secondary_accessions": lambda x: x.split("|"),
                }
        )

        self.database = hmdb_database.where((pd.notnull(hmdb_database)), None)
        self._instance_filename = os.path.join(id_mapping_dir,
                                               'hmdb_instance.p')
        try:
            self.reload()

            print('Loading class data')
        except:
            print('Initializing Chemical mapping')
            self.load()

    def load(self):

        self.hmdb_accession_to_chemical_name = self._to_dict("accession",
                                                             "name")
        self.chemical_name_to_hmdb_accession = self._to_dict("name",
                                                             "accession")
        self.hmdb_accession_to_kegg = self._to_dict("accession", "kegg_id")
        self.drugbank_to_hmdb = self._to_dict("drugbank_id", "accession")
        self.hmdb_acc_to_secondary_acc = self._to_dict('accession',
                                                       'secondary_accessions')
        self.kegg_to_hmdb_accession = self._to_dict("kegg_id", "accession")
        self.hmdb_accession_to_protein = self.convert_to_dict_from_list(
                "accession", "protein_associations")
        self.synonyms_to_hmdb = None
        self.save()

    def merge_accesions(self, accesion):
        secondary_acc = self.hmdb_acc_to_secondary_acc[accesion]
        all_accesions = secondary_acc.append(secondary_acc)
        all_accesions = set(all_accesions)
        return '|'.join(i for i in all_accesions)

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
        return {k: list(v) for k, v in self.database.groupby(key)[value]}

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
        >>> cm.check_synonym_dict(term='dodecene', format_name='accession')
        HMDB59874

        """

        for index, row in self.database.iterrows():
            each = row.synonyms
            if isinstance(each, list):
                if term in each:
                    return self.database.iloc[index][format_name]
            else:
                print(type(each), repr(each))
                each = literal_eval(str(each))
                print(type(each), repr(each))
                quit()
                if term in each:
                    return self.database.iloc[index][format_name]
        return None

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
        with open(self._instance_filename, 'wb') as f:
            f.write(pickle.dumps(self.__dict__, protocol=-1))

    def reload(self):
        """ loads class instance

        Returns
        -------

        """

        with open(self._instance_filename, 'rb') as f:
            data = f.read()

        try:
            self.__dict__ = pickle.loads(data, encoding='utf-8')
        except:
            self.__dict__ = pickle.loads(data)


if __name__ == "__main__":
    cm = ChemicalMapper()
    # print(list(cm.hmdb_accession_to_protein)[:10])
    # print(cm.hmdb_accession_to_protein['HMDB00005'])
    print(cm.check_synonym_dict(term='dodecene', format_name='accession'))
