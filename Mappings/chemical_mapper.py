try:
    import cPickle
except:  # python3 doesnt have cPickle
    import pickle as cPickle
import os

import pandas as pd

directory = os.path.dirname(__file__)


class ChemicalMapper(object):
    """
    converts chemical species
    Datbase was creating by pulling down everything from HMDB
    """
    __slots__ = ('database', 'hmdb_accession_to_chemical_name',
                 'chemical_name_to_hmdb_accession', 'hmdb_accession_to_kegg',
                 'kegg_to_hmdb_accession', 'hmdb_accession_to_protein',
                 'synonyms_to_hmdb', '__dict__'
                 )

    def __init__(self):
        self.database = None
        self.hmdb_accession_to_chemical_name = {}
        self.chemical_name_to_hmdb_accession = {}
        self.hmdb_accession_to_kegg = {}
        self.kegg_to_hmdb_accession = {}
        self.hmdb_accession_to_protein = {}
        self.synonyms_to_hmdb = {}
        try:
            self.reload()
            print('Loading class data')
        except:
            print('Initializing Chemical mapping')
            self.load()

    def load(self):
        filename = os.path.join(directory, 'HMDB_processing',
                                'hmdb_dataframe.csv.gz')
        hmdb_database = pd.read_csv(filename)

        self.database = hmdb_database.where((pd.notnull(hmdb_database)), None)
        self.hmdb_accession_to_chemical_name = self.convert_to_dict(
                "accession", "name")
        self.chemical_name_to_hmdb_accession = self.convert_to_dict(
                "name", "accession")
        self.hmdb_accession_to_kegg = self.convert_to_dict("accession",
                                                           "kegg_id")
        self.kegg_to_hmdb_accession = self.convert_to_dict("kegg_id",
                                                           "accession")
        self.hmdb_accession_to_protein = self.convert_to_dict_from_list(
                "accession", "protein_associations")
        self.synonyms_to_hmdb = None
        self.save()

    def convert_to_dict(self, key, value):
        """ creates a dictionary from hmdb with a list of values for each key

        :param key:
        :param value:
        :return:
        """
        return {k: list(v) for k, v in self.database.groupby(key)[value]}

    def convert_to_dict_from_list(self, key, value):
        """ creates a dictionary from hmdb with a list of values for each key

        :param key:
        :param value:
        :return:
        """
        return {
        k: v.tolist()[0].strip("[]").replace(' ', '').replace('\'', '').split(
            ',') for k, v in self.database.groupby(key)[value]}

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

        print('Number of HMDB accessions = {0}'.format(
                len(self.database['accession'].unique())))
        print('Number of unique KEGG ids = {0}'.format(
                len(self.hmdb_accession_to_kegg.keys())))
        print('Number of HMDB to KEGG mappings = {0}'.format(
                len(self.kegg_to_hmdb_accession.values())))

    def save(self):
        """ save class instance

        :return:
        """
        print('Saving class data')
        with open(os.path.join(directory, 'hmdb_instance.txt'), 'w') as f:
            f.write(cPickle.dumps(self.__dict__))

    def reload(self):
        """ loads class instance

        :return:
        """

        with open(os.path.join(directory, 'hmdb_instance.txt'), 'r') as f:
            data = f.read()
            f.close()
        self.__dict__ = cPickle.loads(data)


if __name__ == '__main__':
    import timeit

    # cm = ChemicalMapper()
    x = timeit.timeit(ChemicalMapper().load)
    print(x)
    quit()
    cm = ChemicalMapper()

    count = 0
    # for k, v in cm.database.groupby('accession')['protein_associations']:
    #     v= v.tolist()[0]
    #     print(1, v,len(v))
    #     v = v.strip("[]").replace(' ', '').replace('\'', '').split(',')
    #     print(2, v, len(v))
    #     print(k,':',v, len(v))
    # if count == 3:
    #     quit()
    # count+=1

    # cm.print_info()
    # print(cm.hmdb_accession_to_kegg.keys())
    print('HMDB43286', cm.hmdb_accession_to_protein['HMDB43286'])
    quit()
    print(cm.hmdb_accession_to_kegg['HMDB06070'])
    print(cm.hmdb_accession_to_kegg['HMDB06525'])
    print(cm.hmdb_accession_to_kegg['HMDB38693'])
    print(cm.hmdb_accession_to_kegg['HMDB06523'])
    for i in cm.hmdb_accession_to_protein.keys():
        print(i, cm.hmdb_accession_to_protein[i])
