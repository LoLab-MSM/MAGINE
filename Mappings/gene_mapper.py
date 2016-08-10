import cPickle as pickle
import os

import pandas as pd

directory = os.path.dirname(__file__)

# mouse genes
mouse_kegg_to_uniprot = pickle.load(open(os.path.join(directory, 'mouse_kegg_mapper.p'), 'rb'))
mouse_uniprot_to_gene_name = pickle.load(open(os.path.join(directory, 'mouse_uniprot_to_gene_mapper.p'), 'rb'))
# human genes
human_uniprot_to_gene_name = pickle.load(open(os.path.join(directory, 'network_based_UPids_to_genes.p'), 'rb'))
human_uniprot_to_gene_name_2 = pickle.load(open(os.path.join(directory, 'kegg_to_uniprot_gene_name.p'), 'rb'))
human_kegg_to_uniprot = pickle.load(open(os.path.join(directory, 'human_kegg_mapper.p'), 'rb'))


class GeneMapper:
    """
    converts chemical species
    Datbase was creating by pulling down everything from HMDB
    """
    categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id', 'biocyc_id', 'synonyms',
                  'pubchem_compound_id', 'iupac_name', 'protein_associations']

    hgnc_categories = ['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name', 'alias_name', 'alias_symbol']

    def __init__(self, species='hsa'):
        self.species = species
        hgnc = pd.read_csv(os.path.join(directory, 'data', 'hngc.gz'))
        self.database = hgnc
        # hgnc
        self.gene_name_to_uniprot = self.convert_to_dict(hgnc, 'symbol', 'uniprot_ids')
        self.gene_name_to_alias_name = self.convert_to_dict(hgnc, 'symbol', 'alias_name')
        self.gene_name_to_ensembl = self.convert_to_dict(hgnc, 'symbol', 'ensembl_gene_id')
        self.gene_name_to_kegg = None

        self.uniprot_to_gene_name = self.convert_to_dict(hgnc, 'uniprot_ids', 'symbol')
        self.uniprot_to_kegg = None
        self.protein_name_to_gene_name = None
        self.protein_name_to_uniprot = None
        self.kegg_to_uniprot = None

    def convert_to_dict(self, data, key, value):
        """ creates a dictionary with a list of values for each key

        :param key:
        :param value:
        :return:
        """
        return {k: list(v) for k, v in data.groupby(key)[value]}

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


gm = GeneMapper()
