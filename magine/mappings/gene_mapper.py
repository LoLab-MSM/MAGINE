try:
    import cPickle as pickle
except:  # python3 doesnt have cPickle
    import pickle
import os

import pandas as pd

from magine.mappings.databases import download_hgnc, download_ncbi, download_uniprot_db

directory = os.path.dirname(__file__)

pd.set_option('display.width', 20000)


def _to_dict(data, key, value):
    """
    creates a dictionary with a list of values for each key

    Parameters
    ----------
    data : pandas.DataFrame
    key : str
    value : str

    Returns
    -------

    """
    return_dict = {}
    for k, v in data.groupby(key)[value]:
        return_dict[k] = list(set(v))
        if None in return_dict[k]:
            return_dict[k].remove(None)
    return return_dict


class GeneMapper(object):
    """
    converts chemical species
    Datbase was creating by pulling down everything from HMDB
    """

    hgnc_categories = ['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name',
                       'alias_name', 'alias_symbol']

    def __init__(self, species='hsa'):
        self.species = species
        self._reload_fname = os.path.join(directory, 'data',
                                          'gene_mapping_instance.p')

        hgnc_name = os.path.join(directory, 'data', 'hgnc.gz')
        if not os.path.exists(hgnc_name):
            download_hgnc()
            assert os.path.exists(hgnc_name)

        # gather data from HGNC
        hgnc = pd.read_csv(hgnc_name)
        self.hgnc = hgnc

        ncbi_name = os.path.join(directory, 'data', 'ncbi.gz')
        if not os.path.exists(ncbi_name):
            download_ncbi()
            assert os.path.exists(ncbi_name)

        ncbi = pd.read_csv(ncbi_name)
        self.ncbi = ncbi

        # gather data from uniprot
        uniprot_path = os.path.join(directory, 'data', 'human_uniprot.csv.gz')
        # check to see if exists, if not create it
        if not os.path.exists(uniprot_path):
            download_uniprot_db()
            assert os.path.exists(uniprot_path)

        uniprot = pd.read_csv(uniprot_path, index_col=0)
        uniprot = pd.pivot_table(uniprot, columns='mapping_type',
                                 index='uniprot', aggfunc='first')
        uniprot.columns = uniprot.columns.droplevel()
        uniprot.reset_index(inplace=True)

        self.uniprot = uniprot

        # All are empty until set in load or reload
        self.gene_name_to_uniprot = {}
        self.gene_name_to_alias_name = {}
        self.gene_name_to_ensembl = {}
        self.gene_name_to_kegg = {}
        self.uniprot_to_gene_name = {}
        self.uniprot_to_kegg = {}
        self.protein_name_to_gene_name = None
        self.protein_name_to_uniprot = None
        self.kegg_to_gene_name = {}
        self.kegg_to_uniprot = {}
        self.ncbi_to_symbol = {}

        try:
            self.reload()
            print('Loading class data')
        except:
            print('Initializing Gene mapping')
            self.load()

    def load(self):

        self.gene_name_to_uniprot = _to_dict(self.hgnc, 'symbol',
                                             'uniprot_ids')
        self.gene_name_to_alias_name = _to_dict(self.hgnc, 'symbol',
                                                'alias_name')
        self.gene_name_to_ensembl = _to_dict(self.hgnc, 'symbol',
                                             'ensembl_gene_id')
        self.gene_name_to_kegg = _to_dict(self.uniprot, 'Gene_Name', 'KEGG')
        self.uniprot_to_kegg = _to_dict(self.uniprot, 'uniprot', 'KEGG')
        self.uniprot_to_gene_name = _to_dict(self.hgnc, 'uniprot_ids',
                                             'symbol')

        self.protein_name_to_gene_name = None
        self.protein_name_to_uniprot = None
        self.kegg_to_gene_name = _to_dict(self.uniprot, 'KEGG', 'Gene_Name')
        self.kegg_to_uniprot = _to_dict(self.uniprot, 'KEGG', 'uniprot')
        self.ncbi_to_symbol = _to_dict(self.ncbi, 'GeneID', 'Symbol')
        self.save()

    def save(self):
        """ save class instance
        """
        print('Saving class data')
        with open(self._reload_fname, 'wb') as f:
            f.write(pickle.dumps(self.__dict__))

    def reload(self):
        """ loads class instance
        """

        with open(self._reload_fname, 'rb') as f:
            data = f.read()
            f.close()
        self.__dict__ = pickle.loads(data)


if __name__ == '__main__':
    gm = GeneMapper()
    gm.load()
    print(gm.ncbi_to_symbol)
