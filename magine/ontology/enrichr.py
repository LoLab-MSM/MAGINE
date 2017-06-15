import json
import requests
import pandas as pd
import os

_path = os.path.dirname(__file__)

_valid_libs = set()
with open(os.path.join(_path, '_valid_enricher_libs.txt'), 'r') as f:
    for n in f.read().split('\n'):
        _valid_libs.add(n)


class Enrichr(object):
    def __init__(self):
        self._url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        self._valid_libs = _valid_libs

    def print_valid_libs(self):
        for lib_name in self._valid_libs:
            print(lib_name)

    def run(self, list_of_genes, gene_set_lib='GO_Biological_Process_2017'):
        """

        Parameters
        ----------
        list_of_genes : list_like
            List of genes using HGNC gene names
        gene_set_lib : str
            Name of gene set library
            To print options use Enrichr.print_valid_libs

        Returns
        -------

        """
        assert isinstance(list_of_genes, list)
        assert gene_set_lib in _valid_libs, \
            "{} not in valid ids {}".format(gene_set_lib, _valid_libs)

        description = 'Example gene list'
        genes_str = '\n'.join(list_of_genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(self._url, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data = json.loads(response.text)
        user_list_id = data['userListId']

        enrichment_url = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        response = requests.get(
            enrichment_url + query_string % (user_list_id, gene_set_lib)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)
        #####
        # ENRICHR return a list of entries with each entry having these terms
        # Rank, Term name, P-value, Z-score, Combined score, Overlapping genes,
        # Adjusted p-value, Old p-value, Old adjusted p-value
        #####
        list_of_dict = []
        for i in data[gene_set_lib]:
            tmp_dict = dict()
            tmp_dict['rank'] = i[0]
            tmp_dict['term_name'] = i[1]
            tmp_dict['p_value'] = i[2]
            tmp_dict['z_score'] = i[3]
            tmp_dict['combined_score'] = i[4]
            tmp_dict['overlapping_genes'] = i[5]
            tmp_dict['adj_p_value'] = i[6]
            tmp_dict['old_p_value'] = i[7]
            tmp_dict['old_adj_p_value'] = i[8]
            list_of_dict.append(tmp_dict)
        cols = ['term_name', 'rank', 'p_value', 'z_score', 'combined_score',
                'adj_p_value', 'overlapping_genes', 'old_adj_p_value',
                'old_p_value']
        df = pd.DataFrame(list_of_dict, columns=cols)

        return df

if __name__ == '__main__':
    e = Enrichr()
    g_list = ['PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2',
              'P2RX7', 'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1',
              'ANAPC16', 'TMCC1', 'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2',
              'PLBD2', 'LARP7', 'TECPR2', 'ZNF302', 'CUX1', 'MOB2', 'CYTH2',
              'SEC22C', 'EIF4E3', 'ROBO2', 'ADAMTS9-AS2', 'CXXC1', 'LINC01314',
              'ATF7', 'ATP5F1']

    df = e.run(g_list, 'GO_Biological_Process_2017')
    print(df.head(10))
