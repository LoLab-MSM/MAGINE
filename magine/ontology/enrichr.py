import json
import os
import re

import numpy as np
import pandas as pd
import requests

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
        """
        Print a list of all available libraries EnrichR has to offer.
        Returns
        -------

        """
        for lib_name in sorted(self._valid_libs):
            print(lib_name)

    def run(self, list_of_genes, gene_set_lib='GO_Biological_Process_2017',
            verbose=False):
        """

        Parameters
        ----------
        list_of_genes : list_like
            List of genes using HGNC gene names
        gene_set_lib : str
            Name of gene set library
            To print options use Enrichr.print_valid_libs
        verbose : bool
            print information
        Returns
        -------

        """
        assert isinstance(list_of_genes, list)
        assert gene_set_lib in _valid_libs, \
            "{} not in valid ids {}".format(gene_set_lib, _valid_libs)

        if verbose:
            print("Running Enrichr with gene set {}".format(gene_set_lib))

        description = 'Example gene list'
        genes_str = '\n'.join(list_of_genes)
        payload = {
            'list':        (None, genes_str),
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
            tmp_dict['overlapping_genes'] = '|'.join(g for g in i[5])
            tmp_dict['adj_p_value'] = i[6]
            list_of_dict.append(tmp_dict)
        cols = ['term_name', 'rank', 'p_value', 'z_score', 'combined_score',
                'adj_p_value', 'overlapping_genes',
                ]
        df = pd.DataFrame(list_of_dict, columns=cols)

        if gene_set_lib.startswith('GO'):
            def get_go_id(row):
                s = row['term_name']
                go_id = re.search(r'\((GO.*?)\)', s).group(1)

                return go_id

            df['GO_id'] = df.apply(get_go_id, axis=1)

            term_names = df['term_name'].replace('(' + df['GO_id'] + ')', '')

            df['term_name'] = term_names
            cols.insert(0, 'GO_id')
        if verbose:
            print("Done calling Enrichr.")
        return df[cols]

    def run_samples(self, sample_lists, sample_ids,
                    gene_set_lib='GO_Biological_Process_2017', save_name=None):
        """

        Parameters
        ----------
        sample_lists : list_like
            List of lists of genes for enrichment analysis
        sample_ids : list
            list of ids for the provided sample list
        gene_set_lib : str
            Type of gene set, refer to Enrichr.print_valid_libs
        save_name : str, optional
            if provided it will save a file as a pivoted table with
            the term_ids vs sample_ids

        Returns
        -------

        """
        assert isinstance(sample_lists, list), "Please provide list of lists"
        assert isinstance(sample_lists[0], list), "Please provide list of lists"
        df_all = []
        for i, j in zip(sample_lists, sample_ids):
            df = self.run(i, gene_set_lib)
            df['sample_id'] = j
            df_all.append(df)
        df_all = pd.concat(df_all)
        index = ['term_name']

        if 'GO_id' in list(df.columns):
            index.insert(0, 'GO_id')

        p_df = pd.pivot_table(df_all, index=index,
                              columns='sample_id',
                              values=['term_name', 'rank', 'p_value', 'z_score',
                                      'combined_score',
                                      'adj_p_value', 'overlapping_genes',
                                      ],
                              aggfunc='first', fill_value=np.nan
                              )

        if save_name:
            p_df.to_excel('{}_enricher.xlsx'.format(save_name),
                          merge_cells=True)
        return p_df


if __name__ == '__main__':
    e = Enrichr()
    g_list = ['PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2',
              'P2RX7', 'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1',
              'ANAPC16', 'TMCC1', 'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2',
              'PLBD2', 'LARP7', 'TECPR2', 'ZNF302', 'CUX1', 'MOB2', 'CYTH2',
              'SEC22C', 'EIF4E3', 'ROBO2', 'ADAMTS9-AS2', 'CXXC1', 'LINC01314',
              'ATF7', 'ATP5F1']

    df = e.run(g_list, 'GO_Biological_Process_2017')
    lists = [['BAX', 'BCL2', 'CASP3'], ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'], save_name='test')
    # print(df.head(10))
