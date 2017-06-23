import io
import os
import sys
import networkx as nx
import pandas as pd
if sys.version_info[0] == 3:
    from urllib.request import urlopen
else:
    from urllib import urlopen

directory = os.path.dirname(__file__)


class BioGridDownload(object):
    def __init__(self):
        self.out_dir = directory
        self.url = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip'
        self.url2 = 'https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-CHEMICALS-LATEST.chemtab.zip'
        self._reverse = {"<-", "?-"}
        self._forward = {"->", "->"}
        self._db_name = 'BioGrid'

    def parse_network(self):
        """
        Parses tab delimited file to networkx.DiGraph


        Returns
        -------

        """

        chemical_int = pd.read_csv(io.BytesIO(urlopen(self.url2).read()),
                                   compression='zip',
                                   delimiter='\t',
                                   error_bad_lines=False,
                                   low_memory=False)
        print(chemical_int.dtypes)
        print(chemical_int.head(10))
        print(chemical_int.shape)
        print(chemical_int['Interaction Type'].unique())

        chemical_int = chemical_int[chemical_int['Organism'] == 'Homo sapiens']
        chem_cols = ['Official Symbol', 'Action',
                     'Chemical Name', 'ATC Codes', 'CAS Number', 'Pubmed ID',
                     'Chemical Type', 'Chemical Source', 'Chemical Source ID']

        chemical_int = chemical_int[chem_cols]
        print(chemical_int.head(10))
        print(chemical_int['Action'].unique())
        print(chemical_int['Chemical Type'].unique())
        chem_table = chemical_int.as_matrix()
        chem_g = nx.DiGraph()
        nodes_added = set()
        for row in chem_table:
            target = row[0]
            source = row[3]
            action = row[2]
            atc_code = row[4]
            cas_num = row[5]
            pubmed = row[6]
            chem_source = row[8]
            chem_source_id = row[9]

            chem_g.add_edge(target, source,
                            attr_dict={'interactionType': action,
                                       'pubmed': pubmed,
                                       })


        table = pd.read_csv(io.BytesIO(urlopen(self.url).read()),
                            compression='zip',
                            delimiter='\t',
                            error_bad_lines=False,
                            low_memory=False)
        print(table.dtypes)
        print(table.dtypes)
        print(table.head(10))
        print(table['Modification'].unique())
        print(table['Phenotypes'].unique())
        print(table['Experimental System Type'].unique())

        protein_cols = ['Official Symbol Interactor A',
                        'Official Symbol Interactor B',
                        'Pubmed ID', 'Score', 'Modification',
                        'Experimental System Type']
        table = table[protein_cols]
        print(table.head(10))
        table = table.as_matrix()
        g = nx.DiGraph()
        added_genes = set()

        for r in table:
            gene1, gene2, inter, prediction, mod, score = self._check_rows2(r)

            if gene1 == gene2:
                continue
            else:
                if gene1 not in added_genes:
                    g.add_node(gene1, sourceDB=self._db_name,
                               speciesType='gene')
                    added_genes.add(gene1)
                if gene2 not in added_genes:
                    g.add_node(gene2, sourceDB=self._db_name,
                               speciesType='gene')
                    added_genes.add(gene2)
                if prediction:
                    g.add_edge(gene1, gene2, interactionType=inter,
                               ptm=mod, score=score, prediction='True',
                               sourceDB=self._db_name)
                else:
                    g.add_edge(gene1, gene2, interactionType=inter,
                               ptm=mod, score=score, sourceDB=self._db_name)
        ge = set(g.nodes())
        g1 = set(table[:, 0])
        g2 = set(table[:, 1])
        g_all = set()
        g_all.update(g1)
        g_all.update(g2)

        for i in g_all:
            if i not in ge:
                print(i)

        print("Reactome network has {} nodes ".format(len(g.nodes())))
        print("Reactome network has {} edges ".format(len(g.edges())))
        nx.write_gml(g, os.path.join(directory, 'biogrid.gml'))

    def return_direction(self, text, gene1, gene2):
        if text in self._reverse:
            return gene2, gene1
        if text in self._forward:
            return gene1, gene2
        else:
            return gene1, gene2

    def _check_rows2(self, row):
        # print(row)

        interaction_type = None
        annotation = row[2]
        g1, g2 = self.return_direction(row[3], row[0], row[1])
        score = row[4]
        modification = ''
        if 'predicted' in annotation:
            prediction = True
        if 'repress' in annotation:
            interaction_type = 'repression'
        elif 'express' in annotation:
            interaction_type = 'expression'
        elif 'reaction' in annotation:
            interaction_type = 'reaction'
        elif 'complex' in annotation:
            interaction_type = 'complex'
        elif 'catalyze' in annotation:
            interaction_type = 'catalyze'
        elif 'deactiv' in annotation:
            interaction_type = 'deactivate'
        elif 'activ' in annotation:
            interaction_type = 'activate'
        elif 'inhib' in annotation:
            interaction_type = 'inhibit'
        elif 'binding' in annotation:
            interaction_type = 'binding'
        elif 'input' in annotation:
            interaction_type = '?'
        elif 'indirect effect' in annotation:
            interaction_type = 'indirect'
        elif 'compound' in annotation:
            interaction_type = 'compound'
        elif 'state change' in annotation:
            interaction_type = 'state change'
        elif 'interaction' in annotation:
            interaction_type = 'binding'
        elif 'PPrel' in annotation:
            interaction_type = 'binding'
        if 'dephosphoryl' in annotation:
            modification = 'dephosphorylation'
            if interaction_type is None:
                interaction_type = 'dephosphorylation'
        elif 'phosphoryl' in annotation:
            modification = 'phosphorylation'
            if interaction_type is None:
                interaction_type = 'phosphorylation'
        elif 'ubiquiti' in annotation:
            modification = 'ubiquitination'
            if interaction_type is None:
                interaction_type = 'ubiquitination'
        if interaction_type is None:
            interaction_type = '?'

        return g1, g2, interaction_type, prediction, modification, score

if __name__ == '__main__':
    import time

    st = time.time()
    rfi = BioGridDownload()
    rfi.parse_network()
    end_time = time.time()
    print(end_time - st)
    # 33.48 seconds with set
