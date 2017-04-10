import networkx as nx
import pandas as pd
import os
import tempfile
import requests
import numpy as np

directory = os.path.dirname(__file__)


class ReactomeFunctionalInteraction(object):
    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = directory
        self.target_file = 'FIsInGene_031516_with_annotations.txt.zip'
        self.url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'
        self.out_path = os.path.join(self.tmp_dir, self.target_file)

    def setup(self):
        r = requests.get(self.url, stream=True)
        response = requests.head(self.url)
        file_size = int(response.headers['content-length'])
        print("Downloading: %s Bytes: %s" % (self.target_file, file_size))
        file_size_dl = 0
        block_sz = 8192
        v = set()
        milestone_markers = range(0, 101, 10)

        with open(self.out_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=block_sz):
                file_size_dl += len(chunk)
                percent_download = int(
                    np.floor(file_size_dl * 100. / file_size))

                if percent_download in milestone_markers:
                    if percent_download not in v:
                        print("{}%".format(percent_download))
                        v.add(percent_download)
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

        print("Downloaded {} and stored {}".format(self.url, self.out_path))

    def parse_network(self):
        table = pd.read_table(self.out_path)
        g = nx.DiGraph()
        added_genes = set()
        for i, row in table.iterrows():
            score = row['Score']
            gene1, gene2, interaction, prediction, mod = self.check_row(row)
            if gene1 == gene2:
                continue
            elif prediction:
                if gene1 not in added_genes:
                    g.add_node(gene1, sourceDB='ReactomeFI', speciesType='gene')
                    added_genes.add(gene1)
                if gene2 not in added_genes:
                    g.add_node(gene2, sourceDB='ReactomeFI', speciesType='gene')
                    added_genes.add(gene2)
                g.add_edge(gene1, gene2, interactionType=interaction, ptm=mod,
                           score=score, prediction='True', sourceDB='ReactomeFI')
            else:
                if gene1 not in added_genes:
                    g.add_node(gene1, sourceDB='ReactomeFI', speciesType='gene')
                    added_genes.add(gene1)
                if gene2 not in added_genes:
                    g.add_node(gene2, sourceDB='ReactomeFI', speciesType='gene')
                    added_genes.add(gene2)

                g.add_edge(gene1, gene2, interactionType=interaction, ptm=mod,
                           score=score, sourceDB='ReactomeFI')
        ge = set(g.nodes())
        g1 = set(table['Gene1'])
        g2 = set(table['Gene2'])
        g_all = set()
        g_all.update(g1)
        g_all.update(g2)
        for i in g_all:
            if i not in ge:
                print(i)

        print("Reactome network has {} nodes ".format(len(g.nodes())))
        print("Reactome network has {} edges ".format(len(g.edges())))
        nx.write_gml(g, os.path.join(directory, 'reactome_fi.gml'))

    @staticmethod
    def return_direction(text, gene1, gene2):
        if text in ("<-", "?-"):
            return gene2, gene1
        if text in ("->", "->"):
            return gene1, gene2
        else:
            return gene1, gene2


    def check_row(self, row):
        prediction = False
        interaction_type = None
        annotation = row['Annotation']
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
        g1, g2 = self.return_direction(row['Direction'], row['Gene1'], row['Gene2'])
        return g1, g2, interaction_type, prediction, modification

if __name__ == '__main__':
    rfi = ReactomeFunctionalInteraction()
    rfi.setup()
    rfi.parse_network()
