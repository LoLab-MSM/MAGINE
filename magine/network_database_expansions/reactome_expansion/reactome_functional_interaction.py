import io
import os
import tempfile
import urllib

import networkx as nx
import numpy as np
import pandas as pd
import requests

directory = os.path.dirname(__file__)


class ReactomeFunctionalInteraction(object):
    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = directory
        self.target_file = 'FIsInGene_031516_with_annotations.txt.zip'
        self.url = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2015/FIsInGene_031516_with_annotations.txt.zip'
        self.out_path = os.path.join(self.tmp_dir, self.target_file)
        self._reverse = {"<-", "?-"}
        self._forward = {"->", "->"}

    def _setup(self):
        """
        Downloads file

        Returns
        -------

        """
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

    # @profile
    def parse_network(self):
        """
        Parses tab delimited file to networkx.DiGraph


        Returns
        -------

        """
        table = pd.read_csv(io.BytesIO(urllib.urlopen(self.url).read()),
                            compression='zip',
                            delimiter='\t',
                            error_bad_lines=False)
        table = table.as_matrix()
        g = nx.DiGraph()
        added_genes = set()

        for r in table:
            gene1, gene2, inter, prediction, mod, score = self._check_rows2(r)

            if gene1 == gene2:
                continue
            else:
                if gene1 not in added_genes:
                    g.add_node(gene1, sourceDB='ReactomeFI',
                               speciesType='gene')
                    added_genes.add(gene1)
                if gene2 not in added_genes:
                    g.add_node(gene2, sourceDB='ReactomeFI',
                               speciesType='gene')
                    added_genes.add(gene2)
                if prediction:
                    g.add_edge(gene1, gene2, interactionType=inter,
                               ptm=mod, score=score, prediction='True',
                               sourceDB='ReactomeFI')
                else:
                    g.add_edge(gene1, gene2, interactionType=inter,
                               ptm=mod, score=score, sourceDB='ReactomeFI')
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
        nx.write_gml(g, os.path.join(directory, 'reactome_fi.gml'))

    def return_direction(self, text, gene1, gene2):
        if text in self._reverse:
            return gene2, gene1
        if text in self._forward:
            return gene1, gene2
        else:
            return gene1, gene2

    def _check_rows2(self, row):
        # print(row)
        prediction = False
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
    rfi = ReactomeFunctionalInteraction()
    rfi.parse_network()
    end_time = time.time()
    print(end_time - st)
    # 33.48 seconds with set
