"""
HMDB database download and processing
"""

import os
import tempfile
import xml.etree.cElementTree as ElementTree
import zipfile

import numpy as np
import pandas as pd
import requests

directory = os.path.dirname(__file__)

categories = ['kegg_id', 'name', 'accession', 'chebi_id', 'chemspider_id',
              'biocyc_id', 'synonyms', 'pubchem_compound_id',
              'protein_associations', 'inchikey', 'iupac_name',
              'ontology',  # 'cellular_location', 'biofunction'
              # 'cellular_location', 'biofunction', 'molecular_framework'
              # 'secondary_accessions',

              # 'normal_concentrations','chemical_formula', 'smiles',
              # 'drugbank_id', 'average_molecular_weight',
              # 'pathways', 'metlin_id',
              ]


class HMDB(object):
    """ Downloads and processes HMDB metabolites database

    """

    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = directory
        self.target_file = 'hmdb_metabolites.zip'
        self.out_name = os.path.join(self.out_dir, '..', 'data',
                                     'hmdb_dataframe.csv.gz')

    def _download_hmdb(self):
        """ Downloads hmdb metabolites xml file

        :return:
        """

        hmdb_db_url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'

        out_path = os.path.join(self.tmp_dir, self.target_file)

        r = requests.get(hmdb_db_url, stream=True)
        response = requests.head(hmdb_db_url)
        file_size = int(response.headers['content-length'])
        print("Downloading: %s Bytes: %s" % (self.target_file, file_size))
        file_size_dl = 0
        block_sz = 1024
        # block_sz = 8024
        v = set()
        milestone_markers = range(0, 101, 10)

        with open(out_path, 'wb') as f:
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

        print("Downloaded {} and stored {}".format(hmdb_db_url, out_path))
        return

    def _unzip_hmdb(self, out_directory):
        """ unzips hmdb metabolites file

        :return:
        """
        hmdb_file = os.path.join(self.tmp_dir, self.target_file)

        if not os.path.exists(hmdb_file):
            print("Downloading metabolites information from HMDB")
            self._download_hmdb()
        print("Unzipping metabolites file")
        zip_ref = zipfile.ZipFile(hmdb_file, 'r')
        zip_ref.extractall(out_directory)
        zip_ref.close()
        print("Done unzipping metabolites file")

    def setup(self):
        """ parse HMDB to Pandas.DataFrame

        """
        # out_dir = b'c:\users\jamesp~1\\appdata\local\\temp\\tmpnvpdle\HMDB'
        out_dir = os.path.join(self.tmp_dir, 'HMDB')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            self._unzip_hmdb(out_dir)

        count = 0

        print("Parsing metabolites information from files")
        tmp_all = [None] * 1000000
        # tmp_all = []
        for i in os.listdir(out_dir):
            print(i)
            filename = os.path.join(out_dir, i)

            # get an iterable
            context = ElementTree.iterparse(filename, events=("start", "end"))

            # turn it into an iterator
            context = iter(context)

            # get the root element
            event, root = next(context)
            for event, elem in context:

                if '}' in elem.tag:
                    elem.tag = elem.tag.split('}', 1)[1]
                if elem.tag == 'metabolite' and event == 'end':
                    # """

                    # if count > 10000:
                    #     break
                    # """
                    template = _create_dict(elem)
                    tmp_all[count] = template

                    count += 1
                    elem.clear()
                    root.clear()

        df = pd.DataFrame(tmp_all[:count], columns=template.keys())
        for i in template.keys():
            df[i] = df[i].str.decode('utf8')
        df.to_csv(self.out_name, compression='gzip', index=False, mode='w',
                  encoding='utf-8')
        print("Done processing HMDB")
        # quit()

    def load_db(self):
        df = pd.read_csv(self.out_name)
        return df


def _create_dict(elem):
    template = {}
    for i in categories:
        n = elem.find(i)
        if i == 'protein_associations':
            output = []
            for pr in n.findall('protein'):
                for gn in pr.findall('gene_name'):
                    gene_name = gn.text
                    if gene_name is not None:
                        output.append(gene_name.encode('utf-8', 'ignore'))
            template[i] = output
        elif i == 'synonyms':
            output = []
            for pr in n.findall('synonym'):
                output.append(pr.text.encode('utf-8', 'ignore'))
            template[i] = output
        elif i == 'ontology':
            bf_all = []
            cl_all = []
            for pr in n.findall('biofunctions'):
                for bf in pr.findall('biofunction'):
                    bf_text = bf.text
                    if bf_text is not None:
                        bf_all.append(bf_text.encode('utf-8', 'ignore'))
            for pr in n.findall('cellular_locations'):
                for cl in pr.findall('cellular_location'):
                    cl_text = cl.text
                    if cl_text is not None:
                        cl_all.append(cl.text.encode('utf-8', 'ignore'))
            template['cellular_locations'] = cl_all
            template['biofunction'] = bf_all
        else:
            output = n.text
            if output is not None:
                output = output.encode('utf-8', 'ignore')
            template[i] = output
    return template

if __name__ == '__main__':
    import time

    # st = time.time()
    hm = HMDB()
    # hm.download_hmdb()
    end = time.time()
    # time1 = end - st
    # print("Download time took {}".format(time1))
    st = time.time()
    hm.setup()
    end = time.time()
    time2 = end - st
    print("Parsing time took {}".format(time2))
    # df = hm.load_db()
    # print(df.head(10))
    # 97.2819998264 lxml
    # 105.569000006 lxml root.clear()
    # 83.878000021 lxml, no root.clear(), no  elem.getprevious()
    # 57.4879999161 cElementTree
