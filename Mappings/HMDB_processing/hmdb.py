"""
HMDB database download and processing
"""
import czipfile as zipfile
import os

import numpy as np
import pandas as pd
import tempfile
import urllib2
import xml.etree.cElementTree as ElementTree

from Mappings.xml_to_dictionary import XmlDictConfig, XmlListConfig

directory = os.path.dirname(__file__)


class HMDB:
    """ Downloads and processes HMDB metabolite database

    """

    def __init__(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.out_dir = directory
        self.target_file = 'hmdb_metabolites.zip'

    def download_hmdb(self):
        """ Downloads hmdb metabolite xml file

        :return:
        """

        hmdb_db_url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
        u = urllib2.urlopen(hmdb_db_url)
        out_path = os.path.join(self.tmp_dir, self.target_file)
        f = open(out_path, 'wb')
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        print("Downloading: %s Bytes: %s" % (self.target_file, file_size))
        file_size_dl = 0
        block_sz = 8192
        v = set()
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            f.write(buffer)
            percent_download = np.floor(file_size_dl * 100. / file_size)
            if percent_download % 10 == 0:
                status = r"%10d  [%3.2f%%]" % (
                    file_size_dl, file_size_dl * 100. / file_size)
                status += chr(8) * (len(status) + 1)
                if percent_download not in v:
                    print(status)
                    v.add(percent_download)
        f.close()
        print("Downloaded {} and stored {}".format(hmdb_db_url, out_path))

    def unzip_hmdb(self, out_directory):
        """ unzips hmdb metabolite file

        :return:
        """
        hmdb_file = os.path.join(self.tmp_dir, self.target_file)
        if not os.path.exists(hmdb_file):
            print("Downloading metabolite information from HMDB")
            self.download_hmdb()
        print("Unzipping metabolite file")
        zip_ref = zipfile.ZipFile(hmdb_file, 'r')
        zip_ref.extractall(out_directory)
        zip_ref.close()

    def parse_hmdb(self):
        """ parse HMDB to Pandas.DataFrame

        """
        out_dir = os.path.join(self.tmp_dir, 'HMDB')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            self.unzip_hmdb(out_dir)

        count = 0
        categories = ['kegg_id', 'name', 'accession', 'chebi_id',
                      'chemspider_id',
                      'biocyc_id', 'synonyms', 'pubchem_compound_id',
                      'protein_associations',
                      # 'secondary_accessions',
                      # 'iupac_name',
                      # 'normal_concentrations','chemical_formula', 'smiles',
                      # 'drugbank_id', 'average_molecular_weight',
                      # 'pathways', 'metlin_id',
                      'inchikey',
                      ]
        print("Parsing metabolite information from files")
        tmp_all = []
        for i in os.listdir(out_dir):
            if i.startswith('H'):
                # if count > 100:
                #     continue
                count += 1
                filename = os.path.join(out_dir, i)
                tree = ElementTree.parse(filename)
                xmldict = XmlDictConfig(tree.getroot())
                tmp = []
                for cat in categories:
                    info = xmldict[cat]
                    if info is None:
                        tmp.append(info)
                    elif cat == 'protein_associations':
                        tmp_list = []
                        if type(info) == str:
                            continue
                        else:
                            if type(info['protein']) == XmlDictConfig:
                                tmp_list.append(info['protein']['gene_name'])

                            elif type(info['protein']) == XmlListConfig:
                                for ii in info['protein']:
                                    tmp_list.append(ii['gene_name'])
                        tmp.append(tmp_list)
                    elif cat == 'synonyms':
                        tmp_list = []
                        if type(info) == str:
                            continue
                        else:
                            if type(info['synonym']) == str:
                                tmp_list.append(info['synonym'])
                            else:
                                for each in info['synonym']:
                                    tmp_list.append(each)
                        tmp.append(tmp_list)
                    else:
                        tmp.append(info.encode('ascii', 'ignore'))
                tmp_all.append(tmp)
        df = pd.DataFrame(tmp_all, columns=categories)
        df.to_csv(os.path.join(self.out_dir, 'hmdb_dataframe.csv.gz'),
                  compression='gzip')
        print("Done processing HMDB")


if __name__ == '__main__':
    hm = HMDB()
    hm.parse_hmdb()
