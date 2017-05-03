import os

import pandas as pd

dir_name = os.path.dirname(__file__)


def create_mouse_dataframe():
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmppinga/by_organism/HUMAN_9606_idmapping.dat.gz',
    headers = ['uniprot', 'mapping_type', 'mapping']
    print("First time initialzing gene_mapper! \n"
          "Downloading from Uniprot. Might take awhile.")

    mouse = pd.read_table(url, delimiter='\t', names=headers,
                          compression='gzip', verbose=True)

    outfile = os.path.join(dir_name, '..', 'data', 'human_uniprot.csv.gz')
    mouse.to_csv(outfile, compression='gzip', columns=headers, header=True)

    return mouse


def create_human_dataframe():
    url2 = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'
    headers = ['uniprot', 'mapping_type', 'mapping']

    print("First time initialzing gene_mapper! \n"
          "Downloading from Uniprot. Might take awhile.")

    human = pd.read_table(url2, delimiter='\t', names=headers,
                          compression='gzip', verbose=True)

    outfile = os.path.join(dir_name, '..', 'data', 'human_uniprot.csv.gz')
    human.to_csv(outfile, compression='gzip', columns=headers, header=True)
    return human


def convert_to_dict(data, key, value):
    """ creates a dictionary with a list of values for each key

    :param key:
    :param value:
    :return:
    """
    return {k: list(v) for k, v in data.groupby(key)[value]}


def uni_to_kegg():
    outfile = os.path.join(dir_name, '..', 'data', 'human_uniprot.csv.gz')
    human = pd.read_csv(outfile, index_col=0)

    data = pd.pivot_table(human, columns='mapping_type', index='uniprot',
                          aggfunc='first')
    pd.set_option('display.width', 20000)
    print(data.head(10))

    data.columns = data.columns.droplevel()
    data.reset_index(inplace=True)

    kegg_only = convert_to_dict(data, 'KEGG', 'uniprot')

    print(kegg_only['hsa:7529'])


if __name__ == '__main__':
    human = create_human_dataframe()
    uni_to_kegg()
    # mouse = create_mouse_dataframe()
