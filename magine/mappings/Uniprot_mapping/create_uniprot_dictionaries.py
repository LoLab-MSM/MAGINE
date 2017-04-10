import pandas as pd
import os

dir_name = os.path.dirname(__file__)

headers = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "RefSeq",
           "GI", "PDB", "GO", "UniRef100", "UniRef90",
           "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM", "UniGene",
           "PubMed", "EMBL", "EMBL-CDS", "Ensembl",
           "Ensembl_TRS", "Ensembl_PRO", "Additional PubMed"]

wanted_headers = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "PDB",
                  "GO", "Ensembl"]


def create_mouse_dataframe():
    mouse = pd.read_table('MOUSE_10090_idmapping_selected.tab.gz',
                          delimiter='\t', names=headers)
    print(mouse.head(10))
    mouse = mouse[wanted_headers]
    mouse.to_csv('../mouse_uniprot.gz', compression='gzip',
                 columns=wanted_headers, header=True)
    mouse = pd.read_csv('../mouse_uniprot.gz')
    return mouse


def create_human_dataframe():
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz'
    url2 = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'
    headers = ['uniprot', 'mapping_type', 'mapping']

    print("First time initialzing gene_mapper! \n"
          "Downloading from Uniprot. Might take awhile.")

    human = pd.read_table(url2, delimiter='\t',
                          names=headers,
                          compression='gzip', verbose=True)
    print(human.dtypes)
    print(human.head(10))
    print(human['mapping_type'].unique())
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
