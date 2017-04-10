import pandas as pd
import os

dirname = os.path.dirname(__file__)
url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt'


def download_hgnc():
    """
    Downloads HGNC and stores it as a pandas.DataFrame
    Returns
    -------

    """

    print("First time initializing gene_mapper! \n"
          "Downloading from HGNC. Might take awhile.")
    hgnc = pd.read_table(url, delimiter='\t', low_memory=False)
    hgnc = hgnc[hgnc['status'] == 'Approved']
    hgnc = hgnc[['symbol',
                 'uniprot_ids',
                 'ensembl_gene_id',
                 'name',
                 'location',
                 'entrez_id',
                 'ucsc_id',
                 'mirbase',
                 'vega_id',
                 'alias_name',
                 'alias_symbol']]
    outfile = os.path.join(dirname, '..', 'data', 'hgnc.gz')
    hgnc.to_csv(outfile, compression='gzip', header=True)
    hgnc = pd.read_csv(outfile)
    return hgnc


if __name__ == '__main__':
    hngc = download_hgnc()
