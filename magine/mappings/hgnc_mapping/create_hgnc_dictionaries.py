import os

import pandas as pd

dirname = os.path.dirname(__file__)


def download_hgnc():
    """
    Downloads HGNC and stores it as a pandas.DataFrame
    Returns
    -------

    """
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt'
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt'
    print("First time initializing gene_mapper! \n"
          "Downloading from HGNC. Might take awhile.")
    hgnc = pd.read_table(url, delimiter='\t', low_memory=False)
    print(hgnc.dtypes)
    print(hgnc['mirbase'].unique())
    print(hgnc.head(10))
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
    print(hgnc.head(10))
    outfile = os.path.join(dirname, '..', 'data', 'hgnc.gz')
    hgnc.to_csv(outfile, compression='gzip', header=True)
    hgnc = pd.read_csv(outfile)
    return hgnc


def download_ncbi():
    url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    ncbi = pd.read_table(url, delimiter='\t', low_memory=False,
                         compression='gzip')
    print(ncbi.dtypes)
    print(ncbi['Symbol'].unique())

    cols = ['GeneID', 'Symbol', 'description']
    ncbi = ncbi[cols]
    print(ncbi.head(10))
    outfile = os.path.join(dirname, '..', 'data', 'ncbi.gz')
    ncbi.to_csv(outfile, compression='gzip', header=True)


if __name__ == '__main__':
    # hngc = download_hgnc()
    download_ncbi()
