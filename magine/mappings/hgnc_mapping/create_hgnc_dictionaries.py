import pandas as pd
import os

dirname = os.path.dirname(__file__)

def create_from_hngc():
    # hngc = pd.read_table('hngc_protein-coding_gene.txt', delimiter='\t', low_memory=False)
    hngc = pd.read_table(
        'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt',
        delimiter='\t', low_memory=False)
    print(hngc.dtypes)
    hngc = hngc[hngc['status'] == 'Approved']
    hngc = hngc[['symbol',
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
    hngc.to_csv(outfile, compression='gzip', header=True)
    hngc = pd.read_csv(outfile)
    return hngc


if __name__ == '__main__':
    hngc = create_from_hngc()
