import pandas as pd

headers = ['hgnc_id', 'symbol', 'name', 'locus_group', 'locus_type', 'status', 'location', 'location_sortable',
           'alias_symbol', 'alias_name', 'prev_symbol', 'prev_name', 'gene_family', 'gene_family_id',
           'date_approved_reserved', 'date_symbol_changed', 'date_name_changed'    'date_modified', 'entrez_id',
           'ensembl_gene_id', 'vega_id', 'ucsc_id', 'ena', 'refseq_accession', 'ccds_id', 'uniprot_ids', 'pubmed_id',
           'mgd_id', 'rgd_id', 'lsdb', 'cosmic', 'omim_id', 'mirbase', 'homeodb', 'snornabase', 'bioparadigms_slc',
           'orphanet', 'pseudogene.org', 'horde_id', 'merops', 'imgt', 'iuphar', 'kznf_gene_catalog', 'mamit-trnadb',
           'cd', 'lncrnadb', 'enzyme_id', 'intermediate_filament_db']

wanted_headers = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "PDB", "GO"]


def create_from_hngc():
    hngc = pd.read_table('hngc_protein-coding_gene.txt', delimiter='\t', low_memory=False)
    print(hngc.dtypes)
    hngc = hngc[hngc['status'] == 'Approved']
    hngc = hngc[['symbol', 'uniprot_ids', 'ensembl_gene_id', 'name', 'alias_name', 'alias_symbol']]
    hngc.to_csv('hngc.gz', compression='gzip', header=True)
    hngc = pd.read_csv('hngc.gz')
    return hngc


hngc = create_from_hngc()
print(hngc['symbol'])


def convert_to_dict(key, value):
    """ creates a dictionary from hmdb with a list of values for each key

    :param key:
    :param value:
    :return:
    """
    return {k: list(v) for k, v in hngc.groupby(key)[value]}


un_to_genename = convert_to_dict('uniprot_ids', 'symbol')
for i in un_to_genename:
    print(i)
