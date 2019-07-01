"""Interface to mapping between species IDs.

We constructed two classes for mappings ids.
1. GeneMapper
2. ChemicalMapper


Currently supported options
---------------------------

1. Genes
    KEGG, HGNC, Uniprot, Entrez, ensembl_gene_id
2. Metabolites/compounds
    kegg_id, name, accession, chebi_id, chemspider_id, biocyc_id, synonyms,
    pubchem_compound_id, protein_associations, inchikey, iupac_name, ontology,
    drugbank_id, chemical_formula, smiles, metlin_id, average_molecular_weight

"""
from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

__all__ = ['ChemicalMapper', 'GeneMapper']

import magine.mappings.databases.download_libraries as dl
import magine.networks.databases as nd


def download_id_mapping():
    dl.download_hgnc()
    dl.download_ncbi()
    dl.download_uniprot()


def download_network_dbs():
    nd.download_reactome_fi()
    nd.download_signor()
    nd.download_biogrid()
    dl.HMDB()


if __name__ == '__main__':
    import time

    st = time.time()
    download_id_mapping()
    download_network_dbs()
    et = time.time()
    print(et - st)
