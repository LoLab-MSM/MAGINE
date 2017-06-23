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

import os
from .chemical_mapper import ChemicalMapper
from .gene_mapper import GeneMapper

__all__ = ['ChemicalMapper', 'GeneMapper']


data_dir = os.path.join(os.path.dirname(__file__), 'data')
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
