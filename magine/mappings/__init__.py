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

