import os
from magine.mappings.chemical_mapper import ChemicalMapper
from magine.mappings.gene_mapper import GeneMapper

__all__ = ['ChemicalMapper', 'GeneMapper']


data_dir = os.path.join(os.path.dirname(__file__), 'data')
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
