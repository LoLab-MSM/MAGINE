from magine.data.datatypes import ExperimentalData
from magine.magine_analysis import Analyzer
from magine.networks.network_generator import build_network
import os

__all__ = ['ExperimentalData', 'Analyzer', 'build_network']

dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'data')
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
