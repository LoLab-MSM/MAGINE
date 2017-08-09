import os
from shutil import copytree, rmtree

from appdirs import user_data_dir

dir_name = user_data_dir('magine')

network_data_dir = os.path.join(dir_name, 'network_data')
id_mapping_dir = os.path.join(dir_name, 'id_data')

# check or create the main data storage directory
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
# check or create the network data storage directory
if os.path.exists(id_mapping_dir):
    rmtree(id_mapping_dir)

# check or create the network data storage directory
if os.path.exists(network_data_dir):
    rmtree(network_data_dir)

_sample_dir = os.path.join(os.path.dirname(__file__), '_sample_databases')
id_dir = os.path.join(_sample_dir, 'id_data')
network_dir = os.path.join(_sample_dir, 'network_data')

copytree(id_dir, id_mapping_dir)
copytree(network_dir, network_data_dir)
