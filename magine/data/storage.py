from appdirs import user_data_dir
import os

dir_name = user_data_dir('magine')

dir_name = os.getenv('MAGINE_DATA', dir_name)

network_data_dir = os.path.join(dir_name, 'network_data')
id_mapping_dir = os.path.join(dir_name, 'id_data')

# check or create the main data storage directory
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

# check or create the id mapping data storage directory
if not os.path.exists(id_mapping_dir):
    os.makedirs(id_mapping_dir)

# check or create the network data storage directory
if not os.path.exists(network_data_dir):
    os.makedirs(network_data_dir)



