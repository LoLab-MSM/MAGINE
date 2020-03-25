import os
from shutil import rmtree

from appdirs import user_data_dir

dir_name = user_data_dir('magine')

dir_name = os.getenv('MAGINE_DATA', dir_name)

network_data_dir = os.path.join(dir_name, 'network_data')
id_mapping_dir = os.path.join(dir_name, 'id_data')


def create_storage_structure():
    # check or create the main data storage directory
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    # check or create the id mapping data storage directory
    if not os.path.exists(id_mapping_dir):
        os.makedirs(id_mapping_dir)

    # check or create the network data storage directory
    if not os.path.exists(network_data_dir):
        os.makedirs(network_data_dir)


def clear_cached_dbs():
    """Remove old database cached downloads"""
    # check or create the main data storage directory
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    # check or create the network data storage directory
    if os.path.exists(id_mapping_dir):
        rmtree(id_mapping_dir)

    # check or create the network data storage directory
    if os.path.exists(network_data_dir):
        rmtree(network_data_dir)


create_storage_structure()
