import os
from shutil import copyfile

from appdirs import user_data_dir

dir_name = user_data_dir('magine')

# check or create the main data storage directory
if not os.path.exists(dir_name):
    os.makedirs(dir_name)


id_mapping_dir = os.path.join(dir_name, 'id_data')
# check or create the network data storage directory
if not os.path.exists(id_mapping_dir):
    os.makedirs(id_mapping_dir)

_sample_dir = os.path.join(os.path.dirname(__file__), '_sample_databases')

hmdb_file = os.path.join(_sample_dir, 'hmdb_dataframe.csv')
hmdb_out = os.path.join(id_mapping_dir, 'hmdb_dataframe.csv')

copyfile(hmdb_file, hmdb_out)
