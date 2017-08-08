import os
from shutil import copyfile

from appdirs import user_data_dir

dir_name = user_data_dir('magine')

id_mapping_dir = os.path.join(dir_name, 'id_data')

_sample_dir = os.path.join(os.path.dirname(__file__), '_sample_databases')

hmdb_file = os.path.join(_sample_dir, 'hmdb_dataframe.csv')
hmdb_out = os.path.join(id_mapping_dir, 'hmdb_dataframe.csv')

copyfile(hmdb_file, hmdb_out)
