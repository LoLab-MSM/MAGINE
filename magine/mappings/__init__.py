import os

dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'data')
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
