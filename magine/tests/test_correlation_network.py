"""
import os

import numpy as np
import pandas as pd
import pathos.multiprocessing as mp

from magine.networks.correlation_networks.calculate_correlation_network import \
    correlation_sampling


def test_correlation():

    dir_name = os.path.join(os.path.dirname(__file__),
                            'Data', 'jak_atra.csv.gz')

    data = pd.read_csv(dir_name)
    data = data[data['species_type'] == 'protein']
    tmp = pd.pivot_table(data, index=['protein'],
                         columns='time', aggfunc='first', fill_value=0.0)
    names = np.array(tmp.index.tolist())

    tmp = tmp['treated_control_fold_change']

    tmp = tmp.values
    good_index = (tmp == 0.0).sum(axis=1) < 3
    tmp = np.array(tmp)
    tmp = tmp[good_index]
    names = names[good_index]

    # tmp = tmp[:100]
    # names = names[:100]

    save_name = 'test'

    correlation_sampling(tmp, names, save_name, create_plots=True, pool=pool)


if __name__ == '__main__':
    pool = mp.Pool(4)
    test_correlation()
"""
