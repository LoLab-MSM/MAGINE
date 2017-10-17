import os

import numpy as np
import pandas as pd

from magine.networks.correlation_networks.calculate_correlation_network import \
    correlation_sampling


def test_correlation():

    dir_name = os.path.join(os.path.dirname(__file__),
                            'Data', 'large_example.csv')
    data = pd.read_csv(dir_name)
    data = data[data['species_type'] == 'protein']
    tmp = pd.pivot_table(data, index=['protein'],
                         columns='time', aggfunc='first')
    names = np.array(tmp.index.tolist())
    tmp = tmp['treated_control_fold_change']
    tmp = tmp.fillna(np.nan)

    tmp = tmp.as_matrix()
    good_index = np.isfinite(tmp).sum(axis=1) > 3
    tmp = np.array(tmp)
    tmp = tmp[good_index]
    names = names[good_index]
    tmp = tmp[:100]
    names = names[:100]
    save_name = 'test'
    correlation_sampling(tmp, names, save_name, create_plots=True)
    correlation_sampling(tmp, names, save_name, create_plots=False,
                         in_parallel=True)

if __name__ == '__main__':
    test_correlation()
