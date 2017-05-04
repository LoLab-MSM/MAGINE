import os

import numpy as np
import pandas as pd


def test_correlation():
    from magine.networks.correlation_networks.calculate_correlation_network import \
        correlation_sampling
    dir_name = os.path.join(os.path.dirname(__file__),
                            'Data', 'large_example.csv')
    data = pd.read_csv(dir_name)

    data = data[data['species_type'] == 'protein']
    # data = data[data['significant_flag']]

    tmp = pd.pivot_table(data, index=['protein'],
                         columns='time', aggfunc='first')
    names = np.array(tmp.index.tolist())
    tmp = tmp['treated_control_fold_change']
    tmp = tmp.fillna(np.nan)

    tmp = tmp.as_matrix()
    good_index = np.isfinite(tmp).sum(axis=1) > 2
    tmp = np.array(tmp)
    tmp = tmp[good_index]
    names = names[good_index]
    print(tmp)
    print(len(tmp))
    print("Total data size = {}".format(np.shape(tmp)))
    total_combinations = len(tmp) * (len(tmp) - 1) / 2.
    print("Data with at least 3 finite numbers = {}".format(np.shape(tmp)))
    print(
    "total number of possible combinations = {}".format(total_combinations))

    save_name = 'test'
    output = correlation_sampling(tmp, names, 1000, save_name,
                                  create_plots=True, sample_all=False)
    output = correlation_sampling(tmp, names, 100, save_name,
                                  create_plots=False, in_parallel=True)

if __name__ == '__main__':
    test_correlation()
