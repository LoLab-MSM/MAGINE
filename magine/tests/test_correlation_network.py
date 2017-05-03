import os

import numpy as np
import pandas as pd


def test_correlation():
    # data = pd.read_csv('../Data/proteins_plus_metabolites_across_time.csv.gz')
    # data = pd.read_csv('../Data/sig_only_proteins_across_time.csv.gz')
    from magine.networks.correlation_networks.calculate_correlation_network import \
        correlation_sampling
    dir_name = os.path.join(os.path.dirname(__file__),
                            'Data', 'large_example.csv')
    data = pd.read_csv(dir_name)

    data = data[data['species_type'] == 'protein']
    data = data[data['significant_flag']]

    tmp = pd.pivot_table(data, index=['protein'],
                         columns='time', aggfunc='first')
    names = np.array(tmp.index.tolist())
    tmp = tmp['treated_control_fold_change']
    tmp = tmp.fillna(0)

    tmp = tmp.as_matrix()

    print("Total data size = {}".format(np.shape(tmp)))
    tmp = np.array(tmp, np.float)
    print(len(tmp))

    total_combinations = len(tmp) * (len(tmp) - 1) / 2.
    print("Data with at least 3 finite numbers = {}".format(np.shape(tmp)))
    print(
    "total number of possible combinations = {}".format(total_combinations))

    save_name = 'test'
    output = correlation_sampling(tmp, names, 100, save_name,
                                  create_plots=True, sample_all=False)


if __name__ == '__main__':
    test_correlation()
