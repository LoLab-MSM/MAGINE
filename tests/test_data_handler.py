import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import pathos.multiprocessing as mp
from magine.data_handler import ExperimentalData

data_directory = os.path.join(os.path.dirname(__file__), 'Data')
# data_directory = os.path.join(os.path.dirname(__file__), 'magine')
exp_data = ExperimentalData('example_apoptosis.csv', data_directory)


def test_table():
    exp_data.create_table_of_data(sig=True, save_name='sig_table')


def test_time_series_volcano():
    exp_data.time_series_volcano('label_free', 'test_label_free',
                                 bh_critera=True)


def test_plot_list():
    l = ['ADORA1', 'PARP1', 'BAX']
    x = list(exp_data.list_proteins)
    processes = []
    import time
    st = time.time()
    for i in x:
        # p = mp.Process(target=exp_data.plot_list_of_genes,
        #                args=(list(x), i, 'proteins', i, False, True))
        processes.append((list(x), 'out', 'proteins', i, False, True))
        # processes.append(p)
        # p.start()
        # exp_data.plot_list_of_genes(list_of_genes=list(x), save_name=i,
        #                             out_dir='proteins', title=i,
        #                             plot_all_x=False, log_scale=False)
    # [x.start() for x in processes]
    # """
    pool = mp.Pool(4)
    # for task in processes:
    #     pool.apply(exp_data.plot_list_of_genes, args=(task,))
    pool.map(exp_data.plot_list_of_genes, processes)
    pool.close()
    pool.join()
    # """
    print(time.time() - st)


def test_html_output():
    exp_data.plot_all_proteins()


# cProfile.run('test_plot_list()', sort=1)
if __name__ == '__main__':
    test_plot_list()
