import numpy as np
import pandas as pd

hmdb = pd.read_pickle('../hmdb_dataframe.p')
hmdb.set_index('name').to_dict()
names = hmdb.keys()
print(names)
def extract_humancyc():
    data = np.loadtxt('HumanCyc_metabolic_pathways.sif', delimiter='\t', dtype=str)


    int_types = ['neighbor-of', 'catalysis-precedes', 'in-complex-with', 'consumption-controlled-by',
                 'controls-production-of', 'used-to-produce', 'reacts-with', 'controls-state-change-of',
                 'controls-transport-of-chemical', 'chemical-affects']


    for i in range(len(data)):
        if data[i, 1] == 'used-to-produce':
            if data[i,0] in names:
                print('hello',hmdb[data[i,0]])
            if data[i, 2] in names:
                print('hello2')
            #print(data[i, 0], data[i, 2])

extract_humancyc()