import pandas as pd
import numpy as np


def subset_data():
    df = pd.read_csv('reformatted_data.csv.gz', low_memory=False)

    print(df.dtypes)
    genes = pd.pivot_table(df, columns='time_points',
                           index=('gene', 'protein'))

    # genes = genes[genes['time'][1].isnull()]
    genes = genes[genes['time'][6].isnull()]
    genes = genes[genes['time'][24].isnull()]
    list_of_index = genes.index.values
    d = np.array(list_of_index)
    list_of_good = []
    letters = ['B', 'D', 'Z', 'X', 'Y', 'N', 'M']
    for i in d:
        gene_name = i[0]
        cont = False
        for i in letters:
            if i in gene_name:
                cont = True
                continue
        if not cont:
            list_of_good.append(gene_name)
    print(len(list_of_good))
    print(df.shape)
    df = df[df['gene'].isin(list_of_good)]
    print(df.shape)
    df.to_csv('large_example.csv', index=False)

    print(genes.head(10))

    compounds = pd.pivot_table(df, columns='time_points',
                               index=('compound', 'name'))

    print(compounds.shape)
    print(compounds.head(10))
