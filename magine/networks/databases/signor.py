import pandas as pd
from magine.networks.standards import edge_standards


def dowload():
    col_names = [
        'ENTITYA', 'TYPEA', 'IDA', 'DATABASEA', 'ENTITYB', 'TYPEB', 'IDB',
        'DATABASEB', 'EFFECT', 'MECHANISM', 'RESIDUE', 'SEQUENCE', 'TAX_ID',
        'CELL_DATA', 'TISSUE_DATA', 'MODULATOR_COMPLEX', 'TARGET_COMPLEX',
        'MODIFICATIONA', 'MODASEQ', 'MODIFICATIONB', 'MODBSEQ', 'PMID',
        'DIRECT', 'NOTES', 'ANNOTATOR', 'SENTENCE', 'SIGNOR_ID', 'nan']

    """
    table = pd.read_csv('https://signor.uniroma2.it/getData.php?organism=9606',
                        names=col_names, delimiter='\t', index_col=None,
                        error_bad_lines=False, encoding='utf-8'
                        )

    table.to_csv('signor.csv.gz', compression='gzip', index=False,
                 encoding='utf-8')
    """

    table = pd.read_csv('signor.csv.gz', compression='gzip', encoding='utf-8')
    print(table.dtypes)
    print(table.head(10))
    print(table.shape)

    # filter out non direct
    table = table[table['DIRECT'] == 't']

    # Filter out non descriptive
    table = table[~table['MECHANISM'].isnull()]

    # Drop SIGNOR edges, these are generally complexes
    table = table[~(table['DATABASEA'] == 'SIGNOR')]
    table = table[~(table['DATABASEB'] == 'SIGNOR')]

    # Not sure what they mean, so will remove. Ideally other DBs have this info
    table = table[~(table['MECHANISM'] == 'post transcriptional regulation')]

    print(table.shape)

    def map_to_activate_inhibit(row):
        effect = ''
        mechanism = row['MECHANISM']
        if 'down-regulates' in row['EFFECT']:
            effect = 'inhibit'
        elif 'up-regulates' in row['EFFECT']:
            effect = 'activate'
        if mechanism in edge_standards:
            mechanism = edge_standards[mechanism]
        elif mechanism == 'transcriptional regulation':
            if effect == 'inhibit':
                mechanism='repression'
            elif effect == 'activate':
                mechanism = 'expression'
        else:
            print(row)
            print(mechanism)
        if effect == '':
            return mechanism
        else:
            return "|".join([effect, mechanism])

    # relabel edge types
    table['edge_type'] = table.apply(map_to_activate_inhibit, axis=1)




if __name__ == '__main__':
    dowload()
