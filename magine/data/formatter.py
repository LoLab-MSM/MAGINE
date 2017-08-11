import os
import re
import warnings

import numpy as np
import pandas as pd

# reload(sys)
# sys.setdefaultencoding('utf-8')

progensis_list_of_attributes = ['compound', 'compound_id', 'name', 'formula',

                                'score', 'top_score', 'retention_time_min',
                                'treated_control_fold_change',
                                'neutral_mass_da', 'm_z', 'charge',
                                'p_value_group_1_and_group_2',
                                'significant_flag',

                                'data_type',
                                'time', 'time_points',
                                'origin_endogenous',
                                'origin_food', 'origin_drug',
                                'description_food', 'description_drug',
                                'species_type']

cols = ['compound', 'compound_id', 'name', 'formula',
        'score', 'top_score', 'treated_control_fold_change',
        'p_value_group_1_and_group_2', 'significant_flag', 'data_type',
        'time', 'time_points', 'species_type'
        # 'origin_endogenous', 'origin_food',
        #  'origin_drug', 'description_food', 'description_drug',
        ]

rename_met_dict = {
    'max_fold_change_treated_vs_control': "treated_control_fold_change",
    "anova_p": "p_value_group_1_and_group_2",
    'phase': 'data_type',
    "experiment_type": "species_type"
}

valid_exp_data_types = ['label_free', 'ph_silac', 'silac', 'rna_seq',
                        'metabolomics']

out_gene_cols = ['gene', 'protein', 'time', 'treated_control_fold_change',
                 'time_points', 'p_value_group_1_and_group_2',
                 'significant_flag', 'species_type', 'data_type']


def load_csv(directory=None, filename=None):
    """ Creates a pandas.Dataframe from omics data files

    Parameters
    ----------
    directory : str, optional
        directory containing file
    filename : str
        name of file

    Returns
    -------
    pandas.Dataframe

    """
    if directory is None:
        file_location = filename
    else:
        file_location = os.path.join(directory, filename)
    print(file_location)
    assert os.path.exists(file_location), "Does not exist"

    if file_location.endswith('csv'):
        data = pd.read_csv(file_location, parse_dates=False, low_memory=False)
    elif file_location.endswith('csv.gz'):
        data = pd.read_csv(file_location, parse_dates=False, low_memory=False,
                           compression='gzip')
    elif file_location.endswith('xlsx'):
        data = pd.read_excel(file_location, parse_dates=False)
    else:
        print("Right now we can use xlsx and csv files one")
        quit()
    if 'time_points' in data:
        # time = convert_time_to_hours(data['time_points'])
        time = data.apply(convert_time_to_hours, axis=1)

        data.loc[:, 'time'] = time

    if 'gene' in data.dtypes:
        data['primary_genes'] = data['gene'].astype(str)
        data = check_data(data, 'primary_genes')
        tmp_sort = np.sort(data['primary_genes'].unique())
        print(list(tmp_sort[0:5]), list(tmp_sort[-5:]))
        return data

    if 'primary_genes' in data.dtypes:
        data['gene'] = data['primary_genes']
        data = check_data(data, 'primary_genes')
        tmp_sort = np.sort(data['primary_genes'].unique())
        print(list(tmp_sort[0:5]), list(tmp_sort[-5:]))
        return data

    else:
        return data


def merge_metabolite_row(row):
    """ compresses metabolite information columns

    For the metabolite data, we want to minimize the number of columns.
    We check to see if the "name" row exists, if it doesn't, we check to see
    if "description" row exists, if not we default to "compound_id".

    This compresses the redundant information stored in "name", "description",
    and "compound_id".

    Parameters
    ----------
    row

    Returns
    -------

    """
    if not isinstance(row['name'], str):
        if isinstance(row['description'], str):
            row['name'] = row['description']
        else:
            row['name'] = row['compound_id']
    return row['name']


def load_directory(dir_name, dtype=None):
    """ loads all files in a provided directory

    Parameters
    ----------
    dir_name: str
        path to directory containing all files
    dtype : str
        type of data
    Returns
    -------
    pandas.Dataframe

    """
    try:
        all_df = []
        for i in os.listdir(dir_name):
            all_df.append(load_csv(dir_name, i))
        df = pd.concat(all_df, ignore_index=True)
        if dtype == 'silac':
            return _process_silac(df)
        if dtype == 'rna_seq':
            return _process_rna_seq(df)
        if dtype == 'metabolites':
            return _process_metabolites(df)
        if dtype == 'label_free':
            return _process_label_free(df)
        if dtype == 'ph_silac':
            return _process_phsilac(df)
        return df
    except OSError:
        warnings.warn('Directory {} does not exist'.format(dir_name))
        raise


def convert_time_to_hours(d):
    """ replaces string of time point with float

    Parameters
    ----------
    d: pandas.Dataframe
        pandas.Dataframe

    Returns
    -------

    """
    t = d['time_points']
    if isinstance(t, float):
        return t
    if 's' in t:
        time = float(t.replace('s', '')) / 3600.
    elif 'min' in t:
        time = float(t.replace('min', '')) / 60.
    elif 'hr' in t:
        time = float(t.replace('hr', ''))
    elif 'h' in t:
        time = float(t.replace('h', ''))
    else:
        print('no time')
        return None
    return time


_names = {}
for i in range(16):
    _names['{0}-Sep'.format(i)] = 'SEPT{0}'.format(i)
    _names['{0}-Mar'.format(i)] = 'MARCH{0}'.format(i)
    _names['{0}-Dec'.format(i)] = 'DEC{0}'.format(i)


def check_data(data, keyword):
    genes = list(data[keyword])
    for n in _names:
        if n in genes:
            data.loc[data[keyword] == n, keyword] = _names[n]
    return data


def process_raptr_folder(directory):
    all_data = []
    for i in os.listdir(directory):
        if 'subcell' in i:
            continue
            df = load_csv(directory, i)
            df = _process_label_free(df, subcell=True)
            all_data.append(df)
        elif '_silac' in i:
            # continue
            df = load_csv(directory, i)
            df = _process_silac(df)
            all_data.append(df)
        elif 'ph-silac' in i:
            # continue
            df = load_csv(directory, i)
            df = _process_phsilac(df)
            all_data.append(df)
        elif 'protalizer_protein' in i:
            # continue
            df = load_csv(directory, i)
            df = _process_label_free(df, subcell=False)
            all_data.append(df)
        elif 'progenesis' in i:
            df = load_csv(directory, i)
            df = _process_metabolites(df)
            all_data.append(df)
        elif 'cuffdiff' in i:
            df = load_csv(directory, i)
            df = _process_rna_seq(df)
            all_data.append(df)
        else:
            print("Dont know {}".format(i))
    return pd.concat(all_data, ignore_index=True)


def process_list_of_dir_and_add_attribute(list_dim2):
    all_df = []
    for entry in list_dim2:
        print(entry)

        if isinstance(entry, dict):
            for i in ['filename', 'directory', 'sample_id', 'data_type']:
                assert i in entry, "Entry missing {}".format(i)
        else:
            assert len(entry) == 4
        filename = entry['filename']
        sample_id = entry['sample_id']
        data_type = entry['data_type']
        directory = entry['directory']

        d = load_csv(directory=directory, filename=filename)
        # assert data_type in valid_exp_data_types, \
        #     'data_type of {} is unknown'.format(data_type)

        if data_type.startswith('label_free'):
            d = _process_label_free(d)
            d['data_type'] = data_type
        elif data_type == 'ph_silac':
            d = _process_phsilac(d)
        elif data_type == 'silac':
            d = _process_silac(d)
        elif data_type == 'rna_seq':
            d = _process_rna_seq(d)
        elif data_type == 'metabolomics':
            d = _process_metabolites(d)

        d['sample_id'] = sample_id
        # d['time'] = d['sample_id'].astype(str) + '_' + d['time'].astype(str)
        all_df.append(d)
    return pd.concat(all_df)


def _process_rna_seq(data):
    """ processes rna_seq to a compact representation

    Parameters
    ----------
    data: pandas.Dataframe

    Returns
    -------
    pandas.Dataframe
    """
    data = data[data['status'] == 'OK']
    # data = data[~data.gene.str.contains(',')]
    data.loc[:, 'p_value_group_1_and_group_2'] = data['q_value']
    data['treated_control_fold_change'] = np.exp2(
        data['log2_fold_change'].astype(float))

    crit_1 = data['treated_control_fold_change'] < 1
    data.loc[crit_1, 'treated_control_fold_change'] = \
        -1. / data[crit_1]['treated_control_fold_change']

    data.loc[:, 'protein'] = data['gene'] + '_rnaseq'
    data.loc[:, 'data_type'] = 'rna_seq'
    data.loc[:, 'species_type'] = 'protein'

    return data[out_gene_cols]


def _process_metabolites(data):
    """ processes metabolite data to a compact representation

    Parameters
    ----------
    data: pandas.Dataframe

    Returns
    -------
    pandas.Dataframe
    """
    data = data.copy()
    # data = data.drop_duplicates(subset=['compound_id'])
    data['name'] = data.apply(merge_metabolite_row, axis=1)
    idx = data.groupby(['compound'])['score'].transform(max) == data['score']

    data.loc[:, 'top_score'] = False
    data.loc[idx, 'top_score'] = True
    data = data.rename(columns=rename_met_dict)
    data['treated_control_fold_change'] = \
        data['treated_control_fold_change'].apply(pd.to_numeric,
                                                  errors='coerce')

    criteria = np.isfinite(data['treated_control_fold_change'])
    data = data[criteria]
    crit_1 = data['treated_control_fold_change'] == np.inf
    crit_2 = data['treated_control_fold_change'] == -np.inf
    data.loc[crit_1, 'treated_control_fold_change'] = 1000.
    data.loc[crit_2, 'treated_control_fold_change'] = -1000.
    data['compound'] = data['compound'].astype(str)
    data['compound'] = data['compound'].str.encode('utf-8')

    data['compound_id'] = data['compound_id'].astype(str)
    # data['compound_id'] = data['compound_id'].str.decode('utf-8', 'replace')
    # data.loc[] = np.nan
    # checks to see if all data types exist
    for i in progensis_list_of_attributes:
        if i not in data.dtypes:
            warnings.warn("{} is not in dtypes. "
                          "Returning entire dataframe".format(i),
                          RuntimeWarning)
            return data

    return data[progensis_list_of_attributes]


def _process_label_free(data, subcell=False):
    label_free = data.copy()
    label_free['data_type'] = 'label_free'
    label_free['gene'] = label_free['primary_genes'].astype(str)
    label_free = label_free[label_free['gene'] != 'nan']

    protein_names = []

    def find_mod(row):
        s = str(row['protein'])
        if s.startswith('Acetylation'):
            change = 'ace'
            residue = s[s.find("(") + 1:s.find(")")]
            loc = re.findall('@\d+', s)[0].strip('@')
            name = "{0}_{1}({2}){3}_lf".format(row['gene'], residue, change,
                                               loc)

        elif s.startswith('Phosphorylation'):
            change = 'ph'
            residue = s[s.find("(") + 1:s.find(")")]
            loc = re.findall('@\d+', s)[0].strip('@')
            name = "{0}_{1}({2}){3}_lf".format(row['gene'], residue, change,
                                               loc, )
            protein_names.append(name)
        else:
            name = row['gene'] + '_lf'
        return name

    label_free['protein'] = label_free.apply(find_mod, axis=1)

    # label_free['significant_flag'] = False
    # label_free.loc[
    #     label_free['Final Significance'] == 1, 'significant_flag'] = True

    label_free['species_type'] = 'protein'
    if subcell:
        label_free['data_type'] = label_free['data_type'] + label_free[
            'sample_component']
        label_free['significant_flag'] = True
        label_free['p_value_group_1_and_group_2'] = 0.049
    else:
        label_free['p_value_group_1_and_group_2'] = \
            label_free['p_value_group_1_and_group_2']

    headers = ['gene', 'treated_control_fold_change', 'protein',
               'p_value_group_1_and_group_2', 'time', 'data_type',
               'time_points', 'significant_flag', 'species_type']

    label_free = label_free[headers]
    # label_free.to_csv('label_free.csv', index=False)
    return label_free


def _process_silac(data):
    data.loc[:, 'data_type'] = 'silac'
    data.loc[:, 'gene'] = data['primary_genes']
    data.loc[:, 'protein'] = data['primary_genes'] + '_silac'

    # data = data[data['n_significant'] == 2]

    data.loc[:, 'treated_control_fold_change'] = \
        data['mean_treated_untreated_fold_change']

    data.loc[:, 'species_type'] = 'protein'
    data.loc[:, 'significant_flag'] = False

    data.loc[:, 'p_value_group_1_and_group_2'] = 1.0

    criteria = (data['tier_level'] == 1) & \
               (data['fold_change_magnitude'] == 2)

    data.loc[criteria, 'p_value_group_1_and_group_2'] = 0.049
    data.loc[criteria, 'significant_flag'] = True

    return data[out_gene_cols]


def _process_phsilac(data):
    data.loc[:, 'phosphorylated_amino_acid'] = \
        data['phosphorylated_amino_acid'].astype(str)

    mod_sites = data[['modified_sequence', 'modified_seq_location',
                      'phosphorylated_amino_acid']]
    protein_names = []
    for i, row in mod_sites.iterrows():
        seq = row['modified_sequence'].strip('_')
        aa = row['phosphorylated_amino_acid']
        loc = row['modified_seq_location']

        if aa == 'nan':
            reg = '_'.join(loc.split(','))
            reg = reg.replace(' ', '')

            out_string = ''
            for word in ['(ca)', '(ox)']:
                if word in seq:
                    out_string += '_' + word
            protein_names.append("{}_{}_phsilac".format(out_string, reg))
            continue
        aa = aa.split(',')

        if loc == 'NA, NA':
            protein_names.append("_{}_phsilac".format(seq))
            continue
        start_loc = int(loc.split(',')[0])
        s = '_'
        for word in ['(ca)', '(ox)']:
            seq = seq.replace(word, '')
        mod = -2
        for n, m in enumerate(re.finditer('(ph)', seq)):
            s += '{0}(ph){1}_'.format(aa[n], m.start() + start_loc + mod)
            mod -= 4
        s += 'phsilac'
        protein_names.append(s)

    data.loc[:, 'protein'] = data['primary_genes'] + protein_names
    data.loc[:, 'gene'] = data['primary_genes']
    data.loc[:, 'treated_control_fold_change'] = data['both_fold_change_mean']
    data.loc[:, 'p_value_group_1_and_group_2'] = 1.0
    data.loc[:, 'significant_flag'] = False

    criteria = (data['overall_significance'].isin(['SIGNIFICANT',
                                                   'significant']))

    data.loc[criteria, 'p_value_group_1_and_group_2'] = 0.049
    data.loc[criteria, 'significant_flag'] = True

    data.loc[:, 'data_type'] = 'ph_silac'
    data.loc[:, 'species_type'] = 'protein'

    return data[out_gene_cols]


def pivot_table_for_export(data, save_name=None):
    """
    creates a pivot table of combined data that is in MAGINE format

    Parameters
    ----------
    data : pd.DataFrame
        output from magine analysis
    save_name : str
        save name for excel merged format

    Returns
    -------

    """
    index = ['GO_id', 'GO_name', 'depth', 'ref', 'aspect']
    tmp = pd.pivot_table(data, index=index, columns='sample_index',
                         aggfunc='first'
                         )
    tmp = tmp[['pvalue', 'enrichment_score', 'genes', 'n_genes', ]]
    if save_name:
        tmp.to_excel('{}.xlsx'.format(save_name),
                     merge_cells=True)
    return tmp


def pivot_raw_gene_data(data, save_name=None):
    headers1 = ['gene', 'treated_control_fold_change', 'protein',
                'p_value_group_1_and_group_2', 'time', 'data_type',
                'time_points', 'significant_flag']
    cols = ['compound', 'compound_id',
            'treated_control_fold_change',
            'p_value_group_1_and_group_2', 'significant_flag',
            'data_type', 'time', 'time_points',
            ]
    prot = data[data['species_type'] == 'protein'][headers1]

    tmp = pd.pivot_table(prot, index=['gene', 'protein', 'data_type'],
                         columns='time', aggfunc='first')

    tmp.to_csv('protein_data_pivot.csv', index=True)
    if save_name:
        tmp.to_excel('{}.xlsx'.format(save_name),
                     merge_cells=True)
    return tmp
    meta = data[data['species_type'] == 'metabolites'][cols]
    tmp = pd.pivot_table(meta, index=['compound', 'compound_id'],
                         columns='time', aggfunc='first')

    tmp.to_csv('metabolomics_data_pivot.csv', index=True)


def pivot_tables_for_export(data, save_name=None):
    headers1 = ['gene', 'treated_control_fold_change', 'protein',
                'p_value_group_1_and_group_2', 'time', 'data_type',
                'significant_flag',  # 'time_points',
                ]
    cols = ['compound', 'compound_id',
            'treated_control_fold_change',
            'p_value_group_1_and_group_2', 'significant_flag',
            'data_type', 'time',  # 'time_points',
            ]
    if 'compound_id' not in data.dtypes:
        cols.remove('compound_id')
    prot = data[data['species_type'] == 'protein'][headers1]
    meta = data[data['species_type'] == 'metabolites'][cols]
    genes = pd.pivot_table(prot, index=['gene', 'protein', 'data_type'],
                           columns='time', aggfunc='first', fill_value=np.nan)

    if save_name:
        genes.to_excel('{}_genes.xlsx'.format(save_name),
                       merge_cells=True)
        meta.to_csv('{}_genes.csv', index=True)
    if 'compound_id' not in data.dtypes:
        meta = pd.pivot_table(meta, index=['compound'],
                              columns='time', aggfunc='first')
    else:
        meta = pd.pivot_table(meta, index=['compound', 'compound_id'],
                              columns='time', aggfunc='first')
    if save_name:
        meta.to_excel('{}_metabolite.xlsx'.format(save_name),
                      merge_cells=True)
        meta.to_csv('{}_metabolites.csv', index=True)
    return genes, meta


def log2_normalize_df(df, fold_change):
    """
    
    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of fold changes
    fold_change : str
        column that contains fold change values

    Returns
    -------

    """
    tmp_df = df.copy()
    greater_than = tmp_df[fold_change] > 0
    less_than = tmp_df[fold_change] < 0
    tmp_df.loc[greater_than, 'log2fc'] = np.log2(
        tmp_df[greater_than][fold_change].astype(np.float64))
    tmp_df.loc[less_than, 'log2fc'] = -np.log2(
        -tmp_df[less_than][fold_change].astype(np.float64))
    return tmp_df
