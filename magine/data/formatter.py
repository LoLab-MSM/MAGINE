import pandas as pd
import numpy as np
import re
import os
import warnings

names = {}
for i in range(16):
    names['{0}-Sep'.format(i)] = 'SEPT{0}'.format(i)
    names['{0}-Mar'.format(i)] = 'MARCH{0}'.format(i)
    names['{0}-Dec'.format(i)] = 'DEC{0}'.format(i)

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

rename_met_dict = {
    'max_fold_change_treated_vs_control': "treated_control_fold_change",
    "anova_p":                            "p_value_group_1_and_group_2",
    'phase':                              'data_type',
    "experiment_type":                    "species_type"
}


def load_csv(directory, filename):
    """ Creates a pandas.Dataframe from omics data files

    Parameters
    ----------
    directory : str
        directory containing file
    filename : str
        name of file

    Returns
    -------
    pandas.Dataframe

    """
    file_location = os.path.join(directory, filename)
    print(file_location)
    if file_location.endswith('csv'):
        data = pd.read_csv(file_location, parse_dates=False, low_memory=False)
    elif file_location.endswith('xlsx'):
        data = pd.read_excel(file_location, parse_dates=False)
    else:
        print("Right now we can use xlsx and csv files one")
        quit()
    if 'time_points' in data:
        time = convert_time_to_hours(data['time_points'])
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

    elif 'compound' in data.dtypes:
        return process_metabolites(data)
    else:
        print("Data does not have the correct headers!")
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


def process_metabolites(data):
    """ processes metabolite data to a compact representation

    Parameters
    ----------
    data: pandas.Dataframe

    Returns
    -------
    pandas.Dataframe
    """
    data = data.copy()
    data = data.drop_duplicates(subset=['compound_id'])
    names = data.apply(merge_metabolite_row, axis=1)
    data['name'] = names
    idx = data.groupby(['compound'])['score'].transform(max) == data['score']

    data.loc[:, 'top_score'] = False
    data.loc[idx, 'top_score'] = True
    data = data.rename(columns=rename_met_dict)

    # checks to see if all data types exist
    for i in progensis_list_of_attributes:
        if i not in data.dtypes:
            warnings.warn("{} is not in dtypes. "
                          "Returning entire dataframe".format(i),
                          RuntimeWarning)
            return data

    return data[progensis_list_of_attributes]


def load_directory(dir_name):
    """ loads all files in a provided directory

    Parameters
    ----------
    dir_name: str
        path to directory containing all files

    Returns
    -------
    pandas.Dataframe

    """
    try:
        all_df = []
        for i in os.listdir(dir_name):
            all_df.append(load_csv(dir_name, i))
        df = pd.concat(all_df, ignore_index=True)
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
    t = d.unique()
    if len(t) != 1:
        print("Non-unique timepoints in this file!")
        print("Must handle separately")

    t = t[0]
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


def load_silac(directory):
    silac = load_directory(directory)
    silac['data_type'] = 'silac'
    silac['gene'] = silac['primary_genes']

    silac = silac.loc[silac['n_significant'] == 2]
    silac['treated_control_fold_change'] = silac[
        'mean_treated_untreated_fold_change']
    silac['protein'] = silac['primary_genes'] + '_silac'
    # silac = silac.dropna(subset=['significant'])
    silac['species_type'] = 'protein'
    silac['significant_flag'] = False
    silac['p_value_group_1_and_group_2'] = 1.0
    critera = (silac['tier_level'] == 1) & (
    silac['fold_change_magnitude'] == 2)
    silac.loc[critera, 'p_value_group_1_and_group_2'] = 0.049
    silac.loc[critera, 'significant_flag'] = True
    silac = silac[['gene', 'protein', 'time', 'treated_control_fold_change',
                   'time_points', 'p_value_group_1_and_group_2',
                   'significant_flag', 'species_type', 'data_type']]
    return silac


def load_rna_data(directory):
    data = load_directory(directory)
    data = data[data['status'] == 'OK']
    print(data.dtypes)
    print(np.shape(data))
    data = data[data.gene.str.contains(',') == False]
    print(np.shape(data))
    data['p_value_group_1_and_group_2'] = data['q_value']

    data['treated_control_fold_change'] = np.exp2(
        data['log2_fold_change'].astype(float))
    data.loc[data['treated_control_fold_change'] < 1,
             'treated_control_fold_change'] = -1. / data[
        data['treated_control_fold_change'] < 1]['treated_control_fold_change']
    data['protein'] = data['gene'] + '_rnaseq'
    data['data_type'] = 'rna_seq'
    data['species_type'] = 'protein'
    data = data[['gene', 'protein', 'time', 'treated_control_fold_change',
                 'time_points', 'p_value_group_1_and_group_2',
                 'significant_flag', 'species_type', 'data_type']]

    return data


def check_data(data, keyword):
    genes = list(data[keyword])
    for n in names:
        if n in genes:
            data.loc[data[keyword] == n, keyword] = names[n]
    return data


def load_metabolite_dir(dirname):
    cols = ['compound', 'compound_id', 'name', 'formula',
            'score', 'top_score', 'treated_control_fold_change',
            'p_value_group_1_and_group_2', 'significant_flag', 'data_type',
            'time', 'time_points', 'species_type'
            # 'origin_endogenous', 'origin_food',
            # 'origin_drug', 'description_food', 'description_drug',
            ]

    return load_directory(dirname)[cols]



def load_label_free(directory):
    data = load_directory(dir_name=directory)
    label_free = process_label_free(data)
    return label_free


def process_label_free(data):
    label_free = data.copy()
    label_free['data_type'] = 'label_free'
    label_free['gene'] = label_free['primary_genes'].astype(str)
    label_free = label_free[label_free['gene'] != 'nan']
    mod_sites = label_free[['gene', 'protein']]
    protein_names = []
    for i, row in mod_sites.iterrows():
        s = str(row['protein'])
        if s.startswith('Acetylation'):
            change = 'ace'
            residue = s[s.find("(") + 1:s.find(")")]
            loc = re.findall('@\d+', s)[0].strip('@')
            name = "{0}_{1}({2}){3}_lf".format(row['gene'], residue, change,
                                               loc, )
            protein_names.append(name)
        elif s.startswith('Phosphorylation'):
            change = 'ph'
            residue = s[s.find("(") + 1:s.find(")")]
            loc = re.findall('@\d+', s)[0].strip('@')
            name = "{0}_{1}({2}){3}_lf".format(row['gene'], residue, change,
                                               loc, )
            protein_names.append(name)
        else:
            protein_names.append(row['gene'] + '_lf')

    print(label_free.dtypes)
    label_free['p_value_group_1_and_group_2'] = label_free[
        'p_value_group_1_and_group_2']
    label_free['protein'] = protein_names
    # label_free['significant_flag'] = False
    # label_free.loc[
    #     label_free['Final Significance'] == 1, 'significant_flag'] = True

    label_free['species_type'] = 'protein'

    headers = ['gene', 'treated_control_fold_change', 'protein',
               'p_value_group_1_and_group_2', 'time', 'data_type',
               'time_points', 'significant_flag', 'species_type']

    label_free = label_free[headers]
    label_free.to_csv('label_free.csv', index=False)
    return label_free


def load_phsilac(directory):
    phsilac = load_directory(dir_name=directory)

    return _process_phsilac(phsilac)


def _process_phsilac(data):
    data['phosphorylated_amino_acid'] = data[
        'phosphorylated_amino_acid'].astype(str)
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

    data['gene'] = data['primary_genes']
    data['protein'] = data['primary_genes'] + protein_names
    data['treated_control_fold_change'] = data['both_fold_change_mean']

    data['p_value_group_1_and_group_2'] = 1.0
    data['significant_flag'] = False

    criteria = (data['overall_significance'].isin(['SIGNIFICANT',
                                                   'significant']))

    data.loc[criteria, 'p_value_group_1_and_group_2'] = 0.049
    data.loc[criteria, 'significant_flag'] = True

    data['data_type'] = 'ph_silac'
    data['species_type'] = 'protein'
    headers = ['gene', 'treated_control_fold_change', 'protein',
               'p_value_group_1_and_group_2', 'time', 'data_type',
               'time_points', 'significant_flag', 'species_type']
    data = data[headers]
    return data


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
    index = ['GO_id', 'GO_name', 'depth', 'ref', 'slim', 'aspect']
    tmp = pd.pivot_table(data, index=index, columns='sample_index',
                         aggfunc='first'
                         )
    if save_name:
        tmp.to_excel('{}.xlsx'.format(save_name),
                     merge_cells=True)
    return tmp


def pivot_raw_gene_data(data):
    headers1 = ['gene', 'treated_control_fold_change', 'protein',
                'p_value_group_1_and_group_2', 'time', 'data_type',
                'time_points', 'significant_flag']
    cols = ['compound', 'compound_id',
            'treated_control_fold_change',
            'p_value_group_1_and_group_2', 'significant_flag',
            'data_type', 'time', 'time_points',
            ]
    prot = data[data['species_type'] == 'protein'][headers1]
    meta = data[data['species_type'] == 'metabolites'][cols]
    tmp = pd.pivot_table(prot, index=['gene', 'protein', 'data_type'],
                         columns='time', aggfunc='first')

    tmp.to_csv('protein_data_pivot.csv', index=True)

    tmp = pd.pivot_table(meta, index=['compound', 'compound_id'],
                         columns='time', aggfunc='first')

    tmp.to_csv('metabolomics_data_pivot.csv', index=True)
