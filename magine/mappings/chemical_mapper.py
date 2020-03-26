from bioservices import UniChem
import pandas as pd
from sortedcontainers import SortedSet, SortedDict

from magine.mappings.databases.download_libraries import HMDB

try:
    import cPickle as pickle
except:  # python3 doesnt have cPickle
    import pickle
try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str
chem = UniChem()


class ChemicalMapper(object):
    """ Convert chemical species across various ids.

    Database was creating using HMDB

    """

    valid_columns = ['kegg_id', 'name', 'accession', 'chebi_id', 'inchikey',
                     'chemspider_id', 'biocyc_id', 'synonyms', 'iupac_name',
                     'pubchem_compound_id', 'protein_associations',
                     'ontology', 'drugbank_id', 'chemical_formula',
                     'smiles', 'metlin_id', 'average_molecular_weight',
                     'secondary_accessions'
                     ]

    def __init__(self, fresh_download=False):
        """

        Parameters
        ----------
        fresh_download: bool
            download new copy of database

        """

        self.database = None
        self._hmdb_to_chem_name = None
        self._chem_name_to_hmdb = None
        self._hmdb_to_kegg = None
        self._kegg_to_hmdb = None
        self.synonyms_to_hmdb = None
        self._drugbank_to_hmdb = None
        self._hmdb_to_protein = None
        self._hmdb_main_to_protein = None
        self._hmdb_accession_to_main = None
        hmdb_database = HMDB().load_db(fresh_download=fresh_download)
        self.database = hmdb_database.where((pd.notnull(hmdb_database)), None)
        self.database['main_accession'] = self.database['accession']
        sub_db = self.database[
            self.database['secondary_accessions'].str.contains('|', na=False)]
        new_df = tidy_split(sub_db, 'secondary_accessions', '|')
        new_df['accession'] = new_df['secondary_accessions']
        self.database = pd.concat([self.database, new_df])
        self.kegg_hmdb = chem.get_mapping("kegg_ligand", "hmdb")
        self.kegg_to_hmdb = self._to_dict("kegg_id", "main_accession")
        self.hmdb_to_chem_name = self._to_dict("main_accession", "name")

    @property
    def hmdb_to_kegg(self):
        if self._hmdb_to_kegg is None:
            self._hmdb_to_kegg = self._to_dict("accession", "kegg_id")
        return self._hmdb_to_kegg

    @property
    def chem_name_to_hmdb(self):
        if self._chem_name_to_hmdb is None:
            self._chem_name_to_hmdb = self._to_dict("name", "main_accession")
        return self._chem_name_to_hmdb

    @property
    def drugbank_to_hmdb(self):
        if self._drugbank_to_hmdb is None:
            self._drugbank_to_hmdb = self._to_dict("drugbank_id",
                                                   "main_accession")
        return self._drugbank_to_hmdb

    @property
    def hmdb_to_protein(self):
        if self._hmdb_to_protein is None:
            self._hmdb_to_protein = self._from_list_dict(
                "accession", "protein_associations"
            )
        return self._hmdb_to_protein

    @property
    def hmdb_main_to_protein(self):
        if self._hmdb_main_to_protein is None:
            self._hmdb_main_to_protein = self._from_list_dict(
                "main_accession", "protein_associations"
            )
        return self._hmdb_main_to_protein

    @property
    def hmdb_accession_to_main(self):
        if self._hmdb_accession_to_main is None:
            self._hmdb_accession_to_main = self._from_list_dict(
                "accession", "main_accession"
            )
        return self._hmdb_accession_to_main

    def _to_dict(self, key, value):
        """ creates a dictionary with a list of values for each key

        Parameters
        ----------
        key : str
        value : str

        Returns
        -------
        dict

        """
        d = self.database[[key, value]].copy()
        d.dropna(how='any', inplace=True)
        return_dict = SortedDict()
        for i, j in d.values:
            i = i.strip()
            if i not in return_dict:
                return_dict[i] = set()
            return_dict[i].add(j)
        return return_dict

    def _from_list_dict(self, key, value):
        d = self.database[[key, value]].copy()
        d.dropna(how='any', inplace=True)
        return_dict = SortedDict()
        for i, j in d.values:
            if i in return_dict:
                return_dict[i].update(j.split('|'))
            else:
                return_dict[i] = SortedSet(j.split('|'))
        return return_dict

    def check_synonym_dict(self, term, format_name):
        """ checks hmdb database for synonyms and returns formatted name

        Parameters
        ----------
        term : str
        format_name : str

        Returns
        -------
        dict


        Examples
        --------
        >>> cm = ChemicalMapper()
        >>> cm.check_synonym_dict(term='dodecene', format_name='main_accession')
        ['HMDB0000933', 'HMDB0059874']

        """
        synonyms = self.database[['synonyms', format_name]].copy()
        synonyms.dropna(how='any', inplace=True)
        synonyms['synonyms'] = synonyms['synonyms'].str.lower()

        hits = synonyms[synonyms['synonyms'].str.contains(term.lower())]
        matches = sorted(set(hits[format_name].values))
        return matches

    def print_info(self):
        """ print information about the dataframe

        Returns
        -------

        """
        print('Number of HMDB accessions = {0}'.format(
            len(self.database['accession'].unique())))
        print('Number of unique KEGG ids = {0}'.format(
            len(self.hmdb_to_kegg.keys())))
        print('Number of HMDB to KEGG mappings = {0}'.format(
            len(self.kegg_to_hmdb.values())))

    def convert_kegg_nodes(self, network):
        """
        Maps network from kegg to gene names

        Parameters
        ----------
        network : nx.DiGraph

        Returns
        -------
        dict

        """

        still_unknown = []
        hits = [i for i in set(network.nodes) if i.startswith('cpd:')]
        net_kegg_names = dict()
        net_chem_names = dict()
        net_cpd_to_hmdb = dict()
        for i in hits:

            name_stripped = i.lstrip('cpd:')
            net_kegg_names[i] = name_stripped

            if name_stripped in self.kegg_to_hmdb:
                mapping = self.kegg_to_hmdb[name_stripped]
                if isinstance(mapping, (list, set, SortedSet)):
                    names = '|'.join(set(mapping))
                    chem_names = set()
                    for name in mapping:
                        try:
                            chem_names.update(self.hmdb_to_chem_name[name])
                        except:
                            continue
                    net_cpd_to_hmdb[i] = names
                    net_chem_names[i] = order_merge(chem_names)

                elif isinstance(mapping, basestring):

                    chem_n = self.hmdb_to_chem_name[mapping]
                    net_cpd_to_hmdb[i] = mapping
                    net_chem_names[i] = '|'.join(chem_n.encode('ascii',
                                                               'ignore'))
                else:
                    print('Returned something else...', mapping)

            elif i in compound_manual:
                loc = compound_manual[i]
                net_cpd_to_hmdb[i] = loc
                if loc in self.hmdb_to_chem_name:
                    net_chem_names[i] = order_merge(
                        self.hmdb_to_chem_name[loc])
            else:
                still_unknown.append(i)
        if len(still_unknown):

            for i in still_unknown:
                name_stripped = i.lstrip('cpd:')
                if name_stripped in self.kegg_hmdb:
                    net_cpd_to_hmdb[i] = self.kegg_hmdb[name_stripped]
                # else:
                #     print("Cannot find a HMDB mapping for %s " % i)
        return net_cpd_to_hmdb, net_kegg_names, net_chem_names


# manually created based on missing in KEGG
compound_manual = {
    'cpd:C07909': 'HMDB0015015',
    'cpd:C16844': 'HMDB0001039',
    'cpd:C00076': 'HMDB0000464',
    'cpd:C00154': 'HMDB0001338',
    'cpd:C01561': 'HMDB0003550',
    'cpd:C04043': 'HMDB0003791',
    'cpd:C01165': 'HMDB0002104',
    'cpd:C00025': 'HMDB0000148',
    'cpd:C00696': 'HMDB0001403',
    'cpd:C00124': 'HMDB0000143',
}


def order_merge(species_set):
    return '|'.join(sorted(species_set))


def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


if __name__ == "__main__":
    cm = ChemicalMapper()
    print(cm.check_synonym_dict(term='dodecene', format_name='main_accession'))
    print(cm.hmdb_accession_to_main['HMDB15015'])
    print(cm.hmdb_accession_to_main['HMDB0015015'])
    # print(cm.hmdb_to_kegg['HMDB0015015'])
    print(cm.kegg_to_hmdb.keys())
    print(cm.kegg_to_hmdb['C07467'])
