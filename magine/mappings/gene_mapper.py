try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str

from bioservices import HGNC, UniProt
import pandas as pd
from sortedcontainers import SortedSet, SortedDict

from magine.mappings.databases import load_hgnc, load_uniprot, load_ncbi

pd.set_option('display.width', 20000)


class GeneMapper(object):
    """
    Mapping class between common gene ids

    Database was creating by pulling down from NCBI, UNIPROT, HGNC

    """

    def __init__(self, species='hsa'):
        self.species = species
        self.hgnc = load_hgnc()
        self.ncbi = load_ncbi()
        self.uniprot = load_uniprot()
        self._gene_name_to_uniprot = None
        self._gene_name_to_alias_name = None
        self._gene_name_to_ensembl = None
        self._gene_name_to_kegg = None
        self._uniprot_to_gene_name = None
        self._uniprot_to_kegg = None
        self._kegg_to_gene_name = None
        self._kegg_to_uniprot = None
        self._ncbi_to_symbol = None
        # hackish way to initilize
        if 'x' in self.kegg_to_gene_name:
            pass
        if 'x' in self.ncbi_to_symbol:
            pass
        if 'x' in self.gene_name_to_uniprot:
            pass

    @property
    def gene_name_to_uniprot(self):
        if self._gene_name_to_uniprot is None:
            self._gene_name_to_uniprot = _dict(self.hgnc, 'symbol',
                                               'uniprot_ids')
        return self._gene_name_to_uniprot

    @property
    def gene_name_to_alias_name(self):
        if self._gene_name_to_alias_name is None:
            self._gene_name_to_alias_name = _dict(self.hgnc, 'symbol',
                                                  'alias_name')
        return self._gene_name_to_alias_name

    @property
    def gene_name_to_ensembl(self):
        if self._gene_name_to_ensembl is None:
            self._gene_name_to_ensembl = _dict(self.hgnc, 'symbol',
                                               'ensembl_gene_id')
        return self._gene_name_to_ensembl

    @property
    def uniprot_to_gene_name(self):
        if self._uniprot_to_gene_name is None:
            self._uniprot_to_gene_name = _dict(self.hgnc, 'uniprot_ids',
                                               'symbol')
        return self._uniprot_to_gene_name

    @property
    def gene_name_to_kegg(self):
        if self._gene_name_to_kegg is None:
            self._gene_name_to_kegg = _dict(self.uniprot, 'Gene_Name', 'KEGG')
        return self._gene_name_to_kegg

    @property
    def uniprot_to_kegg(self):
        if self._uniprot_to_kegg is None:
            self._uniprot_to_kegg = _dict(self.uniprot, 'uniprot', 'KEGG')
        return self._uniprot_to_kegg

    @property
    def kegg_to_gene_name(self):
        if self._kegg_to_gene_name is None:
            self._kegg_to_gene_name = _dict(self.uniprot, 'KEGG', 'Gene_Name')
            self._kegg_to_gene_name.update(manual_dict)
        return self._kegg_to_gene_name

    @property
    def kegg_to_uniprot(self):
        if self._kegg_to_uniprot is None:
            self._kegg_to_uniprot = _dict(self.uniprot, 'KEGG', 'uniprot')
        return self._kegg_to_uniprot

    @property
    def ncbi_to_symbol(self):
        if self._ncbi_to_symbol is None:
            self._ncbi_to_symbol = _dict(self.ncbi, 'GeneID', 'Symbol')
        return self._ncbi_to_symbol

    def check_synonym_dict(self, term, format_name):
        """ checks hmdb database for synonyms and returns formatted name

        Parameters
        ----------
        term : str
        format_name : str

        Returns
        -------
        dict

        """
        index = 'alias_symbol'
        synonyms = self.hgnc.copy()
        synonyms = synonyms.loc[~synonyms[index].isna()]
        synonyms[index] = synonyms[index].str.upper()

        hits = synonyms.loc[synonyms[index].str.contains(term.upper())].copy()
        hits[index] = hits[index].str.split('|')

        for i, row in hits.iterrows():
            if term in row[index]:
                return [row[format_name]]
        matches = sorted(set(hits[format_name].values))
        return matches

    def convert_kegg_nodes(self, network, species='hsa'):
        """ Convert kegg ids to HGNC gene symbol.

        Parameters
        ----------
        network : nx.DiGraph
        species : str {'hsa'}
            Main support for humans only.

        Returns
        -------
        kegg_to_gene_name, kegg_short : dict, dict
        """

        # Create the dictionary to store all conversions to be returned
        missing = set()
        unknown_genes = set()
        hits = {i for i in set(network.nodes) if i.startswith(species)}
        prefix = species + ':'
        kegg_short = {i: i.lstrip(prefix) for i in hits}

        kegg_to_gene_name = {i: self.kegg_to_gene_name[i]
                             for i in hits if i in self.kegg_to_gene_name}
        for gene, gn in kegg_to_gene_name.items():
            if isinstance(gn, basestring):
                kegg_to_gene_name[gene] = gn
            elif len(gn) == 1:
                kegg_to_gene_name[gene] = gn[0]
            else:
                found = False
                for g in gn:
                    if g in self.gene_name_to_uniprot:
                        kegg_to_gene_name[gene] = g
                        found = True
                if not found:
                    print("Found species that cannot be mapped to HGNC.")
                    print(gene, gn)
                    # kegg_to_gene_name[gene] = gn[0]
                    missing.add(gene)
        missing.update(hits.difference(kegg_to_gene_name))

        # check stores dictionaries
        for gene in missing:
            name_stripped = gene.lstrip(prefix)
            if int(name_stripped) in self.ncbi_to_symbol:
                new = self.ncbi_to_symbol[int(name_stripped)][0]
                if not isinstance(new, float):
                    kegg_to_gene_name[gene] = new
                else:
                    unknown_genes.add(gene)
            else:
                unknown_genes.add(gene)

        if len(unknown_genes) > 0:
            add_dict, missed = self.kegg_to_symbol_through_uniprot(
                unknown_genes)
            kegg_to_gene_name.update(add_dict)

            add_dict, final_missed = self.kegg_to_hugo(missed, species)
            kegg_to_gene_name.update(add_dict)

            print("{} mappings not found from kegg to"
                  " gene name".format(len(final_missed)))
            print(final_missed)

        return kegg_to_gene_name, kegg_short

    def kegg_to_symbol_through_uniprot(self, unknown_genes):
        # create string to call uniprot for mapping
        search_string = '\t'.join(unknown_genes)
        kegg_to_gene_name = dict()
        missing = set()
        uniprot = UniProt(verbose=True)
        # This is where it gets tricky. Checking to see if there is a uniprot
        # mapping for the species, if not, trying from KEGG side. Sometimes
        # kegg  links to a different uniprot, or uniprot links to a diff kegg.
        uni_dict = dict(uniprot.mapping("KEGG_ID", "ACC", query=search_string))
        for i in unknown_genes:
            if i in uni_dict:
                for n in uni_dict[i]:
                    x = uniprot.search("accession:{}".format(n),
                                       columns='genes(PREFERRED),reviewed,id',
                                       limit=1)
                    _, data = x.rstrip('\n').split('\n')
                    name, review, entry = data.split('\t')
                    if n != entry:
                        print(i, n, entry, x, "dont match")
                    elif review == 'reviewed':
                        kegg_to_gene_name[i] = name

            else:
                missing.add(i)
        print("{} mappings not found from kegg to"
              " gene name".format(len(missing)))
        print(missing)
        return kegg_to_gene_name, missing

    def kegg_to_hugo(self, genes, species='hsa'):
        """
        Converts all KEGG names to HGNC

        Parameters
        ----------
        genes : list
        species : str

        Returns
        -------
        dict
        """
        prefix = species + ':'
        hugo = HGNC(verbose=True)
        hugo_dict = {}
        not_found = set()
        for i in genes:
            tmp_name = i.lstrip(prefix)
            mapping = hugo.search(tmp_name)
            if 'response' in mapping:
                response = mapping['response']
                if 'numFound' in response:
                    if response['numFound'] == 0:
                        not_found.add(i)
                        continue
                    elif response['numFound'] == 1:
                        docs = response['docs'][0]
                        hugo_dict[i] = docs['symbol']
                        continue
                    else:
                        if 'symbol' in response['docs'][0]:
                            hugo_dict[i] = response['docs'][0]['symbol']
            else:
                not_found.add(i)
        if not_found != 0:
            print("{} not found after HGNC mapping".format(len(not_found)))
            print("{} ".format(not_found))
        return hugo_dict, not_found


def _dict(data, key, value):
    """
    creates a dictionary with a list of values for each key

    Parameters
    ----------
    data : pandas.DataFrame
    key : str
    value : str

    Returns
    -------

    """
    return_dict = SortedDict()

    d = data[[key, value]].copy()
    d.dropna(how='any', inplace=True)

    for i, j in d.values:
        if i in return_dict:
            return_dict[i].add(j)
        else:
            return_dict[i] = SortedSet([j])
    return return_dict


manual_dict = {'hsa:857': 'CAV1',
               'hsa:2250': 'FGF5',
               'hsa:5337': 'PLD1',
               'hsa:4312': 'MMP1',
               'hsa:102723407': 'IGHV4OR15-8',
               'hsa:100132074': 'FOXO6',
               'hsa:728635': 'DHRS4L1',
               'hsa:10411': 'RAPGEF3',
               'hsa:100101267': 'POM121C',
               'hsa:2768': 'GNA12',
               'hsa:2044': 'EPHA5',
               'hsa:100533467': 'BIVM-ERCC5',
               'hsa:7403': 'KDM6A',
               'hsa:1981': 'EIF4G1',
               'hsa:2906': 'GRIN2D',
               'hsa:4088': 'SMAD3',
               'hsa:6776': 'STAT5A',
               'hsa:182': 'JAG1',
               'hsa:3708': 'ITPR1',
               'hsa:1293': 'COL6A3',
               'hsa:93034': 'NT5C1B',
               'hsa:574537': 'UGT2A2',
               'hsa:11044': 'PAPD7',
               'hsa:57292': 'KIR2DL5A',
               'hsa:107984026': 'ACAP3'
               }

if __name__ == '__main__':
    gm = GeneMapper()
