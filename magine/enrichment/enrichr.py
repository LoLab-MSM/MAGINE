import json
import os
import re
import time
# Will be OK in Python 2
try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str

import numpy as np
import pandas as pd
import requests

from magine.plotting.species_plotting import write_table_to_html
from magine.enrichment.enrichment_result import EnrichmentResult

_path = os.path.dirname(__file__)

_valid_libs = set()
with open(os.path.join(_path, '_valid_enricher_libs.txt'), 'r') as f:
    for n in f.read().splitlines():
        _valid_libs.add(n)

gene = 'gene'

db_types = {
    'histone': [
        'Epigenomics_Roadmap_HM_ChIP-seq',
        'ENCODE_Histone_Modifications_2015',
        'ESCAPE',
    ],
    'mrna': [
        'TargetScan_microRNA_2017',
        'miRTarBase_2017',
    ],
    'kinases': [
        'KEA_2015',
        'LINCS_L1000_Kinase_Perturbations_down',
        'LINCS_L1000_Kinase_Perturbations_up',
        'Kinase_Perturbations_from_GEO_down',
        'Kinase_Perturbations_from_GEO_up',
        'Phosphatase_Substrates_from_DEPOD',
        'ARCHS4_Kinases_Coexp',
        'ARCHS4_IDG_Coexp',
    ],
    'transcription': [
        'ChEA_2016',
        'TRANSFAC_and_JASPAR_PWMs',
        'ARCHS4_TFs_Coexp',
        'Enrichr_Submissions_TF-Gene_Coocurrence',
        'Genome_Browser_PWMs',
        'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
        'ENCODE_TF_ChIP-seq_2015',
        'TF-LOF_Expression_from_GEO',
        'Transcription_Factor_PPIs'
    ],
    'pathways': [
        'KEGG_2016',
        'WikiPathways_2016',
        'Reactome_2016',
        'BioCarta_2016',
        'NCI-Nature_2016',
        'Panther_2016',
        'BioPlex_2017',
    ],
    'complex': [
        'NURSA_Human_Endogenous_Complexome',
        'CORUM',
        'PPI_Hub_Proteins',
    ],
    'ontologies': [
        'GO_Cellular_Component_2017b',
        'GO_Biological_Process_2017b',
        'GO_Molecular_Function_2017b',
        'MGI_Mammalian_Phenotype_2017',
        'Jensen_COMPARTMENTS',
    ],
    'drug': [
        'DrugMatrix',
        'Drug_Perturbations_from_GEO_2014',
        'Old_CMAP_up',
        'Old_CMAP_down',
        'LINCS_L1000_Chem_Pert_up',
        'LINCS_L1000_Chem_Pert_down',
        'LINCS_L1000_Ligand_Perturbations_up',
        'LINCS_L1000_Ligand_Perturbations_down',
    ],
    'disease': [
        'OMIM_Disease',
        'OMIM_Expanded',
        'Jensen_DISEASES',
        'Human_Phenotype_Ontology'
    ],
    'cell_type': [

        'ARCHS4_Tissues',
        'ARCHS4_Cell-lines',
        'Allen_Brain_Atlas_up',
        'Allen_Brain_Atlas_down',
        'Cancer_Cell_Line_Encyclopedia',
        'Human_Gene_Atlas',
        'NCI-60_Cancer_Cell_Lines',
        'Jensen_TISSUES',
    ],
    'crowd': [
        'Disease_Perturbations_from_GEO_down',
        'Disease_Perturbations_from_GEO_up',
        'Drug_Perturbations_from_GEO_down',
        'Drug_Perturbations_from_GEO_up',
        'Single_Gene_Perturbations_from_GEO_up',
        'Single_Gene_Perturbations_from_GEO_down',
        'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
        'Aging_Perturbations_from_GEO_up',
        'Aging_Perturbations_from_GEO_down',
        'Ligand_Perturbations_from_GEO_up',
        'Ligand_Perturbations_from_GEO_down',
        'MCF7_Perturbations_from_GEO_up',
        'MCF7_Perturbations_from_GEO_down',
        'Microbe_Perturbations_from_GEO_up',
        'Microbe_Perturbations_from_GEO_down',
        'SysMyo_Muscle_Gene_Sets',
    ]

}


class Enrichr(object):
    query = '{url}/enrich?userListId={list_id}&backgroundType={lib}'

    def __init__(self, verbose=False):
        self._url = 'http://amp.pharm.mssm.edu/Enrichr'
        self._valid_libs = _valid_libs
        self.verbose = verbose

    def print_valid_libs(self):
        """
        Print a list of all available libraries EnrichR has to offer.
        Returns
        -------

        """
        for lib_name in sorted(self._valid_libs):
            print(lib_name)

    def run(self, list_of_genes, gene_set_lib='GO_Biological_Process_2017'):
        """

        Parameters
        ----------
        list_of_genes : list_like
            List of genes using HGNC gene names
        gene_set_lib : str or list
            Name of gene set library
            To print options use Enrichr.print_valid_libs
        Returns
        -------

        """
        assert isinstance(list_of_genes, (list, set))

        if self.verbose:
            print("Running Enrichr with gene set {}".format(gene_set_lib))

        user_list_id = self._add_gene_list(list_of_genes)
        if isinstance(gene_set_lib, str):

            df = self._run_id(user_list_id, gene_set_lib)
        else:
            df = self._run_list_of_dbs(user_list_id, gene_set_lib)

        init_size = len(df)
        df['term_name'] = df.apply(clean_term_names, axis=1)
        after_size = len(df)
        assert init_size == after_size, 'not the same shape {}'.format(
            gene_set_lib)

        if self.verbose:
            print("Done calling Enrichr.")

        return df

    def _run_id(self, list_id, gene_set_lib):
        if gene_set_lib not in _valid_libs:
            print("{} not in valid ids {}".format(gene_set_lib, _valid_libs))
            return EnrichmentResult()

        q = self.query.format(url=self._url, list_id=list_id, lib=gene_set_lib)
        while True:
            try:
                response = requests.get(q)
                break
            except:
                time.sleep(1)

        if not response.ok:
            # quick way to wait for response
            while not response.ok:
                response = requests.get(q)

        data = json.loads(response.text)
        if len(data[gene_set_lib]) == 0:
            return EnrichmentResult()
        #####
        # ENRICHR return a list of entries with each entry having these terms
        # Rank, Term name, P-value, Z-score, Combined score, Overlapping genes,
        # Adjusted p-value, Old p-value, Old adjusted p-value
        #####

        df = EnrichmentResult(
            data[gene_set_lib],
            columns=['rank', 'term_name', 'p_value', 'z_score',
                     'combined_score', 'gene_hits', 'adj_p_value', '_', '_']
        )

        def compress_genes(row):
            return ','.join(g for g in sorted(row['gene_hits']))

        def get_length(row):
            return len(row['gene_hits'])

        df['genes'] = df.apply(compress_genes, axis=1)
        df['n_genes'] = df.apply(get_length, axis=1)

        cols = ['term_name', 'rank', 'p_value', 'z_score', 'combined_score',
                'adj_p_value', 'genes', 'n_genes']

        df = df[~df['term_name'].isnull()]
        df = df[cols]
        df['db'] = gene_set_lib
        return df

    def _add_gene_list(self, gene_list):

        genes_str = '\n'.join(gene_list)

        payload = {
            'list': (None, genes_str),
            'description': (None, 'MAGINE analysis')
        }

        response = requests.post(self._url + '/addList', files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data = json.loads(response.text)
        return data['userListId']

    def _run_list_of_dbs(self, gene_list_id, databases):
        for count, i in enumerate(databases):
            if count == 0:
                data = self._run_id(gene_list_id, i)
            else:
                data.append(self._run_id(gene_list_id, i))
            if self.verbose:
                print('\t\t{}/{} databases'.format(count, len(databases)))
        return data

    def run_samples(self, sample_lists, sample_ids,
                    database='GO_Biological_Process_2017', save_name=None,
                    create_html=False, out_dir=None, run_parallel=False,
                    exp_data=None, pivot=False):
        """

        Parameters
        ----------
        sample_lists : list_like
            List of lists of genes for enrichment analysis
        sample_ids : list
            list of ids for the provided sample list
        database : str, list
            Type of gene set, refer to Enrichr.print_valid_libs
        save_name : str, optional
            if provided it will save a file as a pivoted table with
            the term_ids vs sample_ids
        create_html : bool
            Creates html of output with plots of species across sample
        out_dir : str
            If create_html, it will place all html plots into this directory
        run_parallel : bool
            If create_html, it will create plots using multiprocessing
        exp_data : magine.data.ExperimentalData
            Must be provided if create_html=True
        pivot : bool

        Returns
        -------
        EnrichmentResult
        """
        assert isinstance(sample_lists, list), "List required"
        assert isinstance(sample_lists[0],
                          (list, set)), "List of lists required"

        for count, (i, j) in enumerate(zip(sample_lists, sample_ids)):
            df = self.run(i, database)
            df['sample_id'] = j
            if count == 0:
                df_all = df
            else:
                df_all.append(df)

        df_all = self._filter_sig_across_term(df_all)

        if save_name:
            s_name = '{}_enrichr'.format(save_name)
            if pivot:
                # pivot output
                index = ['term_name', 'db']

                p_df = pd.pivot_table(
                    df_all, index=index, columns='sample_id', aggfunc='first',
                    fill_value=np.nan,
                    values=['term_name', 'rank', 'p_value', 'z_score',
                            'combined_score', 'adj_p_value', 'genes'],
                )
                # save files
                p_df.to_excel('{}.xlsx'.format(s_name), merge_cells=True)
            df_all.to_csv('{}.csv'.format(s_name), index=False)

        if create_html:
            if exp_data is None:
                print("exp_data required to make plots over samples")
                return df_all

            write_table_to_html(data=df_all, save_name=save_name,
                                out_dir=out_dir, run_parallel=run_parallel,
                                exp_data=exp_data)
        return df_all

    @staticmethod
    def _filter_sig_across_term(data, p_value_thresh=0.05):
        """
        removes terms if there are no significant terms
        Parameters
        ----------
        data : EnrichmentResult


        Returns
        -------

        """
        non_sig = []
        for i, g in data.groupby(['term_name', 'db']):
            if len(g[g['adj_p_value'] <= p_value_thresh]) == 0:
                non_sig.append(i[0])

        return data[~data['term_name'].isin(non_sig)]


def clean_term_names(row):
    term_name = row['term_name']
    if not isinstance(term_name, basestring):
        return term_name

    db = row['db']

    if db in ['GO_Biological_Process_2017', 'GO_Biological_Process_2017b',
              'GO_Molecular_Function_2017', 'GO_Molecular_Function_2017b',
              'GO_Cellular_Component_2017', 'GO_Cellular_Component_2017b']:

        if 'GO:' in term_name:
            term_name = term_name.split('(GO:', 1)[0]
        elif 'go:' in term_name:
            term_name = term_name.split('(go:', 1)[0]
    if db == 'Human_Phenotype_Ontology':
        if 'HP:' in term_name:
            term_name = term_name.split('(HP:', 1)[0]

    if db == 'MGI_Mammalian_Phenotype_2017':
        if term_name.startswith('MP:'):
            term_name = term_name.split('_', 1)[1]

    if db == 'DrugMatrix':
        drug_name = re.search(r'^(.*)(-\d*.*\d_)', term_name).group(1)
        direction = re.search(r'-(.{2})$', term_name).group(0)
        term_name = drug_name + direction
    if db in ['LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Chem_Pert_down', ]:
        exp_id, rest = term_name.split('-', 1)
        drug_name, _ = rest.rsplit('-', 1)
        term_name = drug_name
    if db in ['Old_CMAP_down', 'Old_CMAP_up', ]:
        term_name = term_name.rsplit('-', 1)[0]
    if db in ['Ligand_Perturbations_from_GEO_down',
              'Ligand_Perturbations_from_GEO_up']:
        term_name = term_name.split('_', 1)[0]

    term_name = term_name.strip()
    term_name = term_name.lower()
    for i, j in replace_pairs:
        if i in term_name:
            term_name = term_name.replace(i, j)

    return term_name


def clean_tf_names(data):
    def return_single_tf(row):
        """
        Gets TF name only
        """
        tf = row['term_name']
        if '_' in tf:
            tf = tf.split('_')[0]
        return tf.upper()

    tf_dbs = ['ARCHS4_TFs_Coexp',
              'ChEA_2016',
              'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
              'ENCODE_TF_ChIP-seq_2015',
              'Enrichr_Submissions_TF-Gene_Coocurrence',
              'TRANSFAC_and_JASPAR_PWMs',
              'TF-LOF_Expression_from_GEO',
              'Transcription_Factor_PPIs']
    tfs = data[data['db'].isin(tf_dbs)].copy()
    tfs['term_name'] = tfs.apply(return_single_tf, axis=1)
    return tfs


replace_pairs = [
    ('mus musculus', 'mus'),
    ('homo sapiens', 'hsa'),
    ('rattus norvegicus', 'rat'),
    ('oncorhynchus mykiss', 'onc_mykiss'),
    ('macaca fascicularis', 'mfa'),
    ('oryzias latipes', 'ola'),
    ('sus scrofa', 'ssc'),
    ('danio rerio', 'dre'),
    ('bos taurus', 'bta'),
    ('dictyostelium discoideum', 'dicty'),
    ('myzus persicae', 'm.persicae'),
    ('staphylococcus aureus', 's.aureus'),
    ('escherichia coli', 'e.coli'),
    ('pseudomonas aeruginosa', 'p.aeruginosa'),
    ('drosophila melanogaste', 'fly'),
    ('mycobacterium tuberculosis', 'm.tuberculosis'),
    ('hordeum vulgare', 'h.vulgare'),
    (' (mouse)', '_mus'),
    (' (human)', '_hsa'),
    ('Homo sapiens', 'hsa'),
    ('Mus musculus', 'hsa'),
    ('hg19', '_hsa'),
    ('mm9', '_mus'),

]


def run_enrichment_for_project(exp_data, project_name):
    """

    Parameters
    ----------
    exp_data : magine.data.experimental_data.ExprerimentalData
    project_name : str

    Returns
    -------

    """

    local_dbs = []
    for i in ['drug', 'disease', 'ontologies', 'pathways', 'transcription',
              'kinases', 'histone', 'cell_type']:
        local_dbs += db_types[i]
    e = Enrichr(verbose=True)
    exp = exp_data
    all_df = []
    _dir = 'enrichment_output'
    if not os.path.exists(_dir):
        os.mkdir(_dir)

    print("Running {} databases, might take some time".format(len(local_dbs)))

    def _run_new(samples, timepoints, category):
        print("Running {}".format(category))
        for genes, sample_id in zip(samples, timepoints):
            print('\t time point = {}'.format(sample_id))
            current = "{}_{}_{}".format(str(category),
                                        str(sample_id),
                                        project_name)
            name = os.path.join(_dir, current + '.csv.gz')
            try:
                df = pd.read_csv(name, index_col=None, encoding='utf-8')
            except:
                df = e.run(genes, local_dbs)
                df['sample_id'] = sample_id
                df['category'] = category
                df.to_csv(name, index=False, encoding='utf-8',
                          compression='gzip')
            df['sample_id'] = sample_id
            df['category'] = category
            all_df.append(df)

    pt = exp.proteins.sample_ids
    rt = exp.rna.sample_ids

    if len(pt) != 0:
        _run_new(exp.proteins.sig.by_sample, pt, 'proteomics_both')
        _run_new(exp.proteins.sig.up_by_sample, pt, 'proteomics_up')
        _run_new(exp.proteins.sig.down_by_sample, pt, 'proteomics_down')

    if len(rt) != 0:
        _run_new(exp.rna.by_sample, rt, 'rna_both')
        _run_new(exp.rna.sig.down_by_sample, rt, 'rna_down')
        _run_new(exp.rna.up.up_by_sample, rt, 'rna_up')

    final_df = pd.concat(all_df, ignore_index=True)

    final_df = final_df[
        ['term_name', 'rank', 'combined_score', 'adj_p_value', 'genes',
         'n_genes', 'sample_id', 'category', 'db']
    ]
    final_df = final_df[~final_df['term_name'].isnull()]
    final_df.to_csv('{}.csv.gz'.format(project_name), encoding='utf-8',
                    compression='gzip')
    print("Done with enrichment")


def get_background_list(lib_name):
    """
    Return reference list for given gene referecen set

    Parameters
    ----------
    lib_name : str

    Returns
    -------

    """

    # http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Genes_Associated_with_NIH_Grants

    enrichment_url = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary'
    query_string = '?userListId&libraryName=%s'
    response = requests.get(
        enrichment_url + query_string % lib_name
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    results = json.loads(response.text)

    term_to_gene = []
    assert lib_name in results

    for term, genes_dict in results[lib_name]['terms'].items():
        genes = sorted(set(i for i in genes_dict))
        term_to_gene.append(
            dict(term=term.lower(), gene_list=genes, n_genes=len(genes))
        )

    return term_to_gene
