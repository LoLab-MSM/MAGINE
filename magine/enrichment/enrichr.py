import json
import os
import re

import numpy as np
import pandas as pd
import requests

from magine.plotting.species_plotting import write_table_to_html

_path = os.path.dirname(__file__)

_valid_libs = set()
with open(os.path.join(_path, '_valid_enricher_libs.txt'), 'r') as f:
    for n in f.read().splitlines():
        _valid_libs.add(n)

gene = 'gene'

ontologies = [

    'GO_Biological_Process_2017b',
    'GO_Molecular_Function_2017b',
    'GO_Cellular_Component_2017b',

    'GO_Biological_Process_2017',
    'GO_Molecular_Function_2017',
    'GO_Cellular_Component_2017',

    'Human_Phenotype_Ontology',
    'MGI_Mammalian_Phenotype_2017',
    'HMDB_Metabolites',
    'Tissue_Protein_Expression_from_ProteomicsDB',
]

pathways = [
    'KEGG_2016',
    'NCI-Nature_2016',
    'Panther_2016',
    'WikiPathways_2016',
    'BioCarta_2016',
    'Humancyc_2016',
    'Reactome_2016',
]

kinase = [
    'KEA_2015',
    'LINCS_L1000_Kinase_Perturbations_down',
    'LINCS_L1000_Kinase_Perturbations_up',
    'Kinase_Perturbations_from_GEO_down',
    'Kinase_Perturbations_from_GEO_up',

]
transcription_factors = [
    'ChEA_2016',  # should use RNAseq
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
    'ENCODE_TF_ChIP-seq_2015',
    'Transcription_Factor_PPIs',
    'TRANSFAC_and_JASPAR_PWMs',
]

disease_drug = [
    'DrugMatrix',
    'Drug_Perturbations_from_GEO_2014',
    'Drug_Perturbations_from_GEO_down',
    'Drug_Perturbations_from_GEO_up',
    'OMIM_Disease',
    'OMIM_Expanded',
]

all_dbs = [
    'GO_Biological_Process_2017',
    'GO_Molecular_Function_2017',
    'GO_Cellular_Component_2017',
    'KEGG_2016',
    'NCI-Nature_2016',
    'Panther_2016',
    'WikiPathways_2016',
    'BioCarta_2016',
    'Humancyc_2016',
    'Reactome_2016',
    'KEA_2015',
    'ChEA_2016',
    'DrugMatrix',
    'Drug_Perturbations_from_GEO_2014',

]


def clean_term_names(row):
    term_name = row['term_name']

    if not isinstance(term_name, str):
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


class Enrichr(object):
    query = '{url}/enrich?userListId={list_id}&backgroundType={lib}'

    def __init__(self, exp_data=None, verbose=False):
        self._url = 'http://amp.pharm.mssm.edu/Enrichr'
        self._valid_libs = _valid_libs
        self.exp_data = exp_data
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

            df = self.run_id(user_list_id, gene_set_lib)
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

    def run_id(self, list_id, gene_set_lib):
        if gene_set_lib not in _valid_libs:
            print("{} not in valid ids {}".format(gene_set_lib, _valid_libs))
            return pd.DataFrame()

        q = self.query.format(url=self._url, list_id=list_id, lib=gene_set_lib)

        response = requests.get(q)

        if not response.ok:
            # quick way to wait for response
            while not response.ok:
                response = requests.get(q)

        data = json.loads(response.text)
        if len(data[gene_set_lib]) == 0:
            return pd.DataFrame()
        #####
        # ENRICHR return a list of entries with each entry having these terms
        # Rank, Term name, P-value, Z-score, Combined score, Overlapping genes,
        # Adjusted p-value, Old p-value, Old adjusted p-value
        #####

        df = pd.DataFrame(
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

        if df.shape[0] == 0:
            return df[cols]
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
        data = []
        for count, i in enumerate(databases):
            if self.verbose:
                print('\t\t{}/{} databases'.format(count, len(databases)))
            data.append(self.run_id(gene_list_id, i))

        return pd.concat(data)

    def run_samples(self, sample_lists, sample_ids,
                    gene_set_lib='GO_Biological_Process_2017', save_name=None,
                    create_html=False, out_dir=None, run_parallel=False,
                    exp_data=None):
        """

        Parameters
        ----------
        sample_lists : list_like
            List of lists of genes for enrichment analysis
        sample_ids : list
            list of ids for the provided sample list
        gene_set_lib : str
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
        Returns
        -------

        """
        assert isinstance(sample_lists, list), "List required"
        assert isinstance(sample_lists[0], list), "List of lists required"

        df_all = []
        for i, j in zip(sample_lists, sample_ids):
            df = self.run(i, gene_set_lib)
            df['sample_id'] = j
            df['db'] = gene_set_lib
            df_all.append(df)

        df_all = pd.concat(df_all)
        df_all = self._filter_sig_across_term(df_all)

        if df_all.shape[0] == 0:
            print("No significant terms for {}".format(gene_set_lib))
            return None

        index = ['term_name']

        if 'term_id' in list(df.columns):
            index.insert(0, 'term_id')

        p_df = pd.pivot_table(df_all, index=index,
                              columns='sample_id',
                              values=['combined_score', 'adj_p_value',
                                      'rank', 'p_value', 'z_score', 'genes'],
                              aggfunc='first', fill_value=np.nan
                              )

        if save_name:
            s_name = '{}_enrichr'.format(save_name)
            # save files
            p_df.to_excel('{}.xlsx'.format(s_name), merge_cells=True)
            df_all.to_csv('{}.csv'.format(s_name), index=False)

        if create_html:
            if exp_data is None:
                exp_data = self.exp_data

            if exp_data is None:
                print("exp_data required to make plots over samples")
                return p_df

            write_table_to_html(data=df_all, save_name=save_name,
                                out_dir=out_dir, run_parallel=run_parallel,
                                exp_data=exp_data)
        return p_df

    def run_key_dbs(self, list_g, labels, save_name, create_html=False,
                    run_parallel=False):
        """

        Parameters
        ----------
        list_g : list_like
            List of lists of genes
        labels : list_like
            List of lists of labels
        save_name : str
            prefix of names
        create_html : bool
            if you want to create html of all the outputs
        run_parallel : bool
            Create plots in parallel
        Returns
        -------

        """

        print("Working on {}".format(save_name))

        def _run(list_dbs, db_name):
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter('{}_{}.xlsx'.format(save_name, db_name),
                                    engine='xlsxwriter')
            if self.verbose:
                print("\tWorking on {}".format(db_name))
            html_list = []
            for i in list_dbs:
                print("\t\tRunning {}".format(i))
                if create_html:
                    outdir = save_name + '_' + i
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                else:
                    outdir = None

                df = self.run_samples(sample_lists=list_g, sample_ids=labels,
                                      gene_set_lib=i, create_html=create_html,
                                      save_name=outdir, out_dir=outdir,
                                      run_parallel=run_parallel,
                                      )
                if df is None:
                    continue

                html_list.append(dict(database=i,
                                      filename='{}_filter.html'.format(
                                          outdir)))
                if len(i) > 30:
                    i = i[:30]
                df.to_excel(writer, sheet_name=i, freeze_panes=(1, 2))

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()
            return [db_name, html_list, save_name + '_' + db_name]

        h1 = _run(transcription_factors, 'transcription_factors')
        h2 = _run(pathways, 'pathways')
        h3 = _run(kinase, 'kinases')
        h4 = _run(disease_drug, 'drug_and_disease')
        h5 = _run(ontologies, 'ontologies')
        return [h1, h2, h3, h4, h5]

    def run_set_of_dbs(self, list_g, db='drug'):
        """

        Parameters
        ----------
        list_g : list_like
            List of lists of genes
        db : str
        Returns
        -------
        pandas.DataFrame

        """

        if db == 'transcription_factors':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=transcription_factors)

        if db == 'pathways':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=pathways)
        if db == 'kinases':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=kinase)
        if db == 'drug':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=disease_drug)
        if db == 'ontologies':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=ontologies)
        if db == 'all':
            return self.run(list_of_genes=list_g,
                            gene_set_lib=all_dbs)

    def run_sample_set_of_dbs(self, sample_lists, sample_ids, databases,
                              save_name=None, pivot=True):
        """

        Parameters
        ----------
        sample_lists : list_like
            List of lists, where each list contains genes
        sample_ids : list_like
            Names for each sample_id
        save_name : str, optional
            If you want to save the output as csv/xlsx
        databases : list_like
            Database set.

        pivot: bool
            Pivot the table to have additional columns per sample_id

        Returns
        -------

        """
        assert isinstance(sample_lists, list), "Please provide list of lists"
        assert isinstance(sample_lists[0],
                          list), "Please provide list of lists"

        assert isinstance(databases, list), 'please provide database list'

        df_all = []
        for i, j in zip(sample_lists, sample_ids):

            if self.verbose:
                print("{} / {}".format(j, sample_ids))
            tmp_df = self.run(i, databases)
            tmp_df['sample_id'] = j
            df_all.append(tmp_df)

        # stack results
        df_all = pd.concat(df_all, ignore_index=True)

        # filter terms that are not significant in any sample
        df_all = self._filter_sig_across_term(df_all)

        # verify results exist
        if df_all.shape[0] == 0:
            print("No significant terms for {}".format(databases))
            return None

        if pivot:
            # pivot output
            index = ['term_name', 'db']

            if 'term_id' in list(df_all.columns):
                index.insert(0, 'term_id')

            df_all = pd.pivot_table(
                df_all, index=index, columns='sample_id', aggfunc='first',
                fill_value=np.nan,
                values=['term_name', 'rank', 'p_value', 'z_score',
                        'combined_score', 'adj_p_value', 'genes'],
            )

        if save_name:
            df_all.to_excel('{}_enricher.xlsx'.format(save_name),
                            merge_cells=True)
            df_all.to_csv('{}_enrichr.csv'.format(save_name), index=False,
                          encoding='utf8')

        return df_all

    def run_dbs_testing(self, list_g, labels, save_name):
        """

        Parameters
        ----------
        list_g : list_like
            List of lists of genes
        labels : list_like
            List of lists of labels
        save_name : str
            prefix of names
        Returns
        -------

        """

        print("Working on {}".format(save_name))

        names = []

        def _run(list_dbs, db_name):
            # Create a Pandas Excel writer using XlsxWriter as the engine.

            print("\tWorking on {}".format(db_name))
            all_df = []
            for i in list_dbs:
                print("\t\tRunning {}".format(i))

                df = self.run_samples(sample_lists=list_g, sample_ids=labels,
                                      gene_set_lib=i
                                      )
                if df is None:
                    continue
                df['db'] = i
                all_df.append(df)

            all_df = pd.concat(all_df)
            sn = '{}_{}.csv'.format(db_name, save_name)
            all_df.to_csv(sn, encoding='utf8')
            names.append(sn)

        _run(transcription_factors, 'transcription_factors')
        _run(pathways, 'pathways')
        _run(kinase, 'kinases')
        _run(disease_drug, 'drug_and_disease')
        _run(ontologies, 'ontologies')

        writer = pd.ExcelWriter('{}.xlsx'.format(save_name),
                                engine='xlsxwriter')
        for i in names:
            df = pd.read_csv(i)
            df.to_excel(writer, sheet_name=i.split('_')[0], merge_cells=True,
                        freeze_panes=(2, 1), encoding='utf8')
        writer.save()

    @staticmethod
    def _filter_sig_across_term(data, p_value_thresh=0.05):
        """
        removes terms if there are no significant terms
        Parameters
        ----------
        data : pd.DataFrame


        Returns
        -------

        """
        non_sig = []
        for i, g in data.groupby(['term_name', 'db']):
            if len(g[g['adj_p_value'] <= p_value_thresh]) == 0:
                non_sig.append(i[0])

        return data[~data['term_name'].isin(non_sig)]


db_types = {
    'transcription': [
        'ChEA_2016',
        'TRANSFAC_and_JASPAR_PWMs',
        'ARCHS4_TFs_Coexp',
        'Enrichr_Submissions_TF-Gene_Coocurrence',
        'Genome_Browser_PWMs',
        'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
        'Epigenomics_Roadmap_HM_ChIP-seq',
        'TargetScan_microRNA_2017',
        'miRTarBase_2017',
        'ENCODE_TF_ChIP-seq_2015',
        'TF-LOF_Expression_from_GEO',
        'ENCODE_Histone_Modifications_2015',
        'Transcription_Factor_PPIs'
    ],
    'pathways': [
        'KEGG_2016',
        'WikiPathways_2016',
        'ARCHS4_Kinases_Coexp',
        'Reactome_2016',
        'BioCarta_2016',
        'NCI-Nature_2016',
        'Panther_2016',
        'BioPlex_2017',
        'PPI_Hub_Proteins',
        'KEA_2015',
        'LINCS_L1000_Kinase_Perturbations_down',
        'LINCS_L1000_Kinase_Perturbations_up',
        'Kinase_Perturbations_from_GEO_down',
        'Kinase_Perturbations_from_GEO_up',
        'NURSA_Human_Endogenous_Complexome',
        'CORUM',
        'Phosphatase_Substrates_from_DEPOD',
    ],
    'ontologies': [
        'GO_Cellular_Component_2017b',
        'GO_Biological_Process_2017b',
        'GO_Molecular_Function_2017b',
        'MGI_Mammalian_Phenotype_2017',
        'Human_Phenotype_Ontology',
        'Jensen_TISSUES',
        'Jensen_COMPARTMENTS',
        'Jensen_DISEASES',
    ],
    'disease_drug': [
        'LINCS_L1000_Chem_Pert_up',
        'LINCS_L1000_Chem_Pert_down',
        'LINCS_L1000_Ligand_Perturbations_up',
        'ARCHS4_IDG_Coexp',
        'DrugMatrix',
        'Old_CMAP_up',
        'Old_CMAP_down',
        'GeneSigDB',
        'OMIM_Disease',
        'OMIM_Expanded',
        'VirusMINT',
        'MSigDB_Oncogenic_Signatures',
        'Virus_Perturbations_from_GEO_up',
        'Virus_Perturbations_from_GEO_down',
    ],
    'cell_type': [
        'Human_Gene_Atlas',
        'ARCHS4_Tissues',
        'ARCHS4_Cell-lines',
        'Allen_Brain_Atlas_up',
        'Allen_Brain_Atlas_down',
        'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
        'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
        'Cancer_Cell_Line_Encyclopedia',
        'NCI-60_Cancer_Cell_Lines',
        'ESCAPE',
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


def run_enrichment_for_project(exp_data, project_name):

    local_dbs = [

        'KEGG_2016',
        'NCI-Nature_2016',
        'Panther_2016',
        'WikiPathways_2016',
        'BioCarta_2016',
        'Humancyc_2016',
        'Reactome_2016',
        'KEA_2015',


        # ontologies
        'GO_Biological_Process_2017',
        'GO_Molecular_Function_2017',
        'GO_Cellular_Component_2017',

        'GO_Biological_Process_2017b',
        'GO_Molecular_Function_2017b',
        'GO_Cellular_Component_2017b',

        # transcription
        'ChEA_2016',
        'TRANSFAC_and_JASPAR_PWMs',
        'ENCODE_TF_ChIP-seq_2015',

        # cell types
        # 'Jensen_TISSUES',
        'Human_Gene_Atlas',
        'Tissue_Protein_Expression_from_ProteomicsDB',
        'ARCHS4_Cell-lines',
        'ARCHS4_Tissues',
        'Cancer_Cell_Line_Encyclopedia',
        'NCI-60_Cancer_Cell_Lines',

        # pertubations
        'Kinase_Perturbations_from_GEO_down',
        'Kinase_Perturbations_from_GEO_up',
        'LINCS_L1000_Kinase_Perturbations_down',
        'LINCS_L1000_Kinase_Perturbations_up',
        'Ligand_Perturbations_from_GEO_down',
        'Ligand_Perturbations_from_GEO_up',
        'Old_CMAP_down',
        'Old_CMAP_up',

        # phenotypes
        'Human_Phenotype_Ontology',
        'MGI_Mammalian_Phenotype_2017',
        'Jensen_DISEASES',
        'dbGaP',
        'DrugMatrix',
        'Drug_Perturbations_from_GEO_2014',

    ]
    local_dbs = []
    for i, j in db_types.items():
        local_dbs += j
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
            all_df.append(df)

    pt = exp.proteomics_sample_ids
    rt = exp.rna_sample_ids
    if len(pt) != 0:
        _run_new(exp.proteomics_by_sample_id, pt, 'proteomics_both')
        _run_new(exp.proteomics_up_by_sample_id, pt, 'proteomics_down')
        _run_new(exp.proteomics_down_by_sample_id, pt, 'proteomics_up')
    if len(rt) != 0:
        _run_new(exp.rna_down_over_time, rt, 'rna_down')
        _run_new(exp.rna_up_over_time, rt, 'rna_up')
        _run_new(exp.rna_over_time, rt, 'rna_both')
    final_df = pd.concat(all_df, ignore_index=True)
    final_df = final_df[
        ['term_name', 'rank', 'combined_score', 'adj_p_value', 'genes',
         'n_genes', 'sample_id', 'category', 'db']
    ]
    print(final_df.shape)
    final_df = final_df[~final_df['term_name'].isnull()]
    print(final_df.shape)
    final_df['term_name'] = final_df.apply(clean_term_names, axis=1)
    print(final_df.shape)
    final_df = final_df[~final_df['term_name'].isnull()]
    print(final_df.shape)
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
        # term_to_gene[term] = '|'.join(sorted(set(i for i in genes_dict)))
        term_to_gene.append(
            dict(term=term.lower(), gene_list=genes, n_genes=len(genes))
        )

    return term_to_gene


if __name__ == '__main__':
    e = Enrichr()
    g_list = ['PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2',
              'P2RX7', 'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1',
              'ANAPC16', 'TMCC1', 'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2',
              'PLBD2', 'LARP7', 'TECPR2', 'ZNF302', 'CUX1', 'MOB2', 'CYTH2',
              'SEC22C', 'EIF4E3', 'ROBO2', 'ADAMTS9-AS2', 'CXXC1', 'LINC01314',
              'ATF7', 'ATP5F1']

    df2 = e.run(g_list, ['GO_Biological_Process_2017b', 'KEA_2015'])
    print(df2['db'].unique())
    print(df2.head(10))
    print(df2['term_name'].values)
    quit()
    print(df2.head(10))
    print(list(df2['term_name'])[0:2])

    lists = [['BAX', 'BCL2', 'CASP3'], ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'], save_name='test')
    df2 = e.run_sample_set_of_dbs(lists, ['1', '2', '3'], save_name='test')
