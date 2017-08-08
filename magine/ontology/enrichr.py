import json
import os
import re
import time

import numpy as np
import pandas as pd
import pathos.multiprocessing as mp
import requests

import magine.html_templates.html_tools as html_tools
from magine.plotting.species_plotting import plot_list_of_genes

_path = os.path.dirname(__file__)

_valid_libs = set()
with open(os.path.join(_path, '_valid_enricher_libs.txt'), 'r') as f:
    for n in f.read().splitlines():
        _valid_libs.add(n)

gene = 'gene'


class Enrichr(object):
    def __init__(self, exp_data=None):
        self._url = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        self._valid_libs = _valid_libs
        self.exp_data = exp_data

    def print_valid_libs(self):
        """
        Print a list of all available libraries EnrichR has to offer.
        Returns
        -------

        """
        for lib_name in sorted(self._valid_libs):
            print(lib_name)

    def run(self, list_of_genes, gene_set_lib='GO_Biological_Process_2017',
            verbose=False):
        """

        Parameters
        ----------
        list_of_genes : list_like
            List of genes using HGNC gene names
        gene_set_lib : str
            Name of gene set library
            To print options use Enrichr.print_valid_libs
        verbose : bool
            print information
        Returns
        -------

        """
        assert isinstance(list_of_genes, list)
        assert gene_set_lib in _valid_libs, \
            "{} not in valid ids {}".format(gene_set_lib, _valid_libs)

        if verbose:
            print("Running Enrichr with gene set {}".format(gene_set_lib))

        description = 'Example gene list'
        genes_str = '\n'.join(list_of_genes)
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        response = requests.post(self._url, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data = json.loads(response.text)
        user_list_id = data['userListId']

        enrichment_url = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        response = requests.get(
            enrichment_url + query_string % (user_list_id, gene_set_lib)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        data = json.loads(response.text)
        #####
        # ENRICHR return a list of entries with each entry having these terms
        # Rank, Term name, P-value, Z-score, Combined score, Overlapping genes,
        # Adjusted p-value, Old p-value, Old adjusted p-value
        #####
        list_of_dict = []
        for i in data[gene_set_lib]:
            tmp_dict = dict()
            tmp_dict['rank'] = i[0]
            tmp_dict['term_name'] = i[1]
            tmp_dict['p_value'] = i[2]
            tmp_dict['z_score'] = i[3]
            tmp_dict['combined_score'] = i[4]
            tmp_dict['genes'] = ','.join(g for g in sorted(i[5]))
            tmp_dict['n_genes'] = len(i[5])
            tmp_dict['adj_p_value'] = i[6]
            list_of_dict.append(tmp_dict)
        cols = ['term_name', 'rank', 'p_value', 'z_score', 'combined_score',
                'adj_p_value', 'genes', 'n_genes']

        df = pd.DataFrame(list_of_dict, columns=cols)

        if df.shape[0] == 0:
            return df[cols]

        if gene_set_lib.startswith('GO'):
            def get_go_id(row):
                s = row['term_name']
                go_id = re.search(r'\((GO.*?)\)', s).group(1)
                return go_id

            def remove_term_id(row):
                s = row['term_id']
                term = row['term_name'].replace('(' + s + ')', '')
                return term

            df['term_id'] = df.apply(get_go_id, axis=1)

            term_names = df.apply(remove_term_id, axis=1)

            df['term_name'] = term_names
            cols.insert(0, 'term_id')
        if verbose:
            print("Done calling Enrichr.")
        return df[cols]

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
        assert isinstance(sample_lists, list), "Please provide list of lists"
        assert isinstance(sample_lists[0],
                          list), "Please provide list of lists"
        df_all = []
        for i, j in zip(sample_lists, sample_ids):
            print("{} / {}".format(j, sample_ids))
            df = self.run(i, gene_set_lib)
            df['sample_id'] = j
            df_all.append(df)
        df_all = pd.concat(df_all)
        non_sig = []
        for i, g in df_all.groupby('term_name'):
            if len(g[g['adj_p_value'] < 0.05]) == 0:
                non_sig.append(i)
        df_all = df_all[~df_all['term_name'].isin(non_sig)]

        if df_all.shape[0] == 0:
            print("No significant terms for {}".format(gene_set_lib))
            return None
        index = ['term_name']

        if 'term_id' in list(df.columns):
            index.insert(0, 'term_id')

        p_df = pd.pivot_table(df_all, index=index,
                              columns='sample_id',
                              values=['term_name', 'rank', 'p_value',
                                      'z_score',
                                      'combined_score', 'adj_p_value',
                                      'genes'],
                              aggfunc='first', fill_value=np.nan
                              )

        if save_name:
            p_df.to_excel('{}_enricher.xlsx'.format(save_name),
                          merge_cells=True)

            df_all.to_csv('{}_enrichr.csv'.format(save_name), index=False)
        if create_html:
            if exp_data is None:
                exp_data = self.exp_data

            if exp_data is None:
                print("exp_data required to make plots over samples")
                quit()
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

        disease_drug = [
            'DrugMatrix',
            'Drug_Perturbations_from_GEO_2014',
            'LINCS_L1000_Chem_Pert_down',
            'LINCS_L1000_Chem_Pert_up',
            'OMIM_Disease',
            'OMIM_Expanded',
        ]

        ontologies = [

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

        def _run(list_dbs, db_name):
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter('{}_{}.xlsx'.format(save_name, db_name),
                                    engine='xlsxwriter')
            html_list = []
            for i in list_dbs:
                print("Running {}".format(i))
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
                df.to_excel(writer, sheet_name=i, freeze_panes=(2, 2))

            # Close the Pandas Excel writer and output the Excel file.
            # writer.save()
            return [db_name, html_list, save_name + '_' + db_name]

        h1 = _run(transcription_factors, 'transcription_factors')
        h2 = _run(pathways, 'pathways')
        h3 = _run(kinase, 'kinases')
        h4 = _run(disease_drug, 'drug_and_disease')
        h5 = _run(ontologies, 'ontologies')

        return [h1, h2, h3, h4, h5]


def write_table_to_html(data, save_name='index', out_dir=None,
                        run_parallel=False, exp_data=None):
    """
    Creates a table of all the plots of all genes for each GO term.

    Uses last calculated enrichment array.

    Parameters
    ----------
    data : pandas.DataFrame
    save_name : str
        name of html output file
    out_dir : str, optional
        output path for all plots
    run_parallel : bool
        Create plots in parallel

    Returns
    -------

    """

    # print(data.dtypes)
    # tmp = pivot_table_for_export(data)
    # print(tmp.dtypes)
    list_of_terms = list(data['term_name'].unique())
    fig_dict, to_remove = create_gene_plots(data=data,
                                            list_of_terms=list_of_terms,
                                            save_name=save_name,
                                            out_dir=out_dir,
                                            exp_data=exp_data,
                                            run_parallel=run_parallel
                                            )

    for i in fig_dict:
        data.loc[data['term_name'] == i, 'term_name'] = fig_dict[i]

    data = data[~data['term_name'].isin(to_remove)]

    index = ['term_name']
    if 'term_id' in list(data.columns):
        index.insert(0, 'term_id')

    tmp = pd.pivot_table(data, index=index,
                         columns='sample_id',
                         values=['term_name', 'rank', 'p_value', 'z_score',
                                 'combined_score', 'adj_p_value', 'genes'],
                         aggfunc='first', fill_value=np.nan
                         )
    html_out = save_name

    print("Saving to : {}".format(html_out))

    html_tools.write_single_table(tmp, html_out, 'MAGINE GO analysis')
    html_out = save_name + '_filter'

    html_tools.write_filter_table(tmp, html_out, 'MAGINE GO analysis')


def create_gene_plots(data, list_of_terms, save_name, out_dir=None,
                      exp_data=None,
                      run_parallel=False, plot_type='plotly'):
    """ Creates a figure for each GO term in data

    Data should be a result of running calculate_enrichment.
    This function creates a plot of all proteins per term if a term is
    significant and the number of the reference set is larger than 5 and
    the total number of species measured is less than 100.


    Parameters
    ----------
    data : pandas.DataFrame
        previously ran enrichment analysis
    save_name : str
        name to save file
    out_dir : str
        output path for file
    exp_data : magine.ExperimentalData
        data to plot
    run_parallel : bool
        To run in parallel using pathos.multiprocessing
    plot_type : str
        plotly or matplotlib

    Returns
    -------
    out_array : dict
        dict where keys are pointers to figure locations
    """

    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if not os.path.exists(os.path.join(out_dir, 'Figures')):
            os.mkdir(os.path.join(out_dir, 'Figures'))
    data = data.copy()
    figure_locations = {}
    plots_to_create = []
    to_remove = set()
    assert plot_type == ('plotly' or 'matplotlib')
    # filter data by significance and number of references
    if len(list_of_terms) == 0:
        print("No significant GO terms!!!")
        return figure_locations, to_remove
    # here we are going to iterate through all sig GO terms and create
    # a list of plots to create. For the HTML side, we need to point to
    # a location

    # create plot of genes over time
    for n, i in enumerate(list_of_terms):
        local_exp_data = exp_data.data.copy()
        # want to plot all species over time
        index = data['term_name'] == i

        name = data[index]['term_name'].unique()

        if len(name) > 0:
            name = name[0]

        gene_set = set()
        genes = data[index]['genes']
        for g in genes:
            if isinstance(g, list):
                each = g
            else:
                each = g.split(',')

            gene_set = {j for j in each}

        if plot_type == 'matplotlib':
            # too many genes isn't helpful on plots, so skip them
            if len(gene_set) > 100:
                figure_locations[i] = '<a>{0}</a>'.format(name)
                continue
        local_save_name = 'Figures/{0}_{1}'.format(n, save_name)
        if out_dir is not None:
            local_save_name = '{0}/{1}'.format(out_dir, local_save_name)

        local_save_name = local_save_name.replace(':', '')
        out_point = '<a href="{0}.html">{1}</a>'.format(local_save_name, name)
        figure_locations[i] = out_point

        title = "{0} : {1}".format(str(i), name)
        local_df = local_exp_data[
            local_exp_data[gene].isin(list(gene_set))].copy()
        p_input = (local_df, list(gene_set), local_save_name, '.', title,
                   plot_type)

        plots_to_create.append(p_input)

    # return figure_locations, to_remove

    print("Starting to create plots for each GO term")
    # just keeping this code just in case using pathos is a bad idea
    # ultimately, using matplotlib is slow.

    if run_parallel:
        st2 = time.time()
        pool = mp.Pool()
        pool.map_async(plot_list_of_genes, plots_to_create)
        # pool.map(plot_list_of_genes, plots_to_create)
        pool.close()
        pool.join()
        end2 = time.time()
        print("parallel time = {}".format(end2 - st2))
        print("Done creating plots for each GO term")

    else:
        st1 = time.time()
        for i in plots_to_create:
            plot_list_of_genes(i)
        end1 = time.time()
        print("sequential time = {}".format(end1 - st1))

    return figure_locations, to_remove


if __name__ == '__main__':
    e = Enrichr()
    g_list = ['PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2',
              'P2RX7', 'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1',
              'ANAPC16', 'TMCC1', 'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2',
              'PLBD2', 'LARP7', 'TECPR2', 'ZNF302', 'CUX1', 'MOB2', 'CYTH2',
              'SEC22C', 'EIF4E3', 'ROBO2', 'ADAMTS9-AS2', 'CXXC1', 'LINC01314',
              'ATF7', 'ATP5F1']

    df = e.run(g_list, 'GO_Biological_Process_2017')
    lists = [['BAX', 'BCL2', 'CASP3'], ['CASP10', 'CASP8', 'BAK'],
             ['BIM', 'CASP3']]
    df2 = e.run_samples(lists, ['1', '2', '3'], save_name='test')
    # print(df.head(10))
