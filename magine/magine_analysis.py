import os
import warnings

import networkx as nx
import pandas as pd

from magine.data.formatter import pivot_table_for_export
from magine.html_templates.html_tools import write_single_table
from magine.networks.go_network_generator import GoNetworkGenerator
from magine.networks.visualization.cytoscape_view import RenderModel
from magine.ontology.ontology_analysis import GoAnalysis
from magine.plotting.species_plotting import write_table_to_html_with_figures


class Analyzer(object):
    """ MAGINE analyzer

    """
    def __init__(self, experimental_data, network=None, species='hsa',
                 output_directory='tmp', save_name='tmp', build_network=False):
        """Class to perform entire MAGINE pipeline.

        1. Creates a network
        2. Performs gene enrichment analysis
        3. Combines the two
            .. note::
                Requires cytoscape session to be opened if you want the output networks

        Parameters
        ----------
        experimental_data : magine.data.ExperimentalData
        network : nx.DiGraph
        species : str
        output_directory : str
        save_name : str
        build_network : bool
            Build a network based on experimental_data

        """
        self.exp_data = experimental_data
        self.species = species
        self.save_name = save_name
        self.out_dir = output_directory

        self.go = GoAnalysis(species=self.species,
                             output_directory=self.out_dir,
                             experimental_data=self.exp_data)

        self.network = None
        self.go_net_gen = None
        if build_network:
            if network is not None:
                warnings.warn("Warning : Passing build_network=True "
                              "and a network file! ", RuntimeWarning)
            self.generate_network(save_name)
            self.go_net_gen = GoNetworkGenerator(self.species, self.network)
        elif network is not None:
            self.network = network
            self.go_net_gen = GoNetworkGenerator(self.species, self.network)

        self.html_names = []
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        if not os.path.exists(os.path.join(self.out_dir, 'Figures')):
            os.mkdir(os.path.join(self.out_dir, 'Figures'))
        if not os.path.exists(os.path.join(self.out_dir, 'Network_files')):
            os.mkdir(os.path.join(self.out_dir, 'Network_files'))

    def generate_network(self, save_name, use_hmdb=True, use_reactome=True):
        """
        Generates network with default parameters


        Parameters
        ----------
        save_name : str
            name of network to be saved
        use_hmdb : bool
            Use HMDB to create network
        use_reactome : bool
            Use Reactome to build network


        Returns
        -------

        """
        from magine.networks.network_generator import build_network
        proteins = self.exp_data.list_sig_proteins
        metabolites = None
        if len(self.exp_data.list_metabolites) != 0:
            metabolites = self.exp_data.list_metabolites
        self.network = build_network(
            proteins, save_name=save_name, species=self.species,
            metabolite_list=metabolites,
            all_measured_list=self.exp_data.list_species,
            use_reactome=use_reactome, use_hmdb=use_hmdb)

    def run_single_go(self, data_type='proteomics', fold_change='up',
                      save=True):
        """
        Runs enrichment analysis


        Parameters
        ----------
        data_type : {'proteomics','rnaseq'}
            data type to run enrichment analysis
        fold_change: {'up', 'down', 'both'}
        save: bool
            saves enrichment array to file


        Returns
        -------
        dict

        """

        if data_type == 'proteomics':
            if fold_change == 'up':
                data = self.exp_data.proteomics_up_over_time
            elif fold_change == 'down':
                data = self.exp_data.proteomics_down_over_time
            elif fold_change == 'both':
                data = self.exp_data.proteomics_over_time
            labels = list(self.exp_data.proteomics_time_points)
        elif data_type == 'rnaseq':
            if fold_change == 'up':
                data = self.exp_data.rna_up_over_time
            elif fold_change == 'down':
                data = self.exp_data.rna_down_over_time
            elif fold_change == 'both':
                data = self.exp_data.rna_over_time
            labels = list(self.exp_data.rna_time_points)
        else:
            warnings.warn("Must provide proteomics or rnaseq for GO analysis")
            quit()

        print("Creating enrichment array")
        enrich_array = self.go.calculate_enrichment(data, labels=labels)
        # enrich_array = None
        save_name = '_'.join([self.save_name, data_type, fold_change])
        filename = '<a href="{}.html">link</a>'.format(
                self.out_dir + '/' + save_name)
        html_dict = {'DataType':   data_type,
                     'FoldChange': fold_change,
                     'fileName':   filename,
                     'save_name':  save_name
                     }
        if save:
            print("Saving enrichment array to {}".format(save_name))
            enrich_array.to_csv('{0}/{1}_all_data.csv'.format(self.out_dir,
                                                              save_name),
                                index=False)
            pt = pivot_table_for_export(enrich_array)
            pt.to_csv('{0}/{1}_all_data_pivot.csv'.format(self.out_dir,
                                                          save_name))
        return enrich_array, html_dict

    def run_proteomics_go(self):
        """ Updated function to run all data

        Returns
        -------

        """

        list_of_go_dict = []

        for i in ['up', 'down', 'both']:
            # for i in ['up']:
            e_array, att_dict = self.run_single_go(data_type='proteomics',
                                                   fold_change=i)
            list_of_go_dict.append(att_dict)

        if 'rnaseq' in self.exp_data.exp_methods:
            for i in ['up', 'down', 'both']:
                e_array, att_dict = self.run_single_go(data_type='rnaseq',
                                                       fold_change=i)
                list_of_go_dict.append(att_dict)

        return list_of_go_dict

    def run_go_and_create_html(self, html_name, visualize=False,
                               create_figure=True, run_parallel=True):
        """
        Runs GO analysis and creates an HTML output

        Parameters
        ----------
        html_name : str
        visualize : bool
            Visualize with cytoscape, cytoscape must be open!
        create_figure : bool
            Create plots for each GO
        run_parallel : bool
            Create plots in parallel

        Returns
        -------

        """

        list_of_go_dict = self.run_proteomics_go()
        for i in list_of_go_dict:

            file_name = '{}/{}_all_data.csv'.format(self.out_dir,
                                                    i['save_name'])
            save_name = "{}_{}_{}".format(self.save_name, i['DataType'],
                                          i['FoldChange'])
            print("Processing : {}".format(save_name))
            if create_figure:
                write_table_to_html_with_figures(file_name, self.exp_data,
                                                 save_name,
                                                 out_dir=self.out_dir,
                                                 run_parallel=run_parallel
                                                 )
            if visualize:
                self.create_selected_go_network(file_name, save_name,
                                                visualize=True)

        df = pd.DataFrame(list_of_go_dict)
        df.drop('save_name', axis=1, inplace=True)
        write_single_table(df, 'MAGINE Output', html_name)
        df['fileName'] = df['fileName'].str.replace('\.html', '_filter.html')
        write_single_table(df, 'MAGINE Output', html_name + '_filter')

    def create_selected_go_network(self, file_name, save_name, go_ids=None,
                                   visualize=False, slim=True):
        """ creates a GO level network

        Parameters
        ----------
        file_name: str or pandas.Dataframe
            Enrichment array
        save_name: str
            prefix of GO network to save
        go_ids: list, optional
            list of GO terms to include
        visualize: bool
            Visualize the network after creation
        slim : bool
            Slim the ontology

        Returns
        --------
        tall: networkx.DiGraph
            GO network created with GoNetworkGenerator
        data : pd.Dataframe
            enrichment data from GOAnalysis.calculate_enrichment
        """

        if self.network is None:
            warnings.warn("Warning : A molecular network is required to "
                          "generate a GO level network!.\n"
                          "You must provide Analyzer with a network "
                          "or pass the build_network flag", RuntimeWarning)
            quit()
        elif self.go_net_gen is None:
            self.go_net_gen = GoNetworkGenerator(self.species, self.network)
        if slim:
            from magine.ontology.ontology_tools import filter_ontology_df

            data = filter_ontology_df(file_name, n_hits_per_time=5,
                                      trim_nodes=True,
                                      additional_ids_to_include=go_ids)
        else:
            data = pd.read_csv(file_name)
            data = data[data['GO_id'].isin(go_ids)]
        list_all_go = data['GO_id'].unique()

        tall = self.go_net_gen.create_network_from_list(
            list_of_go_terms=list_all_go,
            save_name='{0}_network'.format(save_name),
            threshold=0)
        if visualize:
            self.visualize_go_network(tall, data, save_name)
        return tall, data

    def visualize_go_network(self, go_network, data, save_name,
                             format_only=False):
        """ Renders GO network with py2cytoscape

        Parameters
        ----------
        go_network : networkx.DiGraph
            GO network created with GoNetworkGenerator
        data : pd.Dataframe
            enrichment data from GOAnalysis.calculate_enrichment
        save_name: str
            prefix to save images of GO network
        format_only: boolean
            option to return formatted network and labels for rendering

        Returns
        -------

        """

        if len(go_network.nodes()) == 0:
            print('No nodes')
            quit()

        labels = data['sample_index'].unique()

        score_array = pd.pivot_table(data, index=['GO_id'],
                                     columns='sample_index')

        x = score_array['enrichment_score'].fillna(0)

        for i in go_network.nodes():
            tmp = go_network.node[i]['go']
            tmp = tmp[:2] + ':' + tmp[2:]
            values = x.loc[tmp]
            go_network.node[i]['color'] = 'red'
            for n, time in enumerate(labels):
                t = 'time_{0:04d}'.format(n)
                go_network.node[i][t] = float(values[time])

        if format_only:
            return go_network, data

        savename = os.path.join(self.out_dir, 'Network_files',
                                '{0}_all_colored.graphml'.format(save_name))

        nx.nx.write_graphml(go_network, savename)
        size_of_data = len(labels)
        rm = RenderModel(go_network, layout='force-directed')
        rm.visualize_by_list_of_time(_create_names(size_of_data),
                                     prefix=save_name,
                                     labels=self.exp_data.timepoints,
                                     out_dir=self.out_dir)


def _create_names(n):
    names = []
    for i in range(n):
        names.append('time_{0:04d}'.format(i))
        print('time_{0:04d}'.format(i))
    return names


