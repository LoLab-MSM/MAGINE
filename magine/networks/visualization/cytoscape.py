import numbers
import os
import time

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests

from py2cytoscape.cyrest.layout import layout as cy_layout
from py2cytoscape.data.cyrest_client import CyRestClient

IP = 'localhost'
PORT = 1234
VERSION = 'v1'

SUID_LIST = 'suid'
url = 'http://{}:{}/{}/'.format(IP, PORT, VERSION)


def create_point(value, lesser, equal, greater):
    return [{'value': str(value), 'lesser': lesser,
             'equal': equal, 'greater': greater}]


def create_slope(min_val=0, max_val=10, values=(1, 10)):
    point_1 = create_point(min_val, values[0], values[0], values[0])
    point_2 = create_point(max_val, values[1], values[1], values[1])
    return point_1 + point_2


class LayoutClient(object):
    def __init__(self):
        self.__url = url + 'apply/layouts'

    def get_all(self):
        return requests.get(self.__url).json()

    def get_option(self, name):
        _url = '{}/{}/columntypes'.format(self.__url, name)
        return requests.get(_url).json()

    def get_parameters(self, name):
        _url = '{}/{}/parameters'.format(self.__url, name)
        return requests.get(_url).json()

    def update(self, name, parameters='color'):
        _url = '{}/{}/{}'.format(self.__url, name, str(parameters))
        requests.get(_url)

    def apply(self, name='force-directed', network=None, params=None):
        self._check_net(network)
        _url = '{}/{}/{}'.format(self.__url, name, network.get_id())
        requests.get(_url, params)

    def bundle_edge(self, network=None):
        self._check_net(network)
        _url = '{}/edgebundling/{}'.format(self.__url, network.get_id())
        requests.get(_url)

    def fit(self, network):
        _url = self.__url + 'apply/fit/' + str(network.get_id())
        requests.get(_url)

    @staticmethod
    def _check_net(network):
        if network is None:
            raise ValueError('Target network is required')


class RenderModel(object):
    def __init__(self, graph, layout="attributes-layout", use_cyto_net=False,
                 clear_cytoscape=False):
        """

        Parameters
        ----------
        graph : nx.DiGraph
        layout : str
        use_cyto_net : bool
            Use already loaded network in cytoscape. Will use first network.
        clear_cytoscape : bool
            Clear cytoscape session. Useful when multiple instance of
            RenderModel.
        """

        self.graph = graph
        self.cy = CyRestClient()
        if clear_cytoscape:
            self.cy.session.delete()
        self.cy.layout2 = LayoutClient()

        for n, (i, j) in enumerate(self.graph.edges()):
            self.graph[i][j]['name'] = '{},{}'.format(i, j)
        if use_cyto_net:
            self.g_cy = self.cy.network.create(self.cy.network.get_all()[0])
        else:
            self.g_cy = self.cy.network.create_from_networkx(self.graph)

        # time.sleep(5)
        self.layout = cy_layout(url)
        view_id_list = self.g_cy.get_views()
        self.view1 = self.g_cy.get_view(view_id_list[0], fmt='view')

        # Marquee Directed Simple
        self.style = self.cy.style.create('MAGINE')
        self.style.update_defaults(default_style)
        # self.style.create_passthrough_mapping('name', vp='NODE_LABEL')
        if layout == 'attributes-layout':
            self.cy.layout2.apply(name=layout, network=self.g_cy,
                                  params={'column': 'color'})
        else:
            self.cy.layout2.apply(
                name=layout, network=self.g_cy,
                params={
                    'defaultSpringLength': 50,
                    'defaultNodeMass': 5,
                }
            )
        self.node_name2id = name2suid(self.g_cy, 'node')
        self.edge_name2id = name2suid(self.g_cy, 'edge')

    def apply_layout(self, layout, params):

        self.cy.layout2.apply(
            name=layout, network=self.g_cy,
            params=params
        )

    def print_options(self):
        """ print cytoscape options for network, style, nodes, edges

        :return:
        """
        vps = pd.Series(self.cy.style.vps.get_all())
        print(vps)
        node_vps = self.cy.style.vps.get_node_visual_props()
        for i in node_vps:
            print(i)
        print('\nNode options\n')
        edge_vps = self.cy.style.vps.get_edge_visual_props()
        for i in edge_vps:
            print(i)
        print('\nEdge options\n')
        network_vps = self.cy.style.vps.get_network_visual_props()
        for i in network_vps:
            print(i)
        style_opts = self.cy.style.get_all()
        print('\nStyle options\n')
        for i in style_opts:
            print(i)

    def map_edge_width(self, column='weight', min_value=None, max_value=None,
                       min_width=3, max_width=10):

        if min_value is None and max_value is None:
            edge_width = {}
            for i, j, d in self.graph.edges(data=True):
                edge_width[self.edge_name2id['{},{}'.format(i, j)]] = d[column]

            _min, _max = min(edge_width.values()), max(edge_width.values())
        elif [isinstance(x, numbers.Number) for x in [min_value, max_value]]:
            raise Exception("Must provide min_value and max_value to "
                            "create continuous mapping of edge width")
        else:
            _min, _max = min_value, max_value

        slope = create_slope(
            min_val=_min,
            max_val=_max,
            values=(min_width, max_width)
        )

        self.style.create_continuous_mapping(
            column=column, col_type='Double', vp='EDGE_WIDTH',
            points=slope
        )

    def visualize_by_list_of_time(self, list_of_time, labels=None,
                                  prefix='tmp', out_dir=None, scale=100):
        """
        create sequences of pdfs and svgs based on list of attributes

        list_of_time should point to attributes of the network.
        This attribute will update the network accordingly.

        Parameters
        ----------
        list_of_time : list_like
        labels : list_like
        prefix : str
        out_dir : str

        """
        if out_dir is not None:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            if not os.path.exists(os.path.join(out_dir, 'Figures')):
                os.mkdir(os.path.join(out_dir, 'Figures'))

        self.map_edge_width(column='weight', min_value=None, max_value=None,
                            min_width=3, max_width=10)

        self.style.create_passthrough_mapping('name', vp='NODE_LABEL')
        self.cy.style.apply(style=self.style, network=self.g_cy)

        node_label_colors = {self.node_name2id[i]: 'black' for i in
                             self.graph.nodes}

        node_color_values = {self.node_name2id[i]: d['color'] for i, d in
                             self.graph.nodes(data=True)}
        # do all node changes
        df_vs_node = pd.DataFrame(
            [node_label_colors, node_color_values],
            index=['NODE_LABEL_COLOR', 'NODE_FILL_COLOR'],
        )

        df_vs_node = df_vs_node.T
        df_vs_node['NODE_BORDER_PAINT'] = df_vs_node['NODE_LABEL_COLOR']

        # do all edge changes

        df_edge = pd.DataFrame(index=self.edge_name2id.values())
        # df_edge = df_edge.T
        df_edge['EDGE_STROKE_UNSELECTED_PAINT'] = 'gray'
        df_edge['EDGE_SOURCE_ARROW_UNSELECTED_PAINT'] = 'gray'
        df_edge['EDGE_TARGET_ARROW_UNSELECTED_PAINT'] = 'gray'

        # update views
        self.view1.batch_update_node_views(df_vs_node)
        self.view1.batch_update_edge_views(df_edge)

        self.cy.layout.fit(network=self.g_cy)
        x = self.view1.get_network_view_as_dict()
        self.view1.update_network_view(
            visual_property='NETWORK_BACKGROUND_PAINT', value=u'#FFFFFF')
        self.view1.update_network_view(
            visual_property='NETWORK_SCALE_FACTOR',
            value=0.8 * x['NETWORK_SCALE_FACTOR']
        )

        x = self.view1.get_network_view_as_dict()
        self.create_png('out.png', 2400)
        trim_photo('out.png', 'x')

        all_node_size = []
        for ind, j in enumerate(list_of_time):
            size = np.array([self.graph.node[n]['sample{}'.format(j)]
                             for n in self.graph.nodes])
            all_node_size.append(size)
        all_node_size = np.array(all_node_size).flatten()
        simple_slope = create_slope(min_val=all_node_size.min(),
                                    max_val=all_node_size.max(),
                                    values=(10, scale))

        for j in sorted(list_of_time):
            time.sleep(1)
            self.style.create_continuous_mapping(column='sample{}'.format(j),
                                                 col_type='Double',
                                                 vp='NODE_SIZE',
                                                 points=simple_slope)
            self.cy.style.apply(style=self.style, network=self.g_cy)

            fig_name = '{0}_{1}'.format(prefix, j)
            print("Saving {}".format(fig_name))

            if out_dir is not None:
                out_file = os.path.join(out_dir, 'Figures',
                                        '{}.png'.format(fig_name))
            else:
                out_file = '{}.png'.format(fig_name)

            if os.path.exists(out_file):
                os.remove(out_file)

            x = self.view1.get_network_view_as_dict()
            self.create_png(out_file, int(x['NETWORK_WIDTH']))
            self.create_svg(out_file, int(x['NETWORK_WIDTH']))
            if labels is None:
                trim_photo(out_file, j)
            else:
                trim_photo(out_file, j)

    def update_by_list_of_nodes(self, list_of_time, prefix='tmp', scale=100):
        """
        create sequences of pdfs and svgs based on list of attributes

        list_of_time should point to attributes of the network.
        This attribute will update the network accordingly.

        Parameters
        ----------
        list_of_time : list_like
        prefix : str


        """
        self.map_edge_width(column='weight', min_value=None, max_value=None,
                            min_width=3, max_width=10)

        self.style.create_passthrough_mapping('name', vp='NODE_LABEL')
        self.cy.style.apply(style=self.style, network=self.g_cy)

        node_label_colors = {self.node_name2id[i]: 'black' for i in
                             self.graph.nodes}

        node_color_values = {self.node_name2id[i]: d['color'] for i, d in
                             self.graph.nodes(data=True)}
        # do all node changes
        df_vs_node = pd.DataFrame(
            [node_label_colors, node_color_values],
            index=['NODE_LABEL_COLOR', 'NODE_FILL_COLOR'],
        ).T
        df_vs_node['NODE_BORDER_PAINT'] = df_vs_node['NODE_LABEL_COLOR']
        self.view1.batch_update_node_views(df_vs_node)
        self.pretty_layout()

        all_node_size = []
        for ind, j in enumerate(list_of_time):
            size = np.array([self.graph.node[n]['sample{}'.format(j)]
                             for n in self.graph.nodes])
            all_node_size.append(size)
        all_node_size = np.array(all_node_size).flatten()
        simple_slope = create_slope(min_val=all_node_size.min(),
                                    max_val=all_node_size.max(),
                                    values=(10, scale))
        # Continuous mapping
        points = [
            {
                'value': str(all_node_size.min()),
                'lesser': 'blue',
                'equal': 'blue',
                'greater': 'blue'
            },
            {
                'value': str(all_node_size.max()),
                'lesser': 'red',
                'equal': 'red',
                'greater': 'red'
            }
        ]
        for sample in sorted(list_of_time):
            time.sleep(1)
            self.style.create_continuous_mapping(
                column='sample{}'.format(sample),
                col_type='Double',
                vp='NODE_SIZE',
                points=simple_slope
            )
            self.style.create_continuous_mapping(
                column='sample{}'.format(sample),
                col_type='Double',
                vp='NODE_FILL_COLOR',
                points=points
            )

            self.cy.style.apply(style=self.style, network=self.g_cy)

            fig_name = '{0}_{1}'.format(prefix, sample)
            print("Saving {}".format(fig_name))

            out_file = '{}.png'.format(fig_name)

            if os.path.exists(out_file):
                os.remove(out_file)

            x = self.view1.get_network_view_as_dict()
            self.create_png(out_file, int(x['NETWORK_WIDTH']))
            self.create_svg(out_file, int(x['NETWORK_WIDTH']))
            trim_photo(out_file, j)

    def pretty_layout(self):
        df_edge = pd.DataFrame(index=self.edge_name2id.values())
        # df_edge = df_edge.T
        df_edge['EDGE_STROKE_UNSELECTED_PAINT'] = 'gray'
        df_edge['EDGE_SOURCE_ARROW_UNSELECTED_PAINT'] = 'black'
        df_edge['EDGE_TARGET_ARROW_UNSELECTED_PAINT'] = 'gray'

        self.view1.batch_update_edge_views(df_edge)

        self.fit()
        x = self.view1.get_network_view_as_dict()
        self.view1.update_network_view(
            visual_property='NETWORK_BACKGROUND_PAINT', value=u'#FFFFFF')
        self.view1.update_network_view(
            visual_property='NETWORK_SCALE_FACTOR',
            value=0.8 * x['NETWORK_SCALE_FACTOR']
        )

    def update_node_edge_by_samples(self, sample_list, labels=None,
                                    prefix='tmp', out_dir=None, scale=50):
        """
        Update node size and edge width by attributes in network

        Created for the magine.networks.annotated_set.create_agn
        It generates images for each sample_list

        Parameters
        ----------
        sample_list : list
        labels : list
        prefix : str
            Output prefix
        out_dir : str
            Directory to place output images
        scale : int
            Max node size for scaling.

        """
        if out_dir is not None:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

        all_node_size = []
        for ind, j in enumerate(sample_list):
            size = np.array([self.graph.node[n]['sample{}'.format(j)]
                             for n in self.graph.nodes])
            all_node_size.append(size)
        node_size = np.array(all_node_size).flatten()
        node_slope = create_slope(min_val=node_size.min(),
                                  max_val=node_size.max(),
                                  values=(10, scale))

        edge_width = []
        for label in sample_list:
            for i, j, d in self.graph.edges(data=True):
                edge_width.append(d['weight{}'.format(label)])

        edge_width = np.array(edge_width).flatten()
        edge_slope = create_slope(min_val=edge_width.min(),
                                  max_val=edge_width.max(),
                                  values=(1, 10))

        self.style.create_passthrough_mapping('name', vp='NODE_LABEL')
        # self.apply_layout('circular', {})
        # self.layout.circular(nodeHorizontalSpacing=50, nodeVerticalSpacing=50)
        self.layout.attribute_circle(NodeAttribute='term', spacing='50')
        self.pretty_layout()

        for sample in sorted(sample_list):
            self.style.create_continuous_mapping(
                column='weight{}'.format(sample),
                col_type='Double',
                vp='EDGE_WIDTH',
                points=edge_slope
            )

            self.style.create_continuous_mapping(
                column='sample{}'.format(sample),
                col_type='Double',
                vp='NODE_SIZE',
                points=node_slope
            )

            self.cy.style.apply(style=self.style, network=self.g_cy)

            fig_name = '{0}_{1}'.format(prefix, sample)
            print("Saving {}".format(fig_name))

            if out_dir is not None:
                fig_name = os.path.join(out_dir, fig_name)

            x = self.view1.get_network_view_as_dict()
            self.create_png(fig_name, int(x['NETWORK_WIDTH']))
            self.create_svg(fig_name, int(x['NETWORK_WIDTH']))
            if labels is None:
                trim_photo('{}.png'.format(fig_name), sample)
            else:
                trim_photo('{}.png'.format(fig_name), sample)

    def update_all_node_attribute(self, nodes, attribute, value):

        # Update only valid nodes
        valid_nodes = {i for i in nodes if i in self.graph.nodes()}

        # assume correct datatype
        node_size = {self.node_name2id[i]: value for i in valid_nodes}
        # # do all node changes
        df_vs_node = pd.DataFrame(
            [node_size, node_size],
            index=[attribute]
        ).T
        self.view1.batch_update_node_views(df_vs_node)

    def update_node_color(self, nodes, color, reset_style=True):
        # resets node colors to default
        if reset_style:
            self.reset_style()

        valid_nodes = {i for i in nodes if i in self.graph.nodes()}

        node_color_values = {self.node_name2id[i]: color for i in valid_nodes}

        # # do all node changes
        df_vs_node = pd.DataFrame(
            [node_color_values, ],
            index=['NODE_FILL_COLOR']
        ).T
        self.view1.batch_update_node_views(df_vs_node)

    def update_node_size(self, nodes, size):
        """ Change list of nodes size to provided integer

        Parameters
        ----------
        nodes : list
            List of nodes. Nodes outside of network will be ignored.
        size : int
            Size of nodes in cytoscape session

        Returns
        -------

        """
        valid_nodes = {i for i in nodes if i in self.graph.nodes()}

        node_size = {self.node_name2id[i]: int(size) for i in valid_nodes}
        # # do all node changes
        df_vs_node = pd.DataFrame(
            [node_size, node_size],
            index=['NODE_HEIGHT', 'NODE_WIDTH']
        ).T
        self.view1.batch_update_node_views(df_vs_node)

    def reset_style(self):
        self.cy.style.apply(style=self.style, network=self.g_cy)

    def create_png(self, out_file, width):
        if not out_file.endswith('png'):
            out_file = f'{out_file}.png'
        with open(out_file, 'wb') as f:
            f.write(self.g_cy.get_png(width))

    def create_svg(self, out_file, width):
        if not out_file.endswith('svg'):
            out_file = f'{out_file}.svg'
        with open(out_file, 'wb') as f:
            f.write(self.get_svg(width))

    def get_svg(self, height=2000):
        url = '%sviews/first.svg?h=%d' % (self.g_cy._CyNetwork__url, height)
        return requests.get(url).content

    def fit(self):
        self.cy.layout.fit(network=self.g_cy)
        time.sleep(1)

    def get_png(self, height=2000):

        url = '%sviews/first.png?h=%d' % (self.g_cy._CyNetwork__url, height)
        return requests.get(url).content


def name2suid(network, obj_type='node'):
    if obj_type is 'node':
        table = network.get_node_table()
    elif obj_type is 'edge':
        table = network.get_edge_table()
    else:
        raise ValueError('No such object type: ' + obj_type)
    table.reset_index(level=0, inplace=True)
    name_to_suid = {}
    for suid, name in table[['SUID', 'name']].values:
        name_to_suid[name] = suid
    return name_to_suid


def trim_photo(im_location, title=None):
    """
    Removes whitespace and adds title to image


    Parameters
    ----------
    im_location: str
        location of file, will be used for output as well
    title : str
        title to provide to add to image
        Will be the sample id

    Returns
    -------

    """

    img = mpimg.imread(im_location)
    mask = img[:, :, 0] < 1

    img = img[np.ix_(mask.any(1), mask.any(0))]

    plt.imshow(img, interpolation='none')
    plt.xticks([])
    plt.yticks([])
    if title is not None:
        plt.title(title, fontsize=24)
    plt.axis('off')
    # out = im_location.replace('.png', '_formatted.png')
    plt.savefig(im_location, dpi=2000, bbox_inches='tight', transparent=True)
    plt.close()


def trip_photo_old(im_location, title=None):
    """
    Removes whitespace and adds title to image


    Parameters
    ----------
    im_location: str
        location of file, will be used for output as well
    title : str
        title to provide to add to image
        Will be the sample id

    Returns
    -------

    """

    img = mpimg.imread(im_location)
    x_dim = img.shape[0]
    y_dim = img.shape[1]
    alpha = np.ones((x_dim, y_dim, 1))
    mask = (img[:, :, 0] == 1.) & (img[:, :, 1] == 1.)

    alpha[mask] = 0
    img = np.append(img, alpha, axis=2)
    to_remove_x = set()
    to_remove_y = set()

    for i in range(x_dim):
        if img[i, :, 0].sum() != y_dim:
            to_remove_x.add(i)
    for i in range(y_dim):
        if img[:, i, 0].sum() != x_dim:
            to_remove_y.add(i)
    if len(to_remove_x) > 2:
        _min, _max = min(to_remove_x), max(to_remove_x)
        _scale = int((.05 * x_dim))
        img = img[_min - _scale:_max + _scale, :, :]

    if len(to_remove_y) > 2:
        _min, _max = min(to_remove_y), max(to_remove_y)
        _scale = int((.05 * y_dim))
        img = img[:, _min - _scale:_max + _scale, :]

    plt.imshow(img, interpolation='none')
    plt.xticks([])
    plt.yticks([])
    if title is not None:
        plt.title(title)
    plt.axis('off')
    out = im_location.replace('.png', '_formatted.png')
    out2 = im_location.replace('.png', '_formatted.svg')
    plt.savefig(out, dpi=300, bbox_inches='tight', transparent=True)
    plt.savefig(out2, bbox_inches='tight', transparent=True)
    plt.close()
    # plt.show()


default_style = {
    u'EDGE_BEND'                         : None,
    u'EDGE_CURVED'                       : True,
    u'EDGE_LABEL'                        : u'',
    u'EDGE_LABEL_COLOR'                  : u'#000000',
    u'EDGE_LABEL_FONT_SIZE'              : 10,
    u'EDGE_LABEL_TRANSPARENCY'           : 255,
    u'EDGE_LABEL_WIDTH'                  : 200.0,
    u'EDGE_LINE_TYPE'                    : u'SOLID',
    u'EDGE_PAINT'                        : u'#808080',
    u'EDGE_SELECTED'                     : False,
    u'EDGE_SELECTED_PAINT'               : u'#FF0000',
    u'EDGE_SOURCE_ARROW_SELECTED_PAINT'  : u'#FFFF00',
    u'EDGE_SOURCE_ARROW_SHAPE'           : u'NONE',
    u'EDGE_SOURCE_ARROW_SIZE'            : 6.0,
    u'EDGE_SOURCE_ARROW_UNSELECTED_PAINT': u'#000000',
    u'EDGE_STROKE_SELECTED_PAINT'        : u'#FF0000',
    u'EDGE_STROKE_UNSELECTED_PAINT'      : u'#404040',
    u'EDGE_TARGET_ARROW_SELECTED_PAINT'  : u'#FFFF00',
    u'EDGE_TARGET_ARROW_SHAPE'           : u'Delta',
    u'EDGE_TARGET_ARROW_SIZE'            : 6.0,
    u'EDGE_TARGET_ARROW_UNSELECTED_PAINT': u'#000000',
    u'EDGE_TOOLTIP'                      : u'',
    u'EDGE_TRANSPARENCY'                 : 150,
    u'EDGE_UNSELECTED_PAINT'             : u'#404040',
    u'EDGE_VISIBLE'                      : True,
    u'EDGE_WIDTH'                        : 2.0,
    u'NETWORK_BACKGROUND_PAINT'          : u'#FFFFFF',
    u'NETWORK_CENTER_X_LOCATION'         : 0.0,
    u'NETWORK_CENTER_Y_LOCATION'         : 0.0,
    u'NETWORK_CENTER_Z_LOCATION': 0.0,
    u'NETWORK_DEPTH': 0.0,
    u'NETWORK_EDGE_SELECTION': True,
    u'NETWORK_HEIGHT': 400.0,
    u'NETWORK_NODE_SELECTION': True,
    u'NETWORK_SCALE_FACTOR': 1.0,
    u'NETWORK_SIZE': 550.0,
    u'NETWORK_TITLE': u'',
    u'NETWORK_WIDTH': 550.0,

    u'NODE_BORDER_PAINT': 'black',
    u'NODE_BORDER_STROKE': u'SOLID',
    u'NODE_BORDER_TRANSPARENCY': 255,
    u'NODE_BORDER_WIDTH': 2.0,
    u'NODE_DEPTH': 0.0,
    u'NODE_FILL_COLOR': 'red',
    u'NODE_HEIGHT': 40.0,
    # u'NODE_LABEL': ' ',
    # u'NODE_LABEL': u'<passthroughMapping attributeName="name" ',
    u'NODE_LABEL_COLOR': u'black',
    u'NODE_LABEL_FONT_FACE': u'SansSerif.plain,plain,12',
    u'NODE_LABEL_FONT_SIZE': 12,
    u'NODE_LABEL_POSITION': u'C,C,c,0.50,0.00',
    u'NODE_LABEL_TRANSPARENCY': 255,
    u'NODE_LABEL_WIDTH': 200.0,
    u'NODE_NESTED_NETWORK_IMAGE_VISIBLE': True,
    u'NODE_PAINT': u'#787878',
    u'NODE_SELECTED': False,
    u'NODE_SELECTED_PAINT': u'#FFFF00',
    u'NODE_SHAPE': u'ELLIPSE',
    u'NODE_SIZE': 50.0,
    u'NODE_TOOLTIP': u'',
    u'NODE_TRANSPARENCY': 255,
    u'NODE_VISIBLE'                      : True,
    u'NODE_WIDTH'                        : 60.0,
    u'NODE_X_LOCATION'                   : 0.0,
    u'NODE_Y_LOCATION'                   : 0.0,
    u'NODE_Z_LOCATION'                   : 0.0
}

# default_style = {}
if __name__ == '__main__':
    # ddn = nx.nx.read_graphml('t_all_colored_pvalue_2.graphml')
    ddn = nx.DiGraph()
    ddn.add_node('x')
    ddn.add_edge('x', 'y')
    rm = RenderModel(ddn, style='Marquee')
    rm.print_options()
    quit()
    rm.visualize_by_list_of_time(
        ['time_0', 'time_1', 'time_2', 'time_3', 'time_4', ])
