import os
import time

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.util_network import NetworkUtil as util

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
        requests.put(_url)

    def apply(self, name='force-directed', network=None, params=None):
        self._check_net(network)
        _url = '{}/{}/{}'.format(self.__url, name, network.get_id())
        requests.get(_url, params)

    def bundle_edge(self, network=None):
        self._check_net(network)
        _url = '{}/edgebundling/{}'.format(self.__url, network.get_id())
        requests.get(_url)

    def _check_net(self, network):
        if network is None:
            raise ValueError('Target network is required')


class RenderModel(object):
    def __init__(self, graph, layout="attributes-layout", style='Directed'):
        # name='FromMAGINE'):
        # force-directed
        # attributes-layout

        self.graph = graph
        self.cy = CyRestClient()
        self.cy.session.delete()
        self.cy.layout2 = LayoutClient()
        self.style = None
        self.edge_name2id = None
        self.node_name2id = None
        self.g_cy = None
        self.view1 = None
        for n, (i, j) in enumerate(self.graph.edges()):
            self.graph[i][j]['name'] = '{},{}'.format(i, j)
        self.g_cy = self.cy.network.create_from_networkx(self.graph)

        time.sleep(2)
        view_id_list = self.g_cy.get_views()
        self.view1 = self.g_cy.get_view(view_id_list[0], format='view')

        # Marquee Directed Simple
        self.style = self.cy.style.create(style)
        # self.style = self.cy.style.create('Marquee')

        options = {
            'NODE_LABEL_FONT_SIZE': 24,
            'EDGE_WIDTH': 2,
            'EDGE_TRANSPARENCY': '150',
            # 'NETWORK_HEIGHT'           : '2800',
            # 'NETWORK_WIDTH'            : '2800',
            'NODE_LABEL_COLOR': 'black',
            'NODE_FILL_COLOR': 'red',
            'NETWORK_BACKGROUND_PAINT': '#00FFFFFF',
            'NODE_SIZE': 80,
            # 'NODE_LABEL_POSITION': 'C,C,c,0.00,-60.00',
            'NETWORK_CENTER_X_LOCATION': 0.0,
            'NETWORK_CENTER_Y_LOCATION': 0.0,
        }

        self.style.update_defaults(options)
        if layout == 'attributes-layout':
            self.cy.layout2.apply(name=layout, network=self.g_cy,
                                  params={'column': 'color'})
        else:
            self.cy.layout2.apply(name=layout, network=self.g_cy)
        self.node_name2id = util.name2suid(self.g_cy, 'node')
        self.edge_name2id = util.name2suid(self.g_cy, 'edge')

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

    def visualize_by_list_of_time(self, list_of_time, labels=None,
                                  prefix='tmp', out_dir=None):
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


        Returns
        -------

        """
        if out_dir is not None:
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            if not os.path.exists(os.path.join(out_dir, 'Figures')):
                os.mkdir(os.path.join(out_dir, 'Figures'))

        node_label_values = {self.node_name2id[i]: d['label'] for i, d in
                             self.graph.nodes(data=True)}
        node_label_colors = {self.node_name2id[i]: 'black' for i in
                             self.graph.nodes}
        node_color_values = {self.node_name2id[i]: d['color'] for i, d in
                             self.graph.nodes(data=True)}

        edge_color_values = {}
        for i in self.edge_name2id:
            edge_color_values[self.edge_name2id[i]] = 'gray'
        edge_width = {}
        for i, j, d in self.graph.edges(data=True):
            edge_width[self.edge_name2id['{},{}'.format(i, j)]] = d['weight']

        _min, _max = min(edge_width.values()), max(edge_width.values())
        slope = create_slope(min_val=_min, max_val=_max, values=(3, 10))

        self.style.create_continuous_mapping(column='weight',
                                             col_type='Double',
                                             vp='EDGE_WIDTH',
                                             points=slope)
        self.cy.style.apply(style=self.style, network=self.g_cy)
        self.cy.layout.fit(network=self.g_cy)
        self.create_png('out.png', 2400)
        trip_photo('out.png', 'x')
        x = self.view1.get_network_view_as_dict()
        self.view1.update_network_view(
            visual_property='NETWORK_SCALE_FACTOR',
            value=0.8 * x['NETWORK_SCALE_FACTOR'])

        all_node_size = []
        for ind, j in enumerate(list_of_time):
            size = np.array([self.graph.node[n]['sample{}'.format(j)]
                             for n in self.graph.nodes])
            all_node_size.append(size)
        size = np.array(all_node_size).flatten()
        simple_slope = create_slope(min_val=size.min(), max_val=size.max(),
                                    values=(10, 50))

        for j in list_of_time:

            self.style.create_continuous_mapping(column='sample{}'.format(j),
                                                 col_type='Double',
                                                 vp='NODE_SIZE',
                                                 points=simple_slope)
            self.cy.style.apply(style=self.style, network=self.g_cy)

            self.view1.update_node_views(visual_property='NODE_LABEL',
                                         values=node_label_values)
            self.view1.update_network_view(
                visual_property='NETWORK_BACKGROUND_PAINT',
                value='rgba(0, 0, 0, 0)', )
            self.view1.update_node_views(visual_property='NODE_LABEL_COLOR',
                                         values=node_label_colors)
            self.view1.update_node_views(visual_property='NODE_FILL_COLOR',
                                         values=node_color_values)
            self.view1.update_node_views(visual_property='NODE_BORDER_PAINT',
                                         values=node_color_values)
            self.view1.update_edge_views(visual_property='EDGE_LABEL_COLOR',
                                         values=edge_color_values)
            self.view1.update_edge_views(
                visual_property='EDGE_STROKE_UNSELECTED_PAINT',
                values=edge_color_values)
            self.view1.update_edge_views(
                visual_property='EDGE_SOURCE_ARROW_UNSELECTED_PAINT',
                values=edge_color_values)
            self.view1.update_edge_views(
                visual_property='EDGE_TARGET_ARROW_UNSELECTED_PAINT',
                values=edge_color_values)
            x = self.view1.get_network_view_as_dict()

            fig_name = 'ont_network_{0}_{1}'.format(prefix, j)
            print("Saving {}".format(fig_name))

            if out_dir is not None:
                out_file = os.path.join(out_dir, 'Figures',
                                        '{}.png'.format(fig_name))
            else:
                out_file = '{}.png'.format(fig_name)

            if os.path.exists(out_file):
                os.remove(out_file)

            self.create_png(out_file, int(x['NETWORK_WIDTH']))
            if labels is None:
                trip_photo(out_file, j)
            else:
                trip_photo(out_file, j)

    def update_node_color(self, attribute, save_name):
        self.cy.style.apply(style=self.style, network=self.g_cy)
        node_color_values = {self.node_name2id[i[0]]: i[1][attribute] for i in
                             self.graph.nodes(data=True)}
        self.view1.update_node_views(visual_property='NODE_FILL_COLOR',
                                     values=node_color_values)

        with open('{0}.svg'.format(save_name), 'wb') as f:
            network_pdf = self.g_cy.get_svg()
            f.write(network_pdf)
            f.close()

        with open('{0}.pdf'.format(save_name), 'wb') as f:
            network_pdf = self.g_cy.get_pdf()
            f.write(network_pdf)
            f.close()

        with open('{0}.png'.format(save_name), 'wb') as f:
            network_pdf = self.g_cy.get_png()
            f.write(network_pdf)
            f.close()

    def create_png(self, out_file, width):
        with open(out_file, 'wb') as f:
            network_svg = self.g_cy.get_png(width)
            f.write(network_svg)

    def get_svg(self, height=2000):
        url = '%sviews/first.svg?h=%d' % (self.g_cy._CyNetwork__url, height)
        return requests.get(url).content


def trip_photo(im_location, title=None):
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


if __name__ == '__main__':
    # ddn = nx.nx.read_graphml('t_all_colored_pvalue_2.graphml')
    ddn = nx.DiGraph()
    rm = RenderModel(ddn, style='Marquee')
    rm.visualize_by_list_of_time(
        ['time_0', 'time_1', 'time_2', 'time_3', 'time_4', ])
