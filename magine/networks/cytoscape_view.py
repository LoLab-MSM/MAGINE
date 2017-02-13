import os

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests

from magine.networks.cytoscape_mod_layout import LayoutClient
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
from py2cytoscape.data.util_network import NetworkUtil as util

class RenderModel:
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
        for n, i in enumerate(self.graph.edges()):
            self.graph[i[0]][i[1]]['name'] = str(i[0]) + ',' + str(i[1])
        self.g_cy = self.cy.network.create_from_networkx(
            self.graph, )  # name=name,
        # collection=name)
        import time
        time.sleep(2)
        view_id_list = self.g_cy.get_views()
        self.view1 = self.g_cy.get_view(view_id_list[0], format='view')

        params = [{u'type'       : u'double', u'name': u'spacingx',
                   u'value'      : 150.0,
                   u'description': u'Horizontal spacing between two partitions in a row:'},
                  {u'type'       : u'double', u'name': u'spacingy',
                   u'value'      : 150.0,
                   u'description': u'Vertical spacing between the largest partitions of two rows:'},
                  {u'name': u'maxwidth', u'value': 1000.0},
                  {u'type'       : u'double', u'name': u'minrad',
                   u'value'      : 20.0,
                   u'description': u'Minimum width of a partition:'},
                  {u'type'       : u'double', u'name': u'radmult',
                   u'value'      : 50.0,
                   u'description': u'Scale of the radius of the partition:'}]


        # Marquee Directed Simple
        self.style = self.cy.style.create(style)

        # self.style = self.cy.style.create('Marquee')

        options = {'NODE_LABEL_FONT_SIZE'     : 24,
                   'EDGE_WIDTH'               : 2,
                   'EDGE_TRANSPARENCY'        : '150',
                   # 'NETWORK_HEIGHT'           : '2800',
                   # 'NETWORK_WIDTH'            : '2800',
                   'NODE_FILL_COLOR': 'red',
                   'NETWORK_BACKGROUND_PAINT' : '#ffffff',
                   'NODE_SIZE' : 80,
                   # 'NODE_LABEL_POSITION': 'C,C,c,0.00,-60.00',
                   'NETWORK_CENTER_X_LOCATION': 0.0,
                   'NETWORK_CENTER_Y_LOCATION': 0.0,
                   }

        self.style.update_defaults(options)
        if layout == 'attributes-layout':
            self.cy.layout2.update(name=layout, parameters=params)
            self.cy.layout2.apply(name=layout, network=self.g_cy,
                                  params={'column': 'color'})
        else:
            self.cy.layout2.apply(name=layout, network=self.g_cy)
        self.node_name2id = util.name2suid(self.g_cy, 'node')
        self.edge_name2id = util.name2suid(self.g_cy, 'edge')
        # self.print_options()
        # quit()

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
                                  prefix='tmp',
                                  out_dir='tmp'):
        """ create sequences of pdfs and svgs based on list of attributes
        list_of_time should point to attributes of the network. This attribute will update the network accordingly.

        :param list_of_time:
        :return:
        """
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if not os.path.exists(os.path.join(out_dir, 'Figures')):
            os.mkdir(os.path.join(out_dir, 'Figures'))
        node_label_values = {self.node_name2id[i[0]]: i[1]['label'] for i in
                             self.graph.nodes(data=True)}
        node_color_values = {self.node_name2id[i[0]]: i[1]['color'] for i in
                             self.graph.nodes(data=True)}
        edge_color_values = {}
        for i in self.edge_name2id:
            edge_color_values[self.edge_name2id[i]] = 'gray'
        edge_width_values = {}
        for i in self.graph.edges(data=True):
            edge_width_values[self.edge_name2id[str(i[0]) + ',' + str(i[1])]] \
                = i[2]['weight']
        simple_slope = StyleUtil.create_slope(
                min=min(edge_width_values.values()),
                max=max(edge_width_values.values()),
                values=(3, 10))

        self.style.create_continuous_mapping(column='weight',
                                             col_type='Double',
                                             vp='EDGE_WIDTH',
                                             points=simple_slope)
        # self.style.create_passthrough_mapping(column='weight', col_type='Double', vp='EDGE_WIDTH')
        self.cy.style.apply(style=self.style, network=self.g_cy)
        self.cy.layout.fit(network=self.g_cy)
        self.create_png('out.png', 2400)
        trip_photo('out.png', 'x')
        x = self.view1.get_network_view_as_dict()
        # for i in x:
        #     print(i, x[i])
        self.view1.update_network_view(
            visual_property='NETWORK_SCALE_FACTOR',
            value=0.8 * x['NETWORK_SCALE_FACTOR'])

        # self.create_png('out.png', 2400)
        # trip_photo('out.png', 'x')

        for ind, j in enumerate(list_of_time):

            size = np.array(
                    [self.graph.node[n][j] for n in self.graph.nodes()])
            simple_slope = StyleUtil.create_slope(min=size.min(),
                                                  max=size.max(),
                                                  values=(10, 50))
            self.style.create_continuous_mapping(column=j, col_type='Double',
                                                 vp='NODE_SIZE',
                                                 points=simple_slope)
            self.cy.style.apply(style=self.style, network=self.g_cy)

            self.view1.update_node_views(visual_property='NODE_LABEL',
                                         values=node_label_values)
            self.view1.update_node_views(visual_property='NODE_LABEL_COLOR',
                                         values=node_label_values)
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

            fig_name = 'go_network_{0}_{1}'.format(prefix, j)
            print("Saving {}".format(fig_name))

            out_file = os.path.join(out_dir, 'Figures',
                                    '{}.png'.format(fig_name))
            if os.path.exists(out_file):
                os.remove(out_file)

            self.create_png(out_file, int(x['NETWORK_WIDTH']))
            if labels is None:
                trip_photo(out_file, j)
            else:
                trip_photo(out_file, labels[ind])

    def update_node_color(self, attribute, save_name):
        self.cy.style.apply(style=self.style, network=self.g_cy)
        node_color_values = {self.node_name2id[i[0]]: i[1][attribute] for i in
                             self.graph.nodes(data=True)}
        self.view1.update_node_views(visual_property='NODE_FILL_COLOR',
                                     values=node_color_values)
        network_pdf = self.g_cy.get_svg()
        with open('{0}.svg'.format(save_name), 'wb') as f:
            f.write(network_pdf)
            f.close()
        network_pdf = self.g_cy.get_pdf()
        with open('{0}.pdf'.format(save_name), 'wb') as f:
            f.write(network_pdf)
            f.close()
        network_pdf = self.g_cy.get_png()
        with open('{0}.png'.format(save_name), 'wb') as f:
            f.write(network_pdf)
            f.close()

    def create_png(self, out_file, width):
        with open(out_file, 'wb') as f:
            network_svg = self.g_cy.get_png(width)
            f.write(network_svg)

    def get_svg(self, height=2000):
        url = '%sviews/first.svg?h=%d' % (self.g_cy._CyNetwork__url, height)
        return requests.get(url).content


def trip_photo(im_location, title):
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
    to_remove_x = set()
    to_remove_y = set()
    x_dim = img.shape[0]
    y_dim = img.shape[1]
    for i in range(x_dim):
        if img[i, :, 0].sum() != y_dim:
            to_remove_x.add(i)
    for i in range(y_dim):
        if img[:, i, 0].sum() != x_dim:
            to_remove_y.add(i)
    if len(to_remove_x) > 2:
        img = img[min(to_remove_x):max(to_remove_x), :, :]
    if len(to_remove_y) > 2:
        img = img[:, min(to_remove_y):max(to_remove_y), :]

    plt.imshow(img)
    plt.xticks([])
    plt.yticks([])
    plt.title(title)
    plt.axis('off')
    out = im_location.replace('.png', '_formatted.png')
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    # plt.show()


if __name__ == '__main__':
    ddn = nx.nx.read_graphml('t_all_colored_pvalue_2.graphml')
    rm = RenderModel(ddn, style='Marquee')
    rm.visualize_by_list_of_time(
            ['time_0', 'time_1', 'time_2', 'time_3', 'time_4', ])
