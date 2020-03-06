import os

import networkx as nx

from magine.networks.exporters import nx_to_dot, export_to_dot
from magine.networks.utils import add_color_graphviz_fmt

IPYTHON = False
try:
    from IPython.display import Image, display

    IPYTHON = True
except RuntimeError:
    pass


def draw_graphviz(network, layout='dot', image_format='png', width=800,
                  save_name=None):
    """ Draws an


    Parameters
    ----------
    network : nx.DiGraph
    layout : str
        Which graphviz engine to use to layout graph
    image_format: str
        image format (png, svg, pdf)
    width : int
        Width of image
    save_name : str, optional
        If you want to export to file

    Returns
    -------

    """
    net_copy = network.copy()
    net_copy.graph.update({'K': '1'})
    net_copy.graph.update({'repulsiveforce ': '1'})
    net_copy.graph.update({'overlap ': 'false'})
    net_copy.graph.update({'splines ': 'true'})
    if save_name is not None:
        if save_name.endswith(image_format):
            out_name = save_name
        else:
            out_name = '{}.{}'.format(save_name, image_format)
        nx_to_dot(net_copy).write(out_name, format=image_format,
                                  prog=layout)
        if IPYTHON and run_from_ipython():
            display(Image(out_name, width=width))
    else:
        img = nx_to_dot(net_copy).create(format='png', prog=layout)
        if IPYTHON and run_from_ipython():
            display(Image(img, width=width))


def paint_network_overtime(graph, list_of_lists, color_list, save_name,
                           labels=None, create_gif=False):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    graph : nx.DiGraph
        Network
    list_of_lists : list_like
        List of lists, where the inner list contains the node to add the
        color
    color_list : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved
    labels: list_like
        list of labels to add to graph per sample

    """

    if len(list_of_lists) != len(color_list):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_of_lists):
            print('Length of labels must be equal to len of data')
            return

    string = 'convert -delay 100 '
    tmp_graph = graph.copy()

    for n, i in enumerate(list_of_lists):
        graph2 = add_color_graphviz_fmt(tmp_graph, i, color_list[n])

        if labels is not None:
            graph2.graph['label'] = labels[n]
            graph2.graph['fontsize'] = 13

        s_name = '%s_%04i.png' % (save_name, n)
        draw_graphviz(graph2, layout='dot', save_name=s_name,
                      image_format='png')
        string += s_name + ' '
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name
    if create_gif:
        os.system(string1)
        os.system(string2)


def paint_network_overtime_up_down(graph, list_up, list_down, save_name,
                                   color_up='red', color_down='blue',
                                   labels=None, create_gif=False):
    """
    Adds color attribute to network over time and creates figure.

    Parameters
    ----------
    graph : nx.DiGraph
        Network
    list_up : list_like
        List of lists, where the inner list contains the node to add the
        color
    list_down : list_like
        list of colors for each time point
    color_up : str
        color for first list of species
    color_down : str
        color of second list of species
    save_name : str
        prefix for images to be saved
    labels: list_like
        list of labels to add to graph per sample
    create_gif : bool
        Create a gif from series of images


    """

    if len(list_up) != len(list_down):
        print('Length of list of data must equal len of color list')
        return
    if labels is not None:
        if len(labels) != len(list_down):
            print('Length of labels must be equal to len of data')
            return
    string = 'convert -delay 100 '
    tmp_graph = graph.copy()

    for n, (up, down) in enumerate(zip(list_up, list_down)):
        tmp_graph = add_color_graphviz_fmt(tmp_graph, up, color_up)
        tmp_graph = add_color_graphviz_fmt(tmp_graph, down, color_down)
        both = set(up).intersection(set(down))
        tmp_graph = add_color_graphviz_fmt(tmp_graph, both, 'yellow')

        if labels is not None:
            tmp_graph.graph['label'] = labels[n]
            tmp_graph.graph['fontsize'] = 13

        s_name = '%s_%04i.png' % (save_name, n)

        export_to_dot(tmp_graph, s_name, 'png', 'dot')
        draw_graphviz(tmp_graph, layout='dot', save_name=s_name,
                      image_format='png')
        string += s_name + ' '

    if create_gif:
        os.system('{}  {}.gif'.format(string, save_name))
        os.system('{}  {}.pdf'.format(string, save_name))


def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
