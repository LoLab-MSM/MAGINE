import os

import seaborn as sns

from magine.networks.exporters import nx_to_igraph
from magine.networks.utils import add_attribute_to_network

try:
    basestring
# Allows isinstance(foo, basestring) to work in Python 3
except:
    basestring = str


def draw_igraph(mol_net, save_name=None, layout='auto', title=None,
                positions=None, cluster=False, node_size=50,
                bbox=None, margin=None, inline=False, node_font_size=24,
                font_size=36):
    """

    Parameters
    ----------

    mol_net : networkx.DiGraph
        networkx graph
    save_name : str
        Save name
    layout : str
        Options = "kk", "fr", "drl", "lgl", "tree", "graphopt", "mds", "sugiyama",
        "auto", "grid_fr"
    title : str, optional
    positions : list, optional
    bbox : list
    cluster : bool
        Add a shape around species that are tagged the same.
         g.node[0][term_name]='dna'
         g.node[0][term_name]='dna'
         This will be grouped together.
    margin : list
        Margin of white space for boxes. This should be considered for the
         labels of the nodes
    node_size : int
    node_font_size: int
        Size of node labels
    font_size : int
        Title and group name (if clustering) font size
    inline : bool
        Return plot as svg view ipython

    Returns
    -------

    """
    try:
        import cairo
    except ImportError:
        ImportError("Please install pycairo to use igraph plotting")

    try:
        import igraph
        from igraph.drawing.text import TextDrawer
        from igraph.drawing.colors import color_name_to_rgba
    except ImportError:
        raise ImportError("No igraph, cannot use plotting function")

    if isinstance(mol_net, igraph.Graph):
        g = mol_net
    else:
        g = nx_to_igraph(mol_net)
        if not isinstance(g, igraph.Graph):
            print("Error converting to Igraph")
            return
    if inline:
        if bbox is None:
            bbox = [500, 500]
        if margin is None:
            margin = [50, 50, 50, 50]

    if bbox is None:
        bbox = [2400, 2400]
    if margin is None:
        margin = [50, 50, 50, 50]

    _valid_layouts = {
        "kk", "fr", "drl", "lgl", "tree", "graphopt", "mds", "sugiyama",
        "auto", "grid_fr",
    }
    if layout not in _valid_layouts:
        raise Exception('layout {} not in {}'.format(layout, _valid_layouts))

    mark_groups = None
    membership = None

    if cluster and 'termName' in g.vs.attributes():
        cl = igraph.VertexClustering(g).FromAttribute(g, attribute='termName')
        membership = cl.membership
        if membership is not None:
            gcopy = g.copy()
            # edges = []
            # for edge in g.es():
            #     if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
            #         edges.append(edge)
            # gcopy.delete_edges(edges)
            if positions is None:
                positions = gcopy.layout(layout)
            n_clusters = len(set(membership))
            rgb_colors = sns.color_palette("tab20", n_clusters)
            colors = rgb_colors.as_hex()
            mark_groups = dict()
            for n, color in enumerate(colors):
                mem = tuple(i for i, j in enumerate(membership) if j == n)
                mark_groups[mem] = color

    if positions is None:
        positions = g.layout(layout)

    visual_style = dict()
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["vertex_size"] = node_size
    visual_style["vertex_label_size"] = node_font_size
    visual_style["layout"] = positions
    visual_style["margin"] = margin
    # visual_style["edge_curved"] = True

    if 'color' in g.es.attributes():
        visual_style["edge_color"] = g.es["color"]

    if save_name is not None:
        if not save_name.endswith('.png'):
            save_name += '.png'

    if cluster:
        # add some white space to add labels
        bbox_plus_margin = [bbox[0] + bbox[0] * .45, bbox[1] + bbox[0] * .1]
        margin[2] += bbox[0] * .25
        margin[0] += bbox[0] * .05
    else:
        bbox_plus_margin = bbox

    # create entire surface
    plot = igraph.Plot(save_name,
                       bbox=bbox_plus_margin,
                       background='white')

    # add plot to surface
    plot.add(g, mark_groups=mark_groups, bbox=bbox, **visual_style)

    # have to redraw to add the plot
    plot.redraw()

    # Grab the surface, construct a drawing context and a TextDrawer
    ctx = cairo.Context(plot.surface)
    ctx.set_font_size(font_size)

    if title is not None:
        drawer = TextDrawer(ctx, title, halign=TextDrawer.CENTER)
        drawer.draw_at(0, 40, width=600)

    if membership is not None:
        labels = dict()
        for vertex in g.vs():
            name = vertex['termName'].replace(',', '\n')
            if name in labels:
                continue
            labels[name] = rgb_colors[membership[vertex.index]]
        spacing = bbox[1] / len(labels)
        ctx.set_font_size(24)
        for n, (label, rgb_c) in enumerate(labels.items()):
            text_drawer = TextDrawer(ctx, text=label, halign=TextDrawer.LEFT)
            x_coord = bbox[0] + 25
            y_coord = 100 + (n * spacing)
            text_drawer.draw_at(x=x_coord, y=y_coord, width=300, wrap=True)
            ctx.set_source_rgba(*rgb_c)
            ctx.arc(x_coord - 1.5 * node_size, y_coord, node_size, 0, 6.28)
            ctx.fill()
            ctx.set_source_rgba(0, 0, 0, 1)
    # Save the plot
    if save_name is not None:
        plot.save()
    if inline:
        return plot

    return plot, positions


def paint_network_overtime(graph, exp_data, color_list, save_name,
                           compile_images=False, cluster=False, fig_width=1500,
                           fig_height=1500, layout='auto', inline=False):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    compile_images
    graph : nx.DiGraph
        Network
    exp_data : list_like
        List of lists, where the inner list contains the node to add the
        color
    color_list : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved
    cluster : bool
        Group nodes based on category 'term_name'
    fig_width : int
        Size of fig
    fig_height : int
        Size of fig
    layout : str
        Layout format for igraph
    Returns
    -------

    """

    labels = []
    measured_list = []
    for i, j in zip(exp_data.species.sig.sample_ids,
                    exp_data.species.sig.by_sample):
        measured_list.append(j)
        labels.append(i)
    if isinstance(color_list, basestring):
        color_list = [color_list, ] * len(measured_list)
    if len(measured_list) != len(color_list):
        print('Length of list'
              'of data must equal len of color list')
        return

    string = 'convert -delay 100 '
    tmp_graph = graph.copy()
    pos = None
    for n, i in enumerate(measured_list):
        graph2 = add_attribute_to_network(tmp_graph, i, 'color', color_list[n],
                                          'white')

        s_name = '%s_%04i_igraph.png' % (save_name, n)
        result, pos = draw_igraph(graph2, s_name, positions=pos, node_size=50,
                                  bbox=[fig_width, fig_height],
                                  margin=[100, 100, 100, 100],
                                  cluster=cluster, title=labels[n],
                                  layout=layout)
        if inline:
            from IPython.display import Image, display
            display(Image(s_name))
        string += ' ' + s_name
        # yield result
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name

    if compile_images:
        os.system(string1)
        os.system(string2)


def color_by_list(graph, species_list, labels, colors, save_name,
                  compile_images=False, cluster=False, fig_width=1500,
                  fig_height=1500, layout='auto', inline=False, **kwargs):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    compile_images
    graph : nx.DiGraph
        Network
    species_list : list_like
        List of lists, where the inner list contains the node to add the
        color
    labels : list_like
        List of labels
    colors : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved
    cluster : bool
        Group nodes based on category 'term_name'
    fig_width : int
        Size of fig
    fig_height : int
        Size of fig
    layout : str
        Layout format for igraph
    inline : bool
        Display figures inline (for ipython/jupyter)
    Returns
    -------

    """
    if inline:
        from IPython.display import Image, display
    if isinstance(colors, str):
        colors = [colors, ] * len(species_list)
    if len(species_list) != len(colors):
        print('Length of list'
              'of data must equal len of color list')
        return
    if labels is None:
        labels = [None] * len(species_list)
    string = 'convert -delay 100 '
    tmp_graph = graph.copy()
    pos = None
    for n, i in enumerate(species_list):
        graph2 = add_attribute_to_network(tmp_graph, i, 'color', colors[n],
                                          'white')

        s_name = '%s_%04i_igraph.png' % (save_name, n)
        result, pos = draw_igraph(graph2, s_name, positions=pos,
                                  title=labels[n], **kwargs)
        if inline:
            display(Image(s_name))
        string += ' ' + s_name
        # yield result
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name

    if compile_images:
        os.system(string1)
        os.system(string2)


if __name__ == '__main__':
    import networkx as nx
    import igraph

    # g = nx.read_gpickle('../background_network.p.gz')
    g = nx.DiGraph()
    g.add_edge("a", "b")
    print(type(g))
    draw_igraph(g, save_name='test', layout='lgl')
