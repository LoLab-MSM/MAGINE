import os
import time
from random import randint

import seaborn as sns

import magine.networks.network_tools as tools
from magine.networks.utils import networkx_to_igraph


def create_figure(mol_net, save_name=None):
    """

    Parameters
    ----------
    mol_net : networkx.DiGraph
        networkx graph
    save_name : str
        Save name

    Returns
    -------

    """
    try:
        import cairo
    except ImportError:
        print("Please install pycairo to use igraph plotting")
        return

    try:
        import igraph
        from igraph.drawing.text import TextDrawer
        from igraph.drawing.colors import color_name_to_rgba
    except ImportError:
        print("No igraph, cannot use plotting function")
        return

    g = networkx_to_igraph(mol_net)
    if not isinstance(g, igraph.Graph):
        return
    cl = igraph.VertexClustering(g).FromAttribute(g, attribute='termName')
    membership = cl.membership

    if membership is not None:
        gcopy = g.copy()
        edges = []
        edges_colors = []
        for edge in g.es():
            if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
                # edges.append(edge)
                # edges_colors.append("gray")
                edges_colors.append("black")
            else:
                edges_colors.append("black")
        # gcopy.delete_edges(edges)
        # layout = gcopy.layout("kk")
        # layout = gcopy.layout("tree")
        # layout = gcopy.layout("drl")
        layout = gcopy.layout("graphopt")
        # layout = gcopy.layout("mds")
        # layout = gcopy.layout("sugiyama")
        # layout = gcopy.layout("fr")
        g.es["color"] = edges_colors
    visual_style = dict()
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["edge_color"] = g.es["color"]
    visual_style["vertex_size"] = 50
    visual_style["layout"] = layout
    # visual_style["bbox"] = (4000, 2500)
    # visual_style["bbox"] = (1000, 1000)
    visual_style["bbox"] = (2000, 2000)
    visual_style["margin"] = 100
    # visual_style["node_label"] = g.vs["label"]
    # for vertex in g.vs():
    #     vertex["label"] = vertex['name']

    colors = []
    for i in range(0, max(membership) + 1):
        # colors.append(k_colors[randint(0, len(k_colors))])
        colors.append('#%06X' % randint(0, 0xFFFFFF))
    _groups = dict()
    for vertex in g.vs():
        _groups[vertex.index] = colors[membership[vertex.index]]
        vertex["color"] = colors[membership[vertex.index]]
    _mark_groups = dict()
    # print(_groups)
    for i in _groups:
        co = _groups[i]
        if co in _mark_groups:
            _mark_groups[co].append(i)
        else:
            _mark_groups[co] = [i]
    mark_groups = dict()
    for i in _mark_groups:
        mark_groups[tuple(_mark_groups[i])] = i
    # print(_mark_groups)
    # print(mark_groups)
    mark_groups = dict()
    visual_style["vertex_color"] = g.vs["color"]

    # igraph.plot(cl, "{}.pdf".format(save_name), mark_groups=True, **visual_style)
    if save_name is not None:
        plot = igraph.Plot("{}.png".format(save_name), bbox=(3000, 2200),
                           background="white")
    else:
        plot = igraph.Plot(bbox=(3000, 2200), background="white")
    plot.add(cl,
             # mark_groups=mark_groups,
             **visual_style)
    # plot = igraph.plot(cl, mark_groups=True, **visual_style)
    plot.redraw()
    # Grab the surface, construct a drawing context and a TextDrawer
    ctx = cairo.Context(plot.surface)
    ctx.set_font_size(36)
    # drawer = TextDrawer(ctx, "Test title", halign=TextDrawer.CENTER)
    # drawer.draw_at(0, 40, width=600)
    labels = dict()
    for vertex in g.vs():
        labels[colors[membership[vertex.index]]] = vertex['termName'].replace(
                ',', '\n')

    for n, i in enumerate(colors):
        text_drawer = TextDrawer(ctx, text=labels[i], halign=TextDrawer.LEFT)
        text_drawer.draw_at(x=2100, y=100 + (n * 200), width=900, wrap=True)
        ctx.set_source_rgba(*color_name_to_rgba(i))
        ctx.arc(2000, 100 + (n * 200), 100 / 2, 0, 2 * 3.14)
        ctx.fill()
        ctx.set_source_rgba(0, 0, 0, 1)

    # Make the plot draw itself on the Cairo surface
    # Save the plot
    if save_name is not None:
        plot.save()
    else:
        return plot._repr_svg_()


def create_igraph_figure(mol_net, save_name=None, layout='auto',
                         positions=None,
                         bbox=(2000, 2000), margin=(10, 200, 10, 10),
                         node_size=50,
                         cluster=False, title=None):
    """

    Parameters
    ----------
    mol_net : networkx.DiGraph
        networkx graph
    save_name : str
        Save name

    Returns
    -------

    """
    try:
        import cairo
    except ImportError:
        print("Please install pycairo to use igraph plotting")
        return

    try:
        import igraph
        from igraph.drawing.text import TextDrawer
        from igraph.drawing.colors import color_name_to_rgba
    except ImportError:
        print("No igraph, cannot use plotting function")
        return
    g = networkx_to_igraph(mol_net)
    if not isinstance(g, igraph.Graph):
        return
    _valid_layouts = {
        "kk", "drl", "lgl", "tree", "graphopt", "mds", "sugiyama", "auto"
    }
    assert layout in _valid_layouts, \
        'layout {} not in {}'.format(layout, _valid_layouts)

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

            rgb_colors = sns.color_palette("Set2", len(membership))
            colors = rgb_colors.as_hex()
            mark_groups = dict()
            for n, color in enumerate(colors):
                l = [i for i, j in enumerate(membership) if j == n]
                mark_groups[tuple(l)] = color

    if positions is None:
        positions = g.layout(layout)

    visual_style = dict()
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["vertex_size"] = node_size
    visual_style["layout"] = positions
    # visual_style["bbox"] = bbox
    visual_style["margin"] = margin

    if 'color' in g.es.attributes():
        visual_style["edge_color"] = g.es["color"]

    if 'color' in g.vs.attributes():
        visual_style["vertex_color"] = g.vs["color"]

    if save_name is not None:
        if not save_name.endswith('.png'):
            save_name += '.png'
    if cluster:
        bbox = (bbox[0] + margin[1] * 1.25, bbox[1] + 100)
    plot = igraph.Plot(save_name, bbox=bbox,
                       background='white')
    plot.add(g,
             mark_groups=mark_groups,
             **visual_style)

    plot.redraw()

    if membership is not None:
        # Grab the surface, construct a drawing context and a TextDrawer
        ctx = cairo.Context(plot.surface)
        ctx.set_font_size(36)
        if title is not None:
            drawer = TextDrawer(ctx, title, halign=TextDrawer.CENTER)
            drawer.draw_at(0, 40, width=600)

        labels = dict()
        for vertex in g.vs():
            labels[vertex['termName'].replace(',', '\n')] = rgb_colors[
                membership[vertex.index]]

        for n, i in enumerate(labels):
            text_drawer = TextDrawer(ctx, text=i, halign=TextDrawer.LEFT)
            y_coord = 100 + (n * 200)
            text_drawer.draw_at(x=bbox[0] - margin[1], y=y_coord, width=400,
                                wrap=True)
            ctx.set_source_rgba(*labels[i])
            ctx.arc(bbox[0] - (margin[1] * 1.2), y_coord - 10, node_size, 0,
                    6.28)
            ctx.fill()
            ctx.set_source_rgba(0, 0, 0, 1)
    # Save the plot
    if save_name is not None:
        plot.save()
    return plot, positions


def paint_network_overtime(graph, exp_data, color_list, save_name,
                           compile_images=False):
    """
    Adds color attribute to network over time.

    Parameters
    ----------
    graph : pygraphviz.AGraph
        Network
    exp_data : list_like
        List of lists, where the inner list contains the node to add the
        color
    color_list : list_like
        list of colors for each time point
    save_name : str
        prefix for images to be saved

    Returns
    -------

    """
    from IPython.display import Image, display
    labels = []
    measured_list = []
    for i, j in sorted(exp_data.sig_species_over_time.items()):
        measured_list.append(j)
        labels.append(i)
    if isinstance(color_list, str):
        color_list = [color_list, ] * len(measured_list)
    if len(measured_list) != len(color_list):
        print('Length of list'
              'of data must equal len of color list')
        return

    string = 'convert -delay 100 '
    tmp_graph = graph.copy()
    tmp_graph = tools._format_to_directions(tmp_graph)
    pos = None
    for n, i in enumerate(measured_list):
        graph2 = tools.add_attribute_to_network(tmp_graph, i, 'color',
                                                color_list[n],
                                                'white')

        s_name = '%s_%04i_igraph.png' % (save_name, n)
        result, pos = create_igraph_figure(graph2, s_name, positions=pos,
                                           node_size=25, bbox=(1000, 1000),
                                           margin=(100, 200, 200, 10),
                                           cluster=True, title=labels[n])

        display(Image(s_name))
        time.sleep(3)
        string += ' ' + s_name
    string1 = string + '  %s.gif' % save_name
    string2 = string + '  %s.pdf' % save_name

    if compile_images:
        os.system(string1)
        os.system(string2)


def _is_running_in_ipython():
    """Internal function that determines whether igraph is running inside
    IPython or not."""
    try:
        # get_ipython is injected into the Python builtins by IPython so
        # this should succeed in IPython but throw a NameError otherwise
        get_ipython
        return True
    except NameError:
        return False
