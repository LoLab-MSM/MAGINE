import os
import time

import seaborn as sns

import magine.networks.utils
from magine.networks.exporters import nx_to_igraph


def create_igraph_figure(mol_net, save_name=None, layout='auto', title=None,
                         positions=None, bbox=(2000, 2000), cluster=False,
                         margin=(10, 200, 10, 10), node_size=50):
    """

    Parameters
    ----------

    mol_net : networkx.DiGraph
        networkx graph
    save_name : str
        Save name
    layout : str
    title : str, optional
    positions : list, optional
    bbox : list
    cluster : bool
    margin : list
    node_size : int
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
    g = nx_to_igraph(mol_net)
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
                mem = [i for i, j in enumerate(membership) if j == n]
                mark_groups[tuple(mem)] = color

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
    pos = None
    for n, i in enumerate(measured_list):
        graph2 = magine.networks.utils.add_attribute_to_network(tmp_graph, i,
                                                                'color',
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
