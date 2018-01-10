from random import randint

from magine.networks.network_tools import networkx_to_igraph


def create_figure(mol_net, save_name):
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
    for vertex in g.vs():
        vertex["label"] = vertex['name']

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
    plot = igraph.Plot("{}.png".format(save_name), bbox=(3000, 2200),
                       background="white")
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
    plot.save()


def create_igraph_figure(mol_net, save_name):
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

    layout = g.layout("kk")
    # layout = g.layout("tree")
    # layout = g.layout("drl")
    # layout = g.layout("graphopt")
    # layout = g.layout("mds")
    # layout = g.layout("sugiyama")
    # layout = g.layout("fr")

    visual_style = dict()
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["vertex_size"] = 50
    visual_style["layout"] = layout
    # visual_style["bbox"] = (4000, 2500)
    # visual_style["bbox"] = (1000, 1000)
    visual_style["bbox"] = (2000, 2000)
    visual_style["margin"] = 100
    visual_style["vertex_color"] = g.vs["color"]
    for vertex in g.vs():
        vertex["label"] = vertex['name']
    # print(_mark_groups)
    # print(mark_groups)

    # igraph.plot(cl, "{}.pdf".format(save_name), mark_groups=True, **visual_style)
    plot = igraph.Plot("{}.png".format(save_name), bbox=(3000, 2200),
                       background="white")
    plot.add(g, **visual_style)
    # plot = igraph.plot(cl, mark_groups=True, **visual_style)
    plot.redraw()
    # Save the plot
    plot.save()
