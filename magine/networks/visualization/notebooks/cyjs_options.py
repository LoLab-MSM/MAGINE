layouts = {
    'breadthfirst': {
        'name': 'breadthfirst',
        'directed': 'true',
        'spacingFactor': 1.5
    },
    'cose': {
        'name': 'cose',
        'spacingFactor': 1.5
    },
    'cose-bilkent': {
        'name': 'cose-bilkent',
        'spacingFactor': 1.5,
        'animate': False
    },
    'dagre': {
        'name': 'dagre',
        'rankDir': 'LR',
        'spacingFactor': 1.5
    },
    'concentric': {
        'name': 'concentric',
        'spacingFactor': 1.5

    }
}

styles = {
    'default':
        [
            {"selector": "node",
             "css": {"text-opacity": 1.0,
                     "background-opacity": 1.0,
                     "font-weight": "normal",
                     "color": "rgb(0,153,234)",
                     "border-width": 3.0,
                     "border-color": "rgb(51,51,51)",
                     "font-size": 9,
                     "height": 50,
                     "width": 50,
                     "label": "data(name)",

                     "text-wrap": "wrap",
                     "text-max-width": 50,
                     "shape": "ellipse",
                     "background-color": "data(color)",
                     "text-halign": "center",
                     "font-family": "SansSerif",
                     "text-valign": "center",
                     "border-opacity": 1.0
                     },
             },
            {'selector': "node[speciesType = 'compound']",
             'css': {"text-opacity": 1.0,
                     "background-opacity": 1.0,
                     "font-weight": "normal",
                     "color": "rgb(0,153,234)",
                     "border-width": 3.0,
                     "border-color": "rgb(51,51,51)",
                     "font-size": 9,
                     "height": 35,
                     "width": 35,
                     'label': 'data(chemName)',

                     "text-wrap": "wrap",
                     "text-max-width": 95,
                     "shape": "ellipse",
                     "background-color": "rgb(255,255,255)",
                     "text-halign": "center",
                     "font-family": "SansSerif",
                     "text-valign": "center",
                     "border-opacity": 1.0},
             },

            {"selector": "$node > node",
             "css": {"font-size": 20,
                     "shape": "ellipse",
                     "padding-right": "10px",
                     "padding-bottom": "10px",
                     "padding-top": "10px",
                     "text-valign": "top",
                     "text-halign": "center",
                     "background-color": "#bbb",
                     "padding-left": "10px"},
             },
            {"selector": "node:selected",
             "css": {"background-color": "rgb(255,0,102)"},
             },

            {"selector": "edge",
             "css": {"opacity": 1.0,
                     "text-opacity": 1.0,
                     "font-size": 12,
                     "font-weight": "normal",

                     'width': '2',

                     "color": "rgb(0,0,0)",
                     "line-color": "rgb(51,51,51)",
                     "source-arrow-color": "rgb(0,0,0)",
                     "target-arrow-color": "rgb(51,51,51)",

                     "curve-style": "bezier",
                     "line-style": "solid",
                     "target-arrow-shape": "triangle",
                     "source-arrow-shape": "none",

                     "font-family": "SansSerif"

                     },
             },

            {"selector": "edge[weight>1]",
             "css": {"opacity": 1.0,
                     "text-opacity": 1.0,
                     "font-size": 12,
                     "font-weight": "normal",

                     # 'width': 'mapData(width, 1, 10, 1, 10)',
                     'width': 'data(width)',

                     "color": "rgb(0,0,0)",
                     "line-color": "rgb(51,51,51)",
                     "source-arrow-color": "rgb(0,0,0)",
                     "target-arrow-color": "rgb(51,51,51)",

                     "curve-style": "bezier",
                     "line-style": "solid",
                     "target-arrow-shape": "triangle",
                     "source-arrow-shape": "none",

                     "font-family": "SansSerif"

                     },
             },
            {"selector": "edge[interactionType *= 'inhibit']"
                         ",edge[interactionType *= 'deactivat']",
             "css": {
                 "text-opacity": 1.0,
                 "font-size": 12,
                 "font-weight": "normal",
                 "font-family": "SansSerif",

                 'curve-style': 'bezier',
                 'source-arrow-shape': 'none',
                 'target-arrow-shape': 'tee',

                 "color": "rgb(0,0,0)",
                 "line-color": "rgb(51,51,51)",
                 "source-arrow-color": "rgb(0,0,0)",
                 "target-arrow-color": "rgb(51,51,51)",
                 "line-style": "solid",

             },
             },
            {"selector": "edge:selected",
             "css": {
                 "line-color": "rgb(255,0,0)",
                 "label": "data(interactionType)",
             },
             },
            {'selector': 'edge[weight>1]:selected',
             'css': {
                 'background-color': 'yellow',
                 'line-color': 'red',
                 'label': 'data(label)',
                 'target-arrow-color': 'black',
                 'source-arrow-color': 'black',
                 'width': 'data(width)',
             },
             }
        ]
}
