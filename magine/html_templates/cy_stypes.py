styles = {
    'main': [
        {
            'selector': 'node',
            'css': {
                'content': 'data(id)',
                'text-valign': 'center',
                'text-halign': 'center'
            }
        },
        {
            'selector': 'node[speciesType = \'compound\']',
            'css': {
                'content': 'data(chemName)',
                'text-valign': 'center',
                'text-halign': 'center',
                'shape': 'rectangle'
            }
        },
        {
            'selector': 'node[speciesType = \'gene\']',
            'css': {
                'content': 'data(id)',
                'text-valign': 'center',
                'text-halign': 'center',
                'shape': 'rectangle'
            }
        },
        {
            'selector': '$node > node',
            'css': {
                'padding-top': '10px',
                'padding-left': '10px',
                'padding-bottom': '10px',
                'padding-right': '10px',
                'text-valign': 'top',
                'text-halign': 'center',
                'background-color': 'green'
            }
        },
        {
            'selector': 'edge',
            'css': {
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'target-arrow-color': 'black'

            }
        },
        {
            'selector': "edge[interactionType *= 'inhibit'],"
                        "edge[interactionType *= 'deactivat']",
            'css': {
                'curve-style': 'bezier',
                'target-arrow-shape': 'tee',
                'target-arrow-color': 'black'
            }
        },
        {
            'selector': 'edge:selected',
            'css': {
                'background-color': 'yellow',
                'line-color': 'red',
                'label': 'data(interactionType)',
                'target-arrow-color': 'black',
                'source-arrow-color': 'black'
            }
        }
    ]
    ,
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
                     "height": 35,
                     "width": 35,
                     "content": "data(name)",

                     "text-wrap": "wrap",
                     "text-max-width": 95,
                     "shape": "ellipse",
                     "background-color": "data(color)",
                     "text-halign": "center",
                     "font-family": "SansSerif",
                     "text-valign": "center",
                     "border-opacity": 1.0},
             },
            {'selector': 'node[speciesType = \'compound\']',
             'css': {"text-opacity": 1.0,
                     "background-opacity": 1.0,
                     "font-weight": "normal",
                     "color": "rgb(0,153,234)",
                     "border-width": 3.0,
                     "border-color": "rgb(51,51,51)",
                     "font-size": 9,
                     "height": 35,
                     "width": 35,
                     'content': 'data(chemName)',

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
             "css": {"opacity": 1.0, "text-opacity": 1.0,
                     "font-size": 12, "font-weight": "normal",
                     "target-arrow-shape": "triangle",
                     "source-arrow-shape": "none",
                     "width": "mapData(weight, 1, 10, 1, 5)",
                     "color": "rgb(0,0,0)",
                     "source-arrow-color": "rgb(0,0,0)",
                     "line-color": "rgb(51,51,51)",
                     "content": "data(interaction)",
                     "curve-style": "bezier",
                     "line-style": "solid",
                     "font-family": "SansSerif",
                     "target-arrow-color": "rgb(51,51,51)"},
             },
            {"selector": "edge[interactionType *= 'inhibit']"
                         ",edge[interactionType *= 'deactivat']",
             "css": {'curve-style': 'bezier',
                     'target-arrow-shape': 'tee',
                     },
             },
            {"selector": "edge:selected",
             "css": {"line-color": "rgb(255,0,0)",
                     "label": "data(interactionType)",
                     },
             }
        ]
}
