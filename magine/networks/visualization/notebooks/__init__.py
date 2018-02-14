import os


def init():
    from IPython.core.display import display, Javascript
    JS_LOADER_FILE = "loader.js"
    path = os.path.abspath(os.path.dirname(__file__)) + "/" + JS_LOADER_FILE
    js_loader = open(path).read()
    display(Javascript(js_loader))


# Load Cytoscape.js and dependent libraries
init()
