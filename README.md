[![Coverage Status](https://coveralls.io/repos/github/LoLab-VU/Magine/badge.svg?branch=master)](https://coveralls.io/github/LoLab-VU/Magine?branch=master)

[![Build Status](https://travis-ci.org/LoLab-VU/Magine.svg?branch=master)](https://travis-ci.org/LoLab-VU/Magine)

# Magine- Mechanism of Action Generator involving Network Expansion

Python package that integrates omic style data.

### Requirments
- python 2 or python 3
- Graphviz (http://www.graphviz.org)


# Installation

### Installing using Anaconda with python 2
- Download and install Miniconda https://repo.continuum.io/miniconda/
- Open a terminal
- conda create -n magine_env python=2
- source activate magine_env
- conda config --add channels conda-forge
- conda install orange jinja2 statsmodels networkx graphviz python-igraph

- If you have MacOSX or Linux, you can install pygraphviz using
- - conda install -c pdrops pygraphviz=1.2
- git clone https://github.com/LoLab-VU/magine
- cd magine
- python setup.py install

### Installing with pip
- pip install git+git:https://github.com/LoLab-VU/Magine