[![Coverage Status](https://coveralls.io/repos/github/LoLab-VU/Magine/badge.svg?branch=master)](https://coveralls.io/github/LoLab-VU/Magine?branch=master)

[![Travis Status](https://travis-ci.org/LoLab-VU/Magine.svg?branch=master)](https://travis-ci.org/LoLab-VU/Magine.svg?branch=master)


#Magine- Mechanism of Action Generator involving Network Expansion

Python package that integrates omic style data.

#Requirments
Graphviz (http://www.graphviz.org)

# Installation
## For people new to python enviroment creation
## Anaconda is a great easy way to get you up and running.
## Open a terminal
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh
## follow prompts to install miniconda
## source ~/.bashrc if you added it to your path
conda create -n magine_env python=2
source activate magine_env

- conda config --add channels conda-forge
- conda install orange
- conda install jinja2
- conda install statsmodels
- conda install networkx
- conda install graphviz
- conda install python-igraph
- If you have MacOSX or Linux, you can install pygraphviz using
- - conda install -c pdrops pygraphviz=1.2
- git clone https://github.com/LoLab-VU/magine
- cd magine
- python setup.py install
