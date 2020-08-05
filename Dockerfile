FROM continuumio/miniconda:4.7.12

RUN useradd -ms /bin/bash magine
USER magine
WORKDIR /home/magine

COPY magine/examples /home/magine/examples

RUN conda config --add channels conda-forge
RUN conda install jinja2 statsmodels networkx graphviz pip python=3.7
RUN conda install -c marufr python-igraph
RUN python setup.py install

EXPOSE 8888
ENTRYPOINT ["jupyter", "notebook", "--ip=*", "--no-browser"]
