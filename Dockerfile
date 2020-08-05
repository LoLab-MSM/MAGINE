FROM continuumio/miniconda:4.7.12

RUN conda config --add channels conda-forge
RUN conda install jinja2 statsmodels networkx graphviz pip python=3.7 matplotlib
RUN conda install -c marufr python-igraph

RUN mkdir /tmp/magine-build
COPY setup.py setup.cfg README.rst /tmp/magine-build/
COPY docs /tmp/magine-build/docs/
COPY scripts /tmp/magine-build/scripts/
COPY magine /tmp/magine-build/magine/
RUN cd /tmp/magine-build && python setup.py install
RUN rm -rf /tmp/magine-build

RUN useradd -ms /bin/bash magine
USER magine
WORKDIR /home/magine

COPY magine/examples /home/magine/examples

EXPOSE 8888
ENTRYPOINT ["jupyter", "notebook", "--ip=*", "--no-browser"]
