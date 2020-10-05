FROM continuumio/miniconda:4.7.12
RUN mkdir /tmp/magine-build

RUN apt-get update -y
RUN apt-get install -y graphviz libgts-dev

COPY docs /tmp/magine-build/docs/
COPY scripts /tmp/magine-build/scripts/
COPY magine /tmp/magine-build/magine/
COPY environment.yml setup.py setup.cfg README.rst MANIFEST.in /tmp/magine-build/
RUN conda config --add channels conda-forge
RUN conda env update -f /tmp/magine-build/environment.yml

RUN cd /tmp/magine-build && python setup.py install
RUN rm -rf /tmp/magine-build

RUN useradd -ms /bin/bash magine
USER magine
WORKDIR /home/magine

COPY magine/examples /home/magine/examples

EXPOSE 8888
ENTRYPOINT [ "jupyter", "notebook", "--ip=*", "--no-browser"]