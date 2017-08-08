MAGINE
======

.. image:: https://coveralls.io/repos/github/LoLab-VU/Magine/badge.svg?branch=master
    :target: https://coveralls.io/github/LoLab-VU/Magine?branch=master

.. image:: https://travis-ci.org/LoLab-VU/Magine.svg?branch=master
    :target: https://travis-ci.org/LoLab-VU/Magine

.. image:: https://landscape.io/github/LoLab-VU/Magine/master/landscape.svg?style=flat
   :target: https://landscape.io/github/LoLab-VU/Magine/master

.. image:: https://readthedocs.org/projects/magine/badge/?version=latest
    :target: http://magine.readthedocs.io/en/latest/?badge=latest


Mechanism of Action Generator involving Network Expansion

Python package that integrates omic style data.



Installation
============

1. **Install Anaconda**

   Our recommended approach is to use `Anaconda`_, which is a distribution of
   Python containing most of the numeric and scientific software needed to
   get started. If you are a Mac or Linux user, have used Python before and
   are comfortable using ``pip`` to install software, you may want to skip
   this step and use your existing Python installation.

   Anaconda has a simple graphical installer which can be downloaded from
   https://www.continuum.io/downloads - select your operating system
   and download the **Python 2.7 version**. The default installer options
   are usually appropriate.


        **Windows users:** If you are unsure whether to use the 32-bit or
           64-bit installer, press the Windows Start button, search for “About
           your PC”, and under “System type” it will specify 32-bit operating
           system or 64-bit operating system

2. **Open a terminal**

    We will install most packages with conda::

       $ conda create -n magine_env python=2
       $ source activate magine_env
       $ conda config --add channels conda-forge
       $ conda install numba bottleneck orange jinja2 statsmodels networkx graphviz python-igraph

3. (Linux and MacOS only) **Install pygraphviz**
    ::

     $conda install -c pdrops pygraphviz=1.2


4. (Windows only) **Install pygraphviz**

    Windows users pygraphviz on Windows can be a little troublesome to
    install. Luckily there are binaries that can be downloaded from here
    http://www.lfd.uci.edu/~gohlke/pythonlibs/ . Select the same python version (27).
    One it is downloaded you can install with::

      $ pip install PATH_TO_DOWNlOAD

    where you would type the .whl file that you download.


5. **Install MAGINE**

   The installation is very straightforward with ``pip`` - type the
   following in a terminal::

      $ git clone https://github.com/LoLab-VU/magine
      $ pip install -r requirements.txt
      $ export PYTHONPATH=`pwd`:$PYTHONPATH

     **Mac users:** To open a terminal on a Mac, open Spotlight search
            (press command key and space), type ``terminal`` and press enter.


6. **Install MAGINE (not currently working, please follow step 5)**

   The installation is very straightforward with ``pip`` - type the
   following in a terminal::

      $ pip install git+git:https://github.com/LoLab-VU/Magine

7. **Start Python and MAGINE**

   If you installed Python using `Anaconda`_ on Windows, search for and select
   ``IPython`` from your Start Menu (Windows). Otherwise, open a terminal
   and type ``python`` to get started (or ``ipython``, if installed).

   You will then be at the Python prompt. Type ``import magine`` to try
   loading magine. If no error messages appear and the next Python prompt
   appears, you have succeeded in installing magine!


Documentation
-------------

The manual is available online at http://docs.pysb.org/. You can also
generate the documentation locally by installing Sphinx and running
the following commands::

    $ cd doc
    $ make html

Then open _build/html/index.html in your web browser.