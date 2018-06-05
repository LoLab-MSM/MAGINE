MAGINE
======

.. image:: https://api.codacy.com/project/badge/Grade/cba1091c58a246bfb07f7ed7f86afe24
   :alt: Codacy Badge
   :target: https://app.codacy.com/app/james.c.pino/MAGINE?utm_source=github.com&utm_medium=referral&utm_content=LoLab-VU/MAGINE&utm_campaign=badger

.. image:: https://travis-ci.org/LoLab-VU/MAGINE.svg?branch=master
    :target: https://travis-ci.org/LoLab-VU/MAGINE

.. image:: https://codecov.io/gh/LoLab-VU/MAGINE/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/LoLab-VU/MAGINE

.. image:: https://readthedocs.org/projects/MAGINE/badge/?version=latest
    :target: http://magine.readthedocs.io/en/latest/?badge=latest

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/LoLab-VU/MAGINE/master


Mechanism of Action Generator involving Network Expansion



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
       $ conda install jinja2 statsmodels networkx graphviz 
       $ conda install -c marufr python-igraph

3. **Install MAGINE**

   The installation is very straightforward with ``pip`` - type the
   following in a terminal::

      $ git clone https://github.com/LoLab-VU/magine
      $ cd magine
      $ pip install -r requirements.txt
      $ export PYTHONPATH=`pwd`:$PYTHONPATH

     **Mac users:** To open a terminal on a Mac, open Spotlight search
            (press command key and space), type ``terminal`` and press enter.



4. **Start Python and MAGINE**

   If you installed Python using `Anaconda`_ on Windows, search for and select
   ``IPython`` from your Start Menu (Windows). Otherwise, open a terminal
   and type ``python`` to get started (or ``ipython``, if installed).

   You will then be at the Python prompt. Type ``import magine`` to try
   loading magine. If no error messages appear and the next Python prompt
   appears, you have succeeded in installing magine!


Documentation
-------------

The manual is available online at http://magine.readthedocs.io.
