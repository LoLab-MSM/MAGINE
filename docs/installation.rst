Installation
============


Install MAGINE natively on your computer
----------------------------------------

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

   .. note::
**Windows users:** If you are unsure whether to use the 32-bit or
       64-bit installer, press the Windows Start button, search for “About
       your PC”, and under “System type” it will specify 32-bit operating
       system or 64-bit operating system

2. **Open a terminal**

    We will install most packages with conda

    :command:`conda config --add channels conda-forge`
    :command:`conda install orange jinja2 statsmodels networkx graphviz python-igraph`

3. (Linux and MacOS only) **Install pygraphviz**

    :command:`conda install -c pdrops pygraphviz=1.2`

4. (Windows only) **Install pygraphviz**

    Windows users pygraphviz on Windows can be a little troublesome to
    install. Luckily there are binaries that can be downloaded from here
    http://www.lfd.uci.edu/~gohlke/pythonlibs/ . Select the same python version (27).
    One it is downloaded you can install with
    :command:`pip install PATH_TO_DOWNlOAD`
    where you would type the .whl file that you download.


5. **Install MAGINE**

   The installation is very straightforward with ``pip`` - type the
   following in a terminal:

       :command:`pip install git+git:https://github.com/LoLab-VU/Magine`

   .. note::
**Mac users:** To open a terminal on a Mac, open Spotlight search
       (press command key and space), type ``terminal`` and press enter.

6. **Start Python and MAGINE**

   If you installed Python using `Anaconda`_ on Windows, search for and select
   ``IPython`` from your Start Menu (Windows). Otherwise, open a terminal
   and type ``python`` to get started (or ``ipython``, if installed).

   You will then be at the Python prompt. Type ``import magine`` to try
   loading magine. If no error messages appear and the next Python prompt
   appears, you have succeeded in installing magine!

Recommended additional software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following software is not required for the basic operation of magine, but
provides extra capabilities and features when installed.

* `matplotlib`_

  This Python package allows you to plot the results of your simulations. It
  is not a hard requirement of magine but many of the example scripts use it.
  `matplotlib`_ is included with `Anaconda`_. Otherwise, it can be installed
  with :command:`pip install matplotlib`.

* `pandas`_

  This Python package provides extra capabilities for examining large
  numerical datasets, with statistical summaries and database-like
  manipulation capabilities. It is not a hard requirement of magine, but it is a
  useful addition, particularly with large sets of simulation results.
  `pandas`_ is included with `Anaconda`_. Otherwise, it can be installed with
  :command:`pip install pandas`.

* `IPython`_

  An alternate interactive Python shell, much improved over the standard one.
  `IPython`_ is included with `Anaconda`_. Otherwise, it can be installed
  with :command:`pip install ipython`.




.. _Anaconda: https://www.continuum.io/downloads
.. _Git: http://git-scm.com/
.. _IPython: http://ipython.org/
.. _GraphViz: http://www.graphviz.org/
.. _pandas: http://pandas.pydata.org/
.. _Python: http://www.python.org/
.. _SciPy: http://www.scipy.org/
.. _NumPy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/

