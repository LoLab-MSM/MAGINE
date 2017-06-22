"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# To use a consistent encoding
from codecs import open
from os import path

# Always prefer setuptools over distutils
from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
        name='MAGINE',
        version='0.0.1',
        description='Package to analyze biological data.',
        long_description=long_description,
        url='https://github.com/LoLab-VU/Magine',
        author='James Pino',
        author_email='james.ch.pino@gmail.com',
        license='MIT',

        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            'Development Status :: 3 - Alpha',
        ],

        keywords=['biological networks', 'rnaseq', 'biological pathways',
                  'proteomics', 'metabolomics'],

        # You can just specify the packages manually here if your project is
        # simple. Or you can use find_packages().
        packages=find_packages(exclude=['contrib', 'docs']),

        # List run-time dependencies here.  These will be installed by pip when
        # your project is installed. For an analysis of "install_requires" vs pip's
        # requirements files see:
        # https://packaging.python.org/en/latest/requirements.html
        install_requires=['matplotlib',
                          'jinja2',
                          'networkx',
                          'requests',
                          'py2cytoscape',
                          'orange-bioinformatics==2.6.21',
                          'pygraphviz',
                          'pandas',
                          'xlrd',
                          'matplotlib',
                          'matplotlib-venn',
                          'bioservices',
                          'pathos',
                          'goatools',
                          'plotly'],

        test_suite='nose.collector',
        tests_require=['nose', 'coverage'],

        extras_require={
            'dev':  ['check-manifest'],
            'test': ['coverage'],
        },
        include_package_data=True,
        scripts=['scripts/create_template_project'],

        # If there are data files included in your packages that need to be
        # installed, specify them here.  If using Python 2.6 or less, then these
        # have to be included in MANIFEST.in as well.
        # package_data={
        #     'sample': ['package_data.dat'],
        # },

        # Although 'package_data' is the preferred approach, in some case you may
        # need to place data files outside of your packages. See:
        # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
        # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
        # data_files=[('my_data', ['data/data_file'])],

)
