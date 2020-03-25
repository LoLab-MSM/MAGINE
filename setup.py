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
with open(path.join(here, 'README.rst'), 'r') as f:
    long_description = f.read()

setup(
    name='MAGINE',
    version='0.0.11',
    description='Package to analyze biological data.',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/LoLab-VU/Magine',
    author='James Pino',
    author_email='james.ch.pino@gmail.com',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],

    keywords=[
        'biological networks', 'biological pathways', 'enrichment analysis',
        'network analysis', 'visualization',
        'multi-omics', 'rnaseq', 'proteomics', 'metabolomics'
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['docs']),

    install_requires=[
        'scipy',
        'numpy',
        'pandas',
        'statsmodels',

        'jinja2',

        'requests',

        'sortedcontainers',
        'defusedxml',
        'xlrd',

        'matplotlib',
        'seaborn',
        'matplotlib-venn',
        'plotly==2.7',
        'wordcloud',

        'bioservices',

        'pathos',
        'jupyter',
        'ipywidgets',
        'py2cytoscape',
        'pydot',
        'pydotplus',
        'networkx>=2.1',

    ],

    test_suite='nose.collector',
    tests_require=['nose', 'coverage'],

    extras_require={
        'test': ['coverage'],
    },
    include_package_data=True,
    scripts=['scripts/create_template_project',
             'scripts/download_databases.py'],

)
