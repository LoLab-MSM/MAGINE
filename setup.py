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

import magine

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), 'r') as f:
    long_description = f.read()

setup(
    name='MAGINE',
    version=magine.__version__,
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
        'bioservices',
        'defusedxml',
        'ipywidgets',
        'jinja2',
        'jupyter',
        'matplotlib==3.1.0', 'matplotlib-venn',
        'networkx>=2.1,<2.4',
        'numpy',
        'pandas',
        'pathos',
        'plotly==2.7',
        'py2cytoscape',
        'pycairo',
        'pydot',
        'pydotplus',
        'requests',
        'scipy',
        'seaborn',
        'sortedcontainers',
        'statsmodels',
        'wordcloud',
        'xlrd',
    ],

    test_suite='nose.collector',
    tests_require=['nose', 'coverage'],

    extras_require={
        'test': ['coverage'],
    },
    include_package_data=True,
    scripts=['scripts/create_template_project',
             'scripts/download_databases.py'],
    entry_points={
        'console_scripts': [
            'download_magine_databases=download_databases:run'
        ]
    }
)
