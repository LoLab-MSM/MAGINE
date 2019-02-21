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
    packages=find_packages(exclude=['docs']),

    install_requires=['jinja2',
                      'networkx',
                      'requests',
                      'pandas',
                      'xlrd',
                      'matplotlib',
                      'matplotlib-venn',
                      'bioservices',
                      'pathos',
                      'plotly'],

    test_suite='nose.collector',
    tests_require=['nose', 'coverage'],

    extras_require={
        'test': ['coverage'],
    },
    include_package_data=True,
    scripts=['scripts/create_template_project',
             'scripts/download_databases.py'],

)
