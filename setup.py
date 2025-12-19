#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='ClusterCatcher',
    version='0.1.0',
    description='Single-cell sequencing analysis pipeline for mutation signature detection and cell annotation',
    author='Jake Lehle',
    author_email='',
    url='https://github.com/JakeLehle/ClusterCatcher',
    license='GPL-3.0',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'snakemake_wrapper': ['Snakefile', 'scripts/*', 'envs/*'],
    },
    entry_points={
        'console_scripts': [
            'ClusterCatcher=cli.cli:main',
        ],
    },
    install_requires=[
        'click',
        'snakemake>=7.0',
        'pyyaml',
        'pandas',
        'numpy',
    ],
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='single-cell sequencing mutational-signatures cancer bioinformatics',
)
