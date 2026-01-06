#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClusterCatcher Setup Script
============================

Single-cell sequencing analysis pipeline for mutation signature detection
and cell annotation.

Installation:
    pip install .
    
    # Or for development:
    pip install -e .
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")
else:
    long_description = "ClusterCatcher: Single-cell sequencing analysis pipeline"

# Read version from __init__.py if it exists
version = "1.1.3"

setup(
    name='ClusterCatcher',
    version=version,
    description='Single-cell sequencing analysis pipeline for mutation signature detection and cell annotation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Jake Lehle',
    author_email='',
    url='https://github.com/JakeLehle/ClusterCatcher',
    license='GPL-3.0',
    
    # Package discovery
    packages=find_packages(exclude=['tests', 'tests.*', 'docs', 'docs.*']),
    
    # Include non-Python files
    include_package_data=True,
    package_data={
        '': ['*.yaml', '*.yml', '*.smk', '*.md'],
        'snakemake_wrapper': [
            'Snakefile',
            'scripts/*.py',
            'envs/*.yaml',
            'envs/*.yml',
            'rules/*.smk',
        ],
    },
    
    # Data files to install alongside package
    data_files=[
        ('snakemake_wrapper', ['snakemake_wrapper/Snakefile']),
        ('snakemake_wrapper/scripts', [
            'snakemake_wrapper/scripts/scomatic_mutation_calling.py',
            'snakemake_wrapper/scripts/signature_analysis.py',
            'snakemake_wrapper/scripts/SingleCellGenotype.py',
            'snakemake_wrapper/scripts/generate_summary.py',
        ]),
        ('snakemake_wrapper/envs', [
            'snakemake_wrapper/envs/scomatic.yaml',
            'snakemake_wrapper/envs/signatures.yaml',
        ]),
        ('snakemake_wrapper/rules', [
            'snakemake_wrapper/rules/scomatic_signatures.smk',
        ]),
    ],
    
    # CLI entry points
    entry_points={
        'console_scripts': [
            'ClusterCatcher=cli.cli:main',
            'clustercatcher=cli.cli:main',  # Lowercase alias
        ],
    },
    
    # Core dependencies (minimal set for CLI)
    # Full analysis dependencies are in environment.yml
    install_requires=[
        'click>=8.0',
        'pyyaml>=5.4',
        'pandas>=1.3',
        'numpy>=1.20',
        'rich>=10.0',  # For nice CLI output
    ],
    
    # Optional dependencies for different use cases
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov',
            'black',
            'flake8',
            'mypy',
            'pre-commit',
        ],
        'docs': [
            'sphinx>=4.0',
            'sphinx-rtd-theme',
            'sphinx-click',
            'myst-parser',
        ],
        'full': [
            # Full analysis stack (also available via conda)
            'snakemake>=7.0',
            'scanpy>=1.9',
            'anndata>=0.8',
            'scipy>=1.7',
            'scikit-learn>=1.0',
            'matplotlib>=3.5',
            'seaborn>=0.11',
            'tqdm',
            'h5py',
        ],
    },
    
    # Python version requirement
    python_requires='>=3.8,<3.12',
    
    # Project classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
    ],
    
    # Keywords for PyPI
    keywords=[
        'single-cell',
        'sequencing',
        'mutational-signatures',
        'cancer',
        'bioinformatics',
        'scRNA-seq',
        'cell-annotation',
        'SComatic',
        'COSMIC',
        'snakemake',
        'pipeline',
    ],
    
    # Project URLs
    project_urls={
        'Bug Reports': 'https://github.com/JakeLehle/ClusterCatcher/issues',
        'Source': 'https://github.com/JakeLehle/ClusterCatcher',
        'Documentation': 'https://github.com/JakeLehle/ClusterCatcher#readme',
    },
    
    # Zip safety
    zip_safe=False,
)

