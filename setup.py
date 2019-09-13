#!/usr/bin/env python
#
# setuptools script for dbcan package
# also installs the bundled Hotpep scripts and data files

from glob import glob
from os import listdir
from os.path import isfile
from setuptools import setup, find_packages

long_description = """This is the standalone version of dbCAN annotation tool for automated CAZyme annotation (known as run_dbCAN.py), written by Le Huang and Tanner Yohe.
"""

setup(name='dbcan',
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      description='Standalone version of dbCAN annotation tool for automated CAZyme annotation',
      long_description=long_description,
      author='Le Huang and Tanner Yohe',
      author_email='505103350@qq.com',
      url='https://github.com/linnabrown/run_dbcan',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      packages=find_packages(),
      package_data={
          'Hotpep': ['*.txt'],
      },
      include_package_data=True,
      scripts=[
          'CGCFinder.py',
          'hmmscan-parser.py',
          'run_dbcan.py',
          'Hotpep/add_functions_orf.py',
          'Hotpep/bact_group_many_proteins_many_patterns.py',
          'Hotpep/list_multidomain_proteins.py',
          'Hotpep/parallel_group_many_proteins_many_patterns_noDNA.py',
          'Hotpep/train_many_organisms_many_families.py',
      ],
      license='GPLv3',
      install_requires=[
          'natsort',
          'setuptools',
      ],
      python_requires='>=3',
     )
