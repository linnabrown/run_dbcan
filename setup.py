#!/usr/bin/env python
#
# distutils setup script for dbcan package
# also installs the bundled Hotpep scripts and data files

from glob import glob
from os import listdir
from os.path import isfile
from setuptools import setup, find_packages

long_description = """This is the standalone version of dbCAN annotation tool for automated CAZyme annotation (known as run_dbCAN.py), written by Tanner Yohe and Le Huang.
"""

setup(name='dbcan',
      version="3.0.1",
      description='Standalone version of dbCAN annotation tool for automated CAZyme annotation',
      long_description=long_description,
      author='Tanner Yohe, Le Huang and Qiwei Ge',
      author_email='lehuang@unc.edu',
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
      packages=find_packages(
        exclude=[
          "db",
          "db.*"
        ]
      ) + ['dbcan_cli'],
      include_package_data=True,
      scripts=[
          'dbcan_cli/hmmscan-parser.py'
      ],
      entry_points={
        "console_scripts":[
          "run_dbcan = dbcan_cli.run_dbcan:cli_main",
        ]
      },
      license='GPLv3',
      install_requires=[
          'natsort',
          'setuptools',
          'scipy',
          'psutil',
          'numpy'
      ],
      python_requires='>=3.5',
      zip_safe=False
     )
