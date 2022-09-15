#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from pathlib import Path


this_directory = Path(Path(__file__).parent)
with open(this_directory / 'README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

DESCRIPTION = 'Multiple HMM - search via synteny structures in Python'
LONG_DESCRIPTION = 'Multiple HMM - search via synteny structures in Python'
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'
NAME = 'pynteny'
AUTHOR = "Semidán Robaina Estévez"
AUTHOR_EMAIL = "srobaina@ull.edu.es"
MAINTAINER = "Semidán Robaina Estévez"
MAINTAINER_EMAIL = "srobaina@gmail.com"
DOWNLOAD_URL = 'http://github.com/robaina/Pynteny'
LICENSE = 'Creative Commons Attribution 4.0 International'
VERSION = '0.0.2'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=find_packages(),
      include_package_data=True,
      install_requires=['biopython', 'pyfastx', 'pandas', 'wget'],
      entry_points ={
            'console_scripts': [
                'pynteny = pynteny.cli:main'
            ]
        }
      )