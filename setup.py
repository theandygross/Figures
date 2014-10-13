#!/usr/bin/env python

from distutils.core import setup

setup(name='PandasFigureWrappers',
      version='0.1',
      description='Wrappers for Scipy and R stats using Pandas objects.',
      author='Andrew Gross',
      author_email='the.andrew.gross@gmail.com',
      url='http://andy-gross.flavors.me',
      package_dir = {'': 'src'},
      packages=['Figures'],
     )