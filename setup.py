#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='calla',
      version='0.1',
      description='Calculation for structural engineers',
      url='https://github.com/warmwaver/calla',
      author='warmwaver',
      author_email='warmwaver@126.com',
      packages=find_packages(exclude=['test']),
      long_description="Open implementation of formulas in specifications or codes for structural engineering.",
      license="MIT",
      platforms=["any"],
     )
