#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='calla',
      version='0.4',
      description='Open implementation of formulas in specifications or codes for structural engineering',
      url='https://github.com/warmwaver/calla',
      author='warmwaver',
      author_email='warmwaver@126.com',
      packages=['calla','calla.GB','calla.JTG','calla.CJJ','calla.TB','callex'],
      #packages=find_packages(exclude=['test','docs']),
      long_description="Open implementation of formulas in specifications or codes for structural engineering.",
      license="MIT",
      platforms=["any"],
      python_requires='>=3.10',
     )
