from setuptools import setup
import os, re

long_description = """
Vtnurbs is a geometry constructor and class that is specifically developed for idealized ventricular anatomy. 
It employs Non-Uniform Rational B-Splines (NURBS) to construct a multipatch geometry of either a left or a 
bi-ventricle geometry. The geometry is used as a basis for the CardIGA module (IGA-based cardiac model) and
the multipatch fitting module.
"""
with open(os.path.join('vtnurbs', '__init__.py')) as f:
  version = next(filter(None, map(re.compile("^version = '([a-zA-Z0-9.]+)'$").match, f))).group(1)

setup(name     = 'vtnurbs',
      version  = '2.0',
      author   = 'Robin Willems',
      packages = ['vtnurbs'],
      description      = 'NURBS-based ventricle template geometry generator',
      download_url     = 'https://github.com/CardiacIGA/multipatch-ventricle/tree/main',
      long_description = long_description,
      zip_safe = False)