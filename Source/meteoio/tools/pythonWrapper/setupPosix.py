"""
setupPosix.py: This script will build the python-wrapper for meteoIO on Linux or Mac. The script has to be
             executed with following arguments: 'python.exe setupPosix.py build_ext --inplace'

Author: Thiemo Theile
Date created: 4/1/2019
Python Version: 2.7 or 3.x
"""

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

import sysconfig
from Cython.Distutils import build_ext

from Cython.Build import cythonize
import os
import collectMeteoIOSourceFiles as fileCollector

#clean up (so that the wrapper has to be rebuild):
if os.path.exists("src/meteoioWrapper.cpp"):
    os.remove("src/meteoioWrapper.cpp")


ext_modules = [
    Extension('meteoioWrapper',
              ['src/meteoioWrapper.pyx'],
              language="c++",
              libraries=['meteoio'],
              library_dirs=['../../lib']
              )
]
setup(
    name='meteoioWrapper',
    version="0.0.1",
    author="Thiemo Theile",
    author_email="thiemotheile@gmx.de",
    description=("This python-wrapper uses Cython (https://cython.org/) to make "
                 "the meteoIO-library accessible from python-scripts."),
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    include_dirs=["../../"]
)
