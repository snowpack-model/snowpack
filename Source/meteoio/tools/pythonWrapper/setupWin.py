"""
setupWin.py: This script will build the python-wrapper for meteoIO on Windows. The script has to be
             executed with following arguments: 'python.exe setupWin.py build_ext --inplace'

Important!!! With the variable howToSetup you can choose how you want to setup the wrapper:
                "linkToMeteoioDLL" or "buildEverything"

Author: Thiemo Theile
Date created: 4/1/2019
Python Version: 2.7 or 3.x
"""


#howToSetup = "linkToMeteoioDLL" #with visual studio 2017 and python 3.7 this should work, but it does not. see readme.txt
howToSetup = "buildEverything" #this works with "Microsoft Visual C++ Compiler for Python 2.7" and python 2.7. It seems to
                               # be the only way to get this running on windows, even if not very elegant. for more details
                               # see readme.txt.

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext

from Cython.Build import cythonize
import os
import collectMeteoIOSourceFiles as fileCollector

#clean up (so that the wrapper has to be rebuild):
if os.path.exists("src/meteoioWrapper.cpp"):
    os.remove("src/meteoioWrapper.cpp")

if(howToSetup=="linkToMeteoioDLL"):
    ext_modules = [
        Extension('meteoioWrapper',
                  ['src/meteoioWrapper.pyx'],
                  language="c++",
                  extra_objects=["libmeteoio.lib"],
                  libraries=['libmeteoio'],
                  library_dirs=['../../lib/RelWithDebInfo'] #adjust this path if necessary!!!
                  )
    ]

    setup(
        name="meteoioWrapper",
        version="0.0.1",
        author="Thiemo Theile",
        author_email="thiemotheile@gmx.de",
        description=("This python-wrapper uses Cython (https://cython.org/) to make "
                     "the meteoIO-library accessible from python-scripts."),
        cmdclass = {'build_ext':build_ext},
        ext_modules = ext_modules,
        include_dirs=["../../"]
    )


if(howToSetup=="buildEverything"):
    sourcefiles = ["src/meteoioWrapper.pyx"]+fileCollector.getListOfSourcefiles()

    extensions = [Extension(
        "meteoioWrapper",
        sourcefiles
        #extra_compile_args=["/EHsc"]
    )]

    setup(
        name="meteoioWrapper",
        version="0.0.1",
        author="Thiemo Theile",
        author_email="thiemotheile@gmx.de",
        description=("This python-wrapper uses Cython (https://cython.org/) to make "
                     "the meteoIO-library accessible from python-scripts."),
        ext_modules = cythonize(extensions),
        include_dirs = ["../../"]
    )