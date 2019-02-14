#important!!! choose here how you want to setup the wrapper:

#howToSetup = "linkToMeteoioDLL" #with visual studio 2017 and python 3.7 this shoudl work, but it does not. see readme.txt
howToSetup = "buildEverything" #this works with "Microsoft Visual C++ Compiler for Python 2.7" and python 2.7. It seems to
                               # be the only way to get this running on windows, even if not very elegant. for more details
                               # see readme.txt.

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

#from distutils import sysconfig
import sysconfig
from Cython.Distutils import build_ext #for dll?

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
                  # extra_link_args=extra_compile_args,
                  #extra_compile_args = extra_compile_args
                  )
    ]

    setup(
        name = 'meteoioWrapper',
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
        ext_modules = cythonize(extensions),
        include_dirs = ["../../"]
    )
