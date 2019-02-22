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
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    include_dirs=["../../"]
)
