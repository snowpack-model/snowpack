# Description: Setup file for cy_func_loop.pyx
#
# Building: python firn_utilities_cy_setup.py build_ext --inplace
#
# Author: Christian Steger, April 2016

# Load modules
from distutils.core import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension # for extensions

ext_modules = [Extension("firn_utilities_cy",
               ["firn_utilities_cy.pyx"],
               libraries = ["m"],
               extra_compile_args = ["-O3", "-ffast-math"])]

setup(name = "firn_utilities_cy",
      cmdclass = {"build_ext": build_ext},
      ext_modules = ext_modules)