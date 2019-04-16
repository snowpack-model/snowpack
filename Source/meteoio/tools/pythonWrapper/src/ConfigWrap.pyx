"""ConfigWrap.pyx: This file wraps the c++-class Config (frommeteoio/Config.h) to the python-class PyConfig.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

# distutils: language = c++

from Config cimport Config

cdef class PyConfig:
    cdef Config c_config

    def __cinit__(self, filename):
        self.c_config = Config(filename)

    def get(self, key, section, dflt):
        return self.c_config.get(key,section,dflt)