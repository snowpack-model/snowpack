"""IOUtils.pxd: This file wraps the function getLibVersion (from meteoio/IOUtils.h).
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

from libcpp.string cimport string
from libcpp cimport bool

# Declare the class with cdef
cdef extern from "meteoio/IOUtils.h" namespace "mio":
    string getLibVersion(const bool& short_version)
