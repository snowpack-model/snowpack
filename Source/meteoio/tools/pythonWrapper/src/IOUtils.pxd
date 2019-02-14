from libcpp.string cimport string
from libcpp cimport bool

# Declare the class with cdef
cdef extern from "meteoio/IOUtils.h" namespace "mio":
    string getLibVersion(const bool& short_version)
