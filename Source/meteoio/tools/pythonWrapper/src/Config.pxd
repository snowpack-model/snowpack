from libcpp.string cimport string
from libcpp cimport bool
#from libcpp cimport exception

# Declare the class with cdef
cdef extern from "meteoio/Config.h" namespace "mio":
    cdef cppclass Config:
        Config() except +
        Config(string& filename_in) except +
        string get(const string& key, const string& section, const string& dflt) const