
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector
#from libcpp cimport exception

###########################



###########################

# Declare the class with cdef
cdef extern from "meteoio/dataClasses/MeteoData.h" namespace "mio":
    cdef cppclass MeteoData:
        MeteoData() except +
        const string toString() const
        const string getStationID() const
        
    ctypedef vector[MeteoData] METEO_SET