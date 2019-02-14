from libcpp.string cimport string
from libcpp.vector cimport vector
#from libcpp cimport exception

from Config cimport Config
from Date cimport Date
from MeteoData cimport METEO_SET

# Declare the class with cdef
cdef extern from "meteoio/IOManager.h" namespace "mio":
    cdef cppclass IOManager:
        #IOManager() except +
        #IOManager(const string& filename_in) except +
        IOManager(const Config& i_cfg) except +
    
        size_t getMeteoData(const Date& i_date, METEO_SET& vecMeteo)

        size_t getMeteoData(const Date& dateStart, const Date& dateEnd, vector[ METEO_SET ]& vecVecMeteo)
        
        void writeMeteoData(const vector[ METEO_SET ]& vecMeteo, const string& option)
