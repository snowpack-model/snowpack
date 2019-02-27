# distutils: language = c++

from MeteoData cimport MeteoData
from libcpp.vector cimport vector

cdef class PyMeteoData:
    cdef MeteoData c_meteodata

    def __cinit__(self):
        self.c_meteodata = MeteoData()

    def toString(self):
        return self.c_meteodata.toString()
        
    def getStationID(self):
        return self.c_meteodata.getStationID();
        
        
#ctypedef list(PyMeteoData) PY_METEO_SET
		    