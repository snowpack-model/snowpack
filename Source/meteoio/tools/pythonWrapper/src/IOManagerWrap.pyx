# distutils: language = c++


from IOManager cimport IOManager
from MeteoData cimport METEO_SET
from MeteoData cimport MeteoData

#import ConfigWrap ###muss nicht importiert werden, kennt er schon durch cimport Config im pxd-file
#import CoordsWrap
#import DateWrap

cdef class PyIOManager:
    cdef IOManager*c_iomanager	

    #def __cinit__(self, filename):
    #    self.c_iomanager = new IOManager(filename)
	
    def __cinit__(self, PyConfig config):
        self.c_iomanager = new IOManager(config.c_config)

    def __dealloc__(self):
        del self.c_iomanager
        
    def getMeteoData(self, PyDate date):
        print("collecting meteo data for this date: ")
        print(date.toString())
        cdef METEO_SET vecMeteo
        #getMeteoData aufrufen, dann vecMeteo in liste umwandeln und zurueckgeben
        nData = self.c_iomanager.getMeteoData(date.c_date, vecMeteo)
        #print(nData)
        results = []
        for item in vecMeteo:
            #print (item.toString())
            itemPy = PyMeteoData()
            itemPy.c_meteodata = item
            results.append(itemPy)
        return results

    def getMeteoDataRange(self, PyDate dateStart, PyDate dateEnd):
        print("collecting meteo data for a time range ")
        cdef vector[METEO_SET] vecVecMeteo
        #getMeteoData aufrufen, dann vecMeteo in liste umwandeln und zurueckgeben
        nData = self.c_iomanager.getMeteoData(dateStart.c_date, dateEnd.c_date, vecVecMeteo)
        results = []
        for vecMeteo in vecVecMeteo:
            resultsTemp = []
            for item in vecMeteo:
                itemPy = PyMeteoData()
                itemPy.c_meteodata = item
                resultsTemp.append(itemPy)
            results.append(resultsTemp)
        return results
        
    def writeMeteoData(self, PyVecVecMeteo, option):
        cdef vector[METEO_SET] vecVecMeteo
        cdef METEO_SET vecMeteo
        for PyVecMeteo in PyVecVecMeteo:
            vecMeteo.clear()
            for PyMeteo in PyVecMeteo: 
                h=<PyMeteoData> PyMeteo
                vecMeteo.push_back( h.c_meteodata )
            vecVecMeteo.push_back( vecMeteo )
        self.c_iomanager.writeMeteoData(vecVecMeteo, option)
        