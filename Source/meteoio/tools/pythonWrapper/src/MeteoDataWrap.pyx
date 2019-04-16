"""MeteoDataWrap.pyx: This file wraps the c++-class MeteoData (frommeteoio/dataClasses/MeteoData.h) to the python-class
                      PyMeteoData.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

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

    def getParameterName(self, parindex):
        return self.c_meteodata.getParameterName(parindex)

    def getStaticParameterIndex(self, parname):
        return self.c_meteodata.getStaticParameterIndex(parname)

    def setDate(self, PyDate in_date):
        self.c_meteodata.setDate(in_date.c_date)

    def addParameter(self, i_paramname):
        self.c_meteodata.addParameter(i_paramname)

    def param_exists(self, parname):
        return self.c_meteodata.param_exists(parname)

    def reset(self):
        self.c_meteodata.reset()

    def isResampled(self):
        return self.c_meteodata.isResampled()

    def setResampled(self, in_resampled):
        self.c_meteodata.setResampled(in_resampled)

    def standardizeNodata(self, plugin_nodata):
        self.c_meteodata.standardizeNodata(plugin_nodata)

    #you can get a certain variable this way: temperature = meteoData["TA"]
    def __getitem__(self, parname):
        return self.c_meteodata(parname)

    def getNameForParameter(self, parindex):
        return self.c_meteodata.getNameForParameter(parindex)

    def getParameterIndex(self, parname):
        return self.c_meteodata.getParameterIndex(parname)

    def getNrOfParameters(self):
        return self.c_meteodata.getNrOfParameters()

    def merge(self, PyMeteoData meteo2):
        self.c_meteodata.merge(meteo2.c_meteodata)

########## public member variables: ########################

    property date:
        def __get__(self):
            pyDate = PyDate()
            pyDate.c_date = self.c_meteodata.date
            return pyDate
        def __set__(self, PyDate pyDate):
            self.c_meteodata.date = pyDate.c_date

    property meta:
        def __get__(self):
            pyStationData = PyStationData()
            pyStationData.c_stationData = self.c_meteodata.meta
            return pyStationData
        def __set__(self, PyStationData pyStationData):
            self.c_meteodata.meta = pyStationData.c_stationData


    property nrOfParameters:                   #can use this to access C++ member nrOfParameters
        def __get__(self):
            return self.c_meteodata.nrOfParameters
        #def __set__(self, int i):
        #    self.c_meteodata.nrOfParameters = i



######### useful helper functions: #####################

def getTimeSeriesOfCertainParameter(parname, stationID, vecVecMeteoData):
    values=[]
    for vecMeteoData in vecVecMeteoData:
        for meteoData in vecMeteoData:
            if(meteoData.getStationID() == stationID):
                values.append(meteoData[parname])
    return values

