"""StationDataWrap.pyx: This file wraps the c++-class StationData (frommeteoio/dataClasses/StationData.h) to the
                        python-class PyStationData.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

# distutils: language = c++

from StationData cimport StationData

cdef class PyStationData:
    cdef StationData c_stationData

    def __cinit__(self):
        self.c_stationData = StationData()

    def getStationID(self):
        return self.c_stationData.getStationID()

    def getStationName(self):
        return self.c_stationData.getStationName()

    def getPosition(self):
        pyCoords = PyCoords()
        pyCoords.c_coords = self.c_stationData.getPosition()
        return pyCoords

    def getHash(self):
        return self.c_stationData.getHash()

    def getAltitude(self):
        return self.c_stationData.getAltitude()

    def getSlopeAngle(self):
        return self.c_stationData.getSlopeAngle()

    def getAzimuth(self):
        return self.c_stationData.getAzimuth()

    def setStationData(self, PyCoords i_position, str i_id, str i_name):
        self.c_stationData.setStationData(i_position.c_coords, i_id, i_name)

    def setSlope(self, double in_slope_angle, double in_azimuth):
        self.c_stationData.setSlope(in_slope_angle,in_azimuth)

    def toString(self):
        return self.c_stationData.toString()


