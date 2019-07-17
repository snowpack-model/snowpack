"""CoordsWrap.pyx: This file wraps the c++-class Coords (frommeteoio/dataClasses/Coords.h) to the python-class PyCoords.
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

# distutils: language = c++

from Coords cimport Coords

cdef class PyCoords:
    cdef Coords c_coords

    def stringToType(self, typeString):
        type=0
        if typeString=="DEBUG":
            type=0
        if typeString=="FULL":
            type=1
        if typeString=="LATLON":
            type=2
        if typeString=="XY":
            type=3
        if typeString=="CARTESIAN":
            type=4
        return type
	
    def __cinit__(self):
        self.c_coords = Coords()
		
    def setLatLon(self, in_lat, in_lon, in_alti, in_update):
        self.c_coords.setLatLon(in_lat,in_lon,in_alti,in_update)
			
    def getLat(self):
        return self.c_coords.getLat()

    def getLon(self):
        return self.c_coords.getLon()

    def toString(self,typeString="DEBUG"):
        type = self.stringToType(typeString)
        return self.c_coords.toString(type)


    