"""Coords.pxd: This file wraps the class Coords (from meteoio/dataClasses/Coords.h).
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "meteoio/dataClasses/Coords.h" namespace "mio":
    cdef cppclass Coords:
        #Keywords for selecting the toString formats
        ctypedef enum FORMATS:
            DEBUG, # As much information as possible, useful for debugging
            FULL, # Provide all the usually necessary information
            LATLON, # Simplified, lat/lon only
            XY, # Simplified cartesian, only easting/northing
            CARTESIAN # Compact representation only containing the X/Y and I/J coordinates

        Coords() except +
        void setLatLon(double, double, double, bool)
        double getLat()
        double getLon()
        double getAltitude()
        const string toString(const FORMATS& type) const
