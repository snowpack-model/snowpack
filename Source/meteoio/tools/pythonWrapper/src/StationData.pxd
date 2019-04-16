"""StationData.pxd: This file wraps the class StationData (from meteoio/dataClasses/StationData.h).
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

from libcpp.string cimport string
from libcpp cimport bool

from Coords cimport Coords

cdef extern from "meteoio/dataClasses/StationData.h" namespace "mio":
    cdef cppclass StationData:

        StationData() except +
        string getStationID() const
        string getStationName() const
        Coords getPosition() const
        string getHash() const
        double getAltitude() const
        double getSlopeAngle() const
        double getAzimuth()

        void setStationData(const Coords& i_position, const string& i_id, const string& i_name)
        void setSlope(const double& in_slope_angle, const double& in_azimuth)

        const string toString() const
