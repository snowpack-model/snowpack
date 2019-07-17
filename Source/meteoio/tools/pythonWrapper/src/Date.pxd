"""Date.pxd: This file wraps the class Date (from meteoio/dataClasses/Date.h).
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

from libcpp.string cimport string
from libcpp cimport bool
#from libcpp cimport exception


# Declare the class with cdef
cdef extern from "meteoio/dataClasses/Date.h" namespace "mio":
    cdef cppclass Date:

        #Keywords for selecting the date formats (the subsecond resolution is dropped when not needed)
        ctypedef enum FORMATS:
            ISO, # ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS.sss (fields might be dropped, in the least to the most significant order)
            ISO_TZ, # ISO 8601 format (same as ISO) but with time zone specification
            ISO_Z, # ISO 8601 format, forcing GMT and Zulu (Z) timezone specification
            FULL, # ISO 8601 followed by the julian date (in parenthesis)
            NUM, # ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
            DIN, #DIN5008 format: DD.MM.YYYY HH:MM:SS.sss
            ISO_WEEK, # ISO 8601 week date: YYYY-Www-D (for example: 2014-W41-1)
            ISO_DATE # ISO 8601 date format without the time (ie YYYY-MM-DD)

        Date() except +
        Date(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& in_timezone, const bool& in_dst) except +
        void setDate(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& in_timezone, const bool& in_dst)
        const string toString() const
        const string toString(const FORMATS& type, const bool& gmt) const

        bool operator==(const Date&) const
        bool operator!=(const Date&) const
        bool operator<(const Date&) const
        bool operator<=(const Date&) const
        bool operator>(const Date&) const
        bool operator>=(const Date&) const

        #Intervals arithmetic
        #const Date operator+(const Date&) const
        #const Date operator-(const Date&) const
        const Date operator+(const double&) const
        const Date operator-(const double&) const
        const Date operator*(const double&) const
        const Date operator/(const double&) const


