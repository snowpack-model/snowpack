"""MeteoData.pxd: This file wraps the class MeteoData (from meteoio/dataClasses/MeteoData.h).
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector
#from libcpp cimport exception

from Date cimport Date
from StationData cimport StationData
###########################

# Declare the class with cdef
cdef extern from "meteoio/dataClasses/MeteoData.h" namespace "mio":
    cdef cppclass MeteoData:

        ctypedef enum MERGE_TYPE:
            STRICT_MERGE=0, # Station1 receives data from station2 only for common timestamps
            EXPAND_MERGE=1, # If station2 can provide some data before/after station1, this extra data is added to station1
            FULL_MERGE=2 # All timestamps from station2 are brought into station1 even if the timestamps don't match

        # anchor meteoparam this enum provides indexed access to meteorological fields
        ctypedef enum Parameters:
            firstparam=0,
            P=firstparam, # Air pressure
            TA, # Air temperature
            RH, # Relative humidity
            TSG, # Temperature of the ground surface
            TSS, # Temperature of the snow surface
            HS, # Height of snow
            VW, # Wind velocity
            DW, # Wind direction
            VW_MAX, # Maximum wind velocity
            RSWR, # Reflected short wave radiation
            ISWR, # Incoming short wave radiation
            ILWR, # Incoming long wave radiation (downwelling)
            TAU_CLD, # Cloud transmissivity or ISWR/ISWR_clear_sky
            PSUM, # Water equivalent of precipitations, either solid or liquid
            PSUM_PH, # Precipitation phase: between 0 (fully solid) and 1(fully liquid)
            lastparam=PSUM_PH

        MeteoData() except +
        const string toString() const
        const string getStationID() const

        const string& getParameterName(const size_t& parindex)
        size_t getStaticParameterIndex(const string& parname)

        void setDate(const Date& in_date)
        size_t addParameter(const string& i_paramname)
        bool param_exists(const string& parname) const
        void reset()

        bool isResampled() const
        void setResampled(const bool& in_resampled)

        void standardizeNodata(const double& plugin_nodata)

        #double& operator()(const size_t& parindex)
        #const double& operator()(const size_t& parindex) const
        double& operator()(const string& parname)
        #const double& operator()(const string& parname) const

        const string& getNameForParameter(const size_t& parindex) const
        size_t getParameterIndex(const string& parname) const
        size_t getNrOfParameters() const

        void merge(const MeteoData& meteo2) except +

        #static std::set<std::string> listAvailableParameters(const std::vector<MeteoData>& vecMeteo);


########## public member variables: #################
        Date date #Timestamp of the measurement
        StationData meta #The meta data of the measurement
        const size_t nrOfParameters

    ctypedef vector[MeteoData] METEO_SET