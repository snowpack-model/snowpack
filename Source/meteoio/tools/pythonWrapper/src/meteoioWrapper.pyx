# distutils: language = c++
# cython: language_level=3
# cython: c_string_enconding=utf-8

"""meteoioWrapper.pyx: This file collects all the single wrapper files for the meteoIO-wrapper. This file will be called
                        when the wrapper is built (see setup*.py)
   Author: Thiemo Theile
   Date created: 4/1/2019
"""

include "ConfigWrap.pyx"
include "CoordsWrap.pyx"
include "StationDataWrap.pyx"
include "MeteoDataWrap.pyx"
include "DateWrap.pyx"
include "IOManagerWrap.pyx"
include "IOUtilsWrap.pyx"