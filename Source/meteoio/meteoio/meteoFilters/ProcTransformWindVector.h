/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef PROCTRANSFORMWINDVECTOR_H
#define PROCTRANSFORMWINDVECTOR_H

//#include <meteoio/meteoFilters/WindowedFilter.h> //use this one for filters relying on a data window, for example std_dev
#include <meteoio/meteoFilters/ProcessingBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcTransformWindVector
 * @brief This filter projects wind direction, and/or wind speed components, from WGS84 to a PROJ4 supported coordinate system,
 * defined by an EPSG code (requires PROJ4).
 * @details
 * The filter assumes the wind direction and/or the U and V wind speed components are defined in WGS84 (i.e, north/south is parallel
 * to longitude, and east/west is parallel to latitude). After transformation to a PROJ4 supported coordinate system defined by an
 * EPSG code, the wind direction is defined such that north/south is parallel to northing, and east/west is parallel to easting.
 * The following arguments are supported:
 * - COORDPARAM: provides the target EPSG code for transformation
 *
 * @note
 * - If no COORDPARAM is provided, it is tried to read COORDPARAM from the [Input] section.
 * - If both DW and wind speed components are present and defined, all three variables will be recalculated, even when the filter is
 * only set to act on wind direction or wind speed components. This ensures consistency.
 * - When applying the filter on wind speed components, the filter needs to be specified for only one component (the other one is
 * automatically recalculated, to maintain consistency).
 * - Wind speed components should be called U, VW_U or WIND_U for east/west and V, VW_V or WIND_V for north/south, respectively.
 * - At the North and South Pole, the transform is undefined. In the limit towards the pole (|latitude| > 89.999), the transform becomes
 * inaccurate.
 * 
 * @ingroup processing
 * @author Nander Wever
 * @date   2020-06-10
 *
 * Example using wind direction:
 * @code
 * DW::filter1		= TRANSFORMWINDVECTOR
 * DW::arg1::COORDPARAM	= 3031		; Antarctic Polar Stereographic
 * @endcode
 * Example using wind speed components:
 * @code
 * U::filter1		= TRANSFORMWINDVECTOR
 * U::arg1::COORDPARAM	= 21781		; CH1903 / LV03 -- Swiss CH1903 / LV03
 * @endcode
 */

class ProcTransformWindVector : public ProcessingBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
//class TEMPLATE : public WindowedFilter { //use this one for filters relying on a data window, for example std_dev
	public:
		ProcTransformWindVector(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
#ifdef PROJ4
		static std::string findUComponent(const MeteoData& md);
		static std::string findVComponent(const MeteoData& md);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config &cfg);
#endif
		std::string t_coordparam;
};

} //end namespace

#endif
