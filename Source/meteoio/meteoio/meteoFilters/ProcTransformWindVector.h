// SPDX-License-Identifier: LGPL-3.0-or-later
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

#include <meteoio/meteoFilters/ProcessingBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcTransformWindVector
 * @brief This filter reprojects wind direction, and/or wind speed components, between PROJ supported coordinate systems.
 * @details
 * The filter reprojects the wind direction and/or the U and V wind speed components defined by a PROJ supported coordinate system.
 * For example, when the source coordinate system is WGS84, U and V components indicate flow along latitudes and longitudes, respectively.
 * After using this filter to reproject to for example EPSG:3031 (Antarctic Polar Stereographic), U and V components will be directed
 * along easting and northing, respectively.
 *
 * The following arguments are supported:
 * - COORDPARAM_SRC:	provides the source coordinate system for the transformation. Provide EPSG code or a filename with a
 * projection string. Leave empty, or use "4326" to specify WGS84.
 * - COORDPARAM:	provides the target coordinate system for the transformation. Provide EPSG code or a filename with a
 * projection string. Leave empty, or use "4326" to specify WGS84.
 * - RACMO2:		Default: false. Set to true to properly deal with the rotated polar grid from RACMO2.
 *
 * @note
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
 * Example using wind direction, with source coordinate system described in the file racmo2.prj:
 * @code
 * DW::filter1			= TRANSFORMWINDVECTOR
 * DW::arg1::COORDPARAM_SRC	= racmo2.prj	; file contains: "-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=-161.0 +lon_0=193.0"
 * DW::arg1::COORDPARAM		= 3031		; Antarctic Polar Stereographic
 * @endcode
 * Example using wind speed components:
 * @code
 * U::filter1			= TRANSFORMWINDVECTOR
 * U::arg1::COORDPARAM_SRC		= 4326		; WGS84
 * U::arg1::COORDPARAM		= 21781		; CH1903 / LV03 -- Swiss CH1903 / LV03
 * @endcode
 */
#if defined(PROJ4)
	#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
	#include <proj_api.h>
#elif defined(PROJ)
	#include <proj.h>
#endif

class ProcTransformWindVector : public ProcessingBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
	public:
		ProcTransformWindVector(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		~ProcTransformWindVector();
		ProcTransformWindVector(const ProcTransformWindVector& c);
		ProcTransformWindVector& operator = (const ProcTransformWindVector& c);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
#if defined PROJ4 || defined PROJ
		static std::string findUComponent(const MeteoData& md);
		static std::string findVComponent(const MeteoData& md);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config &cfg);
		void initPROJ(void);
		void TransformCoord(const double& X_in, const double& Y_in, double& X_out, double& Y_out);

#if defined(PROJ4)
		projPJ pj_src, pj_dest;
#elif defined(PROJ)
		PJ_CONTEXT* pj_context;
		PJ* pj_trans;
#endif

		// The declarations below are needed when the copy constructor is called.
		std::vector< std::pair<std::string, std::string> > vecArgs_i;
		std::string name_i;
		Config cfg_i;
#endif
		std::string s_coordparam, t_coordparam;
		bool RACMO2;		// to properly deal with the rotated polar grids from RACMO2
};

} //end namespace

#endif
