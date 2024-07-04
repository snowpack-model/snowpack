// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2024 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef RADPROJ_H
#define RADPROJ_H

#include <meteoio/meteoFilters/ProcessingBlock.h>

#include <vector>
#include <string>

namespace mio {

/**
 * @class RadProj
 * @brief Project incoming short wave radiation from any slope (angle and azimuth) to any other slope.
 * @details
 * This filter is useful to generate horizontal radiation values from measurements performed parallel to the slope
 * or to generate slope-parallel values from horizontal measurements. It can also be used to correct for a sensor that tilted 
 * over time and therefore was not parallel anymore. It takes as arguments at least one set of angle/azimuth for the source and/or the 
 * destination sensor. It can actually project from any inclination to any other inclination. In case of self shading 
 * on the source sensor, the projected values will be set to nodata.
 * 
 * It takes the following arguments:
 *    - SRC_AZI and SRC_ANGLE: the inclination parameters (azimuth in degrees north, angle in degrees 
 *      from the horizontal) from the source sensor (default: horizontal);
 *    - DEST_AZI and DEST_ANGLE: the inclination parameters (azimuth in degrees north, angle in degrees 
 *      from the horizontal) from the destination sensor (default: horizontal);
 * 
 * @note The angles and azimuth that are relevant are the ones measured at the sensor, not necessarily the slope under the sensor!
 * 
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2024-06-14
 * @code
 * ISWR::filter1	= PROJECT_RADIATION
 * ISWR::arg1::DEST_SLOPE = 35
 * ISWR::arg1::DEST_AZI = 180
 * @endcode
 */

class RadProj : public ProcessingBlock {
	public:
		RadProj(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
        double src_azi, dest_azi, src_angle, dest_angle;
};

} //end namespace

#endif
