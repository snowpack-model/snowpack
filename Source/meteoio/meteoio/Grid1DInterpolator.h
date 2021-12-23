// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 MobyGIS Srl, Trento, Italy                                      */
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

#ifndef GRID1DINTERPOLATOR_H
#define GRID1DINTERPOLATOR_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/Config.h>
#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

#include <string>
#include <vector>
#include <map>

namespace mio {

/**
 * @class Grid1DInterpolator
 * @brief A class to temporally resample grid objects.
 * @author Michael Reisecker
 * @date 2021-09
 */
class Grid1DInterpolator {
	public:
		Grid1DInterpolator(const Config& in_cfg);
		Grid1DInterpolator(const Grid1DInterpolator& org) = default;
		~Grid1DInterpolator();
		Grid1DInterpolator& operator=(const Grid1DInterpolator&);
		bool resampleData(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& available_grids, Grid2DObject& resampled_grid);
		double getWindowSize() const { return grid_window_size; }

	private:
		std::string getGridAlgorithmForParameter(const std::string& parname) const;

		static const std::string section_name;
		std::map<std::string, GridResamplingAlgorithm*> algorithm_map; //per parameter interpolation algorithms
		const Config& cfg;
		double grid_window_size = 86400.; ///< in seconds, default is 2 Julian days
		bool enable_grid_resampling = true; ///< easy way to turn grid resampling on/off
};

} //end namespace

#endif
