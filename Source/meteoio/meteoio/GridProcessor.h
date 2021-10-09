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

#ifndef GRIDPROCESSOR_H
#define GRIDPROCESSOR_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/Config.h>
#include <meteoio/Grid1DInterpolator.h>

#include <map>
#include <string>
#include <vector>

namespace mio {

/**
 * @class GridProcessor
 * @brief This class is handled by a GridManager and performs grid filtering and temporal resampling.
 * @author Michael Reisecker
 * @date 2021-09
 */
class GridProcessor {
	public:
		GridProcessor(const Config& cfg);
		bool resample(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid);
		static std::map<Date, Grid2DObject>::const_iterator seek(const Date& date, const std::map<Date, Grid2DObject>& grids, const bool& exact_match = false);
		static std::map<Date, Grid2DObject>::const_iterator seek_before(const Date& date, const std::map<Date, Grid2DObject>& grids);
		static std::map<Date, Grid2DObject>::const_iterator seek_after(const Date& date, const std::map<Date, Grid2DObject>& grids);
		double getWindowSize() const { return gi1d.getWindowSize(); };

	private:
		static std::set<std::string> getParameters(const Config& cfg);

		Grid1DInterpolator gi1d;
		bool enable_grid_filtering = false;
};

} //end namespace

#endif
