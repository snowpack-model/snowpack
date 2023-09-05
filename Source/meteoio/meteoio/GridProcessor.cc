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

#include <meteoio/GridProcessor.h>

#include <algorithm>

namespace mio {

/**
 * @brief Constructor for a grid processor.
 * @details A GridProcessor object manages grid filtering and grid resampling; the latter
 * through a Grid1dInterpolator.
 * @param[in] cfg The current simulation's configuration.
 */
GridProcessor::GridProcessor(const Config& cfg) : gi1d(cfg)
{
	cfg.getValue("ENABLE_GRID_FILTERING", "GridFilters", enable_grid_filtering, IOUtils::nothrow);
	if (enable_grid_filtering)
		std::cout << "[W] Grid filtering is not implemented yet." << std::endl;
}

/**
 * @brief This function forwards interpolation requests.
 * @details Delegate a resampling request to the Grid1dInterpolator (which will in turn forward
 * the call to the actual interpolation routine).
 * @param[in] date Date to resample to.
 * @param[in] parameter Meteo parameter to resample.
 * @param[in] all_grids A list of all grids to use for resampling, including the corresponding dates.
 * @param[out] resampled_grid Grid filled with resampled data.
 */
bool GridProcessor::resample(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	return gi1d.resampleData(date, parameter, all_grids, resampled_grid);
}

/**
 * @brief Search a list of grids for a specific date.
 * @details This function takes a map of date/grid combinations and searches it for a specific date.
 * @param[in] date Date to look for.
 * @param[in] grids A map with the dates as keys, and corresponding grids as values.
 * @param[in] exact_match Does the date have to match exactly?
 * @return Iterator to the found entry (or the grid map end if n/a). If the match does not have to
 * be exact, the grid found right before will be returned (as iterator).
 */
std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek(const Date& date, const std::map<Date, Grid2DObject>& grids, const bool& exact_match)
{
	if (grids.empty())
		return grids.end();
	if (grids.size() == 1)
		return grids.begin(); //if there is only 1 grid, it is closest to any date...
	std::map<Date, Grid2DObject>::const_iterator it = grids.find(date);
	if (it != grids.end() || exact_match)
		return it;
	for (it = ++grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return --it;
	}
	return --grids.end(); //if the date is after all grids, the last one is closest
}

/**
 * @brief Find a grid that is available right after a specific date.
 * @details This function takes a map of date/grid combinations and searches it for the first grid
 * that comes after a specified date.
 * @param[in] date Date to look for.
 * @param[in] grids A map with the dates as keys, and corresponding grids as values.
 * @return Iterator to the found entry (or the grid map end if n/a).
 */
std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek_before(const Date& date, const std::map<Date, Grid2DObject>& grids)
{
	if (grids.empty())
		return grids.end();
	if (date < grids.begin()->first) //there is no element before date
		return grids.end();
	if (grids.size() == 1)
		return grids.begin(); //if there is only 1 grid, it is closest to any date...
	for (auto it = ++grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return --it;
	}
	return grids.end();
}

/**
 * @brief Find a grid that is available before a specific date.
 * @details This function takes a map of date/grid combinations and searches it for the grid
 * that comes right before a specified date.
 * @param[in] date Date to look for.
 * @param[in] grids A map with the dates as keys, and corresponding grids as values.
 * @return Iterator to the found entry (or the grid map end if n/a).
 */
std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek_after(const Date& date, const std::map<Date, Grid2DObject>& grids)
{
	if (grids.empty())
		return grids.end();
	for (auto it = grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return it;
	}
	return grids.end();
}

} //namespace
