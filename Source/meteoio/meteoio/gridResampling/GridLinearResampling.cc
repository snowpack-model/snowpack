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

#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/GridProcessor.h>
#include <meteoio/IOUtils.h>

#include <sstream>
#include <iomanip>

namespace mio {

/**
 * @brief Constructor for a grid resampling algorithm.
 * @details On initialization, a resampling object stores its user settings.
 * @param[in] algoname The current algorithm's semantic name.
 * @param[in] i_parname The current meteo parameter's identifier.
 * @param[in] dflt_window_size The default grid resampling window size.
 * @param[in] vecArgs Vector of arguments (user settings) for this algorithm.
 */
GridLinearResampling::GridLinearResampling(const std::string& algoname, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(algoname, i_parname, dflt_window_size, vecArgs)
{
	//do nothing
}

/**
 * @brief Print this algorithm's properties to a stream.
 * @return Semantic description of the algorithm's setup.
 */
std::string GridLinearResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) <<
		algo << "[ window_size=" << grid_window_size << " ]";
	return ss.str();
}

/**
 * @brief Perform temporal grid resampling.
 * @details This function performs the interpolation routine and returns a grid resampled to
 * the desired date.
 * @param[in] date Date to resample the data to.
 * @param[in] all_grids List of all grids available to this resampling algorithm, as well as
 * their corresponding dates.
 * @param[out] resampled_grid The temporally resampled grid.
 */
void GridLinearResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	auto it_before( GridProcessor::seek_before(date, all_grids) );
	auto it_after( GridProcessor::seek_after(date, all_grids) );
	if (it_before == all_grids.end() || it_after == all_grids.end() || all_grids.size() == 1)
		throw IOException("Grids not loaded to cover linear interpolation date (is your buffer size big enough?)", AT);
	const Grid2DObject& grid_before = it_before->second;
	const Grid2DObject& grid_after = it_after->second;
	resampled_grid.set(grid_before, IOUtils::nodata);

	//solve:
	//(y - y1)/(x - x1) = (y2 - y1)/(x2 - x1)
	//==> y = y1 + (y2 - y1)/(x2 - x1) * (x - x1)
	const double x1 = it_before->first.getJulian();
	const double x2 = it_after->first.getJulian();
	const double xx = date.getJulian();
	if (x1 == x2)
		throw IOException("Equal start and end date for grid linear interpolation", AT);

	for (size_t jj = 0; jj < grid_before.size(); ++jj) {
		const double& y1 = grid_before(jj);
		const double& y2 = grid_after(jj);
		if ((y1 == IOUtils::nodata) || (y2 == IOUtils::nodata))
			continue; //already at nodata
		const double aa = (y2 - y1) / (x2 - x1);
		resampled_grid(jj) = y1 + aa * (xx - x1);
	} //endfor jj
}

} //namespace
