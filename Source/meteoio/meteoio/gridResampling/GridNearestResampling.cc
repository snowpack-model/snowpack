// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2025 SLF, Davos, Switzerland                                         */
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

#include <meteoio/gridResampling/GridNearestResampling.h>
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
 * @param[in] dflt_max_gap_size The default grid resampling window size.
 * @param[in] vecArgs Vector of arguments (user settings) for this algorithm.
 */
GridNearestResampling::GridNearestResampling(const std::string& algoname, const std::string& i_parname,
	const double& dflt_max_gap_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(algoname, i_parname, dflt_max_gap_size, vecArgs) {}

/**
 * @brief Print this algorithm's properties to a stream.
 * @return Semantic description of the algorithm's setup.
 */
std::string GridNearestResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) <<
		algo << "[ window_size=" << max_gap_size << " ]";
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
void GridNearestResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	auto it_before( GridProcessor::seek_before(date, all_grids) );
	auto it_after( GridProcessor::seek_after(date, all_grids) );
	if (it_before == all_grids.end() || it_after == all_grids.end() || all_grids.size() == 1)
		throw IOException("Grids not loaded to cover linear interpolation date (is your buffer size big enough?)", AT);
	const Grid2DObject& grid_before = it_before->second;
	const Grid2DObject& grid_after = it_after->second;
	resampled_grid.clear();

	const double x1 = it_before->first.getJulian();
	const double x2 = it_after->first.getJulian();
	const double xx = date.getJulian();
	if (x1 == x2) throw IOException("Equal start and end date for nearest neighbour interpolation", AT);
	if (x2-x1 > max_gap_size) return;

	resampled_grid.set(grid_before, IOUtils::nodata);
	for (size_t jj = 0; jj < grid_before.size(); ++jj) {
		const double diff1 = xx - x1;
		const double diff2 = x2 - xx;
		const double y1 = grid_before(jj);
		const double y2 = grid_after(jj);
		
		if (y1 == IOUtils::nodata) {
			resampled_grid(jj) = y2; //if y2 is also nodata, then we keep nodata
			continue;
		}
		if (y2 == IOUtils::nodata) {
			resampled_grid(jj) = y1;
			continue;
		}
		
		if (IOUtils::checkEpsilonEquality(diff1, diff2, 0.1/1440.)) { //within 6 seconds
			resampled_grid(jj) = 0.5*(y1+y2);
		} else if (diff1<diff2) {
			resampled_grid(jj) = y1;
		} else {
			resampled_grid(jj) = y2;
		}
	} //endfor jj
}

} //namespace
