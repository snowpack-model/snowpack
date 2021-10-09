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

#include <meteoio/gridResampling/GridTimeseriesResampling.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <meteoio/GridProcessor.h>
#include <meteoio/IOUtils.h>

#include <sstream>

namespace mio {

/**
 * @brief Constructor for a grid resampling algorithm.
 * @details On initialization, a resampling object stores its user settings.
 * @param[in] i_algoname The current algorithm's semantic name.
 * @param[in] i_parname The current meteo parameter's identifier.
 * @param[in] dflt_window_size The default grid resampling window size.
 * @param[in] vecArgs Vector of arguments (user settings) for this algorithm. Note that settings must
 * be given for this algorithm's name (e. g. TA::TIMESERIES::EXTRAPOLATE = T for TA::TIMESERIES::ALGORITHM = LINEAR).
 */
GridTimeseriesResampling::GridTimeseriesResampling(const std::string& i_algoname, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(i_algoname, i_parname, dflt_window_size, vecArgs), vecArgs_(vecArgs),
	base_algorithm_("LINEAR")
{
	for (size_t ii = 0; ii < vecArgs.size(); ++ii) {
		if (vecArgs[ii].first == "ALGORITHM")
			base_algorithm_ = vecArgs[ii].second;
	}
}

/**
 * @brief Print this algorithm's properties to a stream.
 * @return Semantic description of the algorithm's setup.
 */
std::string GridTimeseriesResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
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
void GridTimeseriesResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	//retrieve an algorithm from the time series resampling algorithm factory:
	resampled_grid.set(all_grids.begin()->second, IOUtils::nodata);
	ResamplingAlgorithms* ts_interpolator = ResamplingAlgorithmsFactory::getAlgorithm(base_algorithm_, parname, 86400., vecArgs_);

	StationData point_meta; //fill this object with what's available of grid point metadata

	for (size_t xx = 0; xx < resampled_grid.getNx(); ++xx) {
		for (size_t yy = 0; yy < resampled_grid.getNy(); ++yy) { //iterate over all grid points separately
			Coords point_coords;
			point_coords.setGridIndex((int)xx, (int)yy, IOUtils::inodata); //TODO: get altitude (e. g. for solar resampling)
			resampled_grid.gridify(point_coords); //calculate coordinates for grid point
			point_meta.setStationData(point_coords, "_GRS_ID_", "_GRID_RES_STAT_"); //some unique dummy ID
			std::vector<MeteoData> vecM; //extract time series to this vector

			ResamplingAlgorithms::ResamplingPosition pos = ResamplingAlgorithms::exact_match;
			size_t index = IOUtils::npos;
			size_t counter = 0;

			MeteoData resampled_pt; //point at which to resample
			for (auto it = all_grids.begin(); it != all_grids.end(); ++it) {
				MeteoData md( it->first, point_meta );
				md(parname) = it->second((int)xx, (int)yy);
				if (it->first > date && index == IOUtils::npos) { //put a nodata point at the date to be resampled
					resampled_pt = md; //copy meta data
					resampled_pt.reset();
					resampled_pt.setDate(date);
					vecM.push_back(resampled_pt);
					index = counter; //remember index of nodata point
				}
				vecM.push_back(md);
				counter++;
			}

			if (index == IOUtils::npos) { //requested date is after available time span
				resampled_pt.setDate(date);
				resampled_pt.meta = point_meta;
				vecM.push_back(resampled_pt);
				index = vecM.size() - 1;
				pos = ResamplingAlgorithms::end;
			} else if (index == 0) { //requested date was inserted before available time span
				pos = ResamplingAlgorithms::begin;
			}

			ts_interpolator->resample(point_meta.getHash(), index, pos,
				resampled_pt.getParameterIndex(parname), vecM, resampled_pt);
			resampled_grid(xx, yy) = resampled_pt(parname);
		} //endfor yy
	} //endfor xx

	delete ts_interpolator;
}

} //namespace
