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

#ifndef GRIDTIMESERIESRESAMPLING_H
#define GRIDTIMESERIESRESAMPLING_H

#include <string>
#include <utility>
#include <vector>

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

namespace mio {

/**
 * @class GridTimeseriesResampling
 * @brief This grid resampling class builds time series at all grid positions and sends them
 * to meteo 1d resampling algorithms to perform the interpolations.
 * @details You can use the algorithms available at \ref resampling like follows:
 * @code
 * [GridInterpolations1D]
 * TA::RESAMPLE = TIMESERIES
 * TA::TIMESERIES::ALGORITHM = LINEAR
 * TA::TIMESERIES::EXTRAPOLATE = TRUE
 * @endcode
 * @note Currently the algorithm has no knowledge of the used DEM (for solar resampling).
 * @author Michael Reisecker
 * @date 2021-09
 */
class GridTimeseriesResampling : public GridResamplingAlgorithm {
	public:
		GridTimeseriesResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid);
		std::string toString() const;

	private:
		std::vector< std::pair<std::string, std::string> > vecArgs_;
		std::string base_algorithm_; ///< Name of timeseries resampling algorithm to use
};

} //end namespace mio

#endif
