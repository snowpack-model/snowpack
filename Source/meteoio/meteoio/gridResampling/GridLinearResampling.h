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

#ifndef GRIDLINEARRESAMPLING_H
#define GRIDLINEARRESAMPLING_H

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

namespace mio {

/**
 * @class GridLinearResampling
 * @brief Point-wise linear interpolation between grids.
 * @details This algorithm performs linear interpolation between available grids.
 * The syntax is the same as for timeseries interpolations; for example:
 * @code
 * [GridInterpolations1D]
 * TA::RESAMPLE = LINEAR
 * @endcode
 * @author Michael Reisecker
 * @date 2021-09
 */
class GridLinearResampling : public GridResamplingAlgorithm {
	public:
		GridLinearResampling(const std::string& algoname, const std::string& i_parname, const double& dflt_window_size,
			const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid) override;
		std::string toString() const override;
};

} //end namespace mio

#endif
