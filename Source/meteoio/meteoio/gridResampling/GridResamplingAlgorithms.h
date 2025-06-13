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

#ifndef GRIDRESAMPLINGALGORITHM_H
#define GRIDRESAMPLINGALGORITHM_H

#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <vector>

namespace mio {

/**
 * @class GridResamplingAlgorithm
 * @brief Interface class for grid resampling algorithms.
 * @details This class provides generic functionality to temporal grid resampling algorithms
 * (which will inherit from this class).
 */
class GridResamplingAlgorithm {

	public:
		GridResamplingAlgorithm(const std::string& algorithm, const std::string& i_parname, const double& dflt_max_gap_size, const std::vector< std::pair<std::string, std::string> >& /*vecArgs*/);
		virtual ~GridResamplingAlgorithm() = default;
		
		void setMaxGapSize(const double& i_max_gap_size);
		virtual void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid) = 0;
		virtual std::string toString() const = 0;

	protected:
		const std::string algo, parname;
		double max_gap_size;
};

/**
 * @class GridResamplingAlgorithmsFactory
 * @brief Object factory for temporal grid resampling algorithms.
 */
class GridResamplingAlgorithmsFactory {
	public:
		static GridResamplingAlgorithm* getAlgorithm(const std::string& i_algorithm, const std::string& parname,
			const double& max_gap_size, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg);
};

} //end namespace

#endif
