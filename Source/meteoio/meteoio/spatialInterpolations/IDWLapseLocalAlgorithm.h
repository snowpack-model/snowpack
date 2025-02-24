// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef LOCALIDWLAPSE_ALGORITHM_H
#define LOCALIDWLAPSE_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class LocalIDWLapseAlgorithm
 * @ingroup spatialization
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * @details
 * The closest n stations to each pixel are
 * used to compute the local lapse rate, allowing to project the contributions of these n stations to the
 * local pixel with an inverse distance weight. It therefore takes the following arguments:
 *  - NEIGHBORS: how many neighbouring stations should be used;
 *  - MAX_DISTANCE: maximum allowed distance (in meters) between the stations and grid points to interpolate to;
 *  - SCALE: this is a scaling parameter to smooth the IDW distribution. In effect, this is added to the distance in order
 * to move into the tail of the 1/d distribution (default: 1000m);
 *  - ALPHA: this is an exponent to the 1/d distribution (default: 1);
 *
 * Either NEIGHBORS or MAX_DISTANCE needs to be specified. When both are specified, both restrictions are used. In such cases,
 * only stations less than MAX_DISTANCE away from the grid point to interpolate to are used, up to a maximum of n stations.
 * This can lead to less than n stations to be included in the interpolation.
 *
 * @note Beware, this method sometimes produces very sharp transitions
 * as it spatially moves from one station's area of influence to another one!
 *
 * @code
 * TA::algorithms            = LIDW_LAPSE
 * TA::lidw_lapse::neighbors = 6
 * @endcode
 */
class LocalIDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		LocalIDWLapseAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm);
		virtual double getQualityRating(const Date& i_date) override;
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid) override;
	private:
		Trend trend;
		double scale, alpha; ///<a scale parameter to smooth out the 1/dist and an exponent
		size_t nrOfNeighbors;
		double MaxDistance;
};

} //end namespace mio

#endif
