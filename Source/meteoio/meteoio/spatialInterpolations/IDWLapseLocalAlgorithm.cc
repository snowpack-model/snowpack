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

#include <meteoio/spatialInterpolations/IDWLapseLocalAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

LocalIDWLapseAlgorithm::LocalIDWLapseAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm)
                      : InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), trend(vecArgs, i_algo, i_param), scale(1e3), alpha(1.), nrOfNeighbors(0), MaxDistance(0.)
{
	const std::string where( "Interpolations2D::"+i_param+"::"+i_algo );
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="NEIGHBORS") {
			IOUtils::parseArg(vecArgs[ii], where, nrOfNeighbors);
		} else if (vecArgs[ii].first=="MAX_DISTANCE") {
			IOUtils::parseArg(vecArgs[ii], where, MaxDistance);
		} else if (vecArgs[ii].first=="SCALE") {
			IOUtils::parseArg(vecArgs[ii], where, scale);
		} else if (vecArgs[ii].first=="ALPHA") {
			IOUtils::parseArg(vecArgs[ii], where, alpha);
		}
	}

	if (nrOfNeighbors==0 && MaxDistance==0.) throw InvalidArgumentException("Please provide either the number of nearest neighbors (NEIGHBORS) or maximum distance (MAX_DISTANCE), or both, to use for "+where, AT);
	if (nrOfNeighbors==1) throw InvalidArgumentException("Please provide a value for NEIGHBORS larger than 1 for "+where, AT);
	if (MaxDistance<0) throw InvalidArgumentException("Please provide a positive value for MAX_DISTANCE to use for "+where, AT);
}

double LocalIDWLapseAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.7;
}

void LocalIDWLapseAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	Interpol2D::LocalLapseIDW(vecData, vecMeta, dem, nrOfNeighbors, MaxDistance, grid, scale, alpha);
	if (nrOfNeighbors>0) {
		info << "using nearest " << nrOfNeighbors << " neighbors";
		if (MaxDistance>0.) info << ", and ";
	}
	if (MaxDistance>0.) info << "maximum distance " << MaxDistance << " meters";
}

} //namespace
