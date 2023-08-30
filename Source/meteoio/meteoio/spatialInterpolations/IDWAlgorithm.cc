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

#include <meteoio/spatialInterpolations/IDWAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/MathOptim.h>

namespace mio {

IDWAlgorithm::IDWAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm)
                         : InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), scale(1e3), alpha(1.)
{
	const std::string where( "Interpolations2D::"+i_param+"::"+i_algo );
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="SCALE") {
			IOUtils::parseArg(vecArgs[ii], where, scale);
		} else if (vecArgs[ii].first=="ALPHA") {
			IOUtils::parseArg(vecArgs[ii], where, alpha);
		}
	}
}

double IDWAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		return 0.3;
	} else if (nrOfMeasurments > 1){
		return 0.5;
	}

	return 0.2;
}

void IDWAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	const bool isWindSpeed = (param == "VW" || param == "VW_MAX" || param == "VW_DRIFT");
	const bool isWindDir = (param == "DW");

	if (isWindSpeed || isWindDir) {
		// Get VW and DW at stations that have both
		std::vector<double> vecDataVW;
		std::vector<double> vecDataDW;
		std::vector<double> vecDataU;
		std::vector<double> vecDataV;
		std::vector<StationData> vecMetaC;
		tsmanager.getMeteoData(date, vecMeteo);
		for (size_t ii=0; ii<vecMeteo.size(); ii++){
			size_t VW = MeteoData::VW;	// first guess for the wind speed parameter. In case of interpolating DW, use VW.
			if (param == "VW_DRIFT" || param == "VW_MAX") VW = vecMeteo[ii].getParameterIndex(param);	// override VW when other wind speed parameter is requested
			// Store only data from stations that have both wind speed and direction
			if ((vecMeteo[ii](VW) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::DW) != IOUtils::nodata)){
				vecDataVW.push_back(vecMeteo[ii](VW));
				vecDataDW.push_back(vecMeteo[ii](MeteoData::DW));
				vecMetaC.push_back(vecMeteo[ii].meta);
			}
		}

		// When no station has both VW and DW available, fall back to default IDW
		if (vecMetaC.size() == 0) {
			Interpol2D::IDW(vecData, vecMeta, dem, grid, scale, alpha);
			return;
		} else {
			info << "using wind speed components from " << vecMetaC.size() << " stations";
			for (size_t ii=0; ii<vecMetaC.size(); ii++) {
				vecDataU.push_back(IOUtils::VWDW_TO_U(vecDataVW[ii], vecDataDW[ii]));
				vecDataV.push_back(IOUtils::VWDW_TO_V(vecDataVW[ii], vecDataDW[ii]));
			}
		}

		// Apply IDW on both components individually
		Grid2DObject gridU(dem, 0.);
		Grid2DObject gridV(dem, 0.);
		Interpol2D::IDW(vecDataU, vecMeta, dem, gridU, scale, alpha);
		Interpol2D::IDW(vecDataV, vecMeta, dem, gridV, scale, alpha);

		//recompute VW, DW in each cell
		grid.set(dem, IOUtils::nodata);
		for (size_t ii=0; ii<grid.size(); ii++) {
			if (gridU(ii) != IOUtils::nodata && gridV(ii) != IOUtils::nodata) {
				if (isWindSpeed) {
					// Wind speed components
					grid(ii) = Optim::fastSqrt_Q3(gridU(ii) * gridU(ii) + gridV(ii) * gridV(ii));
				} else {
					// Wind direction
					grid(ii) = IOUtils::UV_TO_DW(gridU(ii), gridV(ii));
				}
			}
		}
	} else {
		Interpol2D::IDW(vecData, vecMeta, dem, grid, scale, alpha);
	}
}

} //namespace
