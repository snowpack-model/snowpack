/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/spatialInterpolations/IDWSlopesAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {
const double IDWSlopesAlgorithm::min_slope = 10.;
const double IDWSlopesAlgorithm::max_slope = 38.;

IDWSlopesAlgorithm::IDWSlopesAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm)
                   : InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), 
                     vecFlat(), vecNorth(), vecEast(), vecSouth(), vecWest(),
                     metaFlat(), metaNorth(), metaEast(), metaSouth(), metaWest(),
                     trend(vecArgs, i_algo, i_param), scale(1e3), alpha(1.)
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

double IDWSlopesAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date;
	tsmanager.getMeteoData(date, vecMeteo);
	
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		const double val = vecMeteo[ii](param);
		if (val==IOUtils::nodata) continue;

		const double slope = vecMeteo[ii].meta.getSlopeAngle();
		if (slope==IOUtils::nodata) continue;
		if (slope<=min_slope) {
			vecFlat.push_back( val );
			metaFlat.push_back( vecMeteo[ii].meta );
			continue;
		}

		const double azimuth = vecMeteo[ii].meta.getAzimuth();
		if (azimuth==IOUtils::nodata) continue;
		if (azimuth<45. || azimuth>315.) {
			vecNorth.push_back( val );
			metaNorth.push_back( vecMeteo[ii].meta );
		} else if (azimuth<135) {
			vecEast.push_back( val );
			metaEast.push_back( vecMeteo[ii].meta );
		} else if (azimuth<225) {
			vecSouth.push_back( val );
			metaSouth.push_back( vecMeteo[ii].meta );
		} else {
			vecWest.push_back( val );
			metaWest.push_back( vecMeteo[ii].meta );
		}
	}

	nrOfMeasurments = vecFlat.size();
	static const size_t minNrStations = 3;
	if (vecFlat.size()<minNrStations || vecNorth.size()<minNrStations || vecEast.size()<minNrStations || vecSouth.size()<minNrStations || vecWest.size()<minNrStations) return 0;

	return 0.8;
}

void IDWSlopesAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	info << vecNorth.size() << " north, " << vecEast.size() << " east, " << vecSouth.size() << " south, " << vecWest.size() << " west stations";

	trend.detrend(metaFlat, vecFlat);
	Interpol2D::IDW(vecFlat, metaFlat, dem, grid, scale, alpha);
	trend.retrend(dem, grid);

	Grid2DObject north;
	trend.detrend(metaNorth, vecNorth);
	Interpol2D::IDW(vecNorth, metaNorth, dem, north, scale, alpha);
	trend.retrend(dem, north);

	Grid2DObject east;
	trend.detrend(metaEast, vecEast);
	Interpol2D::IDW(vecEast, metaEast, dem, east, scale, alpha);
	trend.retrend(dem, east);

	Grid2DObject south;
	trend.detrend(metaSouth, vecSouth);
	Interpol2D::IDW(vecSouth, metaSouth, dem, south, scale, alpha);
	trend.retrend(dem, south);

	Grid2DObject west;
	trend.detrend(metaWest, vecWest);
	Interpol2D::IDW(vecWest, metaWest, dem, west, scale, alpha);
	trend.retrend(dem, west);

	//now we merge the grids together as a weighted average of the aspects and slope angle
	for (size_t ii=0; ii<dem.size(); ++ii) {
		if (dem(ii)==IOUtils::nodata) continue;

		const double slope = dem.slope(ii);
		const double azi = dem.azi(ii);
		if (slope==IOUtils::nodata || azi==IOUtils::nodata) continue;

		if (slope<=min_slope) continue; //we keep the flat value
		const double w_flat = (slope>=max_slope)? 0. : (max_slope - slope)/(max_slope - min_slope);
		
		if (azi<90) {
			const double w_azi = 1. - azi / 90.;
			grid(ii) = w_flat*grid(ii) + (1.-w_flat)*(w_azi*north(ii) + (1.-w_azi)*east(ii));
		} else if (azi<180) {
			const double w_azi = 1. - fabs(azi - 90.) / 90.;
			grid(ii) = w_flat*grid(ii) + (1.-w_flat)*(w_azi*east(ii) + (1.-w_azi)*south(ii));
		} else if (azi<270) {
			const double w_azi = 1. - fabs(azi - 180.) / 90.;
			grid(ii) = w_flat*grid(ii) + (1.-w_flat)*(w_azi*south(ii) + (1.-w_azi)*west(ii));
		} else {
			const double w_azi = 1. - fabs(azi - 270.) / 90.;
			grid(ii) = w_flat*grid(ii) + (1.-w_flat)*(w_azi*west(ii) + (1.-w_azi)*north(ii));
		}
		
	}

}

} //namespace
