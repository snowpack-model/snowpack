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

#include <meteoio/dataGenerators/MeteoIndexGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Sun.h>

namespace mio {

const double MeteoIndex::soil_albedo = .23; //grass
const double MeteoIndex::snow_albedo = .85; //snow
const double MeteoIndex::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

void MeteoIndex::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( section+"::"+algo );
	bool has_type=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );

			if (user_algo=="WINDCHILL") model = WINDCHILL;
			else if (user_algo=="HEATINDEX") model = HEATINDEX;
			else if (user_algo=="WET_BULB") model = WET_BULB;
			//else if (user_algo=="WBGT_INDEX") model = WBGT_INDEX; //curently not accurate enough
			else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for "+where, AT);

			has_type = true;
		}
	}

	if (!has_type) throw InvalidArgumentException("Please provide a TYPE for "+where, AT);
}

bool MeteoIndex::windChill(const size_t& param, MeteoData& md)
{
	const double TA = md(MeteoData::TA);
	const double VW = md(MeteoData::VW);
	
	if (TA!=IOUtils::nodata && VW!=IOUtils::nodata) {
		md(param) = Atmosphere::windChill(TA, VW);
		return true;
	} else
		return false;
}

bool MeteoIndex::heatIndex(const size_t& param, MeteoData& md)
{
	const double TA = md(MeteoData::TA);
	const double RH = md(MeteoData::RH);
	
	if (TA!=IOUtils::nodata && RH!=IOUtils::nodata) {
		md(param) = Atmosphere::heatIndex(TA, RH);
		return true;
	} else
		return false;
}

bool MeteoIndex::wetBulbTemperature(const size_t& param, MeteoData& md)
{
	const double TA = md(MeteoData::TA);
	const double RH = md(MeteoData::RH);
	const double alt = md.meta.getAltitude();
	
	if (TA!=IOUtils::nodata && RH!=IOUtils::nodata && alt!=IOUtils::nodata) {
		md(param) = Atmosphere::wetBulbTemperature(TA, RH, alt);
		return true;
	} else
		return false;
}

bool MeteoIndex::WBGT_index(const size_t& param, MeteoData& md)
{
	const double TA = md(MeteoData::TA);
	const double RH = md(MeteoData::RH);
	const double VW = md(MeteoData::VW);
	const double ISWR = md(MeteoData::ISWR);
	if (TA==IOUtils::nodata || RH==IOUtils::nodata || VW==IOUtils::nodata || ISWR==IOUtils::nodata) return false;
	
	const double lat = md.meta.position.getLat();
	const double lon = md.meta.position.getLon();
	const double alt = md.meta.position.getAltitude();
	if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;
	
	const double julian_gmt = md.date.getJulian(true);
	const double HS = md(MeteoData::HS);
	double albedo = 0.5;
	if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
		albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
	
	sun.setLatLon(lat, lon, alt);
	sun.setDate(julian_gmt, 0.);
	sun.calculateRadiation(TA, RH, albedo);
	const double Md = sun.getSplitting(ISWR);
	
	double azimuth, elevation;
	sun.position.getHorizontalCoordinates(azimuth, elevation);
	const double cos_Z = cos( (90. - elevation)*Cst::to_rad );
	
	md(param) = Atmosphere::WBGT_index(TA, RH, VW, (1.-Md)*ISWR, Md*ISWR, cos_Z, alt);;
	return true;
}

bool MeteoIndex::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		if (model==WINDCHILL) return windChill(param, md);
		if (model==HEATINDEX) return heatIndex(param, md);
		if (model==WET_BULB) return wetBulbTemperature(param, md);
		if (model==WBGT_INDEX) return WBGT_index(param, md);
	}

	return true; //all missing values could be filled
}

bool MeteoIndex::create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool all_filled = true;
	for (size_t ii=ii_min; ii<ii_max; ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if (status==false) all_filled=false;
	}

	return all_filled; //could all missing values be filled?
}

} //namespace
