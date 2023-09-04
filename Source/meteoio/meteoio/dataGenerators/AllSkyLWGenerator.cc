// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013-2021 WSL Institute for Snow and Avalanche Research    SLF-DAVOS  */
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

#include <meteoio/dataGenerators/AllSkyLWGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Sun.h>

#include <string>
#include <utility>

namespace mio {

AllSkyLWGenerator::AllSkyLWGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ, const Config &i_cfg)
                   : TauCLDGenerator(vecArgs, i_algo, i_section, TZ, i_cfg), sun(), model(OMSTEDT)
{ 
	//TauCLDGenerator will do its own arguments parsing, then AllSkyLWGenerator
	//so make sure that we don't use here the same name as an argument to TauCLDGenerator!
	parse_args(vecArgs); 
}

void AllSkyLWGenerator::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( section+"::"+algo );
	bool has_type=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );

			if (user_algo=="CARMONA") model = CARMONA;
			else if (user_algo=="OMSTEDT") model = OMSTEDT;
			else if (user_algo=="KONZELMANN") model = KONZELMANN;
			else if (user_algo=="UNSWORTH") model = UNSWORTH;
			else if (user_algo=="CRAWFORD") {
				model = CRAWFORD;
				if (cloudiness_model==TauCLDGenerator::DEFAULT)
					cloudiness_model = TauCLDGenerator::CLF_CRAWFORD;
			} else if (user_algo=="LHOMME") {
				model = LHOMME;
				if (cloudiness_model==TauCLDGenerator::DEFAULT)
					cloudiness_model = TauCLDGenerator::CLF_LHOMME;
			} else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for "+where, AT);

			has_type = true;
		}
	}
	
	if (!has_type) throw InvalidArgumentException("Please provide a TYPE for "+where, AT);
}

bool AllSkyLWGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), TAU_CLD=md(MeteoData::TAU_CLD);
		const double CLD = (md.param_exists("CLD"))? md("CLD") : IOUtils::nodata;
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;
		double cloudiness = (TAU_CLD!=IOUtils::nodata)? Atmosphere::Kasten_cloudiness( TAU_CLD ) : IOUtils::nodata;
		
		if (CLD!=IOUtils::nodata) {
			//Synop sky obstructed from view -> fully cloudy
			if (CLD>9. || CLD<0.) throw InvalidArgumentException("Cloud cover CLD should be between 0 and 8!", AT);
			cloudiness = std::max(std::min(CLD/8., 1.), 0.1);
		}

		const std::string station_hash( md.meta.stationID + ":" + md.meta.stationName );
		const double julian_gmt = md.date.getJulian(true);
		bool cloudiness_from_cache = false;

		//try to get a cloudiness value
		if (cloudiness==IOUtils::nodata) {
			const double lat = md.meta.position.getLat();
			const double lon = md.meta.position.getLon();
			const double alt = md.meta.position.getAltitude();
			sun.setLatLon(lat, lon, alt);
			sun.setDate(julian_gmt, 0.);

			bool is_night;
			cloudiness = TauCLDGenerator::getCloudiness(md, sun, is_night);
			if (cloudiness==IOUtils::nodata && !is_night) return false;

			if (is_night) { //interpolate the cloudiness over the night
				const auto& cloudiness_point = last_cloudiness.find(station_hash); //we get a pair<julian_date, cloudiness>
				if (cloudiness_point==last_cloudiness.end()) return false;

				cloudiness_from_cache = true;
				const double last_cloudiness_julian = cloudiness_point->second.first;
				const double last_cloudiness_value = cloudiness_point->second.second;
				if ((julian_gmt - last_cloudiness_julian) < 1.) cloudiness = last_cloudiness_value;
				else return false;
			}
		}

		//save the last valid cloudiness
		if (!cloudiness_from_cache)
			last_cloudiness[station_hash] = std::pair<double,double>( julian_gmt, cloudiness );

		//run the ILWR parametrization
		if (model==LHOMME)
			value = Atmosphere::Lhomme_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, cloudiness);
		else if (model==CARMONA)
			value = Atmosphere::Carmona_ilwr(RH, TA, cloudiness);
		else if (model==OMSTEDT)
			value = Atmosphere::Omstedt_ilwr(RH, TA, cloudiness);
		else if (model==KONZELMANN)
			value = Atmosphere::Konzelmann_ilwr(RH, TA, cloudiness);
		else if (model==UNSWORTH)
			value = Atmosphere::Unsworth_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, cloudiness);
		else if (model==CRAWFORD) {
			int year, month, day;
			md.date.getDate(year, month, day);
			value = Atmosphere::Crawford_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, static_cast<unsigned char>(month), cloudiness);
		}
	}

	return true; //all missing values could be filled
}

bool AllSkyLWGenerator::create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool all_filled = true;
	for (size_t ii=ii_min; ii<ii_max; ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if (status==false) all_filled=false;
	}

	return all_filled;
}

} //namespace
