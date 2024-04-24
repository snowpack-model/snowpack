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

bool AllSkyLWGenerator::generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& /*vecMeteo*/)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;
		const double cloudiness = TauCLDGenerator::getCloudiness(md);
		if (cloudiness==IOUtils::nodata) return false;

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
		const bool status = generate(param, vecMeteo[ii], vecMeteo);
		if (status==false) all_filled=false;
	}

	return all_filled;
}

} //namespace
