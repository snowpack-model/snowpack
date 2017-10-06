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

#include <meteoio/dataGenerators/PPHASEGenerator.h>

namespace mio {

void PPhaseGenerator::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "generators::"+algo );
	bool has_type=false, has_snow=false, has_rain=false;
	double snow_thresh=273.15, rain_thresh=273.15; //to silence a warning

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );

			if (user_algo=="THRESH") model = THRESH;
			else if (user_algo=="RANGE") model = RANGE;
			else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for "+where+" generator", AT);

			has_type = true;
		} else if(vecArgs[ii].first=="SNOW") {
			IOUtils::parseArg(vecArgs[ii], where, snow_thresh);
			has_snow = true;
		} else if(vecArgs[ii].first=="RAIN") {
			IOUtils::parseArg(vecArgs[ii], where, rain_thresh);
			has_rain = true;
		}
	}

	if (!has_type) throw InvalidArgumentException("Please provide a TYPE for "+where, AT);
	if (model == THRESH) {
		if (!has_snow) throw InvalidArgumentException("Please provide a snow/rain threshold for "+where, AT);
		fixed_thresh = snow_thresh;
	}
	if (model == RANGE) {
		if (!has_snow || !has_rain) throw InvalidArgumentException("Please provide a a snow and a rain threshold for "+where, AT);
		if (snow_thresh==rain_thresh) throw InvalidArgumentException(where+" : the two provided threshold must be different", AT);
		if (snow_thresh>rain_thresh) std::swap(snow_thresh, rain_thresh);
		range_start = snow_thresh;
		range_norm = 1. / (rain_thresh-snow_thresh);
	}
}

bool PPhaseGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA);
		if (TA==IOUtils::nodata) return false;
		
		if (model==THRESH) {
			value = (TA>=fixed_thresh)? 1. : 0.;
		} else if (model==RANGE) {
			const double tmp_rainfraction = range_norm * (TA - range_start);
			value = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}

	return true; //all missing values could be filled
}

bool PPhaseGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;
	
	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			all_filled = false;
	}

	return all_filled;
}

} //namespace
