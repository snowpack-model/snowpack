// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/FilterRate.h>
#include <cmath>
#include <algorithm>

using namespace std;

namespace mio {

FilterRate::FilterRate(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
           : ProcessingBlock(vecArgs, name, cfg), min_rate_of_change(0.), max_rate_of_change(0.), methodParam(LEFT)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::both; //for the rest: default values
}

size_t FilterRate::findNextPoint(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& start_idx)
{
	for (size_t ii=start_idx; ii<vecM.size(); ii++) {
		if (vecM[ii](param) != IOUtils::nodata) return ii;
	}

	return IOUtils::npos;
}

double FilterRate::getRate(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& idx, const size_t& cmp_idx)
{
	static const double day2sec = (24.*3600.);

	const double curr_value = vecM[idx](param);
	if (curr_value == IOUtils::nodata) return IOUtils::nodata;
	const double curr_time = vecM[idx].date.getJulian();

	const double prev_value = vecM[cmp_idx](param);
	const double prev_time = vecM[cmp_idx].date.getJulian();

	//this work both for left and right, as signs cancel out between the values and the times
	const double local_rate = (curr_value-prev_value) / ((curr_time-prev_time+1e-12)*day2sec); //per seconds

	return local_rate;
}

bool FilterRate::filterOut(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& idx, const size_t& last_good, const size_t& next_good) const
{
	const Date dt_start(2012, 10, 29, 17, 0, 1.);
	const Date dt_end(2012, 10, 29, 18, 30, 1.);

	const double curr_value = vecM[idx](param);
	if (curr_value == IOUtils::nodata) return false; //this should have been tested before we called filterOut()

	if (methodParam == LEFT) {
		const double left_rate = getRate(vecM, param, idx, last_good);
		if (left_rate == IOUtils::nodata) return false;
		const bool filter_left = (left_rate>max_rate_of_change || left_rate<min_rate_of_change);
		return filter_left;

	} else if (methodParam == RIGHT) {
		if (next_good == IOUtils::npos) return false;

		const double right_rate = getRate(vecM, param, idx, next_good);
		if (right_rate == IOUtils::nodata) return false;
		const bool filter_right = (right_rate>max_rate_of_change || right_rate<min_rate_of_change);
		return filter_right;

	} else {
		const double left_rate = getRate(vecM, param, idx, last_good);
		if (left_rate == IOUtils::nodata) return false;
		const bool filter_left = (left_rate>max_rate_of_change || left_rate<min_rate_of_change);

		if (next_good == IOUtils::npos)
			return filter_left;

		const double right_rate = getRate(vecM, param, idx, next_good);
		//if (right_rate == IOUtils::nodata) return false; //not necessary as if curr_value==nodata, it would have been caught be left_rate
		const bool filter_right = (right_rate>max_rate_of_change || right_rate<min_rate_of_change);

		if (methodParam == LEFT_AND_RIGHT) {
			return (filter_left && filter_right);
		} else {
			return (filter_left || filter_right);
		}
	}

	return false;
}

void FilterRate::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	size_t next_good = IOUtils::npos; //next point after the current one that is not nodata

	//last point before the current one that is not nodata
	size_t last_good = findNextPoint(ovec, param, 0);

	if (last_good == IOUtils::npos) //can not find a good point to start
		return;

	for (size_t ii=(last_good+1); ii<ovec.size(); ii++) {
		double& curr_value = ovec[ii](param);
		if (curr_value == IOUtils::nodata) continue;

		//only update next_good when we we'll use it and we've reached it
		if (methodParam != LEFT && (next_good == IOUtils::npos || next_good<=ii))
			next_good = findNextPoint(ovec, param, ii+1);

		const bool filter_point = filterOut(ovec, param, ii, last_good, next_good);
		if (filter_point) {
			curr_value = IOUtils::nodata;
		} else {
			last_good = ii;
		}
	}
}

void FilterRate::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	bool has_max=false, has_min=false;

	for (const auto& arg : vecArgs) {
		if (arg.first=="MIN") {
			IOUtils::parseArg(arg, where, min_rate_of_change);
			has_min = true;
		} else if (arg.first=="MAX") {
			IOUtils::parseArg(arg, where, max_rate_of_change);
			has_max = true;
		} else if (arg.first=="METHOD") {
			const std::string type_str( IOUtils::strToUpper(arg.second) );
			if (type_str=="LEFT") methodParam = LEFT;
			else if (type_str=="RIGHT") methodParam = RIGHT;
			else if (type_str=="LEFT_AND_RIGHT") methodParam = LEFT_AND_RIGHT;
			else if (type_str=="LEFT_OR_RIGHT") methodParam = LEFT_OR_RIGHT;
			else
				throw InvalidArgumentException("Invalid type \""+arg.second+"\" for \""+where+"\". Please use \"LEFT\", \"RIGH\", \"LEFT_AND_RIGHT\" or \"LEFT_OR_RIGHT\".", AT);
		}
	}

	if (!has_max) throw InvalidArgumentException("Please provide a MAX value for "+where, AT);
	if (has_min && min_rate_of_change>max_rate_of_change) throw InvalidArgumentException("MIN can not be larger than MAX for "+where+". Usually MIN is set to a negative value to constrain the descent rate while max constraints the ascent rate.", AT);
	if (!has_min) min_rate_of_change = -max_rate_of_change;
}

}
