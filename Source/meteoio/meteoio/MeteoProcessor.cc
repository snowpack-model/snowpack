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
#include <cstddef>
#include <meteoio/MeteoProcessor.h>
#include <meteoio/meteoFilters/TimeFilters.h>
#include <algorithm>
#include <regex>

using namespace std;

namespace mio {

MeteoProcessor::MeteoProcessor(const Config& cfg, const char& rank, const IOUtils::OperationMode &mode) : mi1d(cfg, rank, mode), processing_stack(), enable_meteo_filtering(true)
{
	//ENABLE_METEO_FILTERING is documented in meteoFilters/ProcessingBlock.cc
	cfg.getValue("ENABLE_METEO_FILTERING", "Filters", enable_meteo_filtering, IOUtils::nothrow);
	
	//Parse [Filters] section, create processing stack for each configured parameter
	const std::set<std::string> set_of_used_parameters( getParameters(cfg) );

	for (std::set<std::string>::const_iterator it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		ProcessingStack* tmp = new ProcessingStack(cfg, *it);
		processing_stack[*it] = tmp;
	}
}

MeteoProcessor::~MeteoProcessor()
{
	//clean up heap memory, processing_stack is std::map<std::string, ProcessingStack*>
	//call delete on the Stack of each parameter
	for (auto const& paramStack : processing_stack) delete paramStack.second;
}

std::set<std::string> MeteoProcessor::getParameters(const Config& cfg)
{
	//regex for selecting the parameter out of the key
	static const std::regex param_regex("^([^:]+)(::?)[^:]+$", std::regex::optimize); //the second group is in order to detect using ':' instead of '::'
	std::smatch param_matches;
	const std::vector<std::string> vec_keys( cfg.getKeys(std::string(), "Filters") );

	std::set<std::string> set_parameters;
	for (const std::string& key : vec_keys) {
		if (std::regex_match(key, param_matches, param_regex)) { //extract the parameter out of the key
			const std::string param( param_matches.str(1) );
			//catching some common syntax errors
			if (param_matches.str(2)!="::") throw InvalidFormatException("Missing ':' in \""+key+"\"", AT);
			//exclude the TIME filters (they are processed earlier)
			if (param==TimeProcStack::timeParamName) continue;

			set_parameters.insert(param);
		}
	}

	return set_parameters;
}

void MeteoProcessor::getWindowSize(ProcessingProperties& o_properties) const
{
	ProcessingProperties tmp;

	for (auto& paramStack : processing_stack) {
		paramStack.second->getWindowSize(tmp);
		compareProperties(tmp, o_properties);
	}

	//Also take the Meteo1DInterpolator into account:
	mi1d.getWindowSize(tmp);
	compareProperties(tmp, o_properties);
}

void MeteoProcessor::compareProperties(const ProcessingProperties& newprop, ProcessingProperties& current)
{
	current.points_before = std::max(current.points_before, newprop.points_before);
	current.points_after = std::max(current.points_after, newprop.points_after);

	if (newprop.time_before > current.time_before)
		current.time_before = newprop.time_before;

	if (newprop.time_after > current.time_after)
		current.time_after = newprop.time_after;
}

void MeteoProcessor::process(std::vector< std::vector<MeteoData> >& ivec,
                             std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	std::swap(ivec, ovec);
	if (processing_stack.empty() || !enable_meteo_filtering) return;

	for (auto& paramStack : processing_stack) {
		std::swap(ovec, ivec);
		paramStack.second->process(ivec, ovec, second_pass);
	}
}

std::set<std::string> MeteoProcessor::initStationSet(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& keyword)
{
	std::set<std::string> results;
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first==keyword) {
			std::istringstream iss(vecArgs[ii].second);
			std::string word;
			while (iss >> word){
				results.insert(word);
			}
		}
	}

	return results;
}

std::vector<DateRange> MeteoProcessor::initTimeRestrictions(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& keyword, const std::string& where, const double& TZ)
{
	std::vector<DateRange> dates_specs;
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first==keyword) {
			std::vector<std::string> vecString;
			const size_t nrElems = IOUtils::readLineToVec(vecArgs[ii].second, vecString, ',');
			
			for (size_t jj=0; jj<nrElems; jj++) {
				Date d1, d2;
				const size_t delim_pos = vecString[jj].find(" - ");
				if (delim_pos==std::string::npos) {
					if (!IOUtils::convertString(d1, vecString[jj], TZ))
						throw InvalidFormatException("Could not process date restriction "+vecString[jj]+" for "+where, AT);
					dates_specs.push_back( DateRange(d1, d1) );
				} else {
					if (!IOUtils::convertString(d1, vecString[jj].substr(0, delim_pos), TZ))
						throw InvalidFormatException("Could not process date restriction "+vecString[jj].substr(0, delim_pos)+" for "+where, AT);
					if (!IOUtils::convertString(d2, vecString[jj].substr(delim_pos+3), TZ))
						throw InvalidFormatException("Could not process date restriction "+vecString[jj].substr(delim_pos+3)+" for "+where, AT);
					dates_specs.push_back( DateRange(d1, d2) );
				}
			}
		}
	}
	
	if (dates_specs.empty()) return dates_specs;
	
	//now sort the vector and merge overlapping ranges
	std::sort(dates_specs.begin(), dates_specs.end()); //in case of identical start dates, the oldest end date comes first
	for (size_t ii=0; ii<(dates_specs.size()-1); ii++) {
		if (dates_specs[ii]==dates_specs[ii+1]) {
			//remove exactly identical ranges
			dates_specs.erase(dates_specs.begin()+static_cast<ptrdiff_t>(ii+1)); //we should have a limited number of elements so this is no problem
			ii--; //we must redo the current element
		} else if (dates_specs[ii].start==dates_specs[ii+1].start || dates_specs[ii].end >= dates_specs[ii+1].start) {
			//remove overlapping ranges. Since the vector is sorted on (start, end), the overlap criteria is very simplified!
			dates_specs[ii].end = std::max(dates_specs[ii].end, dates_specs[ii+1].end);
			dates_specs.erase(dates_specs.begin()+static_cast<ptrdiff_t>(ii+1));
			ii--; //we must redo the current element
		}
	}
	
	return dates_specs;
}

const std::string MeteoProcessor::toString() const {
	std::ostringstream os;
	os << "<MeteoProcessor>\n";
	os << mi1d.toString();
	os << "Processing stacks:\n";
	map<string, ProcessingStack*>::const_iterator it;
	for (it=processing_stack.begin(); it != processing_stack.end(); ++it){
		os << (*it->second).toString();
	}
	os << "</MeteoProcessor>\n";
	return os.str();
}


RestrictionsIdx::RestrictionsIdx(const METEO_SET& vecMeteo, const std::vector<DateRange>& time_restrictions)
                : start(), end(), index(0)
{
	if (time_restrictions.empty()) {
		start.push_back( 0 );
		end.push_back( vecMeteo.size() );
	} else {
		const Date start_dt( vecMeteo.front().date ), end_dt( vecMeteo.back().date );
		
		size_t ts_idx = 0;
		while (ts_idx < time_restrictions.size()) {
			if (time_restrictions[ts_idx].end < start_dt) { //this time_restrictions is before vecMeteo, look for the next one!
				ts_idx++;
				continue;
			}
			if (time_restrictions[ts_idx].start > end_dt) break; //no more time_restrictions applicable to vecMeteo
			
			size_t jj=0;
			//look for time_restriction start
			while (jj<vecMeteo.size() && vecMeteo[jj].date < time_restrictions[ts_idx].start) jj++;
			start.push_back( jj );
			
			//look for time_restriction end
			while (jj<vecMeteo.size() && vecMeteo[jj].date <= time_restrictions[ts_idx].end) jj++;
			end.push_back( jj );
			
			//move to next time restriction period
			ts_idx++;
		}
		
		//no applicable time_restrictions found for vecMeteo
		if (start.empty()) index = IOUtils::npos;
	}
}

size_t RestrictionsIdx::getStart() const
{
	if (index == IOUtils::npos) return IOUtils::npos;
	return start[ index ];
}

size_t RestrictionsIdx::getEnd() const
{
	if (index == IOUtils::npos) return IOUtils::npos;
	return end[ index ];
}

RestrictionsIdx& RestrictionsIdx::operator++()
{
	if (index!=IOUtils::npos) {
		index++;
		if (index >= start.size()) index=IOUtils::npos;
	}
	
	return *this;
}

const std::string RestrictionsIdx::toString() const
{
	std::ostringstream os;
	os << "[ ";
	for (size_t ii=0; ii<start.size(); ii++)
		os << "(" << start[ii] << "," << end[ii] << ") ";
	
	os << "]";
	return os.str();
}


} //namespace
