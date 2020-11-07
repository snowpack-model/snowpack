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
#include <meteoio/MeteoProcessor.h>
#include <meteoio/meteoFilters/TimeFilters.h>
#include <algorithm>

using namespace std;

namespace mio {

MeteoProcessor::MeteoProcessor(const Config& cfg, const char& rank, const IOUtils::OperationMode &mode) : mi1d(cfg, rank, mode), processing_stack()
{
	//Parse [Filters] section, create processing stack for each configured parameter
	const std::set<std::string> set_of_used_parameters( getParameters(cfg) );

	for (std::set<std::string>::const_iterator it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		ProcessingStack* tmp = new ProcessingStack(cfg, *it);
		processing_stack[*it] = tmp;
	}
}

MeteoProcessor::~MeteoProcessor()
{
	//clean up heap memory
	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it)
		delete it->second;
}

std::set<std::string> MeteoProcessor::getParameters(const Config& cfg)
{
	const std::vector<std::string> vec_keys( cfg.getKeys(std::string(), "Filters") );

	std::set<std::string> set_parameters;
	for (size_t ii=0; ii<vec_keys.size(); ++ii){
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			if (vec_keys[ii].length()<=(found+2))
				throw InvalidFormatException("Invalid syntax: \""+vec_keys[ii]+"\"", AT);
			if (vec_keys[ii][found+1]!=':')
				throw InvalidFormatException("Missing ':' in \""+vec_keys[ii]+"\"", AT);
				
			const std::string tmp( vec_keys[ii].substr(0,found) );
			if (tmp==TimeProcStack::timeParamName) continue; //exclude the TIME filters (they are processed earlier)
			set_parameters.insert(tmp);
		}
	}

	return set_parameters;
}

void MeteoProcessor::getWindowSize(ProcessingProperties& o_properties) const
{
	ProcessingProperties tmp;

	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it){
		(*(it->second)).getWindowSize(tmp);
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
	if (processing_stack.empty()) return;
	
	for (std::map<std::string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it) {
		std::swap(ovec, ivec);
		(*(it->second)).process(ivec, ovec, second_pass);
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
				const size_t delim_pos = vecString[jj].find(" - ");
				if (delim_pos==std::string::npos)
					throw InvalidFormatException("Invalid time restriction syntax for " + where + ": the two dates must be separated by ' - '", AT);
				
				Date d1, d2;
				if (!IOUtils::convertString(d1, vecString[jj].substr(0, delim_pos), TZ))
					throw InvalidFormatException("Could not process date "+vecString[jj].substr(0, delim_pos)+" for "+where, AT);
				if (!IOUtils::convertString(d2, vecString[jj].substr(delim_pos+3), TZ))
					throw InvalidFormatException("Could not process date "+vecString[jj].substr(delim_pos+3)+" for "+where, AT);
				dates_specs.push_back( DateRange(d1, d2) );
			}
		}
	}
	
	if (dates_specs.empty()) return dates_specs;
	
	//now sort the vector and merge overlapping ranges
	std::sort(dates_specs.begin(), dates_specs.end()); //in case of identical start dates, the oldest end date comes first
	for (size_t ii=0; ii<(dates_specs.size()-1); ii++) {
		if (dates_specs[ii]==dates_specs[ii+1]) {
			//remove identical ranges
			dates_specs.erase(dates_specs.begin()+ii+1); //we should have a limited number of elements so this is no problem
			ii--; //we must redo the current element
		} else if (dates_specs[ii].start==dates_specs[ii+1].start || dates_specs[ii].end >= dates_specs[ii+1].start) {
			//remove overlapping ranges
			dates_specs[ii].end = dates_specs[ii+1].end;
			dates_specs.erase(dates_specs.begin()+ii+1);
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
		//os << setw(10) << it->first.toString() << "::"; //the processing stack already contains it
		os << (*it->second).toString();
	}
	os << "</MeteoProcessor>\n";
	return os.str();
}

} //namespace
