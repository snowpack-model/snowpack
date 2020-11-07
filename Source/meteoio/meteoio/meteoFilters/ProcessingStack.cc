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
#include <meteoio/meteoFilters/ProcessingStack.h>

#include <algorithm>

using namespace std;

namespace mio {
const char ProcessingStack::NUM[] = "0123456789";
const std::string ProcessingStack::filter_key( "::FILTER" );
const std::string ProcessingStack::arg_key( "::ARG" );

std::vector< std::pair<std::string, std::string> > ProcessingStack::parseArgs(const Config& cfg, const std::string& key, const std::string& parname)
{
	//extract the filter number and perform basic checks on the syntax
	const size_t end_filter = key.find(filter_key); //we know this will be found since it has been matched in cfg.getValues()
	const size_t start_filter_nr = key.find_first_of(NUM, end_filter+filter_key.length());
	const size_t end_filter_nr = key.find_first_not_of(NUM, end_filter+filter_key.length());
	if (start_filter_nr==std::string::npos || end_filter_nr!=std::string::npos) throw InvalidArgumentException("Syntax error: "+key, AT);

	unsigned int filter_nr;
	const std::string filter_nr_str( key.substr(start_filter_nr) );
	if ( !IOUtils::convertString(filter_nr, filter_nr_str) ) InvalidArgumentException("Can not parse filter number in "+key, AT);

	//read the arguments and clean them up (ie remove the {Param}::{args##}:: in front of the argument key itself)
	std::ostringstream arg_str;
	arg_str << parname << arg_key << filter_nr;
	std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getValues(arg_str.str(), "FILTERS") );
	for (size_t jj=0; jj<vecArgs.size(); jj++) {
		const size_t beg_arg_name = vecArgs[jj].first.find_first_not_of(":", arg_str.str().length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+vecArgs[jj].first+"'", AT);
		vecArgs[jj].first = vecArgs[jj].first.substr(beg_arg_name);
	}
	
	return vecArgs;
}

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : filter_stack(), param_name(parname), data_qa_logs(false)
{
	cfg.getValue("DATA_QA_LOGS", "GENERAL", data_qa_logs, IOUtils::nothrow);
	
	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecFilters( cfg.getValues(parname+filter_key, "FILTERS") );
	filter_stack.reserve( vecFilters.size() );
	for (size_t ii=0; ii<vecFilters.size(); ii++) {
		const std::string block_name( IOUtils::strToUpper( vecFilters[ii].second ) );
		if (block_name=="NONE") continue;
		
		const std::vector< std::pair<std::string, std::string> > vecArgs( parseArgs(cfg, vecFilters[ii].first, parname) );
		filter_stack.push_back( BlockFactory::getBlock(block_name, vecArgs, cfg) );
	}
}

void ProcessingStack::getWindowSize(ProcessingProperties& o_properties) const
{
	o_properties.points_before = 0;
	o_properties.points_after = 0;
	o_properties.time_after = Duration(0.0, 0.);
	o_properties.time_before = Duration(0.0, 0.);

	for (size_t jj=0; jj<filter_stack.size(); jj++){
		const ProcessingProperties properties( (*filter_stack[jj]).getProperties() );

		o_properties.points_before = std::max(o_properties.points_before, properties.points_before);
		o_properties.points_after = std::max(o_properties.points_after, properties.points_after);

		if (properties.time_before > o_properties.time_before)
			o_properties.time_before = properties.time_before;

		if (properties.time_after > o_properties.time_after)
			o_properties.time_after = properties.time_after;
	}
}

bool ProcessingStack::applyFilter(const size_t& param, const size_t& jj, const std::vector<MeteoData>& ivec, std::vector<MeteoData> &ovec)
{
	const std::vector<DateRange> time_restrictions( filter_stack[jj]->getTimeRestrictions() );
	
	if (time_restrictions.empty()) {
		filter_stack[jj]->process(static_cast<unsigned int>(param), ivec, ovec);
		return true;
	} else {
		//we know there is at least 1 element
		const Date start( ivec.front().date ), end( ivec.back().date );
		bool filterApplied = false;
		
		//filter for all time restrictions that fit into ivec
		for (size_t ii=0; ii<time_restrictions.size(); ii++) {
			if (time_restrictions[ii].end<start) continue; //sub-range before our data
			if (time_restrictions[ii].start>end) break; //sub-range after our data
			
			const Date sub_start( std::max(time_restrictions[ii].start, start) );
			const Date sub_end( std::min(time_restrictions[ii].end, end) );
			
			//now cut a buffer to this sub-range and filter it
			size_t ivec_start_idx=0; //position of the first sub-range within the full ivec/ovec
			while (ivec[ivec_start_idx].date<sub_start) ivec_start_idx++;
			
			//prepare the sub-range input data from the original ivec
			std::vector<MeteoData> tmp_ivec, tmp_ovec;
			tmp_ivec.reserve( ivec.size()-ivec_start_idx ); //worst case scenario
			for (size_t kk=ivec_start_idx; kk<ivec.size(); kk++) {
				if (ivec[kk].date>sub_end) break;
				tmp_ivec.push_back( ivec[kk] );
			}
			
			filter_stack[jj]->process(static_cast<unsigned int>(param), tmp_ivec, tmp_ovec);
			
			//put back the sub-range filtered data into the full ovec
			for (size_t kk=0; kk<tmp_ovec.size(); kk++) {
				ovec[kk+ivec_start_idx] = tmp_ovec[kk];
			}
			
			filterApplied = true;
		}
		
		return filterApplied;
	}
}

//ivec is passed by value, so it makes an efficient copy
bool ProcessingStack::filterStation(std::vector<MeteoData> ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass, const size_t& stat_idx)
{
	bool filterApplied = false;
	
	//pick one element and check whether the param_name parameter exists
	const size_t param = ivec.front().getParameterIndex(param_name);
	if (param == IOUtils::npos) return filterApplied;
	
	const size_t nr_of_filters = filter_stack.size();
	const std::string statID( ivec.front().meta.getStationID() ); //we know there is at least 1 element (we've already skipped empty vectors)

	//Now call the filters one after another for the current station and parameter
	for (size_t jj=0; jj<nr_of_filters; jj++) {
		if (filter_stack[jj]->skipStation( statID ))
			continue;

		const ProcessingProperties::proc_stage filter_stage( filter_stack[jj]->getProperties().stage );
		if ( second_pass && ((filter_stage==ProcessingProperties::first) || (filter_stage==ProcessingProperties::none)) )
			continue;
		if ( !second_pass && ((filter_stage==ProcessingProperties::second) || (filter_stage==ProcessingProperties::none)) )
			continue;

		//if the filter has not been applied (ie time restriction), move to the next one directly
		if (!applyFilter(param, jj, ivec, ovec[stat_idx])) continue;
		
		filterApplied = true; //at least one filter has been applied in the whole stack
		const size_t output_size = ovec[stat_idx].size();

		if (ivec.size() != output_size) {
			ostringstream ss;
			ss << "The filter \"" << (*filter_stack[jj]).getName() << "\" received " << ivec.size();
			ss << " timestamps and returned " << output_size << " timestamps!";
			throw IndexOutOfBoundsException(ss.str(), AT);
		}

		for (size_t kk=0; kk<output_size; kk++) {
			const double orig = ivec[kk](param);
			const double filtered = ovec[stat_idx][kk](param);
			if (orig!=filtered) {
				ovec[stat_idx][kk].setFiltered(param);
				if (data_qa_logs) {
					const std::string statName( ovec[stat_idx][kk].meta.getStationName() );
					const std::string stat = (!statID.empty())? statID : statName;
					const std::string filtername( (*filter_stack[jj]).getName() );
					cout << "[DATA_QA] Filtering " << stat << "::" << param_name << "::" << filtername << " " << ivec[kk].date.toString(Date::ISO_TZ) << " [" << ivec[kk].date.toString(Date::ISO_WEEK) << "]\n";
				}
			}
		}

		if ((jj+1) != nr_of_filters) {//not necessary after the last filter
			ivec = ovec[stat_idx];
		}
	}

	return filterApplied;
}

//this method applies the whole processing stack for all the stations, all the data points for one meteo param
//(as defined in the constructor)
void ProcessingStack::process(const std::vector< std::vector<MeteoData> >& ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	const size_t nr_stations = ivec.size();
	ovec.resize( nr_stations );

	for (size_t ii=0; ii<nr_stations; ii++) { //for every station
		if ( ivec[ii].empty() ) continue; //no data, nothing to do!
		
		const bool filterApplied = filterStation(ivec[ii], ovec, second_pass, ii);
		//if not even a single filter was applied, just copy input to output
		if (!filterApplied) ovec[ii] = ivec[ii];
	}
}

const std::string ProcessingStack::toString() const
{
	std::ostringstream os;
	//os << "<ProcessingStack>";
	os << setw(10) << param_name << "::";

	for (size_t ii=0; ii<filter_stack.size(); ii++) {
		os << setw(10) << filter_stack[ii]->toString();
	}

	//os << "</ProcessingStack>";
	os << "\n";
	return os.str();
}

} //end namespace
