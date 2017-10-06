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

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : filter_stack(), param_name(parname)
{
	static const  char NUM[] = "0123456789";
	static const std::string filter_key( "::FILTER" );
	static const std::string arg_key( "::ARG" );

	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecFilters( cfg.getValues(parname+filter_key, "FILTERS") );
	for (size_t ii=0; ii<vecFilters.size(); ii++) {
		const std::string block_name( IOUtils::strToUpper( vecFilters[ii].second ) );
		if (block_name=="NONE") continue;

		//extract the filter number and perform basic checks on the syntax
		const std::string key( vecFilters[ii].first );
		const size_t end_filter = key.find(filter_key); //we know this will be found since it has been matched in cfg.getValues()
		const size_t start_filter_nr = key.find_first_of(NUM, end_filter+filter_key.length());
		const size_t end_filter_nr = key.find_first_not_of(NUM, end_filter+filter_key.length());
		if (start_filter_nr==std::string::npos || end_filter_nr!=std::string::npos) throw InvalidArgumentException("Syntax error: "+key, AT);

		unsigned int filter_nr;
		const std::string filter_nr_str( key.substr(start_filter_nr) );
		if ( !IOUtils::convertString(filter_nr, filter_nr_str) ) InvalidArgumentException("Can not parse filter number in "+key, AT);

		//read the arguments and clean them up (ie remove the {Param}::{args##}:: in front of the argument key itself)
		std::ostringstream arg_str;
		arg_str << param_name << arg_key << filter_nr;
		std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getValues(arg_str.str(), "FILTERS") );
		for (size_t jj=0; jj<vecArgs.size(); jj++) {
			const size_t beg_arg_name = vecArgs[jj].first.find_first_not_of(":", arg_str.str().length());
			if (beg_arg_name==std::string::npos)
				throw InvalidFormatException("Wrong argument format for '"+vecArgs[jj].first+"'", AT);
			vecArgs[jj].first = vecArgs[jj].first.substr(beg_arg_name);
		}

		//construct the filter with its name and arguments
		filter_stack.push_back( BlockFactory::getBlock(block_name, vecArgs, cfg) );
	}
}

ProcessingStack::~ProcessingStack()
{
	for (size_t ii=0; ii<filter_stack.size(); ii++)
		delete filter_stack[ii];
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

//this method applies the whole processing stack for all the stations, all the data points for one meteo param
//(as defined in the constructor)
void ProcessingStack::process(const std::vector< std::vector<MeteoData> >& ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	const size_t nr_of_filters = filter_stack.size();
	const size_t nr_stations = ivec.size();
	ovec.resize( nr_stations );

	for (size_t ii=0; ii<nr_stations; ii++) { //for every station
		if ( ivec[ii].empty() ) continue; //no data, nothing to do!

		//pick one element and check whether the param_name parameter exists
		const size_t param = ivec[ii].front().getParameterIndex(param_name);
		if (param != IOUtils::npos) {
			const std::string statID( ivec[ii][0].meta.getStationID() ); //we know there is at least 1 element (see above)
			std::vector<MeteoData> tmp( ivec[ii] );

			//Now call the filters one after another for the current station and parameter
			bool appliedFilter = false;
			for (size_t jj=0; jj<nr_of_filters; jj++) {
				if ((*filter_stack[jj]).skipStation( statID ))
					continue;

				const ProcessingProperties::proc_stage filter_stage( filter_stack[jj]->getProperties().stage );
				if ( second_pass && ((filter_stage==ProcessingProperties::first) || (filter_stage==ProcessingProperties::none)) )
					continue;
				if ( !second_pass && ((filter_stage==ProcessingProperties::second) || (filter_stage==ProcessingProperties::none)) )
					continue;

				appliedFilter = true;
				(*filter_stack[jj]).process(static_cast<unsigned int>(param), tmp, ovec[ii]);

				if (tmp.size() != ovec[ii].size()) {
					ostringstream ss;
					ss << "The filter \"" << (*filter_stack[jj]).getName() << "\" received " << tmp.size();
					ss << " timestamps and returned " << ovec[ii].size() << " timestamps!";
					throw IndexOutOfBoundsException(ss.str(), AT);
				}

				#ifdef DATA_QA
				for (size_t kk=0; kk<ovec[ii].size(); kk++) {
					const double orig = tmp[kk](param);
					const double filtered = ovec[ii][kk](param);
					if (orig!=filtered) {
						const std::string statName( ovec[ii][kk].meta.getStationName() );
						const std::string stat = (!statID.empty())? statID : statName;
						const std::string filtername( (*filter_stack[jj]).getName() );
						cout << "[DATA_QA] Filtering " << stat << "::" << param_name << "::" << filtername << " " << tmp[kk].date.toString(Date::ISO_TZ) << " [" << tmp[kk].date.toString(Date::ISO_WEEK) << "]\n";
					}
				}
				#endif
				if ((jj+1) != nr_of_filters) {//not necessary after the last filter
					for (size_t kk=0; kk<ovec[ii].size(); kk++) {
						tmp[kk](param) = ovec[ii][kk](param);
					}
				}
			}

			if (!appliedFilter) //if not a single filter was applied
				ovec[ii] = ivec[ii]; //just copy input to output
		} else {
			ovec[ii] = ivec[ii]; //just copy input to output
		}
	}
}

const std::string ProcessingStack::toString() const
{
	std::ostringstream os;
	//os << "<ProcessingStack>";
	os << setw(10) << param_name << "::";

	for (size_t ii=0; ii<filter_stack.size(); ii++) {
		os << setw(10) << (*filter_stack[ii]).toString();
	}

	//os << "</ProcessingStack>";
	os << "\n";
	return os.str();
}

} //end namespace
