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
#include <meteoio/meteoFilters/ProcessingStack.h>

#include <algorithm>
#include <iomanip>

using namespace std;

namespace mio {
const std::string ProcessingStack::filter_section( "FILTERS" );
const std::string ProcessingStack::filter_pattern( "::FILTER" );
const std::string ProcessingStack::arg_pattern( "::ARG" );

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : filter_stack(), param_name(parname), data_qa_logs(false)
{
	cfg.getValue("DATA_QA_LOGS", "GENERAL", data_qa_logs, IOUtils::nothrow);

	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecFilters( cfg.getValues(parname+filter_pattern, filter_section) );
	filter_stack.reserve( vecFilters.size() );
	for (size_t ii=0; ii<vecFilters.size(); ii++) {
		const std::string block_name( IOUtils::strToUpper( vecFilters[ii].second ) );
		if (block_name=="NONE") continue;

		const unsigned int cmd_nr = Config::getCommandNr(filter_section, parname+filter_pattern, vecFilters[ii].first);
		const std::vector< std::pair<std::string, std::string> > vecArgs( cfg.parseArgs(filter_section, parname, cmd_nr, arg_pattern) );
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
		//the filters only receive tmp_ovec and expand it to tmp_ivec so ovec may not have been resized...
		ovec = ivec;

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

	std::vector<double> heights_for_parameter(ivec.front().getHeightsForParameter(param_name));
	for (auto& height : heights_for_parameter) {
		std::string full_paramname = param_name;
		if (height != IOUtils::nodata) {
			full_paramname = full_paramname + "@" + MeteoData::convertHeightToString(height);
		} else {
			height = ProcessingBlock::default_height; // change the default to a more user friendly height instead of -999
		}

		//pick one element and check whether the param_name parameter exists
		const size_t param = ivec.front().getParameterIndex(full_paramname);
		if (param == IOUtils::npos) continue;

		const size_t nr_of_filters = filter_stack.size();
		const std::string statID( ivec.front().meta.getStationID() ); //we know there is at least 1 element (we've already skipped empty vectors)

		//Now call the filters one after another for the current station and parameter
		for (size_t jj=0; jj<nr_of_filters; jj++) {
			if (filter_stack[jj]->skipStation( statID ))
				continue;

			if (filter_stack[jj]->skipHeight(height))
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

			if (ivec.size() == output_size) {
				for (size_t kk=0; kk<ivec.size(); kk++) {
					const double orig = ivec[kk](param);
					const double filtered = ovec[stat_idx][kk](param);
					if (orig!=filtered) {
						ovec[stat_idx][kk].setFiltered(param);
						if (data_qa_logs) {
							const std::string statName( ovec[stat_idx][kk].meta.getStationName() );
							const std::string stat = (!statID.empty())? statID : statName;
							const std::string filtername( (*filter_stack[jj]).getName() );
							cout << "[DATA_QA] Filtering " << stat << "::" << full_paramname << "::" << filtername << " " << ivec[kk].date.toString(Date::ISO_TZ) << " [" << ivec[kk].date.toString(Date::ISO_WEEK) << "]\n";
						}
					}
				}
			} else { //filters such as SHIFT might change the number of points
				size_t kk_out=0;
				for (size_t kk=0; kk<ivec.size(); kk++) {
					while (kk_out<output_size && ovec[stat_idx][kk_out].date < ivec[kk].date) { //new points inserted
						ovec[stat_idx][kk_out].setFiltered(param);
						if (data_qa_logs) {
							const std::string statName( ovec[stat_idx][kk_out].meta.getStationName() );
							const std::string stat = (!statID.empty())? statID : statName;
							const std::string filtername( (*filter_stack[jj]).getName() );
							cout << "[DATA_QA] Filtering " << stat << "::" << full_paramname << "::" << filtername << " " << ivec[kk].date.toString(Date::ISO_TZ) << " [" << ivec[kk].date.toString(Date::ISO_WEEK) << "]\n";
						}
						kk_out++;
					}
					if (kk_out==output_size) break;

					if (ovec[stat_idx][kk_out].date == ivec[kk].date) {
						const double orig = ivec[kk](param);
						const double filtered = ovec[stat_idx][kk_out](param);
						if (orig!=filtered) {
							ovec[stat_idx][kk_out].setFiltered(param);
							if (data_qa_logs) {
								const std::string statName( ovec[stat_idx][kk_out].meta.getStationName() );
								const std::string stat = (!statID.empty())? statID : statName;
								const std::string filtername( (*filter_stack[jj]).getName() );
								cout << "[DATA_QA] Filtering " << stat << "::" << full_paramname << "::" << filtername << " " << ivec[kk].date.toString(Date::ISO_TZ) << " [" << ivec[kk].date.toString(Date::ISO_WEEK) << "]\n";
							}
						}
					}
				}
			}

			if ((jj+1) != nr_of_filters) {//not necessary after the last filter
				ivec = ovec[stat_idx];
			}
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
