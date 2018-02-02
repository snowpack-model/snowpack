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
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <algorithm>

#include <meteoio/meteoFilters/TimeFilters.h>
#include <meteoio/meteoFilters/ProcessingStack.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

static inline bool IsUndef (const MeteoData& md) { return md.date.isUndef(); }

TimeSuppr::TimeSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ)
          : ProcessingBlock(vecArgs, name), suppr_dates(), range(IOUtils::nodata)
{
	const std::string where( "Filters::"+block_name );
	properties.stage = ProcessingProperties::first; //for the rest: default values
	const size_t nrArgs = vecArgs.size();

	if (nrArgs!=1)
		throw InvalidArgumentException("Wrong number of arguments for " + where, AT);

	if (vecArgs[0].first=="FRAC") {
		if (!IOUtils::convertString(range, vecArgs[0].second))
			throw InvalidArgumentException("Invalid range \""+vecArgs[0].second+"\" specified for "+where, AT);
		if (range<0. || range>1.)
			throw InvalidArgumentException("Wrong range for " + where + ", it should be between 0 and 1", AT);
	} else if (vecArgs[0].first=="SUPPR") {
		const std::string in_filename( vecArgs[0].second );
		const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
		const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );

		suppr_dates = ProcessingBlock::readDates(block_name, filename, TZ);
	} else
		throw UnknownValueException("Unknown option '"+vecArgs[0].first+"' for "+where, AT);
}

void TimeSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	if (param!=IOUtils::unodata)
		throw InvalidArgumentException("The filter "+block_name+" can only be applied to TIME", AT);
	
	ovec = ivec;
	if (ovec.empty()) return;
	
	if (!suppr_dates.empty()) {
		supprByDates(ovec);
	} else { //only remove a given fraction
		supprFrac(ovec);
	}
}

//this assumes that the DATES_RANGEs in suppr_dates have been sorted by increasing starting dates
void TimeSuppr::supprByDates(std::vector<MeteoData>& ovec) const
{
	const std::string station_ID( ovec[0].meta.stationID ); //we know it is not empty
	const std::map< std::string, std::vector<dates_range> >::const_iterator station_it( suppr_dates.find( station_ID ) );
	if (station_it==suppr_dates.end()) return;

	const std::vector<dates_range> &suppr_specs = station_it->second;
	const size_t Nset = suppr_specs.size();
	size_t curr_idx = 0; //we know there is at least one
	for (size_t ii=0; ii<ovec.size(); ii++) {
		if (ovec[ii].date<suppr_specs[curr_idx].start) continue;

		if (ovec[ii].date<=suppr_specs[curr_idx].end) { //suppress the interval
			ovec[ii].date.setUndef(true); //mark the point to be removed
			continue;
		} else { //look for the next interval
			curr_idx++;
			if (curr_idx>=Nset) break; //all the suppression points have been processed
		}
	}
	
	//now really remove the points from the vector
	ovec.erase( std::remove_if(ovec.begin(), ovec.end(), IsUndef), ovec.end());
}

void TimeSuppr::supprFrac(std::vector<MeteoData>& ovec) const
{
	const size_t set_size = ovec.size();
	const size_t nrRemove = static_cast<size_t>( round( (double)set_size*range ) );

	srand( static_cast<unsigned int>(time(NULL)) );
	size_t ii=1;
	while (ii<nrRemove) {
		const size_t idx = (unsigned)rand() % set_size;
		if (ovec[idx].date.isUndef()) continue; //the point was already removed

		ovec[idx].date.setUndef(true);
		ii++;
	}
	
	//now really remove the points from the vector
	ovec.erase( std::remove_if(ovec.begin(), ovec.end(), IsUndef), ovec.end());
}


TimeUnDST::TimeUnDST(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ)
        : ProcessingBlock(vecArgs, name), dst_changes()
{
	const std::string where( "Filters::"+block_name );
	properties.stage = ProcessingProperties::first; //for the rest: default values
	const size_t nrArgs = vecArgs.size();
	
	if (nrArgs!=1)
		throw InvalidArgumentException("Wrong number of arguments for " + where, AT);

	if (vecArgs[0].first=="CORRECTIONS") {
		const std::string in_filename( vecArgs[0].second );
		const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
		const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );

		dst_changes = ProcessingBlock::readCorrections(block_name, filename, TZ, 2);
		if (dst_changes.empty())
			throw InvalidArgumentException("Please provide at least one DST correction for " + where, AT);
	} else
		throw UnknownValueException("Unknown option '"+vecArgs[0].first+"' for "+where, AT);
}

void TimeUnDST::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{
	if (param!=IOUtils::unodata)
		throw InvalidArgumentException("The filter "+block_name+" can only be applied to TIME", AT);
	
	static const double sec2Jul = 1./(24.*3600.);
	ovec = ivec;
	if (ovec.empty()) return;
	
	const size_t Nset = dst_changes.size(); //we know there is at least one
	size_t next_idx=0;
	double offset = 0.;
	Date prev_date = ovec[0].date - 1.; //so we are sure to be < when checking if timestamps are increasing
	size_t ii=0;
	for (; ii<ovec.size(); ii++) {
		bool apply_change = (ovec[ii].date>=dst_changes[next_idx].date);
		
		//when reverting back to winter time, timestamps are not in increasing order for an overlap period
		if (ovec[ii].date<=prev_date) {
			const double overlap = (dst_changes[next_idx].offset*sec2Jul-offset);
			if (ovec[ii].date>=(dst_changes[next_idx].date - overlap))
				apply_change = true;
		}
		
		if (apply_change) {
			offset = dst_changes[next_idx].offset * sec2Jul;
			next_idx++;
			if (next_idx==Nset) break; //no more new corrections to expect
		}
		
		prev_date = ovec[ii].date;
		if (offset!=0.) ovec[ii].date += offset;
	}
	
	if (offset==0) return; //no more corrections to apply
	
	//if some points remained after the last DST correction date, process them
	for (; ii<ovec.size(); ii++) {
		ovec[ii].date += offset;
	}
}


const std::string TimeProcStack::timeParamName( "TIME" );
TimeProcStack::TimeProcStack(const Config& cfg) : filter_stack()
{
	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecFilters( cfg.getValues(timeParamName+ProcessingStack::filter_key, "FILTERS") );
	for (size_t ii=0; ii<vecFilters.size(); ii++) {
		const std::string block_name( IOUtils::strToUpper( vecFilters[ii].second ) );
		if (block_name=="NONE") continue;
		
		const std::vector< std::pair<std::string, std::string> > vecArgs( ProcessingStack::parseArgs(cfg, vecFilters[ii].first, timeParamName) );
		filter_stack.push_back( BlockFactory::getTimeBlock(block_name, vecArgs, cfg) );
	}
}

//ivec is passed by value, so it makes an efficient copy
void TimeProcStack::process(std::vector< std::vector<MeteoData> >& ivec, const bool& second_pass)
{
	const size_t nr_of_filters = filter_stack.size();
	const size_t nr_stations = ivec.size();

	std::vector<MeteoData> ovec;
	for (size_t ii=0; ii<nr_stations; ii++) { //for every station
		if ( ivec[ii].empty() ) continue; //no data, nothing to do!
		
		const std::string statID( ivec[ii].front().meta.getStationID() ); //we know there is at least 1 element (we've already skipped empty vectors)
		//Now call the filters one after another for the current station and parameter
		for (size_t jj=0; jj<nr_of_filters; jj++) {
			if ((*filter_stack[jj]).skipStation( statID ))
				continue;

			const ProcessingProperties::proc_stage filter_stage( filter_stack[jj]->getProperties().stage );
			if ( second_pass && ((filter_stage==ProcessingProperties::first) || (filter_stage==ProcessingProperties::none)) )
				continue;
			if ( !second_pass && ((filter_stage==ProcessingProperties::second) || (filter_stage==ProcessingProperties::none)) )
				continue;

			(*filter_stack[jj]).process(IOUtils::unodata, ivec[ii], ovec);
			ivec[ii] = ovec;
		}
	}
}

/** 
 * @brief check that timestamps are unique and in increasing order
 * @param[in] vecVecMeteo all the data for all the stations
*/
void TimeProcStack::checkUniqueTimestamps(const std::vector<METEO_SET>& vecVecMeteo)
{
	for (size_t stat_idx=0; stat_idx<vecVecMeteo.size(); ++stat_idx) { //for each station
		const size_t nr_timestamps = vecVecMeteo[stat_idx].size();
		if (nr_timestamps==0) continue;

		Date previous_date( vecVecMeteo[stat_idx].front().date );
		for (size_t ii=1; ii<nr_timestamps; ++ii) {
			const Date current_date( vecVecMeteo[stat_idx][ii].date );
			if (current_date<=previous_date) {
				const StationData& station( vecVecMeteo[stat_idx][ii].meta );
				if (current_date==previous_date)
					throw IOException("Error for station \""+station.stationName+"\" ("+station.stationID+") at time "+current_date.toString(Date::ISO)+": timestamps must be unique!", AT);
				else
					throw IOException("Error for station \""+station.stationName+"\" ("+station.stationID+"): jumping from "+previous_date.toString(Date::ISO)+" to "+current_date.toString(Date::ISO), AT);
			}
			previous_date = current_date;
		}
	}
}

} //end namespace
