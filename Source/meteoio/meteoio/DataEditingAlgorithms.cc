/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research  SLF-DAVOS        */
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
#include <meteoio/DataEditingAlgorithms.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/MeteoData.h> //needed for the merge strategies
#include <meteoio/MeteoProcessor.h> //needed for the time restrictions

#include <algorithm>
#include <fstream>

using namespace std;

namespace mio {
/**
 * @page data_editing Input Data Editing
 * Before any filters, resampling algorithms or data generators are applied, it is possible to edit the original data. There are several
 * edition commands that can be stacked at will, per station ID. This is similar to the way that filters (\ref processing "processing elements") are
 * also stacked together. In order to rectrict any editing to a specific set of time ranges, use the **when** option followed by a comma 
 * delimited list of date intervals (represented by two ISO formatted dates seperated by ' - '), similarly to the 
 * \ref processing "Filters" (but right now, the automerge and merge commands don't support these time restrictions). The general 
 * syntax is ('#' represent a number, so each key remains unique):
 * @code
 * {stationID}::edit# = {command}
 * {stationID}::arg#::{argument name} = {values}
 * 
 * #here is an example
 * WFJ2::edit1 = EXCLUDE
 * WFJ2::arg1::exclude = VW DW ISWR RSWR
 * WFJ2::arg1::when = 2019-12-01T13:00 - 2019-12-25 , 2020-03-05 - 2020-04-15T12:30
 * @endcode
 * 
 * It is also possible to apply a stack of edits to all stations by using the '*' wildcard instead of the station ID 
 * (in such a case, the wildcard stack will be applied before any other stack). Currently, the data creation is applied 
 * after all stacks have been processed but this will change in the near future...
 * 
 * The following Input Data Editing commands are available:
 *     - SWAP: swap two parameters, see EditingSwap
 *     - MOVE: rename one or more parameters into a new name, see EditingMove
 *     - EXCLUDE: delete a list of parameters, see EditingExclude
 *     - KEEP: only keep a list of parameters and reject the others, see EditingKeep
 *     - AUTOMERGE: merge together stations sharing the same station ID, see EditingAutoMerge
 *     - MERGE: merge together one or more stations, see EditingMerge
 *     - COPY: make a copy of a given parameter under a new name, see EditingCopy
 * 
 * @section data_creation Data creation (CREATE)
 * Finally, it is possible to create new data based on some parametrizations. If the requested parameter does not exists, it will be created. Otherwise,
 * any pre-existing data is kept and only missing values in the original data set are filled with the generated values, keeping the original sampling rate. As
 * with all raw data editing, this takes place *before* any filtering/resampling/data generators. As the available algorithms are the same as for the
 * data generators, they are listed in the \ref generators_keywords "data generators section" (but the data creators must be declared in the [InputEditing] section).
 * @code
 * [InputEditing]
 * P::create = STD_PRESS			#the pressure is filled with STD_PRESS if no measured values are available
 * ISWR_POT::create = clearSky_SW		#a new parameter "ISWR_POT" is created and filled with Clear Sky values
 * @endcode
 */

EditingBlock::EditingBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg) 
             : time_restrictions( MeteoProcessor::initTimeRestrictions(vecArgs, "WHEN", "InputEditing::"+name+" for station "+i_stationID, cfg.get("TIME_ZONE", "Input")) ), stationID(i_stationID), block_name(name) {}
             
const std::string EditingBlock::toString() const 
{
	std::ostringstream os;
	os << "[" << stationID << " - " << block_name << "]";
	return os.str();
}

EditingBlock* EditingBlockFactory::getBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
{
	if (name == "SWAP"){
		return new EditingSwap(i_stationID, vecArgs, name, cfg);
	} else if (name == "MOVE"){
		return new EditingMove(i_stationID, vecArgs, name, cfg);
	} else if (name == "EXCLUDE"){
		return new EditingExclude(i_stationID, vecArgs, name, cfg);
	} else if (name == "KEEP"){
		return new EditingKeep(i_stationID, vecArgs, name, cfg);
	} else if (name == "MERGE"){
		return new EditingMerge(i_stationID, vecArgs, name, cfg);
	} else if (name == "AUTOMERGE"){
		return new EditingAutoMerge(i_stationID, vecArgs, name, cfg);
	} if (name == "COPY"){
		return new EditingCopy(i_stationID, vecArgs, name, cfg);
	} else {
		throw IOException("The input data editing block '"+name+"' does not exist! " , AT);
	}
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


////////////////////////////////////////////////// SWAP
EditingSwap::EditingSwap(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), dest_param(), src_param()
{
	parse_args(vecArgs);
}

void EditingSwap::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param ),
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingSwap::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
				
				const size_t src_index = vecMeteo[station][jj].addParameter( src_param ); //either add or just return the proper index
				const double src_value = vecMeteo[station][jj]( src_param );
				
				const size_t dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
				
				vecMeteo[station][jj]( src_index ) = vecMeteo[station][jj]( dest_index );
				vecMeteo[station][jj]( dest_index ) = src_value;
			}
		}
	}
}


////////////////////////////////////////////////// MOVE
EditingMove::EditingMove(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), src_params(), dest_param()
{
	parse_args(vecArgs);
}

void EditingMove::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), src_params);
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src || src_params.empty()) throw InvalidArgumentException("Please provide a valid SRC value for "+where, AT);
}

void EditingMove::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (std::set<std::string>::const_iterator it_set=src_params.begin(); it_set != src_params.end(); ++it_set) { //loop over the parameters to move
			const size_t src_index = vecMeteo[station].front().getParameterIndex( *it_set );
			if (src_index == IOUtils::npos) continue; //no such parameter for this station, skipping
			
			//the next two lines are required to offer time restrictions
			for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
				for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
					
					const size_t dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
					if (vecMeteo[station][jj]( dest_index ) == IOUtils::nodata) {
						vecMeteo[station][jj]( dest_index ) = vecMeteo[station][jj]( src_index );
						vecMeteo[station][jj]( src_index ) = IOUtils::nodata;
					}
				}
			}
		}
	}
}


////////////////////////////////////////////////// EXCLUDE
EditingExclude::EditingExclude(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), exclude_params()
{
	parse_args(vecArgs);
}

void EditingExclude::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_excludes=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="EXCLUDE") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), exclude_params);
			has_excludes = true;
		}
	}

	if (!has_excludes || exclude_params.empty()) throw InvalidArgumentException("Please provide a valid EXCLUDE value for "+where, AT);
}

void EditingExclude::processStation(METEO_SET& vecMeteo, const size_t& startIdx, const size_t& endIdx, const std::set< std::string >& params)
{
	for (size_t jj = startIdx; jj < endIdx; ++jj) { //loop over the timesteps
		for (std::set<std::string>::const_iterator it_set=params.begin(); it_set != params.end(); ++it_set) { //loop over the parameters to exclude
			const std::string param( *it_set );
			if (vecMeteo[jj].param_exists(param))
				vecMeteo[jj](param) = IOUtils::nodata;
		}
	}
}

void EditingExclude::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			processStation(vecMeteo[station], editPeriod.getStart(), editPeriod.getEnd(), exclude_params);
		}
	}
}


////////////////////////////////////////////////// KEEP
EditingKeep::EditingKeep(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), keep_params()
{
	parse_args(vecArgs);
}

void EditingKeep::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_excludes=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="KEEP") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), keep_params);
			has_excludes = true;
		}
	}

	if (!has_excludes || keep_params.empty()) throw InvalidArgumentException("Please provide a valid KEEP value for "+where, AT);
}

void EditingKeep::processStation(METEO_SET& vecMeteo, const size_t& startIdx, const size_t& endIdx, const std::set< std::string >& params)
{
	for (size_t jj=startIdx; jj<endIdx; ++jj) {//loop over the timesteps
		MeteoData& md_ref( vecMeteo[jj] );
		MeteoData md( md_ref );
		md.reset(); //delete all meteo fields

		for (std::set<std::string>::const_iterator it_set=params.begin(); it_set != params.end(); ++it_set) { //loop over the parameters to keep
			const std::string param( *it_set);
			if (!md.param_exists(param)) continue;
			md(param) = md_ref(param);
		}

		//copy back the new object into vecMeteo
		md_ref = md;
	}
}

void EditingKeep::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			processStation(vecMeteo[station], editPeriod.getStart(), editPeriod.getEnd(), keep_params);
		}
	}
}


////////////////////////////////////////////////// AUTOMERGE
EditingAutoMerge::EditingAutoMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), merge_strategy(MeteoData::FULL_MERGE), merge_conflicts(MeteoData::CONFLICTS_PRIORITY)
{
	parse_args(vecArgs);
}

void EditingAutoMerge::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="MERGE_STRATEGY") {
			merge_strategy = MeteoData::getMergeType( vecArgs[ii].second );
		} else if (vecArgs[ii].first=="MERGE_CONFLICTS") {
			merge_conflicts = MeteoData::getMergeConflicts( vecArgs[ii].second );
		}
	}
}

void EditingAutoMerge::mergeStations(const size_t& toStationIdx, STATIONS_SET& vecStation)
{
	const std::string toStationID( vecStation[toStationIdx].stationID );
	
	//stations before toStationIdx are not == stationID both for the "*" station and for any specific stationID
	for (size_t jj=toStationIdx+1; jj<vecStation.size(); jj++) { //loop over the stations
		const std::string fromStationID( IOUtils::strToUpper(vecStation[jj].stationID) );
		if (fromStationID==toStationID) {
			vecStation[toStationIdx].merge( vecStation[jj] );
			std::swap( vecStation[jj], vecStation.back() );
			vecStation.pop_back();
			jj--; //we need to redo the current jj, because it contains another station
		}
	}
}

void EditingAutoMerge::editTimeSeries(STATIONS_SET& vecStation)
{
	if (stationID=="*") {
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			mergeStations(ii, vecStation);
		}
	} else {
		//find our current station in vecStation
		size_t toStationIdx=IOUtils::npos;
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (vecStation[ii].stationID==stationID) {
				toStationIdx = ii;
				break;
			}
		}
		
		if (toStationIdx==IOUtils::npos) return;
		mergeStations(toStationIdx, vecStation);
	}
}

void EditingAutoMerge::mergeMeteo(const size_t& toStationIdx, std::vector<METEO_SET>& vecMeteo) const
{
	const std::string toStationID( vecMeteo[toStationIdx].front().getStationID() );
	size_t nr_conflicts = 0;
	
	for (size_t jj=toStationIdx+1; jj<vecMeteo.size(); jj++) { //loop over the stations
		if (vecMeteo[jj].empty())  continue;
		const std::string fromStationID( IOUtils::strToUpper(vecMeteo[jj].front().getStationID()) );
		if (fromStationID==toStationID) {
			nr_conflicts += MeteoData::mergeTimeSeries(vecMeteo[toStationIdx], vecMeteo[jj], merge_strategy, merge_conflicts); //merge timeseries for the two stations
			std::swap( vecMeteo[jj], vecMeteo.back() );
			vecMeteo.pop_back();
			jj--; //we need to redo the current jj, because it contains another station
		}
	}
	
	if (nr_conflicts>0) std::cerr << "[E] " << nr_conflicts << " automerge conflicts on station " <<  toStationID << "\n";
}

void EditingAutoMerge::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	if (stationID=="*") {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			mergeMeteo(ii, vecMeteo);
		}
	} else {
		//find our current station in vecMeteo
		size_t toStationIdx=IOUtils::npos;
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (vecMeteo[ii].front().getStationID()==stationID) {
				toStationIdx = ii;
				break;
			}
		}
		
		if (toStationIdx==IOUtils::npos) return;
		//stations before toStationIdx are not == stationID, see above
		mergeMeteo(toStationIdx, vecMeteo);
	}
}


////////////////////////////////////////////////// MERGE
EditingMerge::EditingMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), merged_stations(), merged_params(), merge_strategy(MeteoData::EXPAND_MERGE), merge_conflicts(MeteoData::CONFLICTS_PRIORITY)
{
	if (i_stationID=="*")
		throw InvalidArgumentException("It is not possible to do a MERGE on the '*' stationID", AT);
	
	parse_args(vecArgs);
}

void EditingMerge::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="MERGE") {
			IOUtils::readLineToVec( IOUtils::strToUpper(vecArgs[ii].second), merged_stations);
		} else if (vecArgs[ii].first=="PARAMS") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), merged_params);
		} else if (vecArgs[ii].first=="MERGE_STRATEGY") {
			merge_strategy = MeteoData::getMergeType( vecArgs[ii].second );
		} else if (vecArgs[ii].first=="MERGE_CONFLICTS") {
			merge_conflicts = MeteoData::getMergeConflicts( vecArgs[ii].second );
		}
	}
	
	//check that each station ID to merge from is only included once
	const std::set<std::string> tmp(merged_stations.begin(), merged_stations.end());
	if (tmp.size()<merged_stations.size())
		throw InvalidArgumentException("Each station to merge from can only appear once in the list for "+where, AT);
	
	//check that the station does not merge with itself
	if (tmp.count(stationID)>0)
		throw InvalidArgumentException("A station can not merge with itself! Wrong argument in "+where, AT);
	
	if (merged_stations.empty()) throw InvalidArgumentException("Please provide a valid MERGE value for "+where, AT);
}

void EditingMerge::editTimeSeries(STATIONS_SET& vecStation)
{
	//find our current station in vecStation
	size_t toStationIdx=IOUtils::npos;
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		if (vecStation[ii].stationID==stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx==IOUtils::npos) return;
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (vecStation[ii].stationID==fromStationID)
				vecStation[toStationIdx].merge( vecStation[ii] );
		}
	}
}

void EditingMerge::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	//find our current station in vecStation
	size_t toStationIdx = IOUtils::npos;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty()) continue;
		if (vecMeteo[ii].front().getStationID() == stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx == IOUtils::npos) return;
	
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (vecMeteo[ii].front().getStationID() != fromStationID) continue;
			
			if (merged_params.empty()) { //merge all parameters
				MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], vecMeteo[ii], merge_strategy, merge_conflicts );
			} else { //only merge some specific parameters
				//apply a KEEP to a temporary copy of the vector to merge from
				std::vector<MeteoData> tmp_meteo( vecMeteo[ii] );
				EditingKeep::processStation(tmp_meteo, 0, tmp_meteo.size(), merged_params);
				
				MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], tmp_meteo, merge_strategy, merge_conflicts );
			}
		}
	}
}

std::set<std::string> EditingMerge::getDependencies() const
{
	const std::set<std::string> stations_to_purge(merged_stations.begin(), merged_stations.end());
	return stations_to_purge;
}


////////////////////////////////////////////////// COPY
EditingCopy::EditingCopy(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), dest_param(), src_param()
{
	parse_args(vecArgs);
}

void EditingCopy::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param ); //HACK not needed if vecArgs is prepared upper case
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param ),
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingCopy::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
				
				const size_t dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
				
				if (vecMeteo[station][jj].param_exists(src_param))
					vecMeteo[station][jj](dest_index) = vecMeteo[station][jj](src_param);
			}
		}
	}
}

} //end namespace
