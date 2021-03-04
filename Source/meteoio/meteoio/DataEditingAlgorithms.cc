// SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <meteoio/dataGenerators/GeneratorAlgorithms.h> //required for the CREATE editing

#include <algorithm>
#include <fstream>

using namespace std;

namespace mio {
/**
 * @page data_editing Input Data Editing
 * Before any filters, resampling algorithms or data generators are applied, it is possible to edit the original data. There are several
 * edition commands in the [InputEditing] section that can be stacked at will, per station ID. This is similar to the way that 
 * filters (\ref processing "processing elements") are also stacked together. In order to rectrict any editing to a 
 * specific set of time ranges, use the **when** option followed by a comma delimited list of date intervals (represented 
 * by two ISO formatted dates seperated by ' - ', ie with a space on both sides of the dash), similarly to 
 * the \ref processing "Filters". The general syntax is ('#' represent a number, so each key remains unique):
 * @code
 * {stationID}::edit#            = {command}
 * {stationID}::arg#::{argument} = {values}
 * 
 * #here is an example
 * [InputEditing]
 * WFJ2::edit1         = EXCLUDE
 * WFJ2::arg1::params  = VW DW ISWR RSWR
 * WFJ2::arg1::when    = 2019-12-01T13:00 - 2019-12-25 , 2020-03-05 - 2020-04-15T12:30
 * @endcode
 * 
 * It is also possible to apply a stack of edits to all stations by using the '*' wildcard instead of the station ID 
 * (in such a case, <b>the wildcard stack will be applied before any other stack</b>). Please note that all station IDs 
 * and meteo parameters are handled case insensitive in comparisons. On the wildcard stack, it is possible to restrict the 
 * processing to specific station IDs, using the <i>exclude</i> or <i>only</i> options followed by a list of 
 * station IDs (see an example in the \ref generators "Data Generators" exemples). This is supported automatically by 
 * all generators *except* MERGE as it makes little sense to apply an identical merge to all stations...
 * 
 * The following Input Data Editing commands are available:
 *     - NONE: this is used when importing another ini file to overwrite a command with an empty one
 *     - SWAP: swap two parameters, see EditingSwap
 *     - MOVE: rename one or more parameters into a new name, see EditingMove
 *     - EXCLUDE: delete a list of parameters, see EditingExclude
 *     - KEEP: only keep a list of parameters and reject the others, see EditingKeep
 *     - AUTOMERGE: merge together stations sharing the same station ID, see EditingAutoMerge
 *     - MERGE: merge together one or more stations, see EditingMerge
 *     - COPY: make a copy of a given parameter under a new name, see EditingCopy
 *     - CREATE: fill missing values or create new parameters based on basic transformations or parametrizations, see EditingCreate
 *     - METADATA: edit station's metadata, see EditingMetadata
 *
 * @note It is possible to turn off all input editing for timeseries by setting the *Enable_Timeseries_Editing* key to 
 * false in the [InputEditing] section.
 */

static inline bool IsUndef (const MeteoData& md) { return md.date.isUndef(); }

EditingBlock::EditingBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg) 
             : excluded_stations( initStationSet(vecArgs, "EXCLUDE") ), kept_stations( initStationSet(vecArgs, "ONLY") ), time_restrictions( MeteoProcessor::initTimeRestrictions(vecArgs, "WHEN", "InputEditing::"+name+" for station "+i_stationID, cfg.get("TIME_ZONE", "Input")) ), stationID( IOUtils::strToUpper(i_stationID) ), block_name(name) {}
 
std::set<std::string> EditingBlock::initStationSet(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& keyword)
{
	std::set<std::string> results;
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first==keyword) {
			std::istringstream iss(vecArgs[ii].second);
			std::string word;
			while (iss >> word){
				results.insert( IOUtils::strToUpper(word) );
			}
		}
	}

	return results;
}

/**
* @brief Return true if this station ID should be skipped
* @details This is based on the comparison between the current station ID and the one
* this EditingBlock has been configured for as well as the "only" or "exclude" arguments
* provided by the user in the configuration file. It also checks if the mete timeseries is 
* empty or not.
* @param[in] vecMeteo timeseries for the station to inquire for
* @return true (yes, skip this station ID) or false (please, process this station ID)
*/
bool EditingBlock::skipStation(const std::vector<MeteoData>& vecMeteo) const
{
	if (vecMeteo.empty()) return true;
	
	const std::string currentID( IOUtils::strToUpper( vecMeteo.front().meta.stationID ) );
	if (stationID!="*" && stationID!=currentID) return true; //not the station this EditingBlock has been configured for
	if (excluded_stations.count(currentID)!=0) return true; //the station is in the excluded list -> skip
	if (kept_stations.empty()) return false; //there are no kept stations -> do not skip

	return (kept_stations.count(currentID)==0);
}

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
	} else if (name == "CREATE"){
		return new EditingCreate(i_stationID, vecArgs, name, cfg);
	} else if (name == "METADATA"){
		return new EditingMetadata(i_stationID, vecArgs, name, cfg);
	} else {
		throw IOException("The input data editing block '"+name+"' does not exist! " , AT);
	}
}

/**
 * @brief Prepare a station that will be merged in case of time restrictions
 * @details In order to apply time restrictions for merges, the easiest is to filter
 * the "from" stations, one by one, to only contain data within the time restrictions periods. 
 * Then any merge type works as intented and the process is not too expensive. 
 * 
 * @note this method is only supposed to be called when !time_restrictions.empty()
 * @param[in] vecMeteo timeseries of the station that provides the data on the from side of merge
 * @return timeseries of the "from" station only containing data in the time restrictions periods
 */
METEO_SET EditingBlock::timeFilterFromStation(const METEO_SET& vecMeteo) const
{
	if (time_restrictions.empty()) return vecMeteo;
	METEO_SET vecResults;
	
	//the next two lines are required to offer time restrictions
	for (RestrictionsIdx editPeriod(vecMeteo, time_restrictions); editPeriod.isValid(); ++editPeriod) {
		for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
			vecResults.push_back( vecMeteo[jj] );
		}
	}
	
	return vecResults;
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
			//IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param );
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingSwap::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (skipStation(vecMeteo[station])) continue;
		
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
			//IOUtils::toUpper( dest_param );
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
		if (skipStation(vecMeteo[station])) continue;
		
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
            : EditingBlock(i_stationID, vecArgs, name, cfg), exclude_params(), wildcard(false)
{
	parse_args(vecArgs);
}

void EditingExclude::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_excludes=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="PARAMS") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), exclude_params);
			has_excludes = true;
		}
	}

	if (!has_excludes || exclude_params.empty()) throw InvalidArgumentException("Please provide a valid EXCLUDE value for "+where, AT);
	if (exclude_params.size()==1 && *exclude_params.begin()=="*") wildcard = true;
}

void EditingExclude::processStation(METEO_SET& vecMeteo, const size_t& startIdx, const size_t& endIdx, const std::set< std::string >& params) const
{
	if (wildcard) { //special case for the wildcard parameter
		for (size_t jj = startIdx; jj < endIdx; ++jj) { //loop over the timesteps
			vecMeteo[jj].date.setUndef(true);
		}
	} else {
		for (size_t jj = startIdx; jj < endIdx; ++jj) { //loop over the timesteps
			for (std::set<std::string>::const_iterator it_set=params.begin(); it_set != params.end(); ++it_set) { //loop over the parameters to exclude
				const std::string param( *it_set );
				if (vecMeteo[jj].param_exists(param))
					vecMeteo[jj](param) = IOUtils::nodata;
			}
		}
	}
}

void EditingExclude::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (skipStation(vecMeteo[station])) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			processStation(vecMeteo[station], editPeriod.getStart(), editPeriod.getEnd(), exclude_params);
		}
		
		if (wildcard) 
			vecMeteo[station].erase( std::remove_if(vecMeteo[station].begin(), vecMeteo[station].end(), IsUndef), vecMeteo[station].end());
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
		if (vecArgs[ii].first=="PARAMS") {
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
		if (skipStation(vecMeteo[station])) continue;
		
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
	const std::string toStationID( IOUtils::strToUpper(vecStation[toStationIdx].stationID) );
	
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
			if (IOUtils::strToUpper(vecStation[ii].stationID)==stationID) {
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
	const std::string toStationID( IOUtils::strToUpper(vecMeteo[toStationIdx].front().getStationID()) );
	size_t nr_conflicts = 0;
	
	for (size_t jj=toStationIdx+1; jj<vecMeteo.size(); jj++) { //loop over the stations
		if (vecMeteo[jj].empty())  continue;
		
		const std::string fromStationID( IOUtils::strToUpper(vecMeteo[jj].front().getStationID()) );
		if (fromStationID==toStationID) {
			if (time_restrictions.empty()) {
				nr_conflicts += MeteoData::mergeTimeSeries(vecMeteo[toStationIdx], vecMeteo[jj], merge_strategy, merge_conflicts); //merge timeseries for the two stations
			} else {
				std::vector<MeteoData> tmp_meteo( timeFilterFromStation(vecMeteo[jj]) );
				nr_conflicts += MeteoData::mergeTimeSeries(vecMeteo[toStationIdx], tmp_meteo, merge_strategy, merge_conflicts); //merge timeseries for the two stations
			}
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
			if (skipStation(vecMeteo[ii])) continue;
			mergeMeteo(ii, vecMeteo);
		}
	} else {
		//find our current station in vecMeteo
		size_t toStationIdx=IOUtils::npos;
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (IOUtils::strToUpper(vecMeteo[ii].front().getStationID())==stationID) {
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
		if (IOUtils::strToUpper(vecStation[ii].stationID)==stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx==IOUtils::npos) return;
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (IOUtils::strToUpper(vecStation[ii].stationID)==fromStationID)
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
		if (IOUtils::strToUpper(vecMeteo[ii].front().getStationID()) == stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx == IOUtils::npos) return;
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (IOUtils::strToUpper(vecMeteo[ii].front().getStationID()) != fromStationID) continue;
			
			if (merged_params.empty()) { //merge all parameters
				if (time_restrictions.empty()) {
					MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], vecMeteo[ii], merge_strategy, merge_conflicts );
				} else {
					std::vector<MeteoData> tmp_meteo( timeFilterFromStation(vecMeteo[ii]) );
					MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], tmp_meteo, merge_strategy, merge_conflicts );
				}
			} else { //only merge some specific parameters
				std::vector<MeteoData> tmp_meteo = (time_restrictions.empty())? vecMeteo[ii] : timeFilterFromStation(vecMeteo[ii]);
				
				//apply a KEEP to a temporary copy of the vector to merge from
				//std::vector<MeteoData> tmp_meteo( vecMeteo[ii] );
				EditingKeep::processStation(tmp_meteo, 0, tmp_meteo.size(), merged_params);
				
				MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], tmp_meteo, merge_strategy, merge_conflicts );
			}
		}
	}
}

std::set<std::string> EditingMerge::requiredIDs() const
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
			//IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param );
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingCopy::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (skipStation(vecMeteo[station])) continue;
		
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

////////////////////////////////////////////////// CREATE
EditingCreate::EditingCreate(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), cfg_copy(cfg), vecArgs_copy( cleanGeneratorArgs(vecArgs) ), algorithm(), dest_param()
{
	parse_args(vecArgs);
}

//removes the arguments for the InputDataEditing stage before forwarding them to the dataGenerators 
//since they might otherwise complain of an unknown argument...
const std::vector< std::pair<std::string, std::string> > EditingCreate::cleanGeneratorArgs(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	std::vector< std::pair<std::string, std::string> > vecTmp;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		const std::string argument( vecArgs[ii].first );
		if (argument=="ALGORITHM" || argument=="PARAM") continue;
		
		vecTmp.push_back( vecArgs[ii] );
	}
	
	return vecTmp;
}

void EditingCreate::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="ALGORITHM") {
			IOUtils::parseArg(vecArgs[ii], where, algorithm);
		} else if (vecArgs[ii].first=="PARAM") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
		} 
	}

	if (dest_param.empty()) throw InvalidArgumentException("Please provide a PARAM argument for "+where, AT);
	if (algorithm.empty()) throw InvalidArgumentException("Please provide the ALGORITHM argument for "+where, AT);
}

void EditingCreate::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	GeneratorAlgorithm *creator = GeneratorAlgorithmFactory::getAlgorithm(cfg_copy, algorithm, "InputEditing", vecArgs_copy);
	
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (skipStation(vecMeteo[station])) continue;
		
		//the next lines dealing with RestrictionsIdx are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			size_t dest_index = IOUtils::npos;
			for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
				dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
			}
			if (dest_index != IOUtils::npos)
				creator->create(dest_index, editPeriod.getStart(), editPeriod.getEnd(), vecMeteo[station]);
		}
	}
}

////////////////////////////////////////////////// METADATA
EditingMetadata::EditingMetadata(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg)
            : EditingBlock(i_stationID, vecArgs, name, cfg), new_name(), new_id(), lat(IOUtils::nodata), lon(IOUtils::nodata), alt(IOUtils::nodata), slope(IOUtils::nodata), azi(IOUtils::nodata), edit_in_place(true)
{
	parse_args(vecArgs);
}

void EditingMetadata::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		const std::string cmd( IOUtils::strToUpper(vecArgs[ii].first) );
		if (cmd=="NAME") {
			IOUtils::parseArg(vecArgs[ii], where, new_name);
		} else if (cmd=="ID") {
			IOUtils::parseArg(vecArgs[ii], where, new_id);
		} else if (cmd=="LATITUDE") {
			IOUtils::parseArg(vecArgs[ii], where, lat);
		} else if (cmd=="LONGITUDE") {
			IOUtils::parseArg(vecArgs[ii], where, lon);
		} else if (cmd=="ALTITUDE") {
			IOUtils::parseArg(vecArgs[ii], where, alt);
		} else if (cmd=="SLOPE") {
			IOUtils::parseArg(vecArgs[ii], where, slope);
		} else if (cmd=="AZIMUTH") {
			IOUtils::parseArg(vecArgs[ii], where, azi);
		}
	}

	if ( (lat!=IOUtils::nodata && lon==IOUtils::nodata) || (lat==IOUtils::nodata && lon!=IOUtils::nodata))
		throw InvalidArgumentException("Please provide both latitude and longitude for "+where, AT);
	if ( (slope!=IOUtils::nodata && azi==IOUtils::nodata) || (slope==IOUtils::nodata && azi!=IOUtils::nodata))
		throw InvalidArgumentException("Please provide both slope and azimuth for "+where, AT);
	
	edit_in_place = (new_id.empty() || (!new_id.empty() && time_restrictions.empty()));
}

std::set<std::string> EditingMetadata::providedIDs() const
{
	if (new_id.empty()) return std::set<std::string>();
	
	std::set<std::string> tmp_set;
	tmp_set.insert( new_id );
	return tmp_set;
}

/**
* @brief Insert the migrated timestamps into their proper position
* @details 
* If dispatching some parts of one station to another one through an ID rename, then it is necessary to
* move the data to the new station by inserting new elements into the new station and reseting the original
* element (so it does not contain data anymore). This is handled through a merge: a temporary storage (vecTmp)
* contains all the data to migrate and is then merged with the new station (it gets created if necessary).
* 
* This way, renaming multiple stations' subsets to the same ID works smoothly (although if they overlap, 
* a conflcit resolution will kick in and the user might not be aware of this).
* 
* @param vecMeteo MeteoData timeseries for the new station
* @param vecTmp temportary storage
*/
void EditingMetadata::mergeMigratedData(std::vector<METEO_SET>& vecMeteo, const std::vector<METEO_SET>& vecTmp) const
{
	size_t new_station_pos = IOUtils::npos;
	
	//do we already have a station with this ID?
	for (size_t station=0; station<vecMeteo.size(); ++station) {
		if (vecMeteo[station].empty()) continue;
		if (vecMeteo[station].front().getStationID()==new_id) {
			new_station_pos = station;
			break;
		}
	}
	
	//create a new station
	if (new_station_pos == IOUtils::npos) {
		vecMeteo.push_back( std::vector<MeteoData>() );
		new_station_pos = vecMeteo.size()-1;
	}
	
	for (size_t station=0; station<vecTmp.size(); ++station) {
		MeteoData::mergeTimeSeries(vecMeteo[new_station_pos], vecTmp[station], MeteoData::FULL_MERGE, MeteoData::CONFLICTS_PRIORITY);
	}
}


void EditingMetadata::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	const size_t nrStations = vecMeteo.size();
	std::vector<METEO_SET> vecTmp; //in case of !edit_in_place, temporary storage for the migrated data
	
	//now perform the edition on the data, populating a new station if necessary
	for (size_t station=0; station<nrStations; ++station) { //for each station, excluding the potential new one
		if (!edit_in_place) vecTmp.push_back( std::vector<MeteoData>() );
		if (skipStation(vecMeteo[station])) continue;
		
		//the next two lines are required to offer time restrictions
		for (RestrictionsIdx editPeriod(vecMeteo[station], time_restrictions); editPeriod.isValid(); ++editPeriod) {
			for (size_t jj = editPeriod.getStart(); jj < editPeriod.getEnd(); ++jj) { //loop over the timesteps
				if (!edit_in_place) {
					vecTmp[station].push_back( vecMeteo[station][jj] );
					vecMeteo[station][jj].reset();
				}
				StationData &sd = (!edit_in_place)? vecTmp[station].back().meta : vecMeteo[station][jj].meta;
				
				if (!new_name.empty()) sd.stationName = new_name;
				if (!new_id.empty()) sd.stationID = new_id;
				if (lat!=IOUtils::nodata) sd.position.setLatLon(lat, lon, sd.getAltitude());
				if (alt!=IOUtils::nodata) sd.position.setAltitude(alt, false);
				if (slope!=IOUtils::nodata) sd.setSlope(slope, azi);
			}
		}
	}
	
	if (!edit_in_place) //if migrating data between stations (ie editing the station ID)
		mergeMigratedData(vecMeteo, vecTmp);
}

void EditingMetadata::editTimeSeries(STATIONS_SET& vecStation)
{
	//find the source station
	size_t src_index = IOUtils::npos;
	for (size_t station=0; station<vecStation.size(); ++station) {
		if (vecStation[station].stationID==stationID) {
			src_index = station;
			break;
		}
	}
	
	 //the source station could not be found, no changes to vecStation
	if (src_index == IOUtils::npos) return;
	
	//if a new station must be inserted, search where
	size_t new_station_pos = IOUtils::npos;
	if (!edit_in_place) {
		for (size_t station=0; station<vecStation.size(); ++station) {
			if (vecStation[station].getStationID()==new_id) {
				new_station_pos = station;
				break;
			}
		}
		
		if (new_station_pos == IOUtils::npos) {
			vecStation.push_back( vecStation[src_index] );
			new_station_pos = vecStation.size()-1;
		}
	}
	
	//now do the metadata changes
	StationData &sd = (!edit_in_place)? vecStation[new_station_pos] : vecStation[src_index];
	if (!new_name.empty()) sd.stationName = new_name;
	if (!new_id.empty()) sd.stationID = new_id;
	if (lat!=IOUtils::nodata) sd.position.setLatLon(lat, lon, sd.position.getAltitude());
	if (alt!=IOUtils::nodata) sd.position.setAltitude(alt, false);
	if (slope!=IOUtils::nodata) sd.setSlope(slope, azi);
}

} //end namespace
