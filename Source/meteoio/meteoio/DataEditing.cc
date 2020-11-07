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
#include <meteoio/DataEditing.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/MeteoData.h> //needed for the merge strategies

#include <algorithm>
#include <fstream>

using namespace std;

namespace mio {

const char DataEditing::NUM[] = "0123456789";
const std::string DataEditing::command_key( "::EDIT" );
const std::string DataEditing::arg_key( "::ARG" );

DataEditing::DataEditing(const Config& cfgreader)
           : timeproc(cfgreader), dataCreator(cfgreader), editingStack()
{
	const std::set<std::string> editableStations( getEditableStations(cfgreader) );
	
	for (std::set<std::string>::const_iterator it = editableStations.begin(); it != editableStations.end(); ++it) {
		editingStack[ *it ] = buildStack(cfgreader, *it);
	}
}

DataEditing::~DataEditing() 
{
	std::map< std::string, std::vector< EditingBlock* > >::const_iterator it;
	for (it = editingStack.begin(); it != editingStack.end(); ++it) {
		for (size_t ii=0; ii<it->second.size(); ii++)
			delete it->second[ii];
	}
}

DataEditing& DataEditing::operator=(const DataEditing& source) 
{
	if (this != &source) {
		timeproc = source.timeproc;
		dataCreator = source.dataCreator;
		editingStack = source.editingStack;
	}
	
	return *this;
}

/**
 * @brief Build a list of station IDs that will be edited
 * @param[in] cfg Config object to read the configuration from
 * @return set of station IDs
 */
std::set<std::string> DataEditing::getEditableStations(const Config& cfg)
{
	const std::vector<std::string> vec_keys( cfg.getKeys(command_key, "INPUTEDITING", true) );

	std::set<std::string> set_stations;
	for (size_t ii=0; ii<vec_keys.size(); ++ii){
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			if (vec_keys[ii].length()<=(found+2))
				throw InvalidFormatException("Invalid syntax: \""+vec_keys[ii]+"\"", AT);
			if (vec_keys[ii][found+1]!=':')
				throw InvalidFormatException("Missing ':' in \""+vec_keys[ii]+"\"", AT);
				
			const std::string tmp( vec_keys[ii].substr(0,found) );
			set_stations.insert(tmp);
		}
	}

	return set_stations;
}

/**
 * @brief Extract the arguments for a given station ID and Input Data Editing command and store them into a 
 * vector of key / value pairs.
 * @param[in] cfg Config object to read the configuration from
 * @param[in] cmd_key the base command key, such as "SLF2::EDIT1" that will be used to extract the command number
 * @param[in] stationID the station ID to process
 * @return All arguments for this command and this station, as vector of key / value pairs
 */
std::vector< std::pair<std::string, std::string> > DataEditing::parseArgs(const Config& cfg, const std::string& cmd_key, const std::string& stationID)
{
	//extract the cmd number and perform basic checks on the syntax
	const size_t end_cmd = cmd_key.find(command_key); //we know this will be found since it has been matched in cfg.getValues()
	const size_t start_cmd_nr = cmd_key.find_first_of(NUM, end_cmd+command_key.length());
	const size_t end_cmd_nr = cmd_key.find_first_not_of(NUM, end_cmd+command_key.length());
	if (start_cmd_nr==std::string::npos || end_cmd_nr!=std::string::npos) throw InvalidArgumentException("Syntax error: "+cmd_key, AT);

	unsigned int cmd_nr;
	const std::string cmd_nr_str( cmd_key.substr(start_cmd_nr) );
	if ( !IOUtils::convertString(cmd_nr, cmd_nr_str) ) InvalidArgumentException("Can not parse command number in "+cmd_key, AT);

	//read the arguments and clean them up (ie remove the {stationID}::{args##}:: in front of the argument cmd_key itself)
	std::ostringstream arg_str;
	arg_str << stationID << arg_key << cmd_nr;
	std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getValues(arg_str.str(), "INPUTEDITING") );
	for (size_t jj=0; jj<vecArgs.size(); jj++) {
		const size_t beg_arg_name = vecArgs[jj].first.find_first_not_of(":", arg_str.str().length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+vecArgs[jj].first+"'", AT);
		vecArgs[jj].first = vecArgs[jj].first.substr(beg_arg_name);
	}
	
	return vecArgs;
}

/**
 * @brief For a given station ID, build the stack of EditingBlock
 * @param[in] cfg Config object to read the configuration from
 * @param[in] stationID the station ID to process
 * @return vector of EditingBlock* to process in this order
 */
std::vector< EditingBlock* > DataEditing::buildStack(const Config& cfg, const std::string& station_ID)
{
	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecCommands( cfg.getValues(station_ID+command_key, "INPUTEDITING") );
	std::vector< EditingBlock* > cmd_stack;
	cmd_stack.reserve( vecCommands.size() );
	
	for (size_t ii=0; ii<vecCommands.size(); ii++) {
		const std::string cmd_name( IOUtils::strToUpper( vecCommands[ii].second ) );
		if (cmd_name=="NONE") continue;
		
		const std::vector< std::pair<std::string, std::string> > vecArgs( parseArgs(cfg, vecCommands[ii].first, station_ID) );
		cmd_stack.push_back( EditingBlockFactory::getBlock(station_ID, vecArgs, cmd_name, cfg) );
	}
	
	return cmd_stack;
}

/**
 * @brief Build the list of all station IDs each station ID depends on
 * @return for each station IDs, a set of IDs it depends on
 */
std::map< std::string, std::set<std::string> > DataEditing::getDependencies() const
{
	std::map< std::string, std::set<std::string> > dependencies; //stationID -> set of IDs it depends on
	
	//build the list of stations that are merged to and merged from
	std::map< std::string, std::vector< EditingBlock* > >::const_iterator it_blocks;
	for (it_blocks = editingStack.begin(); it_blocks != editingStack.end(); ++it_blocks) {
		const std::string stat_id( it_blocks->first );
		for (size_t jj=0; jj<it_blocks->second.size(); jj++) {
			const std::set<std::string> tmp_set( it_blocks->second[jj]->getDependencies() );
			if (!tmp_set.empty()) {
				dependencies[ stat_id ].insert(tmp_set.begin(), tmp_set.end());
			} else {
				if (dependencies.count(stat_id)==0) 
					dependencies[stat_id] = std::set<std::string>();
			}
		}
	}
	
	return dependencies;
}

/**
 * @brief Build the list of all station IDs that are merged from and therefore will be purged in the end
 * @param[in] dependencies the map of which station ID depends on which ones
 * @return All station IDs that should be purged
 */
std::set<std::string> DataEditing::getMergedFromIDs(const std::map< std::string, std::set<std::string> >& dependencies)
{
	std::set<std::string> mergedFromIDs;
	
	std::map< std::string, std::set<std::string> >::const_iterator it_deps;
	for (it_deps = dependencies.begin(); it_deps != dependencies.end(); ++it_deps) {
		if (it_deps->second.empty()) continue;
		mergedFromIDs.insert( it_deps->second.begin(), it_deps->second.end() );
	}
	
	return mergedFromIDs;
}

/**
 * @brief Resolve the dependencies in order to find out in which order the station IDs should be processed
 * @details
 * For each station ID declared in [InputEditing], we have a list of stations it depends on (dependencies):
 * Ex: station A -> B, C ; station B -> C, E ; station C -> nothing ; station E -> nothing ; processing_order = []
 * each station that has an empty list (ie does not depends on another) is ready to be processed 
 * and is pushed to processing_order vector as well as removed from the dependencies map.
 * Ex: station A -> B, C ; station B -> C, E ; processing_order = [C, E]
 * Each dependency element in a dependency set that does not have dependencies itself is removed from the set of dependencies of
 * any station that has it (as it does not block processing anymore).
 * Ex: station A -> B ; station B -> nothing ; processing_order = [C, E]
 * Then we redo the whole logic until the dependency map is empty (if nothing gets rmeoved in a round, this means we have a circular dependency)
 * Ex: station A -> nothing ; processing_order = [C, E, B]
 * Then empty dependency map ; processing_order = [C, E, B, A]
 * 
 * @param[in] dependencies the map of which station ID depends on which ones
 * @return station IDs in the order they should be processed
 */
std::vector<std::string> DataEditing::getProcessingOrder(std::map< std::string, std::set<std::string> > dependencies)
{
	std::vector<std::string> processing_order;
	
	//find what is the appropriate processing order
	while (!dependencies.empty()) {
		bool dependency_removed = false;
		std::map< std::string, std::set<std::string> >::iterator it_deps;
		//because 'erase' invalidates the iterators, the increment syntax is a little special...
		for (it_deps = dependencies.begin(); it_deps != dependencies.end(); ) {
			const std::string stat_id( it_deps->first );
			for (std::set<std::string>::iterator it = it_deps->second.begin(); it != it_deps->second.end(); ) {
				//if A depends on B and C but C is not listed in dependencies, remove C from A's list
				if (dependencies.count( *it ) == 0) it_deps->second.erase( it++ );
				else ++it;
			}
			
			//if A does not depend on anything anymore: add it to the processing order and remove it from dependencies
			if (it_deps->second.empty()) {
				processing_order.push_back( stat_id );
				dependencies.erase( it_deps++ );
				dependency_removed = true;
			} else {
				++it_deps;
			}
		}
		
		if (!dependency_removed)
			throw InvalidArgumentException("Potential \'circular merge\', this is not supported (see documentation)", AT);
	}
	
	return processing_order;
}

void DataEditing::editTimeSeries(STATIONS_SET& vecStation)
{
	const std::map< std::string, std::set<std::string> > dependencies( getDependencies() );
	const std::set<std::string> mergedFromIDs( getMergedFromIDs(dependencies) );
	const std::vector<std::string> processing_order( getProcessingOrder(dependencies) );
	
	//process widlcard commands first, knowing that '*' merges are prohibited
	if (editingStack.count("*")>0) {
		for (size_t jj=0; jj<editingStack["*"].size(); jj++)
			editingStack["*"][jj]->editTimeSeries(vecStation);
	}
	
	for (size_t ll=0; ll<processing_order.size(); ll++) {
		const std::string current_ID( processing_order[ll] );
		
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (vecStation[ii].getStationID() == current_ID) {
				if (editingStack.count(current_ID) > 0) {
					for (size_t jj=0; jj<editingStack[current_ID].size(); jj++) {
						editingStack[current_ID][jj]->editTimeSeries(vecStation);
					}
				}
			}
		}
	}
	
	//remove the stations that have been merged into other ones, if necessary
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string stationID( IOUtils::strToUpper(vecStation[ii].stationID) );
		if (mergedFromIDs.count( stationID ) >  0) {
			std::swap( vecStation[ii], vecStation.back() );
			vecStation.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

void DataEditing::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	//TODO handle CREATE command
	const std::map< std::string, std::set<std::string> > dependencies( getDependencies() );
	const std::vector<std::string> processing_order( getProcessingOrder(dependencies) );
	
	//process widlcard commands first, knowing that '*' merges are prohibited
	if (editingStack.count("*")>0) {
		for (size_t jj=0; jj<editingStack["*"].size(); jj++) {
			editingStack["*"][jj]->editTimeSeries(vecMeteo);
		}
	}
	
	//process in the order that has been computed above
	for (size_t ll=0; ll<processing_order.size(); ll++) {
		const std::string current_ID( processing_order[ll] );
		
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (vecMeteo[ii].front().getStationID() == current_ID) {
				if (editingStack.count(current_ID) > 0) {
					for (size_t jj=0; jj<editingStack[current_ID].size(); jj++) {
						editingStack[current_ID][jj]->editTimeSeries(vecMeteo);
					}
				}
			}
		}
	}
	
	//remove the stations that have been merged into other ones, if necessary
	const std::set<std::string> mergedFromIDs( getMergedFromIDs(dependencies) );
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty())  continue;
		const std::string stationID( IOUtils::strToUpper(vecMeteo[ii][0].meta.stationID) );
		if (mergedFromIDs.count( stationID ) >  0) {
			std::swap( vecMeteo[ii], vecMeteo.back() );
			vecMeteo.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
	
	//remove trailing pure nodata MeteoData elements (if any)
	purgeTrailingNodata(vecMeteo);

	timeproc.process(vecMeteo);
	TimeProcStack::checkUniqueTimestamps(vecMeteo);

	dataCreator.createParameters(vecMeteo);
}

//some vector of MeteoData might have trailing elements that are purely nodata
void DataEditing::purgeTrailingNodata(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		//purge trailing nodata
		for (size_t jj=vecMeteo[ii].size(); jj>0; jj--) {
			if (!vecMeteo[ii][jj-1].isNodata()) {
				if (jj!=vecMeteo[ii].size()) vecMeteo[ii].resize( jj );
				break;
			}
		}
	}
}

const std::string DataEditing::toString() const
{
	std::ostringstream os;
	os << "<DataEditing>\n";

	os << dataCreator.toString();

	os << "</DataEditing>\n";
	return os.str();
}

} //end namespace
