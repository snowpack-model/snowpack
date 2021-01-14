// SPDX-License-Identifier: LGPL-3.0-or-later
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

#include <meteoio/DataGenerator.h>
#include <meteoio/MeteoProcessor.h> //required to provide RestrictionsIdx

using namespace std;

namespace mio {
const std::string DataGenerator::cmd_section( "GENERATORS" );
const std::string DataGenerator::cmd_pattern( "::GENERATOR" );
const std::string DataGenerator::arg_pattern( "::ARG" );

DataGenerator::DataGenerator(const Config& cfg)
              : mapAlgorithms(), data_qa_logs(false)
{
	cfg.getValue("DATA_QA_LOGS", "GENERAL", data_qa_logs, IOUtils::nothrow);
	
	const std::set<std::string> set_of_used_parameters( getParameters(cfg) );

	std::set<std::string>::const_iterator it;
	for (it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		const std::string parname( *it );
		mapAlgorithms[parname] = buildStack(cfg, parname); //a stack of all generators for this parameter
	}
}

DataGenerator::~DataGenerator()
{ //we have to deallocate the memory allocated by "new GeneratorAlgorithm()"
	std::map< std::string, std::vector<GeneratorAlgorithm*> >::iterator it;

	for (it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		std::vector<GeneratorAlgorithm*> &vec( it->second );
		for (size_t ii=0; ii<vec.size(); ii++)
			delete vec[ii];
	}
}

DataGenerator& DataGenerator::operator=(const DataGenerator& source)
{
	if (this != &source) {
		mapAlgorithms = source.mapAlgorithms;
		data_qa_logs = source.data_qa_logs;
	}
	return *this;
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecMeteo vector containing one point for each station
 */
void DataGenerator::fillMissing(METEO_SET& vecMeteo) const
{
	if (mapAlgorithms.empty()) return; //no generators defined by the end user

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for (it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators( it->second );

		for (size_t station=0; station<vecMeteo.size(); ++station) { //process this parameter on all stations
			const size_t param = vecMeteo[station].getParameterIndex(it->first);
			if (param==IOUtils::npos) continue;

			const std::string statID( vecMeteo[station].meta.getStationID() );
			//these are only required by data_qa_logs
			const double old_val = vecMeteo[station](param);
			const std::string statName( vecMeteo[station].meta.getStationName() );
			const std::string stat = (!statID.empty())? statID : statName;

			bool status = false;
			size_t jj=0;
			while (jj<vecGenerators.size() && status != true) { //loop over the generators
				if (!vecGenerators[jj]->skipStation( statID ) && !vecGenerators[jj]->skipTimeStep( vecMeteo.front().date ) ) {
					status = vecGenerators[jj]->generate(param, vecMeteo[station]);
					if (vecMeteo[station](param) != old_val) {
						vecMeteo[station].setGenerated(param);
						if (data_qa_logs) {
							const std::string parname( it->first );
							const std::string algo_name( vecGenerators[jj]->getAlgo() );
							const Date date( vecMeteo[station].date );
							cout << "[DATA_QA] Generating " << stat << "::" << parname << "::" << algo_name << " " << date.toString(Date::ISO_TZ) << " [" << date.toString(Date::ISO_WEEK) << "]\n";
						}
					} //endif new=old
				}
				jj++;
			}
		}
	}
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecVecMeteo vector containing a timeserie for each station
 */
void DataGenerator::fillMissing(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (mapAlgorithms.empty()) return; //no generators defined by the end user

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for (it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators( it->second );

		for (size_t station=0; station<vecVecMeteo.size(); ++station) { //process this parameter on all stations
			if (vecVecMeteo[station].empty()) continue; //the station does not have any data
			
			const size_t param = vecVecMeteo[station][0].getParameterIndex(it->first);
			if (param==IOUtils::npos) continue;

			const std::string statID( vecVecMeteo[station][0].meta.getStationID() );
			//these are only required by data_qa_logs
			const METEO_SET old_val( vecVecMeteo[station] );
			const std::string statName( old_val[0].meta.getStationName() );
			const std::string stat = (!statID.empty())? statID : statName;

			bool status = false;
			size_t jj=0;
			while (jj<vecGenerators.size() && status != true) { //loop over the generators
				if (!vecGenerators[jj]->skipStation( statID )) {
					
					//loop over time restrictions periods
					status = true; //so if any time restriction period returns false, status will be set to false
					for (RestrictionsIdx editPeriod(vecVecMeteo[station], vecGenerators[jj]->getTimeRestrictions()); editPeriod.isValid(); ++editPeriod) 
						status &= vecGenerators[jj]->create(param, editPeriod.getStart(), editPeriod.getEnd(), vecVecMeteo[station]);
					
					//compare the resulting data with the original copy to see if there are some changes for DATA_QA
					for (size_t kk=0; kk<old_val.size(); kk++) {
						if (old_val[kk](param) != vecVecMeteo[station][kk](param)) {
							vecVecMeteo[station][kk].setGenerated(param);
							if (data_qa_logs) {
								const std::string parname( it->first );
								const std::string algo_name( vecGenerators[jj]->getAlgo() );
								cout << "[DATA_QA] Generating " << stat << "::" << parname << "::" << algo_name << " " << old_val[kk].date.toString(Date::ISO_TZ) << "\n";
							}
						}
					}
				}
				jj++;
			}
		}
	}
}

/**
 * @brief Build a list of station IDs that will be edited
 * @param[in] cfg Config object to read the configuration from
 * @return set of station IDs
 */
std::set<std::string> DataGenerator::getParameters(const Config& cfg)
{
	const std::vector<std::string> vec_keys( cfg.getKeys(cmd_pattern, cmd_section, true) );

	std::set<std::string> set_stations;
	for (size_t ii=0; ii<vec_keys.size(); ++ii){
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			if (vec_keys[ii].length()<=(found+2))
				throw InvalidFormatException("Invalid syntax: \""+vec_keys[ii]+"\"", AT);
			if (vec_keys[ii][found+1]!=':')
				throw InvalidFormatException("Missing ':' in \""+vec_keys[ii]+"\"", AT);
				
			const std::string tmp( vec_keys[ii].substr(0,found) );
			set_stations.insert(tmp); //we keep the case of the parameters
		}
	}

	return set_stations;
}

/**
 * @brief For a given parameter name, build the stack of GeneratorAlgorithm
 * @param[in] cfg Config object to read the configuration from
 * @param[in] parname the parameter to process
 * @return vector of GeneratorAlgorithm* to process in this order
 */
std::vector< GeneratorAlgorithm* > DataGenerator::buildStack(const Config& cfg, const std::string& parname)
{
	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecGenerators( cfg.getValues(parname+cmd_pattern, cmd_section) );
	std::vector< GeneratorAlgorithm* > generators_stack;
	generators_stack.reserve( vecGenerators.size() );
	
	for (size_t ii=0; ii<vecGenerators.size(); ii++) {
		const std::string cmd_name( IOUtils::strToUpper( vecGenerators[ii].second ) );
		if (cmd_name=="NONE") continue;
		
		const unsigned int cmd_nr = Config::getCommandNr(cmd_section, parname+cmd_pattern, vecGenerators[ii].first);
		const std::vector< std::pair<std::string, std::string> > vecArgs( cfg.parseArgs(cmd_section, parname, cmd_nr, arg_pattern) );
		generators_stack.push_back( GeneratorAlgorithmFactory::getAlgorithm( cfg, cmd_name, cmd_section, vecArgs) );
	}
	
	return generators_stack;
}

const std::string DataGenerator::toString() const {
	std::ostringstream os;
	os << "<DataGenerator>\n";
	os << "Generators defined: " << std::boolalpha << !mapAlgorithms.empty() << std::noboolalpha << "\n";
	if (!mapAlgorithms.empty()) {
		os << "User list of generators:\n";
		for (std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
			os << setw(10) << iter->first << " :: ";
			for (size_t jj=0; jj<iter->second.size(); jj++) {
				os << iter->second[jj]->getAlgo() << " ";
			}
			os << "\n";
		}
	}

	os << "</DataGenerator>\n";
	return os.str();
}

} //namespace
