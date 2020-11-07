/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef DATAEDITING_H
#define DATAEDITING_H

#include <meteoio/IOInterface.h>
#include <meteoio/DataCreator.h>
#include <meteoio/DataEditingAlgorithms.h>
#include <meteoio/meteoFilters/TimeFilters.h>

#include <map>
#include <set>
#include <string>

namespace mio {

/**
* @file DataEditing.h
* @class DataEditing
* @brief 
* @details This class handles the whole MeteoIO Data Editing step. It builds for each station ID (including the 
* '*' wildcard) a stack of all EditingBlock to apply, resolve the dependencies (ie in which order all the
* station IDs should be processed) and walks through the stacks to apply the processing.
* 
* In the process, it also extract all the arguments for each EditingBlock and store them into a vector of
* key/value pairs that is then provided to the EditingBlock to be parsed. It also purges the station IDs 
* that have been merged into other ones from the final results.
* @author Mathias Bavay
*/
class DataEditing {
	public:
		DataEditing(const Config&);
		
		DataEditing& operator=(const DataEditing&); ///<Assignement operator
		
		virtual ~DataEditing();

		static void purgeTrailingNodata(std::vector<METEO_SET>& vecMeteo);
		
		void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		void editTimeSeries(STATIONS_SET& vecStation);
		
		const std::string toString() const;

		TimeProcStack timeproc;
		
	private:
		static std::set<std::string> getEditableStations(const Config& cfg);
		static std::vector< std::pair<std::string, std::string> > parseArgs(const Config& cfg, const std::string& cmd_key, const std::string& stationID);
		static std::vector< EditingBlock* > buildStack(const Config& cfg, const std::string& station_ID);
		std::map< std::string, std::set<std::string> > getDependencies() const;
		static std::set<std::string> getMergedFromIDs(const std::map< std::string, std::set<std::string> >& dependencies);
		static std::vector<std::string> getProcessingOrder(std::map< std::string, std::set<std::string> > dependencies);
		
		DataCreator dataCreator;
		std::map< std::string, std::vector< EditingBlock* > > editingStack;
		static const std::string command_key, arg_key;
		static const char NUM[];
};

} //namespace

#endif
