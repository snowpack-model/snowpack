// SPDX-License-Identifier: LGPL-3.0-or-later
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
#ifndef DATAEDITINGALGS_H
#define DATAEDITINGALGS_H

#include <meteoio/IOInterface.h>
#include <meteoio/meteoFilters/TimeFilters.h>

#include <map>
#include <set>
#include <string>

namespace mio {

/** 
 * @class EditingBlock
 * @brief Base class for DataEditing algorithms
 */
class EditingBlock {
	public:
		EditingBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual ~EditingBlock() {}
		
		/**
		 * @brief Apply this editing block
		 * @details This applies the editing block for its station that has been declared in the constructor
		 * on the provided MeteoData timeseries.
		 * @param vecMeteo MeteoData timeseries for all stations
		 */
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo) {(void)vecMeteo;}
		
		/**
		 * @brief Apply this editing block to the StationData
		 * @details This applies the editing block for its station that has been declared in the constructor
		 * on the provided StationData timeseries.
		 * @param vecStation StationData timeseries for all stations
		 */
		virtual void editTimeSeries(STATIONS_SET& vecStation) {(void)vecStation;}
		
		/**
		 * @brief Get the station IDs this editing block depends on for this station
		 * @return a set of station IDs it depends on
		 */
		virtual std::set<std::string> requiredIDs() const {return std::set<std::string>();}
		
		/**
		 * @brief Get the station IDs this editing block provides based on this station
		 * @return a set of station IDs it provides
		 */
		virtual std::set<std::string> providedIDs() const {return std::set<std::string>();}
		
		/**
		 * @brief Get the station IDs to purge after using them for this station ID
		 * @return a set of station IDs to purge after processing
		 */
		virtual std::set<std::string> purgeIDs() const {return requiredIDs();}
		
		const std::string toString() const;
		
	protected:
		std::string getName() const {return block_name;}
		static std::set<std::string> initStationSet(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& keyword);
		bool skipStation(const std::vector<MeteoData>& vecMeteo) const;
		METEO_SET timeFilterFromStation(const METEO_SET& vecMeteo) const; //merge and automerge need this method
		
		const std::set<std::string> excluded_stations, kept_stations;
		const std::vector<DateRange> time_restrictions;
		const std::string stationID, block_name;
};

/** 
 * @class EditingSwap
 * @ingroup processing
 * @brief SWAP input editing command
 * @details
 * It is possible to swap pairs of parameters with the SWAP key. This supports both standard \ref meteoparam "meteorological parameters" as well
 * as non-standard parameters (ie not in the list in the link). If a parameter does not exists, it will be transparently added with a nodata value.
 * 
 * @code
 * [InputEditing]
 * FLU2::edit2      = SWAP
 * FLU2::arg2::dest = ISWR
 * FLU2::arg2::src  = RSWR
 * @endcode
 */
class EditingSwap : public EditingBlock {
	public:
		EditingSwap(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::string dest_param, src_param;
};


/** 
 * @class EditingMove
 * @ingroup processing
 * @brief MOVE input editing command
 * @details
 * It is possible to rename a meteorological parameter thanks to the MOVE key. This key can take 
 * multiple source names that will be processed in the order of declaration. Original names that are not found in the current
 * dataset will silently be ignored, so it is safe to provide a list that contain many possible names:
 * 
 * @code
 * [InputEditing]
 * SLF2::edit1      = MOVE
 * SLF2::arg1::dest = TA
 * SLF2::arg1::src  = air_temp air_temperature temperature_air
 * @endcode
 * 
 * This can be used to rename non-standard parameter names into standard ones. In this example, if TA already had some values, it will keep
 * those and only points not having a value will be filled by either air_temp or air_temperature or temperature_air (the first one in
 * the list to have a value has the priority).
 */
class EditingMove : public EditingBlock {
	public:
		EditingMove(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::set< std::string > src_params;
		std::string dest_param;
};

/** 
 * @class EditingExclude
 * @ingroup processing
 * @brief EXCLUDE input editing command
 * @details
 * It is possible to exclude specific parameters with the "exclude" command. This is done by providing a space delimited list of 
 * \ref meteoparam "meteorological parameters" to exclude for the station with the EXCLUDE argument (the exact opposite can also be done, excluding 
 * ALL parameters except the ones declared with the "keep" command). Giving as parameters the wildcard character \em * deletes the whole timestamp.
 *
 * @code
 * [InputEditing]
 * FLU2::edit3         = EXCLUDE
 * FLU2::arg3::params  = TA RH TSS TSG
 * 
 * SLF2::edit1         = EXCLUDE
 * SLF2::arg3::params  = *
 * SLF2::arg3::when    = 2020-09-01 - 2020-09-03
 * @endcode
 */
class EditingExclude : public EditingBlock {
	public:
		EditingExclude(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void processStation(METEO_SET& vecMeteo, const size_t& startIdx, const size_t& endIdx, const std::set< std::string >& params) const;
		std::set< std::string > exclude_params;
		bool wildcard;
};

/** 
 * @class EditingKeep
 * @ingroup processing
 * @brief KEEP input editing command
 * @details
 * It is possible to exclude ALL parameters except the ones declared with the "keep" command (it is the exact opposite of the
 * EditingExclude command). If such a command has been used for the wildcard station ID '*', it will be applied first so 
 * it is not possible to "keep" some parameters that have already been excluded with the wildcard station ID. 
 * Here below an example relying on wildcards:
 * @code
 * [InputEditing]
 * *::edit1        = KEEP           ;all stations will keep TA, RH and HS and reject the other parameters
 * *::arg1::params = TA RH HS
 * 
 * WFJ2::edit1        = KEEP        ;WFJ2 will only keep HS since ISWR has been removed by the '*' station above
 * WFJ2::arg1::params = HS ISWR
 * @endcode
 */
class EditingKeep : public EditingBlock {
	public:
		EditingKeep(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		static void processStation(METEO_SET& vecMeteo, const size_t& startIdx, const size_t& endIdx, const std::set< std::string >& params); //for use in DataEditingAlgorithms
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::set< std::string > keep_params;
};

/** 
 * @class EditingMerge
 * @ingroup processing
 * @brief MERGE input editing command
 * @details
 * It is possible to merge different data sets together, with the MERGE command. This is useful, for example, to 
 * provide measurements from different stations that actually share the same measurement location or to build 
 * "composite" station from multiple real stations (in this case, using EXCLUDE and/or KEEP commands to fine tune 
 * how the composite station(s) is/are built). It is also possible to restrict which parameters are to be merged
 * with the PARAMS keyword (providing a space-delimited list of parameters).
 * 
 * Please note that the order of declaration defines the priority (ie the first station that has a value for a given 
 * parameter has priority). Please also note that which timestamps will be merged depends on the chosen merge 
 * strategy with the MERGE_STRATEGY option (see MeteoData::Merge_Type, by default it is MeteoData::EXPAND_MERGE). The handling 
 * of merge conflicts can be configured with the MERGE_CONFLICTS optional argument (see MeteoData::Merge_Conflicts, 
 * by default it is MeteoData::CONFLICTS_PRIORITY). Furthermore, a station can be merged into multiple other stations, 
 * but circular dependencies are prohibited (and checked for).
 *
 * @code
 * [Input]
 * METEO = SMET
 * METEOPATH = ./input
 * STATION1  = STB
 * STATION2  = WFJ2
 * STATION3  = WFJ1
 * STATION4  = DAV1
 * [...]
 *
 * [InputEditing]
 * STB::edit1        = EXCLUDE
 * STB::arg1::params = ILWR PSUM
 * 
 * WFJ2::edit1        = KEEP
 * WFJ2::arg1::params = PSUM ILWR RSWR
 *
 * STB::edit2       = MERGE
 * STB::arg2::merge = WFJ2 WFJ1
 * STB::arg2::merge_strategy = FULL_MERGE
 * 
 * DAV1::edit1        = MERGE
 * DAV1::arg1::merge  = WFJ2
 * DAV1::arg1::params = HS RSWR PSUM
 * @endcode
 * 
 * @note One limitation when handling "extra" parameters (ie parameters that are not in the default \ref meteoparam) is that these extra
 * parameters must be known from the beginning. So if station2 appears later in time with extra parameters, make sure that the buffer size
 * is large enough to reach all the way to this new station (by setting General::BUFFER_SIZE at least to the number of days from
 * the start of the first station to the start of the second station)
 */
class EditingMerge : public EditingBlock {
	public:
		EditingMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		virtual void editTimeSeries(STATIONS_SET& vecStation);
		
		std::set<std::string> requiredIDs() const;
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::vector< std::string > merged_stations;
		std::set< std::string > merged_params;
		MeteoData::Merge_Type merge_strategy;
		MeteoData::MERGE_CONFLICTS merge_conflicts;
};

/** 
 * @class EditingAutoMerge
 * @ingroup processing
 * @brief AUTOMERGE input editing command
 * @details
 * This is a special case of merge: only station's have the exact same ID will get merge together. This is useful when reading data
 * for the same station from multiple source in order to rebuild a consistent dataset. If merge conflicts are encountered (such as 
 * identical fields having different values at the same timestamp), warnings will be printed out and the chosen 
 * conflict resolution (provided by the MERGE_CONFLICTS option) will be used (default: MeteoData::CONFLICTS_AVERAGE). 
 * By default, it does a MeteoData::FULL_MERGE but it is possible to provide a different type of merge with 
 * the MERGE_STRATEGY option.
 * 
 * @code
 * [InputEditing]
 * *::edit1 = AUTOMERGE                        ;all stations having the same ID will be merged together
 * @endcode
 */
class EditingAutoMerge : public EditingBlock {
	public:
		EditingAutoMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		virtual void editTimeSeries(STATIONS_SET& vecStation);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		static void mergeStations(const size_t& toStationIdx, STATIONS_SET& vecStation);
		void mergeMeteo(const size_t& toStationIdx, std::vector<METEO_SET>& vecMeteo) const;
		MeteoData::Merge_Type merge_strategy;
		MeteoData::MERGE_CONFLICTS merge_conflicts;
};

/** 
 * @class EditingCopy
 * @ingroup processing
 * @brief COPY input editing command
 * @details
 * It is also possible to duplicate a meteorological parameter as another meteorological parameter. This is done with the COPY command, 
 * such as:
 * 
 * @code
 * [InputEditing]
 * DAV::edit1      = COPY
 * DAV::arg1::dest = TA_copy
 * DAV::arg1::src  = TA
 * @endcode
 * 
 * This creates a new parameter TA_copy that starts as an exact copy of the raw data of TA, for the DAV station. This newly created parameter is
 * then processed as any other meteorological parameter (thus going through filtering, generic processing, spatial interpolations). 
 */
class EditingCopy : public EditingBlock {
	public:
		EditingCopy(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::string dest_param, src_param;
};

/** 
 * @class EditingCreate
 * @ingroup processing
 * @brief CREATE input editing command
 * @details
 * By calling a choice of algorithms, it is possible to convert a parameter into another one (for example, a dew point
 * temperature to a relative humidity), to generate a parameter thanks to a parametrization based on other parameters
 * (for example ILWR based on TA, RH and ISWR) or to generate synthetic values (for example, a purely yearly sinusoidal
 * variation for TA). This input editing command takes two fixed options as well as an undefined number of options 
 * depending on the \ref generators "generator algorithm" that has been chosen:
 *     - ALGORITHM: specify which \ref generators "generator algorithm" to use;
 *     - PARAM: provides the meteorological parameter to generate values for 
 * (either a MeteoData::Parameters or any other, non-standard name).
 * 
 * Then the arguments of the chosen data generator algorithm must also be provided (according to its documentation).
 * 
 * If the destination parameter does not exists, it will be created. Otherwise, any pre-existing data is kept and only 
 * missing values in the original data set are filled with the generated values, keeping the original sampling rate. As the 
 * available algorithms are the same as for the data generators, they are listed in the 
 * \ref generators_keywords "data generators section" (but the data creators must be declared here in the [InputEditing] section
 * as part of the Input Data Editing stack).
 * 
 * @code
 * [InputEditing]
 * DAV::edit1           = CREATE
 * DAV::arg1::algorithm = CST	;use a constant data generator
 * DAV::arg1::param     = RH	;generate data for the Relative Humidity (RH)
 * DAV::arg1::value     = 0.7	;generate a constant value of 0.7 whenever there is no other data for RH
 * @endcode
 */
class EditingCreate : public EditingBlock {
	public:
		EditingCreate(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		static const std::vector< std::pair<std::string, std::string> > cleanGeneratorArgs(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		
		const Config &cfg_copy;
		const std::vector< std::pair<std::string, std::string> > vecArgs_copy;
		std::string algorithm, dest_param;
};

/** 
 * @class EditingMetadata
 * @ingroup processing
 * @brief METADATA input editing command
 * @details
 * This Input Data Editing algorithm allows to manipulate the metadata of a given station. it takes the following arguments:
 *     - NAME: set a new station name;
 *     - ID: set a new station ID;
 *     - LATITUDE: set a new latitude (then the longitude argument MUST also be provided);
 *     - LONGITUDE: set a new longitude (then the latitude argument MUST also be provided);
 *     - ALTITUDE: set a new altitude;
 *     - SLOPE: set a new slope (then the aziumth argument MUST also be provided);
 *     - AZIMUTH: set a new azimuth (then the slope argument MUST also be provided).
 * 
 * Please note that setting a new ID in effect creates a new station. If there are no time restrictions, the old station will
 * disappear and be replaced by the new one. If there are time restrictions, both stations will remain side by side, although with
 * an obviously different time coverage. Since there is a dependency resolution, you can declare some Input Data Editing on the
 * new station ID, it will be applied properly and any editing declared *after* the ID editing command will only apply for data
 * outside of time restrictions (if any) for the renaming command. When creating a new station, it makes sense to also provide the
 * geographic coordinates for it...
 * 
 * @code
 * SLF2::edit1 = METADATA
 * SLF2::arg1::altitude = 1560     ;set the altitude to 1560m
 * 
 * 
 * FLU2::edit1 = METADATA
 * FLU2::arg1::id = TST2
 * FLU2::arg1::WHEN = 2020-01-01 - 2020-02-01	;all data in this time range will move to the new TST2 station
 * 
 * FLU2::edit2 = EXCLUDE            ;this will only be applied to data outside the above-defined range
 * FLU2::arg2::params = TA
 * 
 * TST2::edit1 = SWAP               ; data in the above-defined range will come to this new station
 * TST2::arg1::dest = TA1           ; and then this SWAP will be applied
 * TST2::arg1::src = TA2
 * @endcode
 */
class EditingMetadata : public EditingBlock {
	public:
		EditingMetadata(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
		virtual void editTimeSeries(STATIONS_SET& vecStation);
		
		std::set<std::string> providedIDs() const;
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void mergeMigratedData(std::vector<METEO_SET>& vecMeteo, const std::vector<METEO_SET>& vecTmp) const;
		
		std::string new_name, new_id;
		double lat, lon, alt, slope, azi;
		bool edit_in_place;
};

class EditingBlockFactory {
	public:
		static EditingBlock* getBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &cfg);
};

} //namespace

#endif
