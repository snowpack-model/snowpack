// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/SyntheticIO.h>
#include <meteoio/meteoLaws/Atmosphere.h>

#include <regex>

using namespace std;

namespace mio {
/**
 * @page synthio SynthIO
 * This plugin is quite special since it does not read any data but generates <a href="https://en.wikipedia.org/wiki/Synthetic_data">synthetic data</a>. 
 * It is designed for numerical experiments where controlled conditions are applied to the numerical setup. It can either generate timestamps
 * with no associated meteorological parameters or also generate values for any number of meteorological parameters (per station) with 
 * some basic functions. Finer control can be achieved by combining this plugin with \ref data_editing "input data editing" (specially
 * with the \ref EditingCreate "data creators" restricted to date intervals).
 * 
 * @section synthio_keywords Keywords
 * This plugin uses the following keywords, all in the [Input] section:
 * - providing the stations' metadata:
 *     - COORDSYS: coordinate system (see Coords);
 *     - COORDPARAM: extra coordinates parameters (see Coords);
 *     - STATION#: coordinates of the station (mandatory, see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax);
 *     - ID#: the (short) station id to use (optional but recommended, default is "ID_#")
 *     - NAME#: a descriptive station name to use (optional, default is "STATION_#");
 * - controlling the timestamps generation:
 *     - TIME_ZONE: the time zone of any dates that are provided; 
 *     - SYNTH_START: when to start generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_END: when to stop generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_SAMPLING: sampling rate in seconds, starting with SYNTH_START or the requested date if no SYNTH_START is provided (mandatory);
 * - controlling the meteorological parameters generation (optional. If not set up, then only the timestamps will be generated with no meteorological
 * parameters associated, thus relying on a data creator or data generator to define some values):
 *     - STATION#::{parameter}::%TYPE where {parameter} is either coming from \ref meteoparam "MeteoData::meteoparam" or some free text. The TYPE arguments gives the synthesizing function to use:
 *           - CST: <a href="https://en.wikipedia.org/wiki/Constant_function">constant function</a> with the additional key: VALUE;
 *           - STEP: <a href="https://en.wikipedia.org/wiki/Step_function">step function</a> with the additional keys: STEP_DATE (as ISO formatted date), VALUE_BEFORE, VALUE_AFTER;
 *           - RECTANGLE: <a href="https://en.wikipedia.org/wiki/Rectangular_function">rectangle function</a> with the additional keys: VALUE, STEP_START (as ISO formatted date), STEP_STOP (as ISO formatted date), VALUE_STEP;
 *           - STDPRESS: constant, standard atmospheric pressure as a function of the station's altitude (no arguments).
 *
 * @section synthio_examples Example
 * Example of use: create a station TST1 with half hourly sampling rate and a parameter TA that is constant at 270 K until 2022-09-26T12:00:00 
 * when it gets constant at 300 K.
 * @code
 * [INPUT]
 * COORDSYS = CH1903
 * TIME_ZONE = 1.00
 * 
 * METEO = SYNTH
 * SYNTH_SAMPLING = 1800.00
 * STATION1 = latlon (46.75, 9.80, 2200)
 * NAME1 = Davos::test
 * ID1 = TST1
 * STATION1::TA::TYPE = STEP
 * STATION1::TA::VALUE_BEFORE = 270.000000
 * STATION1::TA::STEP_DATE = 2022-09-26T12:00:00
 * STATION1::TA::VALUE_AFTER = 300.000000
 * @endcode
 */

SynthIO::SynthIO(const std::string& configfile) : cfg(configfile), mapSynthesizers(), vecStations(), dt_start(), dt_end(), dt_step(0.), TZ(0.)
{
	init();
}

SynthIO::SynthIO(const Config& cfgreader) : cfg(cfgreader), mapSynthesizers(), vecStations(), dt_start(), dt_end(), dt_step(0.), TZ(0.)
{
	init();
}

SynthIO::~SynthIO()
{
	//deallocate the Synthesizer* pointer in the mapSynthesizers map
	for (auto& station_map : mapSynthesizers) { // loop over the stations
		for (auto& item : station_map.second) { //loop over the parameters
			delete item.second;
		}
	}
}

void SynthIO::init()
{
	//read start / end time as well as sampling rate
	cfg.getValue("TIME_ZONE", "INPUT", TZ);
	cfg.getValue("SYNTH_SAMPLING", "INPUT", dt_step);
	const std::string dt_start_spec = cfg.get("SYNTH_START", "INPUT", "");
	if (!dt_start_spec.empty() && !IOUtils::convertString(dt_start, dt_start_spec, TZ))
		throw InvalidFormatException("Could not process start date "+dt_start_spec+" for the SYNTH plugin", AT);
	const std::string dt_end_spec = cfg.get("SYNTH_END", "INPUT", "");
	if (!dt_end_spec.empty() && !IOUtils::convertString(dt_end, dt_end_spec, TZ))
		throw InvalidFormatException("Could not process start date "+dt_end_spec+" for the SYNTH plugin", AT);
	
	//read the stations' basic metadata
	std::string coordin, coordinparam; //projection parameters
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	std::vector<std::string> vecIDs, vecNames;
	
	const std::vector< std::pair<std::string, std::string> > coords_specs( cfg.getValuesRegex("STATION[0-9]+", "INPUT") );
	//cfg.getValues("STATION", "INPUT", coords_specs);
	cfg.getValues("ID", "INPUT", vecIDs);
	cfg.getValues("NAME", "INPUT", vecNames);
	const bool has_ids = !vecIDs.empty();
	const bool has_names = !vecNames.empty();
	
	//check for consistency
	if (coords_specs.empty()) throw InvalidArgumentException("Please provide at least one STATION# key containing coordinates for the SYNTH plugin", AT);
	if (has_ids && vecIDs.size()!=coords_specs.size()) throw InvalidArgumentException("Please either provide no IDS or the exact same number as STATION for the plugin SYNTH", AT);
	if (has_names && vecNames.size()!=coords_specs.size()) throw InvalidArgumentException("Please either provide no NAMES or the exact same number as STATION for the plugin SYNTH", AT);
	
	//build the stations' metadata
	for (size_t ii=0; ii<coords_specs.size(); ii++) {
		const Coords loc(coordin, coordinparam, coords_specs[ii].second);
		const std::string id = (has_ids)? vecIDs[ii] : "ID_"+IOUtils::toString( ii+1 );
		const std::string name = (has_names)? vecNames[ii] : "STATION_"+IOUtils::toString( ii+1 );
		const StationData sd(loc, id, name);
		
		vecStations.push_back( sd );
		mapSynthesizers[ id ] = getSynthesizer( coords_specs[ii].first, sd );
	}
}

//get all necessary generators for the current station ID identified by its STATION# user-provided key, for all declared MeteoParameters
std::map< std::string, Synthesizer* > SynthIO::getSynthesizer(const std::string& stationRoot, const StationData &sd) const
{
	//extract the meteo parameter name and its subkey
	//we are matching something like "STATION1::TA::VALUE_BEFORE = 270" and want to extract "TA", "VALUE_BEFORE" and "270"
	//to populate mapArgs and construct a map of (parname, Synthesizer*)
	const std::regex parname_regex(stationRoot+"::([^:]+)::([^:]+)");
	const std::vector<std::string> vec_keys( cfg.getKeys(stationRoot, "Input") );
	
	//map of (parname, vector of ( pair(argument, value) )) such as mapArgs["TA"] = <("VALUE_BEFORE","270"), ("VALUE_AFTER","300")>
	std::map< std::string, std::vector< std::pair<std::string, std::string> > > mapArgs;
	//map of (parname, synthetic function name), such as mapTypes["TA"] = "CST"
	std::map< std::string, std::string > mapTypes;
	
	for (auto& key : vec_keys) {
		std::smatch index_matches;
		if (std::regex_match(key, index_matches, parname_regex)) { //retrieve the parname index
			const std::string parname( index_matches.str(1) );
			const std::string subkey( index_matches.str(2) );
			const std::string value( cfg.get(key, "INPUT", "") );
			
			if (subkey=="TYPE") {
				mapTypes[parname] = value;
			} else {
				mapArgs[parname].push_back( std::make_pair(subkey, value) );
			}
		}
	}
	
	std::map< std::string, Synthesizer* > resultsMap;
	for (const auto& item : mapTypes) { //looping over types since not all Synthesizer have arguments
		const std::string parname( item.first );
		resultsMap[parname] = SynthFactory::getSynth( item.second, stationRoot, sd, parname, mapArgs[parname], TZ);
	}
	
	
	return resultsMap;
}

void SynthIO::readStationData(const Date& date, std::vector<StationData>& vecStation)
{
	if (!dt_start.isUndef() && dt_start>date) return;	//requested date is before the user-defined range
	if (!dt_end.isUndef() && dt_end<date) return;	//requested date is after the user-defined range
	
	vecStation = vecStations;
}

void SynthIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	static const double days_to_sec = 24.*3600.;
	const double dt_step_days = dt_step / days_to_sec;
	vecMeteo.clear();
	
	//check if there is anything to deliver or not
	if (!dt_start.isUndef() && dt_start>dateEnd) return;	//requested period is before the user-defined range
	if (!dt_end.isUndef() && dt_end<dateStart) return;	//requested period is after the user-defined range
	
	const Date true_end = (dt_end.isUndef())? dateEnd : std::min(dateEnd, dt_end);
	Date true_start( dateStart );
	
	const double julian_start = (dt_start.isUndef())? dateStart.getJulian(true) : std::max(dateStart.getJulian(true), dt_start.getJulian(true));
	if (!dt_start.isUndef()) {
		if (dt_start<dateStart) {
			const int nr_steps_to_start = static_cast<int>( std::ceil( (julian_start - dt_start.getJulian(true)) / dt_step_days ) );	//watch out, dt_step is in seconds
			const double julian_true_start = dt_start.getJulian(true) + nr_steps_to_start * dt_step_days;
			true_start.setDate( julian_true_start, 0. ); //the julian calculations were done in gmt
		} else if (dt_start>dateStart) {
			true_start.setDate( dt_start );
		}
	}
	
	for (auto &station : vecStations) {
		std::map< std::string, Synthesizer* > *station_synth = &mapSynthesizers[ station.getStationID() ];
		
		std::vector<MeteoData> vecM;
		for (Date dt=true_start; dt<=true_end; dt+=dt_step_days) {
			MeteoData md( dt, station );
			
			for (const auto& item : (*station_synth) ) {
				md( item.first ) = item.second->generate( dt );
			}
			
			vecM.push_back( md );
		}
		vecMeteo.push_back( vecM );
	}
}


///////////////////////////////////////////////////////
//below the constructors (and argument parsing) and generate()  methods for all Synthesizers

CST_Synth::CST_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs) 
          : value(IOUtils::nodata)
{
	const std::string where( "SYNTH::CST, " + station + "::" + parname );
	bool has_value = false;
	
	//parse the arguments
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="VALUE") {
			IOUtils::parseArg(vecArgs[ii], where, value);
			has_value=true;
		}
	}
	
	if (!has_value) throw InvalidArgumentException("Please provide the VALUE argument for the "+where, AT);
}

double CST_Synth::generate(const Date& /*dt*/) const
{
	return value;
}

STEP_Synth::STEP_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ) 
          : dt_step(), value_before(IOUtils::nodata), value_after(IOUtils::nodata)
{
	const std::string where( "SYNTH STEP, " + station + "::" + parname );
	bool has_step_date = false;
	bool has_value_before = false, has_value_after = false;
	
	//parse the arguments
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="STEP_DATE") {
			if (!IOUtils::convertString(dt_step, vecArgs[ii].second, TZ))
				throw InvalidArgumentException("Can not parse argument '"+vecArgs[ii].first+"' for " + where, AT);
			has_step_date=true;
		} else if (vecArgs[ii].first=="VALUE_BEFORE") {
			IOUtils::parseArg(vecArgs[ii], where, value_before);
			has_value_before=true;
		} else if (vecArgs[ii].first=="VALUE_AFTER") {
			IOUtils::parseArg(vecArgs[ii], where, value_after);
			has_value_after=true;
		}
	}
	
	if (!has_value_before) throw InvalidArgumentException("Please provide the VALUE_BEFORE argument for the "+where, AT);
	if (!has_value_after) throw InvalidArgumentException("Please provide the VALUE_AFTER argument for the "+where, AT);
	if (!has_step_date) throw InvalidArgumentException("Please provide the STEP_DATE argument for the "+where, AT);
}

double STEP_Synth::generate(const Date& dt) const
{
	if (dt<dt_step)
		return value_before;
	else
		return value_after;
}

RECT_Synth::RECT_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ) 
          : step_start(), step_stop(), value(IOUtils::nodata), value_step(IOUtils::nodata)
{
	const std::string where( "SYNTH RECTANGLE, " + station + "::" + parname );
	bool has_step_start = false, has_step_stop = false;
	bool has_value = false, has_value_step = false;
	
	//parse the arguments
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="STEP_START") {
			if (!IOUtils::convertString(step_start, vecArgs[ii].second, TZ))
				throw InvalidArgumentException("Can not parse argument '"+vecArgs[ii].first+"' for " + where, AT);
			has_step_start=true;
		} else if (vecArgs[ii].first=="STEP_STOP") {
			if (!IOUtils::convertString(step_stop, vecArgs[ii].second, TZ))
				throw InvalidArgumentException("Can not parse argument '"+vecArgs[ii].first+"' for " + where, AT);
			has_step_stop=true;
		} if (vecArgs[ii].first=="VALUE") {
			IOUtils::parseArg(vecArgs[ii], where, value);
			has_value=true;
		} else if (vecArgs[ii].first=="VALUE_STEP") {
			IOUtils::parseArg(vecArgs[ii], where, value_step);
			has_value_step=true;
		}
	}
	
	if (!has_value) throw InvalidArgumentException("Please provide the VALUE argument for the "+where, AT);
	if (!has_value_step) throw InvalidArgumentException("Please provide the VALUE_STEP argument for the "+where, AT);
	if (!has_step_start) throw InvalidArgumentException("Please provide the STEP_START argument for the "+where, AT);
	if (!has_step_stop) throw InvalidArgumentException("Please provide the STEP_STOP argument for the "+where, AT);
}

double RECT_Synth::generate(const Date& dt) const
{
	if (dt<step_start || dt>step_stop)
		return value;
	else
		return value_step;
}

STDPRESS_Synth::STDPRESS_Synth(const StationData& sd) 
          : altitude( sd.getAltitude() ) {}

double STDPRESS_Synth::generate(const Date& /*dt*/) const
{
	if (altitude==IOUtils::nodata) return IOUtils::nodata;
	return Atmosphere::stdAirPressure(altitude);
}

///////////////////////////////////////////////////////
// Object factory
Synthesizer* SynthFactory::getSynth(std::string type, const std::string& station, const StationData& sd, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ)
{
	IOUtils::toUpper( type );
	
	if (type == "CST"){
		return new CST_Synth(station, parname, vecArgs);
	} else if (type == "STEP"){
		return new STEP_Synth(station, parname, vecArgs, TZ);
	} else if (type == "RECTANGLE"){
		return new RECT_Synth(station, parname, vecArgs, TZ);
	} else if (type == "STDPRESS"){
		return new STDPRESS_Synth(sd);
	} else {
		throw IOException("The Synthesizer '"+type+"' does not exist!" , AT);
	}
}

} //namespace
