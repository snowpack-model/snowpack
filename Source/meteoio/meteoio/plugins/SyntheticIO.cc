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

using namespace std;

namespace mio {
/**
 * @page synthio SynthIO
 * This plugin is quite special since it does not read any data but generates synthetic data. It is designed for
 * numerical experiments where controlled conditions are applied to the numerical setup.
 * 
 * @section synthio_keywords Keywords
 * This plugin uses the following keywords, all in the [Input] section:
 * - controlling the timestamps generation:
 *     - TIME_ZONE: the time zone of any dates that are provided; 
 *     - SYNTH_START: when to start generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_END: when to stop generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_SAMPLING: sampling rate in seconds, starting with SYNTH_START or the requested date if no SYNTH_START is provided (mandatory);
 * - providing the stations' metadata:
 *     - COORDSYS: coordinate system (see Coords);
 *     - COORDPARAM: extra coordinates parameters (see Coords);
 *     - STATION#: coordinates of the station (mandatory, see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax);
 *     - ID#: the (short) station id to use (optional but recommended, default is "ID_#")
 *     - NAME#: a descriptive station name to use (optional, default is "STATION_#");
 * 
 * @section synthio_examples
 * Example of use:
 * @code
 * [INPUT]
 * COORDSYS = CH1903
 * TIME_ZONE = 1.00
 * 
 * METEO = SYNTH
 * SYNTH_SAMPLING = 1800
 * SYNTH_START = 2022-09-02T12:40
 * STATION = latlon (46.8, 9.81, 1500)
 * @endcode
 */

SynthIO::SynthIO(const std::string& configfile) : cfg(configfile), vecStations(), dt_start(), dt_end(), dt_step(0.)
{
	init();
}

SynthIO::SynthIO(const Config& cfgreader) : cfg(cfgreader), vecStations(), dt_start(), dt_end(), dt_step(0.)
{
	init();
}

void SynthIO::init()
{
	//read start / end time as well as sampling rate
	const double TZ = cfg.get("TIME_ZONE", "INPUT");
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
	std::vector<std::string> coords_specs, vecIDs, vecNames;
	cfg.getValues("STATION", "INPUT", coords_specs);
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
		const Coords loc(coordin, coordinparam, coords_specs[ii]);
		const std::string id = (has_ids)? vecIDs[ii] : "ID_"+IOUtils::toString( ii+1 );
		const std::string name = (has_names)? vecNames[ii] : "STATION_"+IOUtils::toString( ii+1 );
		const StationData sd(loc, id, name);
		
		vecStations.push_back( sd );
	}
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
	
	/*std::cout << "dt_start.isUndef()=" << dt_start.isUndef() << "\n";
	std::cout << "dt_start=" << dt_start.toString(Date::ISO);
	std::cout << " dateStart=" << dateStart.toString(Date::ISO) << "\n";
	std::cout << "Starting from " << true_start.toString(Date::ISO) << " with dt_step_days=" << dt_step_days << "\n";*/
	
	for (auto &station : vecStations) {
		std::vector<MeteoData> vecM;
		for (Date dt=true_start; dt<=true_end; dt+=dt_step_days) {
			MeteoData md( dt, station );
			md(MeteoData::TA) = 280.0;
			vecM.push_back( md );
		}
		vecMeteo.push_back( vecM );
	}
	//std::cout << "vecMeteo[0].size()=" << vecMeteo[0].size() << "\n";
}

} //namespace
