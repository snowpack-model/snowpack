// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 SLF                                                                                                                                */
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
#include <cstddef>
#include <meteoio/plugins/DBO.h>
#include <meteoio/plugins/JsonWrapper.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteoLaws/Meteoconst.h>

#include <regex>

namespace mio {
/**
 * @page dbo DBO
 * @section dbo_format Format
 * This plugin reads meteorological data from DBO
 * via the RESTful web service. To compile the plugin you need to have the <a href="http://curl.haxx.se/">CURL library</a> with its headers present.
 * \warning This plugin is for SLF's internal use only! You can also access the data through a web interface at <a href="https://measurements.slf.ch/">measurements.slf.ch</a>.
 *
 * You can have a look at the stations that are available through this web service on <a href="https://map.geo.admin.ch/?zoom=5&lang=en&topic=ech&bgLayer=ch.swisstopo.pixelkarte-farbe&layers=KML%7C%7Chttps:%2F%2Fstationdocu.slf.ch%2Fkml%2Fnetwork-map.kml&E=2782095.39&N=1179586.56">this map</a>.
 * Please keep in mind that some stations might be overlaid on top of each other and will require you to zoom in quite a lot in order to differentiate them!
 *
 * @section dbo_keywords Keywords
 * This plugin uses the following keywords:
 * - DBO_URL: The URL of the RESTful web service (default: https://pgdata.int.slf.ch)
 * - DBO_PROXY: The URL of a <A HREF="https://linuxize.com/post/how-to-setup-ssh-socks-tunnel-for-private-browsing/">SOCKS5 proxy</A> for the connection to go through (optional, specified as {host}:{port} such as *localhost:8080*, see <A HREF="https://stackoverflow.com/questions/51579063/curl-https-via-an-ssh-proxy">this</A> for more)
 * - STATION#: station code for the given station, prefixed by the network it belongs to (for example: IMIS::SLF2, by default the network is assumed to be IMIS). Valid networks are IMIS, SMN, BEOB, IMIS_RELAIS, VIRTUAL).
 * - DBO_TIMEOUT: timeout (in seconds) for the connection to the server (default: 60s)
 * - DBO_COVERAGE_RESTRICT: only request data from within the provided DateRange (see DateRange::setRange) (this is useful when merging data from several sources, for example to only get the new data from DBO and otherwise use the old data from files)
 * - DBO_DEBUG: print the full requests/answers from the server when something does not work as expected (default: false)
 *
 * @code
 * METEO	= DBO
 * STATION1	= WFJ2
 * STATION2	= SMN::*WFJ1
 * @endcode
 *
 * @section dbo_dependencies Picojson
 * This plugin relies on an embedded version of <A HREF="https://github.com/kazuho/picojson/">picojson</A> for reading and parsing
 * <A HREF="https://en.wikipedia.org/wiki/JSON">JSON</A> data. Picojson is released under a
 * <A HREF="https://opensource.org/licenses/BSD-2-Clause">2-Clause BSD License</A>. Please find here below
 * the full license agreement for picojson:
 *
 * @code
 * Copyright 2009-2010 Cybozu Labs, Inc.
 * Copyright 2011-2014 Kazuho Oku
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * @endcode
 */

static const double dbo_tz = 0.; //assuming GMT

//the pair <date, value> is contained in an array of size 2
inline bool parseTsPoint(const picojson::value& v, Date& datum, double& value)
{
	if (!v.is<picojson::array>()) return false;

	const picojson::array& array = v.get<picojson::array>();
	if (array.size()!=2) return false;

	//reading the date
	if (array[0].is<std::string>())
		IOUtils::convertString(datum, array[0].get<std::string>(), dbo_tz);
	else
		return false;

	//reading the value
	if (array[1].is<double>())
		value = array[1].get<double>();
	else {
		if (!array[1].is<picojson::null>()) return false;
		value = IOUtils::nodata;
		return true;
	}

	return true;
}

/*************************************************************************************************/
//example metadata query: https://pgdata.int.slf.ch/data/stations/IMIS/WFJ2
const int DBO::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string DBO::endpoint_default = "https://pgdata.int.slf.ch";
const std::string DBO::metadata_api = "/data/stations/";
const std::string DBO::data_api = "/data/timeseries/";

DBO::DBO(const std::string& configfile)
      : vecStationName(), vecMeta(), vecTsMeta(), coordin(), coordinparam(),
        endpoint(endpoint_default), coverageRestrict(), json( new JsonWrapper() ), 
        dbo_debug(false)
{
	const Config cfg( configfile );
	initDBOConnection(cfg);
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

DBO::DBO(const Config& cfgreader)
      : vecStationName(), vecMeta(), vecTsMeta(), coordin(), coordinparam(),
        endpoint(endpoint_default), coverageRestrict(), json( new JsonWrapper() ),
        dbo_debug(false)
{
	initDBOConnection(cfgreader);
	IOUtils::getProjectionParameters(cfgreader, coordin, coordinparam);
	cfgreader.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

DBO::DBO(const DBO& c)
    : vecStationName(c.vecStationName), vecMeta(c.vecMeta), vecTsMeta(c.vecTsMeta), coordin(c.coordin), coordinparam(c.coordinparam),
    endpoint(c.endpoint), coverageRestrict(c.coverageRestrict), json( new JsonWrapper() ), dbo_debug(c.dbo_debug) {}

DBO::~DBO()
{
	delete json;
}

DBO& DBO::operator=(const mio::DBO& c)
{
	if (this != &c) {
		vecStationName = c.vecStationName;
		vecMeta = c.vecMeta;
		vecTsMeta = c.vecTsMeta;
		coordin = c.coordin;
		coordinparam = c.coordinparam;
		endpoint = c.endpoint;
		coverageRestrict = c.coverageRestrict;
		json = new JsonWrapper();
		dbo_debug = c.dbo_debug;
	}
	
	return *this;
}

void DBO::initDBOConnection(const Config& cfg)
{
	cfg.getValue("DBO_URL", "Input", endpoint, IOUtils::nothrow);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	std::cerr << "[I] Using DBO URL: " << endpoint << std::endl;
	
	int http_timeout = http_timeout_dflt;
	cfg.getValue("DBO_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("DBO_DEBUG", "INPUT", dbo_debug, IOUtils::nothrow);
	const std::string proxy = cfg.get("DBO_PROXY", "INPUT", "");
	json->setConnectionParams(proxy, http_timeout, dbo_debug);
	
	std::string dateRangeHint;
	cfg.getValue("DBO_COVERAGE_RESTRICT", "INPUT", dateRangeHint, IOUtils::nothrow);
	if (!dateRangeHint.empty()) coverageRestrict.setRange( dateRangeHint, dbo_tz );
}

void DBO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	if (vecMeta.empty()) fillStationMeta();
	vecStation = vecMeta;
}

void DBO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	if (vecMeta.empty()) fillStationMeta();
	if (!coverageRestrict.isUndef() && (coverageRestrict.start>dateEnd || coverageRestrict.end<dateStart)) return;
	
	const Date trueStart = (coverageRestrict.isUndef())? dateStart : std::max(dateStart, coverageRestrict.start);
	const Date trueEnd = (coverageRestrict.isUndef())? dateEnd : std::min(dateEnd, coverageRestrict.end);
	
	vecMeteo.resize(vecMeta.size());
	for(size_t ii=0; ii<vecMeta.size(); ii++)
		readData(trueStart, trueEnd, vecMeteo[ii], ii);
}

/**
* @brief Read and cache the stations' metadata
*/
void DBO::fillStationMeta()
{
	static const std::regex stat_id_regex("([^:]+)::([^:]+)", std::regex::optimize);
	std::smatch stat_id_matches;

	vecMeta.clear();
	vecMeta.resize( vecStationName.size() );
	vecTsMeta.resize( vecStationName.size() );

	for(size_t ii=0; ii<vecStationName.size(); ii++) {
		std::string station_id( IOUtils::strToUpper(vecStationName[ii]) );
		std::string network = "IMIS";

		if (std::regex_match(station_id, stat_id_matches, stat_id_regex)) {
			network = stat_id_matches.str(1);
			station_id = stat_id_matches.str(2);
		}

		const std::string request( endpoint + metadata_api + network + "/" + station_id );
		json->readAndParse(request, station_id);
		const std::string error_msg( json->getString("$.error") );
		if (!error_msg.empty()) 
			throw UnknownValueException("Error with station '"+network+"::"+station_id+"': "+error_msg, AT);
		
		Coords position(coordin, coordinparam);
		position.setLatLon(json->getDouble("$.location.lat"), json->getDouble("$.location.lon"), json->getDouble("$.location.elevation"));
		const StationData sd(position, json->getString("$.code"), json->getString("$.label"));
		vecMeta[ii] = sd;
		
		//from the json, parse and store the time series belonging to this station
		vecTsMeta[ii] = getTsProperties();

		if (dbo_debug) {
			std::cout << "<Station " << station_id << ">\n";
			for (const auto& ts : vecTsMeta[ii]) std::cout << ts.toString() << "\n";
			std::cout << "</Station " << station_id << ">\n";
		}
	}
}

//get the properties of all timeseries belonging to the current station
std::vector<DBO::tsMeta> DBO::getTsProperties() const
{
	std::vector<tsMeta> tsVec;
	const std::vector<picojson::value> results( json->JSONQuery("$.timeseries") );

	for (size_t ii=0; ii<results.size(); ii++) {
		if (!results[ii].is<picojson::array>()) continue;

		for (const auto& ts_obj : results[ii].get<picojson::array>()) { //loop over all provided timeseries
			if (ts_obj.is<picojson::null>()) continue;
			
			std::string code, device_code, agg_type;
			double id = -1., seqNr = 1;
			double interval=0, ts_offset=0;
			Date since, until;

			for (const auto& keyValue : ts_obj.get<picojson::object>()) { //loop over all key/values of a given timeseries
				const std::string key( keyValue.first );

				if (key=="sequenceNumber" && keyValue.second.is<double>()) seqNr = keyValue.second.get<double>();
				if (key=="id" && keyValue.second.is<double>()) id = keyValue.second.get<double>();
				if (key=="measurandCode" && keyValue.second.is<std::string>()) code = keyValue.second.get<std::string>();
				if (key=="deviceCode" && keyValue.second.is<std::string>()) device_code = keyValue.second.get<std::string>();
				if (key=="since" && keyValue.second.is<std::string>()) IOUtils::convertString(since, keyValue.second.get<std::string>(), 0.);
				if (key=="until" && keyValue.second.is<std::string>()) IOUtils::convertString(until, keyValue.second.get<std::string>(), 0.);
				if (key=="aggregationType" && keyValue.second.is<std::string>()) agg_type = keyValue.second.get<std::string>();
				if (key=="measureIntervalInMinutes" && keyValue.second.is<double>()) interval = keyValue.second.get<double>(); //TODO check that it can be cast
				if (key=="measureIntervalOffsetInMinutes" && keyValue.second.is<double>()) ts_offset = keyValue.second.get<double>();
			}

			//reject some timeseries
			if (seqNr!=1) continue; //HACK per WIS, only consider seq number 1
			if (device_code=="BATTERY" || device_code=="LOGGER") continue;
			if (device_code=="MODEL_SNOWPACK") continue; //TODO add an option to use Snowpack computed parameters
			if (agg_type=="SD") continue; //we don't care about standard deviation anyway
			if (id==-1.) continue; //no id was provided

			const std::string param_dbo( IOUtils::strToUpper( code.substr(0, code.find('_')) ) );
			const std::string parname( getParameter(param_dbo, agg_type) );
			if (parname.empty()) continue;
			
			tsMeta tmp(param_dbo, since, until, agg_type, static_cast<size_t>(id), static_cast<unsigned int>(interval*60.), static_cast<unsigned int>(ts_offset*60.), static_cast<unsigned int>(seqNr));
			tmp.parname = parname;
			setUnitsConversion(tmp);
			tsVec.push_back( tmp );
		}
	}

	return tsVec;
}

/**
* @brief Identify the relevant MeteoData::Parameters from DBO provided information
* @param[in] param_str DBO string representation of the meteo parameter
* @param[in] agg_type DBO aggregation type
* @return standardized parameter name or empty string
*/
std::string DBO::getParameter(const std::string& param_str, const std::string& agg_type) const
{
	//not mapped yet: wet bulb temperature TPSY

	if (param_str=="P") return "P";
	else if (param_str=="PASL") return "P_SEA";
	else if (param_str=="TA") return "TA";
	else if (param_str=="RH") return "RH";
	else if (param_str=="TDP") return "TD";
	else if (param_str=="TS0") return "TSG";
	else if (param_str=="TSS") return "TSS";
	else if (param_str=="HS") return "HS";
	else if (param_str=="HN") return "HN";
	else if (param_str=="HNW") return "HN_SWE";
	else if (param_str=="VW" && agg_type=="MAX") return "VW_MAX";
	else if (param_str=="VW") return "VW";
	else if (param_str=="DW") return "DW";
	else if (param_str=="RSWR") return "RSWR";
	else if (param_str=="ISWR") return "ISWR";
	else if (param_str=="ILWR") return "ILWR";
	else if (param_str=="RR") return "PSUM";
	else if (param_str=="TS25") return "TS@25";
	else if (param_str=="TS50") return "TS@50";
	else if (param_str=="TS100") return "TS@100";
	else if (param_str=="TG10") return "TSOIL@10";
	else if (param_str=="TG30") return "TSOIL@30";
	else if (param_str=="TG50") return "TSOIL@50";

	if (dbo_debug) return param_str;
	return "";
}

/**
* @brief Provide the way to convert the DBO units into standardized units (SI).
* It is assume that we can first multiply by a factor, then add an offset.
* @param[in] ts DBO timeseries properties
*/
void DBO::setUnitsConversion(tsMeta& ts)
{
	//compute the conversion parameters (C to K, cm to m, % to [0-1], PINT to PSUM
	if (ts.parname=="TA" || ts.parname=="TD" || ts.parname=="TSG" || ts.parname=="TSS") {
		ts.units_offset = Cst::t_water_freezing_pt;
	} else if(ts.parname=="RH" || ts.parname=="HS" || ts.parname=="HN"|| ts.parname=="HN_SWE") {
		ts.units_factor = 0.01;
	} else if(ts.parname=="PSUM") {
		ts.units_factor = 3600. / ts.interval;
	} else if(ts.parname.find("TS@", 0)!=std::string::npos) {
		ts.units_offset = Cst::t_water_freezing_pt;
	} else if(ts.parname.find("TSOIL@", 0)!=std::string::npos) {
		ts.units_offset = Cst::t_water_freezing_pt;
	}
	
	return;
}

std::vector<DBO::tsData> DBO::getTimeSerie(const size_t& tsID, const double& factor, const double& offset) const
{
	const picojson::value ts( json->goToJSONPath("$.measurements") );
	if (!ts.is<picojson::array>())
		throw InvalidFormatException("Could not parse timeseries "+IOUtils::toString(tsID), AT);

	const picojson::array vecRaw( ts.get<picojson::array>() );
	if (vecRaw.empty()) return std::vector<tsData>();

	std::vector<tsData> vecData( vecRaw.size() );
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		double value;
		Date datum;
		if (!parseTsPoint(vecRaw[ii], datum, value))  {
			JsonWrapper::printJSON(vecRaw[ii], 0);
			throw InvalidFormatException("Error parsing element "+IOUtils::toString(ii)+" of timeseries "+IOUtils::toString(tsID), AT);
		}

		if (value!=IOUtils::nodata) value = value * factor + offset;
		vecData[ii] = tsData(datum, value);
	}

	return vecData;
}

//read all data for the given station
void DBO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	const std::string Start( dateStart.toString(Date::ISO_Z) );
	const std::string End( dateEnd.toString(Date::ISO_Z) );
	const StationData sd( vecMeta[stationindex] );

	//now get the data
	for (const tsMeta &ts : vecTsMeta[stationindex]) {
		
		//check if the current ts contains our period of interest
		const Date tsStart(ts.since), tsEnd(ts.until);
		if ((!tsStart.isUndef() && tsStart>dateEnd) || (!tsEnd.isUndef() && tsEnd<dateStart)) continue;

		const std::string tsID( IOUtils::toString(ts.id) );
		const std::string request( endpoint + data_api + tsID + "?from=" + Start + "&until=" + End );
		json->readAndParse(request, tsID);
		const std::string error_msg( json->getString("$.error") );
		if (!error_msg.empty()) throw UnknownValueException("Error while parsing JSON for '"+vecMeta[stationindex].getStationID()+"': "+error_msg, AT);
		
		const std::vector<tsData> vecData( getTimeSerie(ts.id, ts.units_factor, ts.units_offset) );
		if (vecData.empty()) {
			if (dbo_debug) 
				std::cout << vecMeta[stationindex].getStationID() << " has no data for " << ts.parname << " in the requested period\n";
			continue;
		}

		MeteoData md_pattern = (vecMeteo.empty())? MeteoData(Date(), sd) : vecMeteo.front(); //This assumes that the station is not moving!
		if (!md_pattern.param_exists(ts.parname)) {
			md_pattern.addParameter( ts.parname );
			for (size_t jj=0; jj<vecMeteo.size(); jj++) vecMeteo[jj].addParameter( ts.parname ); //TODO rewrite addParameter to create and attribute a value
		}
		
		const size_t parindex = md_pattern.getParameterIndex( ts.parname );
		mergeTimeSeries(md_pattern, parindex, vecData, vecMeteo);
	}
}

/**
* @brief Merge a newly read timeseries into vecMeteo
* @param[in] md_pattern pattern MeteoData to be used to insert new elements
* @param[in] param index of the current meteo parameter
* @param[in] vecData the raw (but parsed) data
* @param vecMeteo the vector that will receive the new values (as well as inserts if necessary)
*/
void DBO::mergeTimeSeries(const MeteoData& md_pattern, const size_t& param, const std::vector<DBO::tsData>& vecData, std::vector<MeteoData>& vecMeteo) const
{
	if (vecData.empty()) return;

	if (vecMeteo.empty()) { //easy case: the initial vector is empty
		vecMeteo.resize( vecData.size() );
		for (size_t ii=0; ii<vecData.size(); ii++) {
			MeteoData md( md_pattern );
			md.date = vecData[ii].date;
			md(param) = vecData[ii].val;
			vecMeteo[ii] = md;
		}
	} else {
		size_t vecM_start = 0; //the index in vecRaw that matches the original start of vecMeteo
		size_t vecM_end = 0; //the index in vecRaw that matches the original end of vecMeteo

		//filling data before vecMeteo
		if (vecData.front().date<vecMeteo.front().date) {
			const Date start_date( vecMeteo.front().date );
			vecM_start = vecData.size(); //if no overlap is found, take all vecData
			for(size_t ii=0; ii<vecData.size(); ii++) { //find the range of elements to add
				if (vecData[ii].date>=start_date) {
					vecM_start = ii;
					break;
				}
			}

			vecMeteo.insert(vecMeteo.begin(), vecM_start, md_pattern);
			for (size_t ii=0; ii<vecM_start; ii++) {
				vecMeteo[ii].date = vecData[ii].date;
				vecMeteo[ii](param) = vecData[ii].val;
			}
		}

		//general case: merge one timestamp at a time
		std::vector<MeteoData> tmp;
		tmp.reserve( vecMeteo.size() + (vecData.size() - vecM_start)); //"worst case" scenario: all elements will be added

		size_t idx2 = vecM_start; //all previous elements were handled before
		size_t last_vM = vecM_start; //last element from vecMeteo that will have to be invalidated
		for(size_t ii=vecM_start; ii<vecMeteo.size(); ii++) {
			const Date curr_date( vecMeteo[ii].date );
			while ((idx2<vecData.size()) && (curr_date>vecData[idx2].date)) {
				tmp.push_back( md_pattern );
				tmp.back().date = vecData[idx2].date;
				tmp.back()(param) = vecData[idx2].val;
				idx2++;
			}
			if (idx2==vecData.size())  break; //nothing left to merge

			if (curr_date==vecData[idx2].date) {
				vecMeteo[ii](param) = vecData[idx2].val;
				idx2++;
			}
			tmp.push_back( vecMeteo[ii] );
			last_vM = ii;
		}

		const size_t new_count = last_vM - vecM_start + 1;
		if (new_count<tmp.size())
			vecMeteo.insert( vecMeteo.begin() + static_cast<ptrdiff_t>(vecM_start), tmp.size()-new_count, tmp.front()); //so room for the extra params is allocated

		for(size_t ii=0; ii<tmp.size(); ii++)
			vecMeteo[vecM_start+ii] = tmp[ii];

		vecM_end = idx2;

		//filling data after vecMeteo
		if (vecMeteo.back().date<vecData.back().date) {
			if (vecM_end!=vecData.size()) {
				for (size_t ii=vecM_end; ii<vecData.size(); ii++) {
					vecMeteo.push_back( md_pattern );
					vecMeteo.back().date = vecData[ii].date;
					vecMeteo.back()(param) = vecData[ii].val;
				}
			}
		}
	}
}

} //namespace
