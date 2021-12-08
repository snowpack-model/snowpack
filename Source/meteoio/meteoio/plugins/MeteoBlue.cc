// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 SLF                                                                                                                                */
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
#include <meteoio/plugins/MeteoBlue.h>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteoLaws/Meteoconst.h>

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstring>

#include <curl/curl.h>
#include <meteoio/thirdParty/picojson.h>

using namespace std;

namespace mio {
/**
 * @page meteoblue MeteoBlue
 * @section meteoblue_format Format
 * This plugin reads meteorological data from <a href="https://www.meteoblue.com">Meteoblue</a>'s 
 * <a href="https://docs.meteoblue.com/en/weather-apis/packages-api/introduction">packages API</a> as a RESTful web service. 
 * To compile the plugin you need to have the <a href="http://curl.haxx.se/">CURL library</a> with its headers present.
 * 
 * The meteorological parameters are defined by Meteoblue in <a href="https://content.meteoblue.com/fr/specifications/variables-meteo/meteoblue_weather_variables-documentation_EN_v06.pdf">this document</a>.
 *
 * @section meteoblue_keywords Keywords
 * This plugin uses the following keywords:
 * - METEOBLUE_URL: the API endpoint to connect to (default: http://my.meteoblue.com/);
 * - METEOBLUE_APIKEY: the key purchased from Meteoblue that gives access to its API (mandatory);
 * - METEOBLUE_PACKAGES: a space delimited list of <a href="https://docs.meteoblue.com/en/weather-apis/packages-api/forecast-data">packages</a> to read data from. Please note that you must have purchased an API key that gives access to all the packages that you list here! (mandatory);
 * - METEOBLUE_TIMEOUT: timeout (in seconds) for the connection to the server (default: 60s);
 * - METEOBLUE_DEBUG: print more information in order to better understand when something does not work as expected (default: false);
 * - STATION\#: provide the lat, lon and altitude or easting, northing and altitude for a station to get the data from (see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax). If your don't have the elevation of your station of interest, you can use <a href="https://latlongdata.com/elevation/">latlongdata.com</a> (mandatory);
 *    - STATION\#_ID: provide an ID for the declared station (default: "STAT\#");
 *    - STATION\#_NAME: provide a name for the declared station (default: "STATION\#");
 * 
 * @note You can only get (at most) data from NOW-3days up to NOW+8days with the non-historical packages and the packages API as used here.
 * Depending on the forecast length, various sources of data are used in the numerical models (see this 
 * <a href="https://docs.meteoblue.com/en/weather-apis/packages-api/introduction#forecast-length">introduction</a>). In this plugin,
 * the past 3 days are always retrieved (as they include satellite observations and ground measurements) up to the next 7-8 days. 
 * This means that the plugin will ignore any time period specifications when called but MeteoIO will later truncate the extracted
 * dataset to the requested period. 
 * 
 * @code
 * METEO	= METEOBLUE
 * METEOBLUE_APIKEY = xxxxxxxxxxx
 * METEOBLUE_PACKAGES = basic-1h clouds-1h solar-1h
 * STATION1 = latlon (47.56, 7.57, 262)
 * STATION1_ID = BSL1
 * STATION1_NAME = Basel
 * @endcode
 *
 * @section meteoblue_installing Installation
 * On Linux, simply install libcurl for your distribution and cmake should be able to find it automatically. On Windows, please install libcurl 
 * (you can find <a href="https://wiki.openssl.org/index.php/Binaries">pre-compiled binaries</a>). Then you will most probably have to 
 * manually provide the file and path to libcurl.dll.a as well as to the "include" directory of libcurl in cmake in order to be able 
 * to compile the meteoblue plugin.
 * 
 * @section meteoblue_dependencies Picojson
 * This plugin relies on <A HREF="https://github.com/kazuho/picojson/">picojson</A> for reading and parsing
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

const int MeteoBlue::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string MeteoBlue::dflt_endpoint = "https://my.meteoblue.com/";
std::map< std::string, MeteoBlue::meteoParam > MeteoBlue::params_map;
const bool MeteoBlue::__init = MeteoBlue::initStaticData();

bool MeteoBlue::initStaticData()
{
	params_map[ "precipitation" ]       = meteoParam(MeteoGrids::PSUM);
	params_map[ "snowfraction" ]        = meteoParam(MeteoGrids::PSUM_PH, -1., 1.); //the opposite of MeteoGrids
	params_map[ "temperature" ]         = meteoParam(MeteoGrids::TA, 1., 273.15);
	params_map[ "relativehumidity" ]    = meteoParam(MeteoGrids::RH, 0.01, 0.);
	params_map[ "windspeed" ]           = meteoParam(MeteoGrids::VW);
	params_map[ "winddirection" ]       = meteoParam(MeteoGrids::DW);
	params_map[ "sealevelpressure" ]    = meteoParam(MeteoGrids::P_SEA, 100., 0.);
	params_map[ "skintemperature" ]     = meteoParam(MeteoGrids::TSS, 1., 273.15);
	params_map[ "ghi_total" ]           = meteoParam(MeteoGrids::ISWR);
	params_map[ "ghi_instant" ]         = meteoParam(MeteoGrids::ISWR);
	params_map[ "dif_total" ]           = meteoParam(MeteoGrids::ISWR_DIFF);
	params_map[ "dif_instant" ]         = meteoParam(MeteoGrids::ISWR_DIFF);
	params_map[ "surfaceairpressure" ]  = meteoParam(MeteoGrids::P, 100., 0.);
	params_map[ "gust" ]                = meteoParam(MeteoGrids::VW_MAX);
	params_map[ "totalcloudcover" ]     = meteoParam(MeteoGrids::CLD, 8./100., 0); //TODO we have to convert this to a transmissivity!
	
	return true;
}

MeteoBlue::MeteoBlue(const std::string& configfile)
      : cfg(configfile), vecMeta(),
        coordin(), coordinparam(),
        endpoint(dflt_endpoint), apikey(), packages(),
        http_timeout(http_timeout_dflt), debug(false)
{
	init();
}

MeteoBlue::MeteoBlue(const Config& cfgreader)
      : cfg(cfgreader), vecMeta(),
        coordin(), coordinparam(),
        endpoint(dflt_endpoint), apikey(), packages(),
        http_timeout(http_timeout_dflt), debug(false)
{
	init();
}

void MeteoBlue::init()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);

	cfg.getValue("METEOBLUE_DEBUG", "INPUT", debug, IOUtils::nothrow);
	cfg.getValue("METEOBLUE_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("METEOBLUE_APIKEY", "Input", apikey);
	cfg.getValue("METEOBLUE_URL", "Input", endpoint, IOUtils::nothrow);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	
	//building the packages arguments
	std::vector<std::string> vecStr;
	cfg.getValue("METEOBLUE_PACKAGES", "Input", vecStr);
	packages = vecStr[0];
	for (size_t ii=1; ii<vecStr.size(); ii++) packages.append( "_"+vecStr[ii] );
	
	//reading the stations' coordinates to retrieve
	const std::vector< std::pair<std::string, std::string> > vecStationSpecs( cfg.getValues("STATION", "Input") );
	for (size_t ii=0; ii<vecStationSpecs.size(); ii++) {
		if (vecStationSpecs[ii].first.find('_') != std::string::npos) continue; //so we skip the other station#_xxx parameters
		
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		const Coords curr_point(coordin, coordinparam, vecStationSpecs[ii].second);
		if (curr_point.isNodata()) continue;
		
		const std::string id_num( vecStationSpecs[ii].first.substr(std::string("STATION").length()) );
		const bool has_id = cfg.keyExists("STATION"+id_num+"_ID", "Input");
		const std::string stat_id = (has_id)? cfg.get("STATION"+id_num+"_ID", "Input"): "STAT"+id_num;
		const bool has_name = cfg.keyExists("STATION"+id_num+"_NAME", "Input");
		const std::string stat_name = (has_name)? cfg.get("STATION"+id_num+"_NAME", "Input"): "STATION"+id_num;
		
		vecMeta.push_back( StationData(curr_point, stat_id, stat_name) );
	}
	if (vecMeta.empty())
		throw InvalidArgumentException("Please provide stations' coordinates for the METEOBLUE plugin!", AT);
	
	curl_global_init(CURL_GLOBAL_ALL);
}

void MeteoBlue::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation = vecMeta;
}

//with the packages API that we use, we can not choose the start and end dates
void MeteoBlue::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	
	vecMeteo.resize(vecMeta.size());
	for(size_t ii=0; ii<vecMeta.size(); ii++) {
		readData(vecMeta[ii], vecMeteo[ii]);
	}
}

picojson::value MeteoBlue::goToJSONPath(const std::string& path, const picojson::value& v)
{
	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	if (v.is<picojson::object>()) {
		const picojson::value::object& obj = v.get<picojson::object>();
		for (std::map<std::string,picojson::value>::const_iterator it = obj.begin(); it != obj.end(); ++it) {
			if (it->first==local_path) {
				if (!remaining_path.empty())
					goToJSONPath(remaining_path, it->second);
				else
					return it->second;
			}
		}
	}

	return picojson::value();
}

void MeteoBlue::readTime(const picojson::value &v, const StationData& sd, std::vector<MeteoData> &vecMeteo) const
{
	const picojson::value& time( goToJSONPath("$.time", v) );
	if (!time.is<picojson::array>()) {
		throw InvalidFormatException("Could not find a time array in the data section '"+time.to_str()+"'", AT);
	}
	const picojson::array vecRaw( time.get<picojson::array>() );
	
	vecMeteo.reserve( vecRaw.size() );
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		if (!vecRaw[ii].is<std::string>()) {
			std::ostringstream ss; ss << "Could not parse date '" << vecRaw[ii].to_str() << "'";
			throw InvalidFormatException(ss.str(), AT);
		}
		Date datum;
		IOUtils::convertString(datum, vecRaw[ii].get<std::string>(), 0.); //hard-coding GMT as this is what is delivered in ISO8601
		vecMeteo.push_back( MeteoData(datum, sd) );
	}
}

void MeteoBlue::readParameter(const picojson::value &v, const std::string& paramname, std::vector<MeteoData> &vecMeteo) const
{
	const std::map< std::string, meteoParam >::const_iterator it = params_map.find( paramname );
	if (it==params_map.end()) {
		if (debug) std::cout << "in MeteoBlue, skipping unknown parameter '" << paramname << "'\n";
		return;
	}
	
	//retrieve the raw data from the JSON
	const std::string mio_parname( it->second.getParameterName() );
	const picojson::value& param( goToJSONPath("$."+paramname, v) );
	if (!param.is<picojson::array>()) {
		throw InvalidFormatException("Could not parse the data section '"+param.to_str()+"'", AT);
	}
	const picojson::array vecRaw( param.get<picojson::array>() );
	if (vecMeteo.size()!=vecRaw.size()) 
		throw InvalidFormatException("Trying to insert "+IOUtils::toString(vecRaw.size())+" values into vecMeteo of size "+IOUtils::toString(vecMeteo.size()), AT);
	
	//convert and insert the raw data into vecMeteo
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		//create new parameters in any case, even if nodata so all MeteoData objects have the same parameters
		const size_t parindex = vecMeteo[ii].addParameter( mio_parname ); //already existing params just return their index
		
		if (vecRaw[ii].is<picojson::null>()) continue; //keep nodata
		if (!vecRaw[ii].is<double>())
			throw InvalidFormatException("Could not parse '"+vecRaw[ii].to_str()+"' as double", AT);
		
		vecMeteo[ii]( parindex ) = it->second.convertValue( vecRaw[ii].get<double>() );
	}
}

void MeteoBlue::readData(const StationData& sd, std::vector<MeteoData> &vecMeteo)
{
	//build the query URL
	std::ostringstream url;
	url << endpoint << "packages/" << packages << "?lat=" << sd.position.getLat() << "&lon=" << sd.position.getLon();
	const double altitude( sd.position.getAltitude() );
	if (altitude!=IOUtils::nodata) url << "&asl=" << altitude;
	url << "&timeformat=iso8601" << "&history_days=3" << "&apikey=" << apikey;
	
	if (debug) std::cout << "MeteoBlue, using the following URL: " << url.str() << "\n";
	
	std::stringstream ss;
	if (curl_read(url.str(), ss)) { //retrieve the page from the formed URL
		if (ss.str().empty()) throw UnknownValueException("No data returned for query '"+url.str()+"'", AT);
		
		//JSON error checking
		picojson::value v;
		const std::string err( picojson::parse(v, ss.str()) );
		if (!err.empty()) throw IOException("Error while parsing JSON: "+ss.str(), AT);
		if (v.contains("error")) {
			const std::string error_msg( v.get<picojson::object>()["error_message"].to_str() );
			throw InvalidArgumentException(error_msg, AT);
		}
		
		//loop over the datasets
		const picojson::value::object& datasets = v.get<picojson::object>();
		for (picojson::value::object::const_iterator it1 = datasets.begin(); it1 != datasets.end(); ++it1) {
			const std::string datasetname( it1->first );
			if (datasetname=="metadata" || datasetname=="units") continue;
			
			const picojson::value dataset( goToJSONPath("$."+datasetname, v) );
			//if multiple datasets create duplicated timestamps, a dataEditing will have to be setup
			readTime(dataset, sd, vecMeteo);
			
			//loop over the parameters (the time has been handled separately)
			const picojson::value::object& params = dataset.get<picojson::object>();
			for (picojson::value::object::const_iterator it2 = params.begin(); it2 != params.end(); ++it2) {
				const std::string paramname( it2->first );
				if (paramname!="time") readParameter(dataset, paramname, vecMeteo);
			}
		}
	}
}

size_t MeteoBlue::data_write(void* buf, const size_t size, const size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os( *static_cast<ostream*>(userp) );
		const std::streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool MeteoBlue::curl_read(const std::string& url, std::ostream& os) const
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	if (curl) {
		if (CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_ACCEPT_ENCODING, ""))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, http_timeout))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str())))
		{
			code = curl_easy_perform(curl);
		}
		curl_easy_cleanup(curl);
	}

	if (code!=CURLE_OK) {
		if (debug)
			std::cout << "****\nRequest: " << url << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	return (code==CURLE_OK);
}

} //namespace
