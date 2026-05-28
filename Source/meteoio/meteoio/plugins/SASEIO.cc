// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2014 Snow and Avalanche Study Establishment    SASE-CHANDIGARH       */
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
#include <meteoio/plugins/SASEIO.h>

#ifdef _WIN32
	#include <winsock.h>
#endif // _WIN32

#include <mysql.h>
#include <stdio.h>
#include <cstring>

#include <algorithm>

using namespace std;

namespace mio {
/**
* @page sase SASEIO
* @section SASE_format Format
* This is the plugin required to get meteorological data from the SASE MySQL database.
*
* @section SASE_units Units
* The units are assumed to be the following:
* - __temperatures__ in celsius
* - __relative humidity__ in %
* - __wind speed__ in m/s
* - __precipitations__ in mm/h
* - __radiation__ in W/m²
*
* @section SASE_keywords Keywords
* This plugin uses the following keywords:
* - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
* - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
* - SASE_HOST: MySQL Host Name (e.g. localhost or 191.168.145.20); [Input] section
* - SASE_DB: MySQL Database (e.g. snowpack); [Input] section
* - SASE_USER: MySQL User Name (e.g. root); [Input] section
* - SASE_PASS: MySQL password; [Input] section
* - TIME_ZONE: For [Input] and [Output] sections
* - STATION#: station code for the given number #; [Input] section
*/

const double SASEIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
// SQL queries use ? placeholders for parameterized queries
const string SASEIO::MySQLQueryStationMetaData = "SELECT StationName, Latitude, Longitude, Altitude FROM metadata WHERE StationID=? AND StationIDInt=?"; ///< Snow station meta data
const string SASEIO::MySQLQueryMeteoData = "SELECT TimeStamp, TA, RH, VW, HS, TSS, OSWR, ISWR, ILWR, HNW, P FROM station_meteoio_data_aws WHERE StationID=? AND StationName=? AND TimeStamp>=? AND TimeStamp<=? ORDER BY TimeStamp ASC"; ///< METEO Data query using MYSQL Timestamp

//const string SASEIO::MySQLQueryMeteoData = "SELECT timestamp, TA, RH, VW, HS, TSS, OSWR FROM data WHERE StationID=? AND timestamp>=? AND timestamp<=? ORDER BY timestamp ASC"; ///< METEO Data query Using SNPK Timestamp

SASEIO::SASEIO(const std::string& configfile)
        : cfg(configfile), vecStationIDs(), vecStationMetaData(),
          mysqlhost(), mysqldb(), mysqluser(), mysqlpass(),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          in_dflt_TZ(5.5), out_dflt_TZ(5.5)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	readConfig();
}

SASEIO::SASEIO(const Config& cfgreader)
        : cfg(cfgreader), vecStationIDs(), vecStationMetaData(),
          mysqlhost(), mysqldb(), mysqluser(), mysqlpass(),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          in_dflt_TZ(5.5), out_dflt_TZ(5.5)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	readConfig();
}

void SASEIO::readConfig()
{
	cfg.getValue("SASE_HOST", "Input", mysqlhost);
	cfg.getValue("SASE_DB", "Input", mysqldb);
	cfg.getValue("SASE_USER", "Input", mysqluser);
	cfg.getValue("SASE_PASS", "Input", mysqlpass);
	cfg.getValue("TIME_ZONE","Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE","Output", out_dflt_TZ, IOUtils::nothrow);
}

void SASEIO::readStationIDs(std::vector<std::string>& vecStationID) const
{
	vecStationID.clear();
	cfg.getValues("STATION", "INPUT", vecStationID);

	if (vecStationID.empty()) {
		cerr << "\tNo stations specified for SASEIO... is this what you want?\n";
	}
}

void SASEIO::readStationMetaData()
{
	vecStationMetaData.clear();
	std::vector<std::string> vecStationID;
	readStationIDs(vecStationID);

	for (size_t ii=0; ii<vecStationID.size(); ii++) {
		// Retrieve the station IDs - this only needs to be done once per instance
		std::string stat_abk, stao_nr;
		parseStationID(vecStationID[ii], stat_abk, stao_nr);

		// Retrieve the station meta data - this only needs to be done once per instance
		std::vector<std::string> stationMetaData;
		getStationMetaData(stat_abk, stao_nr, MySQLQueryStationMetaData, stationMetaData);

		std::string stao_name;
		double lat, longi, alt;
		if (!IOUtils::convertString(stao_name, stationMetaData.at(0)))
			throw ConversionFailedException("Invalid station name for station "+vecStationID[ii]+": "+stationMetaData.at(0), AT);
		if (!IOUtils::convertString(lat, stationMetaData.at(1), std::dec))
			throw ConversionFailedException("Invalid latitude for station "+vecStationID[ii]+": "+stationMetaData.at(1), AT);
		if (!IOUtils::convertString(longi, stationMetaData.at(2), std::dec))
			throw ConversionFailedException("Invalid longitude for station "+vecStationID[ii]+": "+stationMetaData.at(2), AT);
		if (!IOUtils::convertString(alt, stationMetaData.at(3), std::dec))
			throw ConversionFailedException("Invalid altitude for station "+vecStationID[ii]+": "+stationMetaData.at(3), AT);

		Coords location(coordin,coordinparam);
		location.setLatLon(lat,longi,alt);
		const std::string station_name = (!stao_name.empty())? vecStationID[ii] + ":" + stao_name : vecStationID[ii];

		vecStationMetaData.push_back(StationData(location, vecStationID[ii], station_name));
	}

}

void SASEIO::getStationMetaData(const std::string& stat_abk, const std::string& stao_nr,
                                 const std::string& sqlQuery, std::vector<std::string>& vecMetaData)
{
	vecMetaData.clear();

	MYSQL *conn = mysql_init(nullptr);
	if (!mysql_real_connect(conn, mysqlhost.c_str(), mysqluser.c_str(), mysqlpass.c_str(), mysqldb.c_str(), 0, nullptr, 0)) {
		throw IOException("Could not initiate connection to Mysql server "+mysqlhost, AT);
	}

	// Use prepared statement
	MYSQL_STMT *stmt = mysql_stmt_init(conn);
	if (!stmt) {
		mysql_close(conn);
		throw IOException("Failed to initialize MySQL statement", AT);
	}

	if (mysql_stmt_prepare(stmt, sqlQuery.c_str(), sqlQuery.length())) {
		mysql_stmt_close(stmt);
		mysql_close(conn);
		throw IOException("Failed to prepare SQL statement: " + std::string(mysql_stmt_error(stmt)), AT);
	}

	// Bind parameters
	MYSQL_BIND bind[2];
	memset(bind, 0, sizeof(bind));

	bind[0].buffer_type = MYSQL_TYPE_STRING;
	bind[0].buffer = (char *)stat_abk.c_str();
	bind[0].buffer_length = stat_abk.length();

	bind[1].buffer_type = MYSQL_TYPE_STRING;
	bind[1].buffer = (char *)stao_nr.c_str();
	bind[1].buffer_length = stao_nr.length();

	if (mysql_stmt_bind_param(stmt, bind)) {
		mysql_stmt_close(stmt);
		mysql_close(conn);
		throw IOException("Failed to bind parameters", AT);
	}

	if (mysql_stmt_execute(stmt)) {
		mysql_stmt_close(stmt);
		mysql_close(conn);
		throw IOException("Query execution failed: " + std::string(mysql_stmt_error(stmt)), AT);
	}

	MYSQL_RES *res = mysql_stmt_result_metadata(stmt);
	if (!res) {
		mysql_stmt_close(stmt);
		mysql_close(conn);
		throw IOException("No result set returned", AT);
	}

	mysql_stmt_store_result(stmt);

	const unsigned int column_no = mysql_num_fields(res);
	const unsigned int column_no = mysql_num_fields(res);
	std::string tmp_str;
	MYSQL_ROW row;
	while ( ( row = mysql_fetch_row(res) ) != nullptr ) {
		for (unsigned int ii=0; ii<column_no; ii++) {
			IOUtils::convertString(tmp_str, row[ii]);
			vecMetaData.push_back(tmp_str);
		}
	}
	mysql_free_result(res);
	mysql_stmt_close(stmt);
	mysql_close(conn);

	const size_t nr_metadata = vecMetaData.size();
	if (nr_metadata==0)
		throw NoDataException("Station " + stat_abk+stao_nr + " not found in the database", AT);
	if (nr_metadata<4)
		throw ConversionFailedException("Error while converting station meta data for station "+stat_abk+stao_nr, AT);
}

void SASEIO::parseStationID(const std::string& stationID, std::string& stat_abk, std::string& stao_nr)
{
	stat_abk = stationID.substr(0, stationID.length()-1); //The station name: e.g. KLO
	stao_nr = stationID.substr(stationID.length()-1, 1); //The station number: e.g. 2
	if (!std::isdigit(stao_nr[0])) {
		//the station is one of these non-imis stations that don't contain a number...
		stat_abk = stationID;
		stao_nr = "0";
	}
}

void SASEIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	readStationMetaData(); //reads all the station meta data into the vecStationMetaData (member vector)
	vecStation = vecStationMetaData; //vecStationMetaData is a global vector holding all meta data
}

void SASEIO::readMeteoData(const Date& dateStart , const Date& dateEnd,
                            std::vector<std::vector<MeteoData> >& vecMeteo)
{
	if (vecStationMetaData.empty()) readStationMetaData();

	vecMeteo.clear();
	vecMeteo.insert(vecMeteo.begin(), vecStationMetaData.size(), vector<MeteoData>());

	for (size_t ii=0; ii<vecStationMetaData.size(); ii++) { //loop through relevant stations
		readData(dateStart, dateEnd, vecMeteo, ii, vecStationMetaData);
	}
}

//read meteo data for one station
void SASEIO::readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
                       const size_t& stationindex, const std::vector<StationData>& vecMeta) const
{
	vecMeteo.at(stationindex).clear();

	Date dateS(dateStart), dateE(dateEnd);
	dateS.setTimeZone(in_dflt_TZ);
	dateE.setTimeZone(in_dflt_TZ);

	std::string stat_abk, stao_nr;
	std::vector<std::string> vecHTS1;
	std::vector< std::vector<std::string> > vecResult;
	parseStationID(vecMeta.at(stationindex).getStationID(), stat_abk, stao_nr);

	getStationData(stat_abk, stao_nr, dateS, dateE, vecHTS1, vecResult);
	MeteoData tmpmd;
	tmpmd.meta = vecMeta.at(stationindex);

	for (size_t ii=0; ii<vecResult.size(); ii++){
		parseDataSet(vecResult[ii], tmpmd);
		convertUnits(tmpmd);
		vecMeteo.at(stationindex).push_back( tmpmd ); //Now insert tmpmd
	}
}

void SASEIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	ta = IOUtils::C_TO_K(ta);

	double& tss = meteo(MeteoData::TSS);
	tss = IOUtils::C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
	rh /= 100.;

	/*double& hs = meteo(MeteoData::HS);
	if (hs != IOUtils::nodata)
	hs /= 100.0;

	double& rswr = meteo(MeteoData::RSWR);
	if (rswr != IOUtils::nodata)
	rswr /= 100.0;*/
}

void SASEIO::parseDataSet(const std::vector<std::string>& i_meteo, MeteoData& md) const
{
	const std::string statID( md.meta.getStationID() );

	if (!IOUtils::convertString(md.date, i_meteo.at(0), in_dflt_TZ))
		throw ConversionFailedException("Invalid timestamp for station "+statID+": "+i_meteo.at(0), AT);
	if (!IOUtils::convertString(md(MeteoData::TA), i_meteo.at(1)))
		throw ConversionFailedException("Invalid TA for station "+statID+": "+i_meteo.at(1), AT);
	if (!IOUtils::convertString(md(MeteoData::RH), i_meteo.at(2)))
		throw ConversionFailedException("Invalid RH for station "+statID+": "+i_meteo.at(2), AT);
	if (!IOUtils::convertString(md(MeteoData::VW), i_meteo.at(3)))
		throw ConversionFailedException("Invalid VW for station "+statID+": "+i_meteo.at(3), AT);
	if (!IOUtils::convertString(md(MeteoData::HS), i_meteo.at(4)))
		throw ConversionFailedException("Invalid HS for station "+statID+": "+i_meteo.at(4), AT);
	if (!IOUtils::convertString(md(MeteoData::TSS), i_meteo.at(5)))
		throw ConversionFailedException("Invalid TSS for station "+statID+": "+i_meteo.at(5), AT);
	if (!IOUtils::convertString(md(MeteoData::RSWR), i_meteo.at(6)))
		throw ConversionFailedException("Invalid RSWR for station "+statID+": "+i_meteo.at(6), AT);
	if (!IOUtils::convertString(md(MeteoData::ISWR), i_meteo.at(7)))
		throw ConversionFailedException("Invalid ISWR for station "+statID+": "+i_meteo.at(7), AT);
	if (!IOUtils::convertString(md(MeteoData::ILWR), i_meteo.at(8)))
		throw ConversionFailedException("Invalid ILWR for station "+statID+": "+i_meteo.at(8), AT);
	if (!IOUtils::convertString(md(MeteoData::PSUM), i_meteo.at(9)))
		throw ConversionFailedException("Invalid PSUM for station "+statID+": "+i_meteo.at(9), AT);
	if (!IOUtils::convertString(md(MeteoData::P), i_meteo.at(10)))
		throw ConversionFailedException("Invalid P for station "+statID+": "+i_meteo.at(10), AT);
}

bool SASEIO::getStationData(const std::string& stat_abk, const std::string& stao_nr,
                             const Date& dateS, const Date& dateE,
                             const std::vector<std::string>& vecHTS1,
                             std::vector< std::vector<std::string> >& vecMeteoData) const
{
	vecMeteoData.clear();
	bool fullStation = true;

	// Creating MySQL TimeStamps
	std::string sDate( dateS.toString(Date::ISO) );
	std::replace( sDate.begin(), sDate.end(), 'T', ' ');
	std::string eDate( dateE.toString(Date::ISO) );
	std::replace( eDate.begin(), eDate.end(), 'T', ' ');

	// Use prepared statement
	MYSQL *conn2 = mysql_init(nullptr);
	if (!mysql_real_connect(conn2, mysqlhost.c_str(), mysqluser.c_str(), mysqlpass.c_str(), mysqldb.c_str(), 0, nullptr, 0)) {
		throw IOException("Could not initiate connection to Mysql server "+mysqlhost, AT);
	}

	MYSQL_STMT *stmt = mysql_stmt_init(conn2);
	if (!stmt) {
		mysql_close(conn2);
		throw IOException("Failed to initialize MySQL statement", AT);
	}

	if (mysql_stmt_prepare(stmt, MySQLQueryMeteoData.c_str(), MySQLQueryMeteoData.length())) {
		mysql_stmt_close(stmt);
		mysql_close(conn2);
		throw IOException("Failed to prepare SQL statement: " + std::string(mysql_stmt_error(stmt)), AT);
	}

	// Bind parameters
	MYSQL_BIND bind[4];
	memset(bind, 0, sizeof(bind));

	bind[0].buffer_type = MYSQL_TYPE_STRING;
	bind[0].buffer = (char *)stao_nr.c_str();
	bind[0].buffer_length = stao_nr.length();

	bind[1].buffer_type = MYSQL_TYPE_STRING;
	bind[1].buffer = (char *)stat_abk.c_str();
	bind[1].buffer_length = stat_abk.length();

	bind[2].buffer_type = MYSQL_TYPE_STRING;
	bind[2].buffer = (char *)sDate.c_str();
	bind[2].buffer_length = sDate.length();

	bind[3].buffer_type = MYSQL_TYPE_STRING;
	bind[3].buffer = (char *)eDate.c_str();
	bind[3].buffer_length = eDate.length();

	if (mysql_stmt_bind_param(stmt, bind)) {
		mysql_stmt_close(stmt);
		mysql_close(conn2);
		throw IOException("Failed to bind parameters", AT);
	}

	if (mysql_stmt_execute(stmt)) {
		mysql_stmt_close(stmt);
		mysql_close(conn2);
		throw IOException("Query execution failed: " + std::string(mysql_stmt_error(stmt)), AT);
	}

	MYSQL_RES *res2 = mysql_stmt_result_metadata(stmt);
	if (!res2) {
		mysql_stmt_close(stmt);
		mysql_close(conn2);
		throw IOException("No result set returned", AT);
	}

	mysql_stmt_store_result(stmt);

	const unsigned int column_no2 = mysql_num_fields(res2);
	while ( ( row2 = mysql_fetch_row(res2) ) != nullptr ) {
		std::vector<std::string> vecData;
		for (unsigned int ii=0; ii<column_no2; ii++) {
			std::string row_02;
				if (!row2[ii]){
					IOUtils::convertString(row_02,"-999.0");// HARD CODED :(
				}else{
					IOUtils::convertString(row_02, row2[ii]);
				 }
			vecData.push_back(row_02);
		}
		if (fullStation) {
			for (unsigned int ii=0; ii<static_cast<unsigned int>(vecHTS1.size()); ii++) {
				vecData.push_back(vecHTS1.at(ii));
			}
		}
		vecMeteoData.push_back(vecData);
	}

	mysql_free_result(res2);
	mysql_stmt_close(stmt);
	mysql_close(conn2);
	return fullStation;
}

} //namespace
