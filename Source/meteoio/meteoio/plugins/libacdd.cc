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
#include <meteoio/plugins/libacdd.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOExceptions.h>

#include <fstream>
#include <cerrno>
#include <cstring>

using namespace std;

namespace mio {

/**
* @brief Read all config keys from the selected section and apply some special processing for some keys.
* @details This is used as some sort of caching, only keeping the section of interest.
* @param[in] cfg Config object to read the configuration keys from
* @param[in] section section to read the keys from (all keys from the section will be read at this point)
* @param[in] allow_multi_line should multi-line content be supported?
*/
void ACDD::setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line)
{
	for (size_t ii=0; ii<name.size(); ii++) {
		cfg.getValue(cfg_key[ii], section, value[ii], mio::IOUtils::nothrow);

		if (cfg_key[ii]=="ACDD_SUMMARY") { //overwrite with the content of summary_file if available
			const std::string summary_file = cfg.get("ACDD_SUMMARY_FILE", section, "");
			if (!summary_file.empty()) {
				std::string buffer;
				std::ifstream fin( summary_file.c_str() );
				if (fin.fail())
					throw mio::AccessException("Error opening ACDD_SUMMARY_FILE \""+summary_file+"\", possible reason: "+std::strerror(errno), AT);

				const char eoln = mio::FileUtils::getEoln(fin); //get the end of line character for the file
				try {
					do {
						std::string line;
						getline(fin, line, eoln); //read complete line
						if (allow_multi_line) buffer.append(line+"\n");
						else buffer.append(line+" ");
					} while (!fin.eof());
					fin.close();
				} catch (const std::exception&){
					if (fin.is_open()) fin.close();
					throw;
				}
				
				value[ii] = buffer;
			}
		}
	}
}

void ACDD::defaultInit()
{
	mio::Date now; 
	now.setFromSys();
	addAttribute("date_created", now.toString(mio::Date::ISO_DATE));
	addAttribute("creator_name", mio::IOUtils::getLogName(), "ACDD_CREATOR");
	addAttribute("creator_email", "", "ACDD_EMAIL");
	addAttribute("source", "MeteoIO-" + mio::getLibVersion(true));
	addAttribute("history", now.toString(mio::Date::ISO_Z) + ", " + mio::IOUtils::getLogName() + "@" + mio::IOUtils::getHostName() + ", MeteoIO-" + mio::getLibVersion(true));
	addAttribute("keywords_vocabulary", "AGU Index Terms");
	addAttribute("keywords", "Cryosphere, Mass Balance, Energy Balance, Atmosphere, Land/atmosphere interactions, Climatology", "ACDD_KEYWORDS");
	addAttribute("title", "", "ACDD_TITLE");
	addAttribute("institution", mio::IOUtils::getDomainName(), "ACDD_INSTITUTION");
	addAttribute("project", "", "ACDD_PROJECT");
	addAttribute("program", "", "ACDD_PROGRAM");
	addAttribute("id", "", "ACDD_ID");
	addAttribute("naming_authority", "", "ACDD_NAMING_AUTHORITY");
	addAttribute("processing_level", "", "ACDD_PROCESSING_LEVEL");
	addAttribute("summary", "", "ACDD_SUMMARY"); //special handling, see setUserConfig()
	addAttribute("comment", "", "ACDD_COMMENT");
	addAttribute("acknowledgement", "", "ACDD_ACKNOWLEDGEMENT");
	addAttribute("metadata_link", "", "ACDD_METADATA_LINK");
	addAttribute("license", "", "ACDD_LICENSE");
	addAttribute("product_version", "1.0", "ACDD_PRODUCT_VERSION");
}

/**
* @brief Add an attribute and its content to the internal list
* @details This allows to create or edit attributes. For the MERGE or APPEND modes, if the attribute name is not found, it will be created.
* @param[in] att_name attribute name
* @param[in] att_value attribute value
* @param[in] att_cfg_key associated configuration key (to read user provided values from a mio::Config object)
* @param[in] mode write mode: MERGE (currently empty values will be replaced by the given arguments), APPEND (the value content will be expanded by
* what is provided in att_value, separated by ", ", REPLACE (the current attribute will be fully replaced by the provided arguments)
*/
void ACDD::addAttribute(const std::string& att_name, const std::string& att_value, const std::string& att_cfg_key, Mode mode)
{
	if (att_name.empty())
		throw mio::InvalidFormatException("The attribute name must be provided", AT);
	
	if (mode==MERGE) {
		const size_t pos = find( att_name );
		if (pos==mio::IOUtils::npos) {
			mode = REPLACE;
		} else {
			if (!att_value.empty()) value[pos] = att_value;
			if (!att_cfg_key.empty()) cfg_key[pos] = att_cfg_key;
			return;
		}
	} else if (mode==APPEND) {
		const size_t pos = find( att_name );
		if (pos==mio::IOUtils::npos) {
			mode = REPLACE;
		} else {
			value[pos] = value[pos] + ", " + att_value;
			return;
		}
	}
	
	if (mode==REPLACE) {
		name.push_back( att_name );
		value.push_back( att_value );
		cfg_key.push_back( att_cfg_key );
		return;
	}
	
	//we should not have come here -> throw
	throw mio::InvalidFormatException("The specified write mode does not exists", AT);
}

void ACDD::addAttribute(const std::string& att_name, const double& att_value, const std::string& att_cfg_key, const Mode& mode)
{
	std::ostringstream os;
	os << att_value;
	addAttribute(att_name, os.str(), att_cfg_key, mode);
}

void ACDD::getAttribute(const size_t ii, std::string &att_name, std::string & att_value) const
{
	if (ii<name.size()) {
		att_name=name[ii];
		att_value=value[ii];
	} else {
		att_name="";
		att_value="";
	}
}

/**
* @brief Given an attribute name, return its associated index (or IOUtils::npos if it does not exists)
* @param[in] search_name attribute name to get the index for
* @return attribute index or IOUtils::npos
*/
size_t ACDD::find(const std::string& search_name) const
{
	for (size_t ii=0; ii<name.size(); ii++) {
		if (name[ii]==search_name) return ii;
	}
	
	return mio::IOUtils::npos;
}

void ACDD::setGeometry(const mio::Grid2DObject& grid, const bool& isLatLon)
{
	mio::Coords urcorner(grid.llcorner);
	urcorner.moveByXY(static_cast<double>(grid.getNx())*grid.cellsize, static_cast<double>(grid.getNy())*grid.cellsize);
	
	std::string epsg_str = "4326";
	std::string geometry;
	if (isLatLon) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid.llcorner.getLon() << " " << grid.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << grid.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << urcorner.getLat() << ", ";
		ss << grid.llcorner.getLon() << " " << urcorner.getLat();
		geometry = ss.str();
	}else {
		std::ostringstream os;
		os << grid.llcorner.getEPSG();
		epsg_str = os.str();
		
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid.llcorner.getEasting() << " " << grid.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << grid.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << urcorner.getNorthing() << ", ";
		ss << grid.llcorner.getEasting() << " " << urcorner.getNorthing();
		geometry = ss.str();
	}
	
	addAttribute("geospatial_bounds_crs", "EPSG:"+epsg_str);
	addAttribute("geospatial_bounds", "Polygon (("+geometry+"))");
}

void ACDD::setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon)
{
	if (vecMeteo.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	bool found = false;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty()) continue;

		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << vecMeteo[ii].front().meta.position.getLon() << " " << vecMeteo[ii].front().meta.position.getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << vecMeteo[ii].front().meta.position.getEasting() << " " << vecMeteo[ii].front().meta.position.getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : vecMeteo[ii].front().meta.position.getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=vecMeteo[ii].front().meta.position.getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = vecMeteo[ii].front().meta.position.getLat();
		const double curr_lon = vecMeteo[ii].front().meta.position.getLon();
		found = true;
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
	}
	if (!found) return;
	
	if (epsg>0) { //ie there is at least one valid point and all further points use the same epsg
		std::ostringstream os;
		os << epsg;
		addAttribute("geospatial_bounds_crs", "EPSG:"+os.str());
		addAttribute("geospatial_bounds", "MultiPoint ("+multiPts+")");
	}
	addAttribute("geospatial_lat_min", lat_min);
	addAttribute("geospatial_lat_max", lat_max);
	addAttribute("geospatial_lon_min", lon_min);
	addAttribute("geospatial_lon_max", lon_max);
}

void ACDD::setGeometry(const std::vector< mio::Coords >& vecLocation, const bool& isLatLon)
{
	if (vecLocation.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	bool found = false;
	for (size_t ii=0; ii<vecLocation.size(); ii++) {
		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << vecLocation[ii].getLon() << " " << vecLocation[ii].getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << vecLocation[ii].getEasting() << " " << vecLocation[ii].getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : vecLocation[ii].getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=vecLocation[ii].getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = vecLocation[ii].getLat();
		const double curr_lon = vecLocation[ii].getLon();
		found = true;
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
	}
	if (!found) return;
	
	const bool singlePoint = (lat_min==lat_max && lon_min==lon_max);
	
	if (epsg>0) { //ie there is at least one valid point and all further points use the same epsg
		std::ostringstream os;
		os << epsg;
		addAttribute("geospatial_bounds_crs", "EPSG:"+os.str());
		
		if (singlePoint)
			addAttribute("geospatial_bounds", "Point ("+multiPts+")");
		else
			addAttribute("geospatial_bounds", "MultiPoint ("+multiPts+")");
	}
	if (!singlePoint) {
		addAttribute("geospatial_lat_min", lat_min);
		addAttribute("geospatial_lat_max", lat_max);
		addAttribute("geospatial_lon_min", lon_min);
		addAttribute("geospatial_lon_max", lon_max);
	}
}

void ACDD::setGeometry(const mio::Coords& location, const bool& isLatLon)
{
	std::string epsg_str = "4326";
	std::string geometry;
	if (isLatLon) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << location.getLon() << " " << location.getLat();
		geometry = ss.str();
	}else {
		std::ostringstream os;
		os << location.getEPSG();
		epsg_str = os.str();
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(0) << location.getEasting() << " " << location.getNorthing();
	}
	addAttribute("geospatial_bounds_crs", "EPSG:"+epsg_str);
	addAttribute("geospatial_bounds", "Point ("+geometry+")");
}

void ACDD::setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo)
{
	if (vecMeteo.empty()) return;
	
	mio::Date set_start( vecMeteo[0].front().date );
	mio::Date set_end( vecMeteo[0].back().date );
	int sampling_period = -1;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) { //we must redo station 0 in order to get sampling_period
		if (vecMeteo[ii].empty()) continue;
		const mio::Date curr_start = vecMeteo[ii].front().date;
		const mio::Date curr_end = vecMeteo[ii].back().date;
		if (set_start>curr_start) set_start = curr_start;
		if (set_end<curr_end) set_end = curr_end;
		
		const size_t npts = vecMeteo[ii].size();
		if (npts>1) {
			const int curr_sampling = static_cast<int>( (curr_end.getJulian() - curr_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
			if (sampling_period<=0 || sampling_period>curr_sampling) sampling_period = curr_sampling;
		}
	}
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	if (sampling_period>0) {
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

void ACDD::setTimeCoverage(const std::vector<mio::MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return;
	
	const mio::Date set_start = vecMeteo.front().date;
	const mio::Date set_end = vecMeteo.back().date;
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	const size_t npts = vecMeteo.size();
	if (npts>1) {
		const int sampling_period = static_cast<int>( (set_end.getJulian() - set_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

void ACDD::setTimeCoverage(const std::vector<std::string>& vec_timestamp)
{
	if (vec_timestamp.empty()) return;
	
	mio::Date set_start, set_end;
	mio::IOUtils::convertString(set_start, vec_timestamp.front(), 0.); //in GMT
	mio::IOUtils::convertString(set_end, vec_timestamp.back(), 0.); //in GMT
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	const size_t npts = vec_timestamp.size();
	if (npts>1) {
		const int sampling_period = static_cast<int>( (set_end.getJulian() - set_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

} //namespace

