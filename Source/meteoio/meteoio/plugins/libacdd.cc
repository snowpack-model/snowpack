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
#include <algorithm>


using namespace std;

namespace mio {

/**
* @brief Set the available ACDD attributes, their matching INI config key and if possible a default value
* @return a map of <acdd_attribute_name, acdd_attrs structure>
*/
std::map<std::string, ACDD::acdd_attrs> ACDD::initAttributes()
{
	std::map<std::string, acdd_attrs> tmp;
	mio::Date now; 
	now.setFromSys();
	
	tmp["date_created"] = ACDD_ATTR("date_created", "", now.toString(mio::Date::ISO_DATE));
	tmp["institution"] = ACDD_ATTR("institution", "ACDD_INSTITUTION", mio::IOUtils::getDomainName());
	
	tmp["creator_name"] = ACDD_ATTR("creator_name", "ACDD_CREATOR", mio::IOUtils::getLogName());
	tmp["creator_email"] = ACDD_ATTR("creator_email", "ACDD_CREATOR_EMAIL");
	tmp["creator_institution"] = ACDD_ATTR("creator_institution", "ACDD_CREATOR_INSTITUTION", mio::IOUtils::getDomainName());
	tmp["creator_url"] = ACDD_ATTR("creator_url", "ACDD_CREATOR_URL", mio::IOUtils::getDomainName());
	tmp["creator_type"] = ACDD_ATTR("creator_type", "ACDD_CREATOR_TYPE", "person");
	tmp["contributor_name"] = ACDD_ATTR("contributor_name", "ACDD_CONTRIBUTOR");
	tmp["contributor_role"] = ACDD_ATTR("contributor_role", "ACDD_CONTRIBUTOR_ROLE");
	tmp["publisher_name"] = ACDD_ATTR("publisher_name", "ACDD_PUBLISHER");
	tmp["publisher_email"] = ACDD_ATTR("publisher_email", "ACDD_PUBLISHER_EMAIL");
	tmp["publisher_url"] = ACDD_ATTR("publisher_url", "ACDD_PUBLISHER_URL");
	tmp["publisher_type"] = ACDD_ATTR("publisher_type", "ACDD_PUBLISHER_TYPE");
	
	tmp["source"] = ACDD_ATTR("source", "ACDD_SOURCE", "MeteoIO-" + mio::getLibVersion(true));
	tmp["history"] = ACDD_ATTR("history", "", now.toString(mio::Date::ISO_Z) + ", " + mio::IOUtils::getLogName() + "@" + mio::IOUtils::getHostName() + ", MeteoIO-" + mio::getLibVersion(true));
	tmp["keywords_vocabulary"] = ACDD_ATTR("keywords_vocabulary", "ACDD_KEYWORDS_VOCABULARY", "AGU Index Terms");
	tmp["keywords"] = ACDD_ATTR("keywords", "ACDD_KEYWORDS", "Cryosphere, Mass Balance, Energy Balance, Atmosphere, Land/atmosphere interactions, Climatology");
	
	tmp["title"] = ACDD_ATTR("title", "ACDD_TITLE");
	tmp["project"] = ACDD_ATTR("project", "ACDD_PROJECT");
	tmp["program"] = ACDD_ATTR("program", "ACDD_PROGRAM");
	tmp["id"] = ACDD_ATTR("id", "ACDD_ID");
	tmp["references"] = ACDD_ATTR("references", "ACDD_REFERENCES");
	tmp["naming_authority"] = ACDD_ATTR("naming_authority", "ACDD_NAMING_AUTHORITY");
	tmp["processing_level"] = ACDD_ATTR("processing_level", "ACDD_PROCESSING_LEVEL");
	tmp["summary"] = ACDD_ATTR("summary", "ACDD_SUMMARY"); //special handling, see setUserConfig()
	tmp["comment"] = ACDD_ATTR("comment", "ACDD_COMMENT");
	tmp["acknowledgement"] = ACDD_ATTR("acknowledgement", "ACDD_ACKNOWLEDGEMENT");
	tmp["metadata_link"] = ACDD_ATTR("metadata_link", "ACDD_METADATA_LINK");
	tmp["license"] = ACDD_ATTR("license", "ACDD_LICENSE");
	tmp["product_version"] = ACDD_ATTR("product_version", "ACDD_PRODUCT_VERSION", "1.0");
	tmp["activity_type"] = ACDD_ATTR("activity_type", "ACDD_ACTIVITY_TYPE");
	tmp["operational_status"] = ACDD_ATTR("operational_status", "ACDD_OPERATIONAL_STATUS");
	tmp["wmo__wsi"] = ACDD_ATTR("wmo__wsi", "WIGOS_ID");
	
	return tmp;
}

/**
* @brief Declare the attributes that support multiple, comma delimited valuesd and in this case must have the same number of elements
* @return a set of pairs of <group_name, acdd_attribute_name>
*/
std::set< std::pair< std::string, std::set<std::string> > > ACDD::initLinks()
{
	static const std::set<std::string> creators( {"creator_name", "creator_email", "creator_institution", "creator_url", "creator_type"} );
	static const std::set<std::string> publishers( {"publisher_name", "publisher_email", "publisher_url", "publisher_type"} );
	static const std::set<std::string> contributors( {"contributor_name", "contributor_role"} );
	
	static const std::set< std::pair< std::string, std::set<std::string> > > tmp( {std::make_pair("CREATOR", creators), std::make_pair("PUBLISHER", publishers), std::make_pair("CONTRIBUTOR", contributors)} );
	return tmp;
}

void ACDD::acdd_attrs::readFromFile(std::string& value, const mio::Config& cfg, const std::string& cfg_key, const std::string& section, const bool& allow_multi_line)
{
	const std::string input_file = cfg.get(cfg_key, section, "");
	std::string buffer;
	
	if (!input_file.empty()) {
		std::ifstream fin( input_file.c_str() );
		if (fin.fail())
			throw mio::AccessException("Error opening "+cfg_key+" \""+input_file+"\", possible reason: "+std::strerror(errno), AT);

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
		
		value = buffer;
	}
}

/**
* @brief Set the value of the attribute, either from a config file key, or from an environment variable of the same name or
* some hard-coded default value (in this order of priorities). The key "ACDD_SUMMARY" is handled in a special way: if
* an "ACDD_SUMMARY_FILE" key is found, the file pointed to by the later will be read and copied into ACDD_SUMMARY. Line
* breaks are either kept as such or converted to a single white space, depending on the parameter <i>allow_multi_line</i>
* @param[in] cfg The configuration file to read the keys from (for the attributes that can be set from such a config key)
* @param[in] section Section in the configuration file to read the keys from
* @param[in] allow_multi_line If an "ACDD_SUMMARY_FILE" is set, either keep line breaks in the file it points to or replace them by spaces.
*/
void ACDD::acdd_attrs::setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line)
{
	if (!cfg_key.empty()) {
		//first priority: read from cfg file
		cfg.getValue(cfg_key, section, value, mio::IOUtils::nothrow);
		if (cfg_key=="ACDD_SUMMARY") { //overwrite with the content of summary_file if available
			readFromFile(value, cfg, "ACDD_SUMMARY_FILE", section, allow_multi_line);
		}
		
		//second priority: read from env. var
		if (value.empty()) {
			char *tmp = getenv( cfg_key.c_str() );
			if (tmp!=nullptr) value = std::string(tmp);
		}
	}
	
	//last priority: set from default
	if (value.empty() && !default_value.empty()) {
		value = default_value;
		Default = true;
	} else {
		Default = false;
	}
}

void ACDD::acdd_attrs::setValue(const std::string& i_value, const Mode& mode)
{
	if (mode==MERGE) {
		if (value.empty()) value = i_value;
	} else if (mode==APPEND) {
		if (value.empty()) value = i_value;
		else value = value + ", " + i_value;
	} else if (mode==REPLACE) {
		value = i_value;
	}
	
	Default = false;
}

/**
* @brief Read all config keys from the selected section and apply some special processing for some keys.
* @details This is used as some sort of caching, only keeping the section of interest.
* @param[in] cfg Config object to read the configuration keys from
* @param[in] section section to read the keys from (all keys from the section will be read at this point)
* @param[in] allow_multi_line should multi-line content be supported?
*/
void ACDD::setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line)
{
	for (auto& acdd_attribute : attributes) {
		acdd_attribute.second.setUserConfig(cfg, section, allow_multi_line);
	}
	
	if (attributes.empty()) attributes = initAttributes();
	if (linked_attributes.empty()) {
		linked_attributes = initLinks();
		checkLinkedAttributes();
	}
}

void ACDD::setEnabled(const bool& i_enable)
{
	enabled = i_enable; 
	if (attributes.empty()) attributes = initAttributes();
	if (linked_attributes.empty()) {
		linked_attributes = initLinks();
		checkLinkedAttributes();
	}
}

size_t ACDD::countCommas(const std::string& str)
{
	const size_t count = std::count_if( str.begin(), str.end(), []( const char& c ){return c ==',';});
	return count;
}

void ACDD::checkLinkedAttributes()
{
	for (const std::pair< std::string, std::set<std::string> >& attr_group : linked_attributes) {
		std::set< std::string > tweakDefaults;
		
		//check if all non-default attributes have the same number of sub-elements
		size_t min_count = IOUtils::npos, max_count = IOUtils::npos;
		
		for (const auto& attr : attr_group.second) {
			if (attributes.count(attr)==0) continue;
			
			const auto& attribute( attributes[attr] );
			const std::string& value( attribute.getValue() );
			if (value.empty()) continue;
			if (attribute.isDefault()) {
				tweakDefaults.insert( attribute.getName() );
				continue;
			}
			
			const size_t num_elems = countCommas( value ) + 1;
			if (max_count==IOUtils::npos || num_elems>max_count) max_count = num_elems;
			if (min_count==IOUtils::npos || num_elems<min_count) min_count = num_elems;
		}
		
		if (min_count!=max_count) throw mio::InvalidFormatException("Please configure the same number of fields for each comma-delimited ACDD fields of type '"+attr_group.first+"'", AT);
		
		//copy more default values if necessary
		for (const std::string& attr : tweakDefaults) {
			const std::string& value( attributes[attr].getValue() );
			std::string tmp( value );
			for (size_t ii=1; ii<max_count; ii++) tmp += ", "+value;
			attributes[attr].setValue( tmp, REPLACE );
		}
	}
}

std::string ACDD::toString() const
{
	std::ostringstream os;
	os << "<ACDD attributes>\n";
	for (auto acdd_attribute : attributes) {
		os << "[" << acdd_attribute.first << " -> " << acdd_attribute.second.getValue() << "]\n";
	}
	os << "</ACDD attributes>\n";

	return os.str();
}

/**
* @brief Add an attribute and its content to the internal list
* @details This allows to create or edit attributes. For the MERGE or APPEND modes, if the attribute name is not found, it will be created.
* @param[in] att_name attribute name
* @param[in] att_value attribute value
* @param[in] mode write mode: MERGE (currently empty values will be replaced by the given arguments), APPEND (the value content will be expanded by
* what is provided in att_value, separated by ", ", REPLACE (the current attribute will be fully replaced by the provided arguments)
*/
void ACDD::addAttribute(const std::string& att_name, const std::string& att_value, const Mode& mode)
{
	if (att_name.empty())
		throw mio::InvalidFormatException("The attribute name must be provided", AT);

	if (mode==MERGE || mode==APPEND) {
		if (attributes.count(att_name)==0) {
			attributes[att_name] = acdd_attrs(att_name, att_value, "", "");
		} else {
			attributes[att_name].setValue(att_value, mode);
		}
	} else if (mode==REPLACE) {
		attributes[att_name] = acdd_attrs(att_name, att_value, "", "");
	}
}

void ACDD::addAttribute(const std::string& att_name, const double& att_value, const Mode& mode)
{
	std::ostringstream os;
	os << att_value;
	addAttribute(att_name, os.str(), mode);
}

/**
* @brief Given an attribute name, return its associated value (or an empty string if it does not exists)
* @param[in] att_name attribute name to get the value for
* @return attribute value or empty string
*/
std::string ACDD::getAttribute(std::string &att_name) const
{
	const auto& it = attributes.find( att_name );
	if (it==attributes.end()) return "";

	return it->second.getValue();
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

	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	const double min_val = grid.getMin();
	if (min_val!=IOUtils::nodata) {
		//if there is a min_val, there is also a max_val
		addAttribute("geospatial_vertical_min", min_val);
		const double max_val = grid.getMax();
		addAttribute("geospatial_vertical_max", max_val);
	}


}

void ACDD::setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon)
{
	if (vecMeteo.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	double alt_min=std::numeric_limits<double>::max(), alt_max=-std::numeric_limits<double>::max();
	bool found = false;
	for (const std::vector<mio::MeteoData>& timeseries : vecMeteo) {
		if (timeseries.empty()) continue;

		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << timeseries.front().meta.position.getLon() << " " << timeseries.front().meta.position.getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << timeseries.front().meta.position.getEasting() << " " << timeseries.front().meta.position.getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : timeseries.front().meta.position.getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=timeseries.front().meta.position.getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = timeseries.front().meta.position.getLat();
		const double curr_lon = timeseries.front().meta.position.getLon();
		const double curr_alt = timeseries.front().meta.position.getAltitude();
		found = true;
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
		if (alt_min>curr_alt) alt_min = curr_alt;
		if (alt_max<curr_alt) alt_max = curr_alt;
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

	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	if (alt_min!=IOUtils::nodata) {
		addAttribute("geospatial_vertical_min", alt_min);
		addAttribute("geospatial_vertical_max", alt_max);
	}
}

void ACDD::setGeometry(const std::vector< mio::Coords >& vecLocation, const bool& isLatLon)
{
	if (vecLocation.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	double alt_min=std::numeric_limits<double>::max(), alt_max=-std::numeric_limits<double>::max();
	
	for (const mio::Coords& location : vecLocation) {
		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << location.getLon() << " " << location.getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << location.getEasting() << " " << location.getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : location.getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=location.getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = location.getLat();
		const double curr_lon = location.getLon();
		const double curr_alt = location.getAltitude();
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
		if (alt_min>curr_alt) alt_min = curr_alt;
		if (alt_max<curr_alt) alt_max = curr_alt;
	}
	
	if (epsg>0) { //ie there is at least one valid point and all further points use the same epsg
		std::ostringstream os;
		os << epsg;
		addAttribute("geospatial_bounds_crs", "EPSG:"+os.str());
		
		const bool singlePoint = (lat_min==lat_max && lon_min==lon_max);
		if (singlePoint)
			addAttribute("geospatial_bounds", "Point ("+multiPts+")");
		else
			addAttribute("geospatial_bounds", "MultiPoint ("+multiPts+")");
	}
	
	addAttribute("geospatial_lat_min", lat_min);
	addAttribute("geospatial_lat_max", lat_max);
	addAttribute("geospatial_lon_min", lon_min);
	addAttribute("geospatial_lon_max", lon_max);
	
	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	if (alt_min!=IOUtils::nodata) {
		addAttribute("geospatial_vertical_min", alt_min);
		addAttribute("geospatial_vertical_max", alt_max);
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
	} else {
		std::ostringstream os;
		os << location.getEPSG();
		epsg_str = os.str();
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(0) << location.getEasting() << " " << location.getNorthing();
	}
	addAttribute("geospatial_bounds_crs", "EPSG:"+epsg_str);
	addAttribute("geospatial_bounds", "Point ("+geometry+")");
	addAttribute("geospatial_lat_min", location.getLat());
	addAttribute("geospatial_lat_max", location.getLat());
	addAttribute("geospatial_lon_min", location.getLon());
	addAttribute("geospatial_lon_max", location.getLon());
	
	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	addAttribute("geospatial_vertical_min", location.getAltitude());
	addAttribute("geospatial_vertical_max", location.getAltitude());
}

void ACDD::setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo)
{
	if (vecMeteo.empty()) return;
	
	mio::Date set_start( vecMeteo[0].front().date );
	mio::Date set_end( vecMeteo[0].back().date );
	int sampling_period = -1;
	for (const std::vector<mio::MeteoData>& timeseries : vecMeteo) { //we must redo station 0 in order to get sampling_period
		if (timeseries.empty()) continue;
		const mio::Date curr_start( timeseries.front().date );
		const mio::Date curr_end( timeseries.back().date );
		if (set_start>curr_start) set_start = curr_start;
		if (set_end<curr_end) set_end = curr_end;
		
		const size_t npts = timeseries.size();
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
	
	const mio::Date set_start( vecMeteo.front().date );
	const mio::Date set_end( vecMeteo.back().date );
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

void ACDD::setTimeCoverage(const std::vector<std::string>& vec_timestamp, const double& TZ)
{
	if (vec_timestamp.empty()) return;
	
	mio::Date set_start, set_end;
	mio::IOUtils::convertString(set_start, vec_timestamp.front(), TZ);
	mio::IOUtils::convertString(set_end, vec_timestamp.back(), TZ);
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

