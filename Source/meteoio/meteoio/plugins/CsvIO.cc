/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/CsvIO.h>

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <utility>
#include <errno.h>

using namespace std;

namespace mio {
/**
 * @page csvio CsvIO
 * @section csvio_format Format
 * This plugins offers a flexible way to read Comma Separated Values (<A HREF="https://en.wikipedia.org/wiki/Comma-separated_values">CSV</A>) files. 
 * It is however assumed that:
 *     - each line contains a data record
 *     - each line contains the same number of fields;
 *     - a single character is consistently used through the file as field delimiter (to split each record into fields);
 *     - missing data are represented by an empty value, so two delimiters follow directly each other or by a special value (see NODATA in 
 * \ref csvio_metadata_extraction "Metadata extraction");
 *     - the file may contain a header that may contain additional information (metadata), see below.
 * 
 * In order to reduce the amount of manual configuration, it is possible to extract metadata from the headers or the filename, 
 * such as the station name, ID, coordinates, etc
 *
 * @section csvio_units Units
 * The final units MUST be SI. If not, the conversion offsets/factors must be provided to convert the data back to SI (see required keywords below)
 * or the units declared (in the headers) and supported by this plugin.
 *
 * @section csvio_keywords Keywords
 * This plugin uses the following keywords, in the [Input] section:
 * - COORDSYS: coordinate system (see Coords);
 * - COORDPARAM: extra coordinates parameters (see Coords);
 * - TIME_ZONE: the timezone that should be used to interpret the dates/times (default: 0);
 * - METEOPATH: the directory where the data files are available (mandatory);
 * - STATION\#: input filename (in METEOPATH). As many meteofiles as needed may be specified. If nothing is specified, the METEOPATH directory will be scanned for files with the extension specified in CSV_FILE_EXTENSION;
 * - METEOPATH_RECURSIVE: if set to true, the scanning of METEOPATH is performed recursively (default: false);
 * - POSITION\#: coordinates of the station (default: reading key "POSITION"; see (see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax));
 * - CSV_FILE_EXTENSION: When scanning the whole directory, look for these files (default: .csv). Note that this matching isn't restricted to the end of the file name so if you had files stat1_jan.csv, stat1_feb.csv and stat2_jan.csv you could select January's data by putting "_jan" here;
 * - CSV_SILENT_ERRORS: if set to true, lines that can not be read will be silently ignored (default: false, has priority over CSV_ERRORS_TO_NODATA);
 * - CSV_ERRORS_TO_NODATA: if true, unparseable fields (like text fields) are set to nodata, but the rest of the line is kept (default: false).
 * 
 * The following keys may either be prefixed by "CSV_" (ie as default for all stations) or by "CSV#_" (as only for the current station):
 * - CSV\#_DELIMITER: field delimiter to use (default: ',');
 * - CSV\#_NR_HEADERS: how many lines should be treated as headers? (default: 1);
 * - CSV\#_COLUMNS_HEADERS: header line to interpret as columns headers (default: 1);
 * - CSV\#_FIELDS: columns headers (if they don't exist in the file or to overwrite them); optional
 * - CSV\#_FILENAME_SPEC: pattern to parse the filename and extract metadata out of it; optional
 * - CSV\#_UNITS_HEADERS: header line providing the measurements units (the subset of recognized units is small, please inform us if one is missing for you); optional
 * - CSV\#_UNITS_OFFSET: offset to add to each value in order to convert it to SI; optional
 * - CSV\#_UNITS_MULTIPLIER: factor to multiply each value by, in order to convert it to SI; optional
 * - CSV\#_DATETIME_SPEC: mixed date and time format specification (defaultis ISO_8601: YYYY-MM-DDTHH24:MI:SS);
 * - CSV\#_DATE_SPEC: date format specification (default: YYYY_MM_DD);
 * - CSV\#_TIME_SPEC: time format specification (default: HH24:MI:SS);
 * - CSV\#_SPECIAL_HEADERS: description of how to extract more metadata out of the headers; optional
 * - CSV\#_NODATA: a value that should be interpreted as *nodata* (default: NAN);
 * - CSV\#_NAME: the station name to use (if provided, has priority over the special headers);
 * - CSV\#_ID: the station id to use (if provided, has priority over the special headers);
 * 
 * If no ID has been provided, an automatic station ID will be generated as "ID{n}" where *n* is the current station's index. Regarding the units handling, 
 * it is only performed through either the CSV_UNITS_OFFSET key or the CSV_UNITS_OFFSET / CSV_UNITS_MULTIPLIER keys. These keys expect a value for each
 * column of the file, including the date and time.
 * 
 * @note Since most parameter won't have names that are recognized by MeteoIO, it is advised to map them to \ref meteoparam "MeteoIO's internal names". 
 * This is done either by using the CSV_FIELDS key or using the \ref data_move "data renaming" feature of the 
 * \ref data_manipulations "Raw Data Editing" stage.
 * 
 * @section csvio_date_specs Date and time specification
 * In order to be able to read any date and time format, the format has to be provided in the configuration file. This is provided as a string containing
 * the following special markers:
 * - YYYY: the 4 digits year;
 * - MM: the two digits month;
 * - DD: the two digits day;
 * - HH24: the two digits hour of the day (0-24);
 * - MI: the two digits minutes (0-59);
 * - SS: the two digts seconds (0-59);
 * - TZ: the numerical timezone as offset to GMT (see note below).
 *
 * Any other character is interpreted as itself, present in the string. It is possible to either provide a combined datetime field (so date and time are combined into
 * one single field) or date and time as two different fields. For example:
 * - YYYY-MM-DDTHH24:MI:SS described an <A HREF="https://en.wikipedia.org/wiki/ISO_8601">ISO 8601</A> datetime field;
 * - MM/DD/YYYY described an anglo-saxon date;
 * - DD.MM.YYYY HH24:MI:SS is for a Swiss formatted datetime.
 * 
 * @note When providing a timezone field, it \em must appear at the end of the string. it can either be numerical (such as "+1.") or an abbreviation
 * such as "CET" (see https://en.wikipedia.org/wiki/List_of_time_zone_abbreviations).
 * 
 * @section csvio_metadata_extraction Metadata extraction
 * Since there is no unified way of providing metadata (such as the location, station name, etc) in CSV files, this information has to
 * be either provided in the configuration file (see \ref csvio_keywords "Configuration keywords") or extracted out of either the file
 * name or the file headers. A specific syntax allows to describe where to find which metadata field type.
 * 
 * @subsection csvio_metadata_field_types Metadata fields types
 * The following field types are supported:
 * - NAME;
 * - ID (this will the used as a handle for the station);
 * - ALT (for the altitude);
 * - LON (for the longitude);
 * - LAT (for the latitude);
 * - SLOPE (in degrees);
 * - AZI (for the slope azimuth, in degree as read from a compass);
 * - NODATA (string to interpret as nodata);
 * - PARAM (to identify the content of a file that only contains the time information and one meteorological parameter);
 * - SKIP (skip this field).
 *
 * @subsection csvio_special_headers Header metadata extraction
 * This is performed with the "CSV#_SPECIAL_HEADERS" configuration key. This key is followed by as many metadata 
 * specifications as necessary, of the form {field}:{line}:{column}.
 *
 * Therefore, if the station name is available on line 1, column 3 and the station id on line 2, column 5, the configuration would be:
 * @code
 * CSV_SPECIAL_HEADERS = name:1:3 id:2:5
 * @endcode
 * 
 * @subsection csvio_filename_parsing Filename metadata extraction
 * This is performed with the "CSV#_FILENAME_SPEC" configuration key. This key is followed by the metadata specification 
 * that will be applied to identify the information to extract as well as substrings that are used as "markers" delimiting 
 * the different fields (enclosed within {}).
 * 
 * For example, to parse the filename "H0118_Generoso-Calmasino_-_Precipitation.DAT" use (please note that the extension is NOT provided):
 * @code
 * CSV_FILENAME_SPEC = {ID}_{NAME}-{SKIP}_-_{PARAM}
 * @endcode
 * 
 * If the CSV_FIELDS key is also present, it will have priority. Therefore, it is possible to define one CSV_FILENAME_SPEC for several files and 
 * only define CSV\#_FIELDS for the files that would require a different handling (for example because their parameter would not be recognized).
 * Moreover, it is possible to set \em "AUTOMERGE" to "true" in the input section, so all files leading to the same station ID will be merged together
 * into one single station.
 * 
 * @note Obviously, the {PARAM} metadata field type can only be used for files that contain the time information (either as datetime or seperate date and time) and one
 * meteorological parameter.
 * 
 * @section csvio_examples Examples
 * In order to read a CSV file produced by a Campbell data logger with Swiss-formatted timestamps, you need the following configuration:
 * @code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_NR_HEADERS = 4
 * CSV_COLUMNS_HEADERS = 2
 * CSV_UNITS_HEADERS = 3
 * CSV_DATETIME_SPEC = DD.MM.YYYY HH24:MI:SS
 * CSV_SPECIAL_HEADERS = name:1:2 id:1:4
 * 
 * STATION1 = DisMa_DisEx.dat
 * POSITION1 = latlon 46.810325 9.806657 2060
 * CSV1_ID = DIS4
 * @endcode
 *
 * In order to read a set of files containing each only one parameter and merge them together (see \ref data_manipulations "raw data editing" for more
 * on the merge feature), extracting the station ID, name and meteorological parameter from the filename:
 *@code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_DELIMITER = ;
 * CSV_HEADER_LINES = 1
 * CSV_DATE_SPEC = DD/MM/YYYY
 * CSV_TIME_SPEC = HH24:MI
 * POSITION = latlon (46.8, 9.80, 1700)
 * CSV_FILENAME_SPEC = {ID}_{NAME}_-_{SKIP}-{PARAM}
 * CSV_COLUMNS_HEADERS = 1
 *
 * STATION1 = H0118_Generoso_-_Calmasino_precipitation.DAT
 *
 * STATION2 = H0118_Generoso_-_Calmasino_temperature.DAT    #the parameter name is ambiguous, it will not be recognized
 * CSV2_FIELDS = DATE TIME TA                               #so we define the parameter manually
 * CSV2_UNITS_OFFSET = 0 0 273.15
 *
 * STATION3 = H0118_Generoso_-_Calmasino_reflected_solar_radiation.DAT
 * STATION4 = H0118_Generoso_-_Calmasino_relative_humidity.DAT
 * STATION5 = H0118_Generoso_-_Calmasino_wind_velocity.DAT
 *
 * AUTOMERGE = true
 * @endcode
 * 
 * In order to read a set of files and merge them together (see \ref data_manipulations "raw data editing" for more
 * on the merge feature):
 *@code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_DELIMITER = ;
 * CSV_HEADER_LINES = 1
 * CSV_DATE_SPEC = DD/MM/YYYY
 * CSV_TIME_SPEC = HH24:MI
 * POSITION = latlon (46.8, 9.80, 1700)
 * CSV_NAME = Generoso
 *
 * CSV1_FIELDS = DATE TIME PSUM HS
 * STATION1 = H0118_lg23456.DAT
 *
 * CSV2_FIELDS = DATE TIME TA
 * STATION2 = H0118_lg7850.DAT
 * CSV2_UNITS_OFFSET = 0 0 273.15
 *
 * CSV3_FIELDS = DATE TIME RSWR ISWR
 * STATION3 = H0118_lg64520.DAT
 *
 * CSV4_FIELDS = DATE TIME RH
 * STATION4 = H0118_lg45302.DAT
 * CSV4_UNITS_MULTIPLIER = 1 1 0.01
 *
 * CSV5_FIELDS = DATE TIME VW
 * STATION5 = H0118_wind_velocity.DAT
 *
 * ID1::MERGE = ID2 ID3 ID4 ID5
 * @endcode
 */

//small helper functions for removal/replacement in strings
inline bool isQuote(const char& c) { if (c=='"' || c=='\'') return true; return false;}
inline bool InvalidChar(const char& c) { return (c<32 || c>126); }
//helper function to sort the static keys used for specifying the date/time formats
inline bool sort_dateKeys(const std::pair<size_t,size_t> &left, const std::pair<size_t,size_t> &right) { return left.first < right.first;}

//parse the user provided special headers specification. It is stored as <line_nr, <column, field_type>> in a multimap
//(since there can be multiple keys on the same line)
std::multimap< size_t, std::pair<size_t, std::string> > CsvParameters::parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec)
{
	std::multimap< size_t, std::pair<size_t, std::string> > meta_spec;
	for (size_t ii=0; ii<vecMetaSpec.size(); ii++) {
		std::vector<std::string> vecArgs;
		if (IOUtils::readLineToVec(vecMetaSpec[ii], vecArgs, ':') !=3)
			throw InvalidFormatException("Wrong format for Metadata specification '"+vecMetaSpec[ii]+"'", AT);
		const int linenr = atoi( vecArgs[1].c_str() );
		const int colnr = atoi( vecArgs[2].c_str() );
		if (linenr<=0 || colnr<=0)
			throw InvalidFormatException("Line numbers and column numbers must be >0 in Metadata specification '"+vecMetaSpec[ii]+"'", AT);
		
		meta_spec.insert( make_pair( linenr, make_pair( colnr, vecArgs[0]) ) );
	}
	
	return meta_spec;
}

//Given a provided field_type, attribute the value to the proper metadata variable
void CsvParameters::assignMetadataVariable(const std::string& field_type, const std::string& field_val, double &lat, double &lon)
{
	if (field_type=="ID") {
		if (id.empty()) id = field_val;
	} else if (field_type=="NAME") {
		if (name.empty()) name = field_val;
	} else if (field_type=="NODATA") {
			nodata = field_val;
	} else if (field_type=="SKIP") {
		return;
	} else if (field_type=="PARAM") {
		std::string param( IOUtils::strToUpper( field_val ) );
		if (MeteoData::getStaticParameterIndex( param )!=IOUtils::npos) {
			single_field = param;
			return;
		}
		
		param.erase(std::remove_if(param.begin(), param.end(), &InvalidChar), param.end()); //remove accentuated characters, etc
		//try to map non-standard names to mio's names
		if (param.compare(0, 15, "TEMPERATURE_AIR")==0 || param.compare(0, 16, "TEMPERATURA_ARIA")==0) param="TA";
		else if (param.compare(0, 13, "PRECIPITATION")==0 || param.compare(0, 14, "PRECIPITAZIONE")==0) param="PSUM";
		else if (param.compare(0, 19, "REFLECTED_RADIATION")==0 || param.compare(0, 26, "RADIAZIONE_SOLARE_RIFLESSA")==0) param="RSWR";
		else if (param.compare(0, 18, "INCOMING_RADIATION")==0 || param.compare(0, 27, "RADIAZIONE_SOLARE_INCIDENTE")==0) param="RSWR";
		else if (param.compare(0, 14, "WIND_DIRECTION")==0 || param.compare(0, 15, "DIREZIONE_VENTO")==0) param="DW";
		else if (param.compare(0, 17, "RELATIVE_HUMIDITY")==0 || param.compare(0, 15, "UMIDIT_RELATIVA")==0) param="RH";
		else if (param.compare(0, 13, "WIND_VELOCITY")==0 || param.compare(0, 13, "VELOCIT_VENTO")==0) param="VW";
		
		single_field = param;
	} else {
		if (field_type=="ALT" || field_type=="LON" || field_type=="LAT" || field_type=="SLOPE" || field_type=="AZI") {
			double tmp;
			if (!IOUtils::convertString(tmp, field_val))
				throw InvalidArgumentException("Could not extract metadata '"+field_type+"' for "+file_and_path, AT);
			
			if (field_type=="ALT") location.setAltitude( tmp, false);
			if (field_type=="LON") lon = tmp;
			if (field_type=="LAT") lat = tmp;
			if (field_type=="SLOPE") slope = tmp;
			if (field_type=="AZI") azi = tmp;
		} else 
			throw InvalidFormatException("Unknown parsing key '"+field_type+"' when extracting metadata", AT);
	}
}

//Using the special headers parsed specification (done in parseHeadersSpecs), some metadata is extracted from the headers
void CsvParameters::parseSpecialHeaders(const std::string& line, const size_t& linenr, const std::multimap< size_t, std::pair<size_t, std::string> >& meta_spec, double &lat, double &lon)
{
	std::vector<std::string> vecStr;
	IOUtils::readLineToVec(line, vecStr, csv_delim);
	
	std::multimap<size_t, std::pair<size_t, std::string> >::const_iterator it;
	for (it=meta_spec.equal_range(linenr).first; it!=meta_spec.equal_range(linenr).second; ++it) {
		const size_t colnr = (*it).second.first;
		const std::string field_type( IOUtils::strToUpper( (*it).second.second) );
		if (colnr>vecStr.size() || colnr==0)
			throw InvalidArgumentException("Metadata specification for '"+field_type+"' refers to a non-existent field for file (either 0 or too large) "+file_and_path, AT);
		
		//remove the quotes from the field
		std::string field_val( vecStr[colnr-1] );
		field_val.erase(std::remove_if(field_val.begin(), field_val.end(), &isQuote), field_val.end());
		
		assignMetadataVariable(field_type, field_val, lat, lon);
	}
}

//Extract metadata from the filename, according to a user-provided specification.
//filename_spec in the form of {ID}_{NAME}-{PARAM}_-_{SKIP} where {ID} is a variable and '_-_' is a constant pattern
void CsvParameters::parseFileName(std::string filename, const std::string& filename_spec, double &lat, double &lon)
{
	filename = FileUtils::removeExtension( FileUtils::getFilename(filename) );
	size_t pos_fn = 0, pos_mt = 0; //current position in the filename and in the filename_spec
	if (filename_spec[0]!='{') { //there is a constant pattern at the begining, getting rid of it
		const size_t start_var = filename_spec.find('{');
		if (start_var==std::string::npos) throw InvalidFormatException("No variables defined for filename parsing", AT);
		const std::string pattern( filename_spec.substr(0, start_var) );
		if (filename.substr(0, start_var)!=pattern) throw InvalidFormatException("The filename pattern does not match with the given filename", AT);
		pos_mt = start_var;
		pos_fn = start_var;
	}
	
	//we now assume that we start with a variable
	do {
		//the start of the next constant pattern defines the end of the current variable
		const size_t start_pattern = filename_spec.find('}', pos_mt);
		const size_t end_pattern = filename_spec.find('{', pos_mt+1);
		if (start_pattern==std::string::npos) {
			if (end_pattern!=std::string::npos) throw InvalidFormatException("Unclosed variable delimiter '}' in filename parsing", AT);
			break; //no more variables to read
		}
		const size_t pattern_len = end_pattern - start_pattern - 1;
		
		size_t len_var = std::string::npos; //until end of string
		if (end_pattern!=std::string::npos) {
			const std::string pattern = filename_spec.substr(start_pattern+1, pattern_len); //skip } and {
			const size_t pos_pattern_fn = filename.find(pattern, pos_fn);
			if (pos_pattern_fn==std::string::npos) throw InvalidFormatException("The filename pattern does not match with the given filename", AT);
			len_var = pos_pattern_fn - pos_fn;
		}

		//read the variable type and value
		const std::string field_type( IOUtils::strToUpper(filename_spec.substr(pos_mt+1, start_pattern-pos_mt-1)) ); //skip { and }
		const std::string value( filename.substr(pos_fn, len_var) );
		assignMetadataVariable(field_type, value, lat, lon);

		if (end_pattern==std::string::npos) break; //nothing more to parse
		pos_mt = end_pattern;
		pos_fn = pos_fn + len_var + pattern_len;
	} while (true);

}

//user provided field names are in fieldNames, header field names are in headerFields
//and user provided fields have priority.
void CsvParameters::parseFields(const std::vector<std::string>& headerFields, std::vector<std::string>& fieldNames, size_t &dt_col, size_t &tm_col)
{
	const bool user_provided_field_names = (!fieldNames.empty());
	if (headerFields.empty() && !user_provided_field_names) 
		throw InvalidArgumentException("No columns names could be found. Please either provide CSV_COLUMNS_HEADERS or CSV_FIELDS", AT);
	
	if (fieldNames.empty()) fieldNames = headerFields;
	for (size_t ii=0; ii<fieldNames.size(); ii++) {
		std::string &tmp = fieldNames[ii];
		IOUtils::toUpper( tmp );
		tmp.erase(std::remove_if(tmp.begin(), tmp.end(), &isQuote), tmp.end());
		if (tmp.empty()) continue;
		
		if (tmp.compare("TIMESTAMP")==0) {
			dt_col = tm_col = ii;
		} else if (tmp.compare("DATE")==0 || tmp.compare("GIORNO")==0) {
			dt_col = ii;
		} else if (tmp.compare("TIME")==0 || tmp.compare("ORA")==0) {
			tm_col = ii;
		}
	}

	//if necessary, set the format to the appropriate defaults
	if (dt_col==tm_col) {
		if (datetime_idx.empty())
			setDateTimeSpec("YYYY-MM-DDTHH24:MI:SS");
		if (fieldNames.size()==2 && !single_field.empty() && !user_provided_field_names)
			fieldNames[1] = single_field; //we have a better name from the filename
	} else {
		if (datetime_idx.empty())
			setDateTimeSpec("YYYY-MM-DD");
		if (time_idx.empty())
			setTimeSpec("HH24:MI:SS");
		if (fieldNames.size()==3 && !single_field.empty() && !user_provided_field_names)
			fieldNames[2] = single_field; //we have a better name from the filename
	}
}

//very basic units parsing: a few hard-coded units are recognized and provide the necessary
//offset and multiplier to convert the values back to SI
void CsvParameters::parseUnits(const std::string& line)
{
	static const std::string stdUnits[9] = {"TS", "RN", "W/M2", "M/S", "K", "M", "V", "VOLT", "DEG"};
	static std::set<std::string> noConvUnits( stdUnits, stdUnits+9 );
	
	std::vector<std::string> units;
	IOUtils::readLineToVec(line, units, csv_delim);
	units_offset.resize(units.size(), 0.);
	units_multiplier.resize(units.size(), 1.);
	
	for (size_t ii=0; ii<units.size(); ii++) {
		std::string tmp( units[ii] );
		IOUtils::toUpper( tmp );
		tmp.erase(std::remove_if(tmp.begin(), tmp.end(), &isQuote), tmp.end());
		if (tmp.empty()) continue; //empty unit
		if (noConvUnits.count(tmp)>0) continue; //this unit does not need conversion
		
		if (tmp=="%" || tmp=="CM") units_multiplier[ii] = 0.01;
		else if (tmp=="C") units_offset[ii] = Cst::t_water_freezing_pt;
		else if (tmp=="MM") units_multiplier[ii] = 1e-3;
		else if (tmp=="IN") units_multiplier[ii] = 0.0254;
		else if (tmp=="FT") units_multiplier[ii] = 0.3048;
		else if (tmp=="F") { units_multiplier[ii] = 5./9.; units_offset[ii] = -32.*5./9.;}
		else if (tmp=="KM/H") units_multiplier[ii] = 1./3.6;
		else if (tmp=="MPH") units_multiplier[ii] = 0.44704;
		else if (tmp=="KT") units_multiplier[ii] = 0.5144444444445;
		else {
			throw UnknownValueException("Can not parse unit '"+tmp+"', please inform the MeteoIO developers", AT);
		}
	}
}

//read and parse the file's headers in order to extract all possible information (including how to interpret the date/time information)
void CsvParameters::setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec, const std::string& filename_spec, const std::string& station_idx)
{
	file_and_path = i_file_and_path;
	const std::multimap< size_t, std::pair<size_t, std::string> > meta_spec( parseHeadersSpecs(vecMetaSpec) );
	double lat = IOUtils::nodata;
	double lon = IOUtils::nodata;
	
	//read and parse the file's headers
	if (!FileUtils::fileExists(file_and_path)) throw AccessException("File "+file_and_path+" does not exists", AT); //prevent invalid filenames
	errno = 0;
	if (!filename_spec.empty()) parseFileName( file_and_path, filename_spec, lat, lon);
	std::ifstream fin(file_and_path.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file " << file_and_path << " for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	
	const bool read_units = (units_headers!=IOUtils::npos && units_offset.empty() && units_multiplier.empty());
	size_t linenr=0;
	std::string line;
	std::vector<std::string> headerFields;
	try {
		eoln = FileUtils::getEoln(fin);
		for (size_t ii=0; ii<header_lines; ii++) {
			getline(fin, line, eoln); //read complete signature line
			if (fin.eof()) {
				std::ostringstream ss;
				ss << "Declaring " << header_lines << " header line(s) for file " << file_and_path << ", but it only contains " << linenr << " lines";
				throw InvalidArgumentException(ss.str(), AT);
			}
			linenr++;
			if (line.empty()) continue;
			if (*line.rbegin()=='\r') line.erase(line.end()-1); //getline() skipped \n, so \r comes in last position
			
			if (meta_spec.count(linenr)>0) 
				parseSpecialHeaders(line, linenr, meta_spec, lat, lon);
			if (linenr==columns_headers) //so user provided csv_fields have priority. If columns_headers==npos, this will also never be true
				IOUtils::readLineToVec(line, headerFields, csv_delim);
			if (read_units && linenr==units_headers)
				parseUnits(line);
		}
	} catch (...) {
		fin.close();
		throw;
	}
	fin.close();
	
	if (lat!=IOUtils::nodata || lon!=IOUtils::nodata) {
		const double alt = location.getAltitude(); //so we don't change previously set altitude
		location.setLatLon(lat, lon, alt); //we let Coords handle possible missing data / wrong values, etc
	}
	
	parseFields(headerFields, csv_fields, date_col, time_col);
	
	if (name.empty()) name = FileUtils::removeExtension( FileUtils::getFilename(i_file_and_path) ); //fallback if nothing else could be find
	if (id.empty()) {
		if (station_idx.empty()) 
			id = name; //really nothing, copy "name"
		else
			id = "ID"+station_idx; //automatic numbering of default IDs
	}
}

//check that the format is usable (and prevent parameters injection / buffer overflows)
void CsvParameters::checkSpecString(const std::string& spec_string, const size_t& nr_params)
{
	const size_t nr_percent = (unsigned)std::count(spec_string.begin(), spec_string.end(), '%');
	const size_t nr_placeholders0 = IOUtils::count(spec_string, "%d");
	const size_t nr_placeholders2 = IOUtils::count(spec_string, "%2d");
	const size_t nr_placeholders4 = IOUtils::count(spec_string, "%4d");
	const size_t nr_placeholders5 = IOUtils::count(spec_string, "%32s");
	size_t nr_placeholders = (nr_placeholders0!=std::string::npos)? nr_placeholders0 : 0;
	nr_placeholders += (nr_placeholders2!=std::string::npos)? nr_placeholders2 : 0;
	nr_placeholders += (nr_placeholders4!=std::string::npos)? nr_placeholders4 : 0;
	nr_placeholders += (nr_placeholders5!=std::string::npos)? nr_placeholders5 : 0;
	const size_t pos_pc_pc = spec_string.find("%%");
	if (nr_percent!=nr_params || nr_percent!=nr_placeholders || pos_pc_pc!=std::string::npos)
		throw InvalidFormatException("Badly formatted date/time specification '"+spec_string+"': argument appearing twice or using '%%'", AT);
}

//from a SPEC string such as "DD.MM.YYYY HH24:MIN:SS", build the format string for scanf as well as the parameters indices
//the indices are based on ISO timestamp, so year=0, month=1, ..., ss=5 and tz is handled separately
void CsvParameters::setDateTimeSpec(const std::string& datetime_spec)
{
	static const char* keys[6] = {"YYYY", "MM", "DD", "HH24", "MI", "SS"};
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (size_t ii=0; ii<6; ii++) {
		const size_t key_pos = datetime_spec.find( keys[ii] );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, ii) );
	}
	
	//fill datetime_idx as a vector of [0-5] indices (for ISO fields) in the order they appear in the user-provided format string
	std::sort(sorting_vector.begin(), sorting_vector.end(), &sort_dateKeys);
	for (size_t ii=0; ii<sorting_vector.size(); ii++)
		datetime_idx.push_back( sorting_vector[ii].second );
	
	datetime_format = datetime_spec;
	const size_t tz_pos = datetime_format.find("TZ");
	if (tz_pos!=std::string::npos) {
		if (tz_pos!=(datetime_format.length()-2))
			throw InvalidFormatException("When providing TZ in a date/time format, it must be at the very end of the string", AT);
		has_tz = true;
		datetime_format.replace(tz_pos, 2, "%32s");
	}
	IOUtils::replace_all(datetime_format, "DD", "%2d");
	IOUtils::replace_all(datetime_format, "MM", "%2d");
	IOUtils::replace_all(datetime_format, "YYYY", "%4d");
	IOUtils::replace_all(datetime_format, "HH24", "%2d");
	IOUtils::replace_all(datetime_format, "MI", "%2d");
	IOUtils::replace_all(datetime_format, "SS", "%d");
	
	const size_t nr_params_check = (has_tz)? datetime_idx.size()+1 : datetime_idx.size();
	checkSpecString(datetime_format, nr_params_check);
}

void CsvParameters::setTimeSpec(const std::string& time_spec)
{
	if (time_spec.empty()) return;
	static const char* keys[3] = {"HH24", "MI", "SS"};
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (size_t ii=0; ii<3; ii++) {
		const size_t key_pos = time_spec.find( keys[ii] );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, ii) );
	}

	//fill time_idx as a vector of [0-3] indices (for ISO fields) in the order they appear in the user-provided format string
	std::sort(sorting_vector.begin(), sorting_vector.end(), &sort_dateKeys);
	for (size_t ii=0; ii<sorting_vector.size(); ii++)
		time_idx.push_back( sorting_vector[ii].second );

	time_format = time_spec;
	const size_t tz_pos = time_format.find("TZ");
	if (tz_pos!=std::string::npos) {
		if (tz_pos!=(time_format.length()-2))
			throw InvalidFormatException("When providing TZ in a date/time format, it must be at the very end of the string", AT);
		has_tz = true;
		time_format.replace(tz_pos, 2, "%32s");
	}
	IOUtils::replace_all(time_format, "HH24", "%2d");
	IOUtils::replace_all(time_format, "MI", "%2d");
	IOUtils::replace_all(time_format, "SS", "%d");

	const size_t nr_params_check = (has_tz)? time_idx.size()+1 : time_idx.size();
	checkSpecString(time_format, nr_params_check);
}

Date CsvParameters::parseDate(const std::string& date_str, const std::string& time_str) const
{
	int args[6] = {0, 0, 0, 0, 0, 0};
	char rest[32] = "";
	
	bool status = false;
	switch( datetime_idx.size() ) {
		case 6:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], &args[ datetime_idx[5] ], rest)>=6);
			break;
		case 5:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], rest)>=5);
			break;
		case 4:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], rest)>=4);
			break;
		case 3:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], rest)>=3);
			break;
	}
	if (!status) return Date(); //we MUST have read successfuly at least the date part

	if (!time_idx.empty()) {
		//there is a +3 offset because the first 3 positions are used by the date part
		switch( time_idx.size() ) {
			case 3:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ], &args[ time_idx[1]+3 ], &args[ time_idx[2]+3 ], rest)>=3);
				break;
			case 2:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ], &args[ time_idx[1]+3 ], rest)>=2);
				break;
			case 1:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ], rest)>=1);
				break;
		}
	}

	if (!status) return Date();
	const double tz = (has_tz)? Date::parseTimeZone(rest) : csv_tz;
	return Date((int)args[0], (int)args[1], (int)args[2], (int)args[3], (int)args[4], (int)args[5], tz);
}

StationData CsvParameters::getStation() const 
{
	StationData sd(location, id, name);
	if (slope!=IOUtils::nodata && azi!=IOUtils::nodata)
		sd.setSlope(slope, azi);
	return sd;
}


///////////////////////////////////////////////////// Now the real CsvIO class starts //////////////////////////////////////////
const size_t CsvIO::streampos_every_n_lines = 2000; //save streampos every 2000 lines of data

CsvIO::CsvIO(const std::string& configfile) 
      : cfg(configfile), indexer_map(), csvparam(), vecStations(),
        coordin(), coordinparam(), silent_errors(false), errors_to_nodata(false) { parseInputOutputSection(); }

CsvIO::CsvIO(const Config& cfgreader)
      : cfg(cfgreader), indexer_map(), csvparam(), vecStations(),
        coordin(), coordinparam(), silent_errors(false), errors_to_nodata(false) { parseInputOutputSection(); }

void CsvIO::parseInputOutputSection()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	
	cfg.getValue("CSV_SILENT_ERRORS", "Input", silent_errors, IOUtils::nothrow);
	cfg.getValue("CSV_ERRORS_TO_NODATA", "Input", errors_to_nodata, IOUtils::nothrow);

	const double in_TZ = cfg.get("TIME_ZONE", "Input");
	const std::string meteopath = cfg.get("METEOPATH", "Input");
	std::vector< std::pair<std::string, std::string> > vecFilenames( cfg.getValues("STATION", "INPUT") );

	if (vecFilenames.empty()) { //scan all of the data path for a given file extension if no stations are specified
		bool is_recursive = false;
		std::string csvext(".csv");
		cfg.getValue("METEOPATH_RECURSIVE", "Input", is_recursive, IOUtils::nothrow);
		cfg.getValue("CSV_FILE_EXTENSION", "Input", csvext, IOUtils::nothrow); 
		
		std::list<std::string> dirlist( FileUtils::readDirectory(meteopath, csvext, is_recursive) );
		dirlist.sort();

		size_t hit = 0;	//human readable iterator
		for (std::list<std::string>::iterator it=dirlist.begin(); it!=dirlist.end(); ++it) {
			hit++;
			std::stringstream ss;
			ss << "STATION" << hit; //assign alphabetically ordered ID to station
			
			const std::pair<std::string, std::string> stat_id_and_name(ss.str(), *it);
			vecFilenames.push_back(stat_id_and_name);
		}
	} 
	
	for (size_t ii=0; ii<vecFilenames.size(); ii++) {
		const std::string idx( vecFilenames[ii].first.substr(string("STATION").length()) );
		static const std::string dflt("CSV_"); //the prefix for a key for ALL stations
		const std::string pre( "CSV"+idx+"_" ); //the prefix for the current station only
		
		CsvParameters tmp_csv(in_TZ);
		std::string coords_specs;
		if (cfg.keyExists("POSITION"+idx, "INPUT")) cfg.getValue("POSITION"+idx, "INPUT", coords_specs);
		else cfg.getValue("POSITION", "INPUT", coords_specs);
		const Coords loc(coordin, coordinparam, coords_specs);
		
		std::string name;
		if (cfg.keyExists(pre+"NAME", "Input")) cfg.getValue(pre+"NAME", "Input", name);
		else cfg.getValue(dflt+"NAME", "Input", name, IOUtils::nothrow);
		
		std::string id;
		if (cfg.keyExists(pre+"ID", "Input")) cfg.getValue(pre+"ID", "Input", id);
		else cfg.getValue(dflt+"ID", "Input", id, IOUtils::nothrow);
		tmp_csv.setLocation(loc, name, id);
		
		if (cfg.keyExists(pre+"NODATA", "Input")) cfg.getValue(pre+"NODATA", "Input", tmp_csv.nodata);
		else cfg.getValue(dflt+"NODATA", "Input", tmp_csv.nodata, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"DELIMITER", "Input")) cfg.getValue(pre+"DELIMITER", "Input", tmp_csv.csv_delim);
		else cfg.getValue(dflt+"DELIMITER", "Input", tmp_csv.csv_delim, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"NR_HEADERS", "Input")) cfg.getValue(pre+"NR_HEADERS", "Input", tmp_csv.header_lines);
		else cfg.getValue(dflt+"NR_HEADERS", "Input", tmp_csv.header_lines, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"COLUMNS_HEADERS", "Input")) cfg.getValue(pre+"COLUMNS_HEADERS", "Input", tmp_csv.columns_headers);
		else cfg.getValue(dflt+"COLUMNS_HEADERS", "Input", tmp_csv.columns_headers, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"FIELDS", "Input")) cfg.getValue(pre+"FIELDS", "Input", tmp_csv.csv_fields);
		else cfg.getValue(dflt+"FIELDS", "Input", tmp_csv.csv_fields, IOUtils::nothrow);
		
		if (tmp_csv.columns_headers==IOUtils::npos && tmp_csv.csv_fields.empty())
			throw InvalidArgumentException("Please provide either CSV_COLUMNS_HEADERS or CSV_FIELDS", AT);
		
		if (cfg.keyExists(pre+"UNITS_HEADERS", "Input")) cfg.getValue(pre+"UNITS_HEADERS", "Input", tmp_csv.units_headers);
		else cfg.getValue(dflt+"UNITS_HEADERS", "Input", tmp_csv.units_headers, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"UNITS_OFFSET", "Input")) cfg.getValue(pre+"UNITS_OFFSET", "Input", tmp_csv.units_offset);
		else cfg.getValue(dflt+"UNITS_OFFSET", "Input", tmp_csv.units_offset, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"UNITS_MULTIPLIER", "Input")) cfg.getValue(pre+"UNITS_MULTIPLIER", "Input", tmp_csv.units_multiplier);
		else cfg.getValue(dflt+"UNITS_MULTIPLIER", "Input", tmp_csv.units_multiplier, IOUtils::nothrow);
		
		//Date and time formats. The defaults will be set when parsing the column names (so they are appropriate for the available columns)
		std::string datetime_spec;
		if (cfg.keyExists(pre+"DATETIME_SPEC", "Input")) cfg.getValue(pre+"DATETIME_SPEC", "Input", datetime_spec);
		else cfg.getValue(dflt+"DATETIME_SPEC", "Input", datetime_spec, IOUtils::nothrow);
		
		std::string date_spec;
		if (cfg.keyExists(pre+"DATE_SPEC", "Input")) cfg.getValue(pre+"DATE_SPEC", "Input", date_spec);
		else cfg.getValue(dflt+"DATE_SPEC", "Input", date_spec, IOUtils::nothrow);
		
		std::string time_spec;
		if (cfg.keyExists(pre+"TIME_SPEC", "Input")) cfg.getValue(pre+"TIME_SPEC", "Input", time_spec);
		else cfg.getValue(dflt+"TIME_SPEC", "Input", time_spec, IOUtils::nothrow);
		
		if (!datetime_spec.empty())
			tmp_csv.setDateTimeSpec(datetime_spec);
		else {
			if (!date_spec.empty())
				tmp_csv.setDateTimeSpec(date_spec);
			if (!time_spec.empty()) 
				tmp_csv.setTimeSpec(time_spec);
		}
		
		std::vector<std::string> vecMetaSpec;
		if (cfg.keyExists(pre+"SPECIAL_HEADERS", "Input")) cfg.getValue(pre+"SPECIAL_HEADERS", "Input", vecMetaSpec);
		else cfg.getValue(dflt+"SPECIAL_HEADERS", "Input", vecMetaSpec, IOUtils::nothrow);
		
		std::string filename_spec;
		if (cfg.keyExists(pre+"FILENAME_SPEC", "Input")) cfg.getValue(pre+"FILENAME_SPEC", "Input", filename_spec);
		else cfg.getValue(dflt+"FILENAME_SPEC", "Input", filename_spec, IOUtils::nothrow);
		
		tmp_csv.setFile(meteopath + "/" + vecFilenames[ii].second, vecMetaSpec, filename_spec, idx);
		csvparam.push_back( tmp_csv );
	}
}

void CsvIO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	for (size_t ii=0; ii<csvparam.size(); ii++)
		vecStation.push_back( csvparam[ii].getStation() );
}

std::vector<MeteoData> CsvIO::readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd)
{
	size_t nr_of_data_fields = params.csv_fields.size(); //this has been checked by CsvParameters
	const bool use_offset = !params.units_offset.empty();
	const bool use_multiplier = !params.units_multiplier.empty();
	if ((use_offset && params.units_offset.size()!=nr_of_data_fields) || (use_multiplier && params.units_multiplier.size()!=nr_of_data_fields))
		throw InvalidFormatException("The declared units_offset / units_multiplier must match the number of columns in the file!", AT);

	//build MeteoData template
	MeteoData template_md( Date(0., 0.), params.getStation() );
	for (size_t ii=0; ii<nr_of_data_fields; ii++)
		template_md.addParameter( params.csv_fields[ii] );

	//now open the file
	const std::string filename( params.getFilename() );
	if (!FileUtils::fileExists(filename)) throw AccessException("File '"+filename+"' does not exists", AT); //prevent invalid filenames
	errno = 0;
	std::ifstream fin(filename.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << filename << "\" for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	
	std::string line;
	size_t linenr=0;
	streampos fpointer = indexer_map[filename].getIndex(dateStart);
	if (fpointer!=static_cast<streampos>(-1))
		fin.seekg(fpointer); //a previous pointer was found, jump to it
	else {
		//skip the headers (they have been read already, so we know this works)
		while (!fin.eof() && linenr<params.header_lines) {
			line.clear();
			getline(fin, line, params.eoln);
			linenr++;
		}
	}
	
	//and now, read the data and fill the vector vecMeteo
	std::vector<MeteoData> vecMeteo;
	std::vector<std::string> tmp_vec;
	const std::string nodata( params.nodata );
	const std::string nodata_with_quotes( "\""+params.nodata+"\"" );
	const std::string nodata_with_single_quotes( "\'"+params.nodata+"\'" );
	Date prev_dt;
	size_t count_asc=0, count_dsc=0; //count how many ascending/descending timestamps are present
	while (!fin.eof()){
		getline(fin, line, params.eoln);
		linenr++;
		IOUtils::trim( line );
		if (line.empty()) continue; //Pure comment lines and empty lines are ignored
		
		const size_t nr_curr_data_fields = IOUtils::readLineToVec(line, tmp_vec, params.csv_delim);
		if (nr_of_data_fields==0) nr_of_data_fields = nr_curr_data_fields;
		if (nr_curr_data_fields!=nr_of_data_fields) {
			std::ostringstream ss;
			ss << "File \'" << filename << "\' declares (either as first data line or columns headers or units offset/multiplier) " << nr_of_data_fields << " columns ";
			ss << "but this does not match line " << linenr << ":\n'" << line << "'\n";
			if (silent_errors) {
				std::cerr << ss.str();
				continue;
			} else throw InvalidFormatException(ss.str(), AT);
		}
		
		const Date dt( params.parseDate(tmp_vec[params.date_col], tmp_vec[params.time_col]) );
		if (dt.isUndef()) {
			const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
			const std::string err_msg( "Date or time could not be read in file \'"+filename+"' at line "+linenr_str );
			if (silent_errors) {
				std::cerr << err_msg << "\n";
				continue;
			} else throw InvalidFormatException(err_msg, AT);
		}
		if (!prev_dt.isUndef()) {
			const bool asc = (dt>prev_dt);
			if (asc) count_asc++;
			else count_dsc++;
		}
		prev_dt = dt;

		if (linenr % streampos_every_n_lines == 0) {
			fpointer = fin.tellg();
			if (fpointer != static_cast<streampos>(-1)) indexer_map[filename].setIndex(dt, fpointer);
		}
		if (dt<dateStart) continue;
		if (dt>dateEnd) break;
		
		MeteoData md(template_md);
		md.setDate(dt);
		bool no_errors = true;
		for (size_t ii=0; ii<tmp_vec.size(); ii++){
			if (ii==params.date_col || ii==params.time_col) continue;
			if (tmp_vec[ii].empty() || tmp_vec[ii]==nodata || tmp_vec[ii]==nodata_with_quotes || tmp_vec[ii]==nodata_with_single_quotes) //treat empty value as nodata, try nodata marker w/o quotes
				continue;
			
			double tmp;
			if (!IOUtils::convertString(tmp, tmp_vec[ii])) {
				const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
				const std::string err_msg( "Could not parse field '"+tmp_vec[ii]+"' in file \'"+filename+"' at line "+linenr_str );
				if (silent_errors) {
					std::cerr << err_msg << "\n";
					no_errors = false;
					continue;
				} else if (errors_to_nodata) {
					tmp = IOUtils::nodata;
				} else throw InvalidFormatException(err_msg, AT);
			}
			if (use_multiplier && tmp!=IOUtils::nodata) tmp *= params.units_multiplier[ii];
			if (use_offset && tmp!=IOUtils::nodata) tmp += params.units_offset[ii];
			md( params.csv_fields[ii] ) = tmp;
		}
		if (no_errors) vecMeteo.push_back( md );
	}
	
	if (count_dsc>count_asc) std::reverse(vecMeteo.begin(), vecMeteo.end()); //since we might have DST, we might have localy asc/dsc order...
	
	return vecMeteo;
}

void CsvIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	for (size_t ii=0; ii<csvparam.size(); ii++) {
		vecMeteo.push_back( readCSVFile(csvparam[ii], dateStart, dateEnd) );
	}
}

} //namespace
