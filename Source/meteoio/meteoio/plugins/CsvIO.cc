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
 *     - a single character is consistently used through the file as field delimiter (to split each record into fields);
 *     - each line contains the same number of fields;
 *     - missing data are represented by an empty value, so two delimiters follow directly each other;
 *     - the file may contain a header that may contain additional information (metadata), see below.
 * 
 * In order to reduce the amount of manual configuration, it is possible to extract metadata from the headers, such as the station name, ID, coordinates, etc
 *
 * @section csvio_units Units
 * The final units MUST be SI. If not, the conversion offsets/factors must be provided to convert the data back to SI (see required keywords below).
 *
 * @section csvio_keywords Keywords
 * This plugin uses the following keywords, in the [Input] section:
 * - COORDSYS: coordinate system (see Coords);
 * - COORDPARAM: extra coordinates parameters (see Coords);
 * - TIME_ZONE: the timezone that should be used to interpret the dates/times (default: 0);
 * - METEOPATH: the directory where the data files are available (mandatory);
 * - STATION\#: input filename (in METEOPATH). As many meteofiles as needed may be specified;
 * - POSITION\#: coordinates of the station (default: reading key "POSITION"; see (see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax));
 * 
 * The following keys may either be prefixed by "CSV_" (ie as default for all stations) or by "CSV#_" (as only for the current station):
 * - CSV\#_DELIMITER: field delimiter to use (default: ',');
 * - CSV\#_NR_HEADERS: how many lines should be treated as headers? (default: 1);
 * - CSV\#_COLUMNS_HEADERS: header line to interpret as columns headers (default: 1);
 * - CSV\#_FIELDS: columns headers (if they don't exist in the file or to overwrite them); optional
 * - CSV\#_UNITS_OFFSET: offset to add to each value in order to convert it to SI; optional
 * - CSV\#_UNITS_MULTIPLIER: factor to multiply each value by, in order to convert it to SI; optional
 * - CSV\#_DATETIME_SPEC: mixed date and time format specification (defaultis ISO_8601: YYYY-MM-DDTHH24:MI:SS);
 * - CSV\#_DATE_SPEC: date format specification (default: YYYY_MM_DD);
 * - CSV\#_TIME_SPEC: time format specification (default: HH24:MI:SS);
 * - CSV\#_SPECIAL_HEADERS: description of how to extract more metadata out of the headers; optional
 * 
 * @section csvio_date_specs Date and time specification
 * In order to be able to read any date and time format, the format has to be provided in the configuration file. This is provided as a string containing
 * the following special markers:
 * - YYYY: the 4 digits year;
 * - MM: the two digits month;
 * - DD: the two digits day;
 * - HH24: the two digits hour of the day (0-24);
 * - MI: the two digits minutes (0-59);
 * - SS: the two digts seconds (0-59).
 *
 * Any other character is interpreted as itself, present in the string. It is possible to either provide a combined datetime field (so date and time are combined into
 * one single field) or date and time as two different fields. For example:
 * - YYYY-MM-DDTHH24:MI:SS described an <A HREF="https://en.wikipedia.org/wiki/ISO_8601">ISO 8601</A> datetime field;
 * - MM/DD/YYYY described an anglo-saxon date;
 * - DD.MM.YYYY HH24:MI:SS is for a Swiss formatted datetime.
 *
 * @section csvio_special_headers Header metadata extraction
 * Some valuable metadata might be written into the headers and it is worthwhile to be able to extract it. This is performed with the
 * "CSV#_SPECIAL_HEADERS" configuration key. This key is followed by as many metadata specifications as necessary, of the form {field}:{line}:{column}.
 * The {field} can be any of:
 * - name;
 * - id;
 * - alt (for the altitude);
 * - lon (for the longitude);
 * - lat (for the latitude);
 * - slope (in degrees);
 * - azi (for the slope azimuth, in degree as read from a compass).
 *
 * Therefore, if the station name is available on line 1, column 3 and the station id on line 2, column 5, the configuration would be:
 * @code
 * CSV_SPECIAL_HEADERS = name:1:3 id:2:5
 * @endcode
 * 
 * @section csvio_examples Examples
 * In order to read a CSV file produced by a Campbell data logger with Swiss-formatted timestamps, you need the following configuration:
 * @code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_NR_HEADERS = 4
 * CSV_COLUMNS_HEADERS = 2
 * CSV_DATETIME_SPEC = DD.MM.YYYY HH24:MI:SS
 * STATION1 = DisMa_DisEx.dat
 * POSITION1 = latlon 46.810325 9.806657 2060
 * CSV_SPECIAL_HEADERS = name:1:2 id:1:4
 * CSV_UNITS_MULTIPLIER = 1 1 1 0.01 1 1 1 0.01 1
 * @endcode
 *
 * In order to read a set of files containing each only one parameter and merge them together (see \ref data_manipulations "raw data editing" for more
 * on the merge feature):
 *@code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_DELIMITER = ;
 * CSV_HEADER_LINES = 1
 * CSV_DATE_SPEC = DD/MM/YYYY
 * CSV_TIME_SPEC = HH24:MI
 * POSITION = latlon (46.8, 9.80, 1700)
 *
 * CSV1_FIELDS = DATE TIME PSUM
 * STATION1 = H0118_precipitation.DAT
 *
 * CSV2_FIELDS = DATE TIME TA
 * STATION2 = H0118_tempeartures.DAT
 *
 * CSV3_FIELDS = DATE TIME RSWR
 * STATION3 = H0118_reflected_solar_radiation.DAT
 *
 * CSV4_FIELDS = DATE TIME RH
 * STATION4 = H0118_relaitve_humidity.DAT
 *
 * CSV5_FIELDS = DATE TIME VW
 * STATION5 = H0118_wind_velocity.DAT
 *
 * ID1::MERGE = ID2 ID3 ID4 ID5
 * @endcode
 */

//parse the user provided special fields specification
std::multimap< size_t, std::pair<size_t, std::string> > CsvParameters::parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec) const
{
	std::multimap< size_t, std::pair<size_t, std::string> > meta_spec;
	for (size_t ii=0; ii<vecMetaSpec.size(); ii++) {
		std::vector<std::string> vecArgs;
		if (IOUtils::readLineToVec(vecMetaSpec[ii], vecArgs, ':') !=3)
			throw InvalidFormatException("Wrong format for Metadata specification '"+vecMetaSpec[ii]+"'", AT);
		const int linenr = atoi( vecArgs[1].c_str() );
		const int colnr = atoi( vecArgs[2].c_str() );
		if (linenr<=0 || colnr<=0)
			throw InvalidFormatException("Line numbers and column number must be >0 in Metadata specification '"+vecMetaSpec[ii]+"'", AT);
		
		meta_spec.insert( make_pair( linenr, make_pair( colnr, vecArgs[0]) ) );
	}
	
	return meta_spec;
}

inline bool isQuote(const char& c) { if (c=='"' || c=='\'') return true; return false;}

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
		
		if (field_type=="NAME") {
			name = field_val;
		} else if (field_type=="ID") {
			id = field_val;
		} else {
			double tmp;
			if (!IOUtils::convertString(tmp, field_val))
				throw InvalidArgumentException("Could not parse Metadata specification '"+field_type+"' for "+file_and_path, AT);
			
			if (field_type=="ALT") location.setAltitude( tmp, false);
			if (field_type=="LON") lon=tmp;
			if (field_type=="LAT") lat=tmp;
			if (field_type=="SLOPE") slope=tmp;
			if (field_type=="AZI") azi=tmp;
		}
	}
}

void CsvParameters::parseFields(std::vector<std::string>& fieldNames, size_t &dt_col, size_t &tm_col)
{
	for (size_t ii=0; ii<fieldNames.size(); ii++) {
		std::string &tmp = fieldNames[ii];
		IOUtils::toUpper( tmp );
		tmp.erase(std::remove_if(tmp.begin(), tmp.end(), &isQuote), tmp.end());
		if (tmp.empty()) continue;
		
		if (tmp.compare("TIMESTAMP")==0) {
			dt_col = tm_col = ii;
		} else if (tmp.compare("DATE")==0) {
			dt_col = ii;
		} else if (tmp.compare("TIME")==0) {
			tm_col = ii;
		}
	}
	
	//if necessary, set the format to the appropriate defaults
	if (dt_col==tm_col) {
		if (datetime_idx.empty())
			setDateTimeSpec("YYYY-MM-DDTHH24:MI:SS");
	} else {
		if (datetime_idx.empty())
			setDateTimeSpec("YYYY-MM-DD");
		if (time_idx.empty())
			setTimeSpec("HH24:MI:SS");
	}
}

//read and parse the file's headers in order to extract all possible information
void CsvParameters::setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec)
{
	file_and_path = i_file_and_path;
	const std::multimap< size_t, std::pair<size_t, std::string> > meta_spec( parseHeadersSpecs(vecMetaSpec) );
	
	//read and parse the file's headers
	if (!FileUtils::fileExists(file_and_path)) throw AccessException("File "+file_and_path+" does not exists", AT); //prevent invalid filenames
	errno = 0;
	std::ifstream fin(file_and_path.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file " << file_and_path << " for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	
	size_t linenr=0;
	std::string line;
	double lat = IOUtils::nodata;
	double lon = IOUtils::nodata;
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
			
			if (meta_spec.count(linenr)>0) 
				parseSpecialHeaders(line, linenr, meta_spec, lat, lon);
			if (linenr==columns_headers && csv_fields.empty()) //so user provided csv_fields have priority. If columns_headers==npos, this will also never be true
				IOUtils::readLineToVec(line, csv_fields, csv_delim);
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
	
	//cleanup potential '\r' char at the end of the line
	if (csv_fields.empty())
		throw InvalidArgumentException("No columns names could be retrieved, please provide them through the configuration file", AT);
	std::string &tmp = csv_fields.back();
	if (*tmp.rbegin()=='\r') tmp.erase(tmp.end()-1); //getline() skipped \n, so \r comes in last position
	
	parseFields(csv_fields, date_col, time_col);
}

struct sort_pred {
	bool operator()(const std::pair<size_t,size_t> &left, const std::pair<size_t,size_t> &right) {
		return left.first < right.first;
	}
};
	
//from a SPEC string such as "DD.MM.YYYY HH24:MIN:SS", build the format string for scanf as well as the parameters indices
//the indices are based on ISO timestamp, so year=0, month=1, etc
void CsvParameters::setDateTimeSpec(const std::string& datetime_spec)
{
	static const char* keys[] = {"YYYY", "MM", "DD", "HH24", "MI", "SS"};
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (size_t ii=0; ii<6; ii++) {
		const size_t key_pos = datetime_spec.find( keys[ii] );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, ii) );
	}
	
	std::sort(sorting_vector.begin(), sorting_vector.end(), sort_pred());
	for (size_t ii=0; ii<sorting_vector.size(); ii++)
		datetime_idx.push_back( sorting_vector[ii].second );
	
	datetime_format = datetime_spec;
	IOUtils::replace_all(datetime_format, "DD", "%d");
	IOUtils::replace_all(datetime_format, "MM", "%d");
	IOUtils::replace_all(datetime_format, "YYYY", "%d");
	IOUtils::replace_all(datetime_format, "HH24", "%d");
	IOUtils::replace_all(datetime_format, "MI", "%d");
	IOUtils::replace_all(datetime_format, "SS", "%d");
	
	//check that the format is usable (and prevent parameters injection / buffer overflows)
	const size_t nr_percent = (unsigned)std::count(datetime_format.begin(), datetime_format.end(), '%');
	const size_t nr_placeholders = IOUtils::count(datetime_format, "%d");
	const size_t pos_pc_pc = datetime_format.find("%%");
	if (nr_percent!=datetime_idx.size() || nr_percent!=nr_placeholders || pos_pc_pc!=std::string::npos)
		throw InvalidFormatException("Badly formatted date/time specification '"+datetime_spec+"': argument appearing twice or using '%%'", AT);
}

void CsvParameters::setTimeSpec(const std::string& time_spec)
{
	if (time_spec.empty()) return;
	static const char* keys[] = {"HH24", "MI", "SS"};
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (size_t ii=0; ii<3; ii++) {
		const size_t key_pos = time_spec.find( keys[ii] );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, ii) );
	}

	std::sort(sorting_vector.begin(), sorting_vector.end(), sort_pred());
	for (size_t ii=0; ii<sorting_vector.size(); ii++)
		time_idx.push_back( sorting_vector[ii].second );

	time_format = time_spec;
	IOUtils::replace_all(time_format, "HH24", "%d");
	IOUtils::replace_all(time_format, "MI", "%d");
	IOUtils::replace_all(time_format, "SS", "%d");

	//check that the format is usable (and prevent parameters injection / buffer overflows)
	const size_t nr_percent = (unsigned)std::count(time_format.begin(), time_format.end(), '%');
	const size_t nr_placeholders = IOUtils::count(time_format, "%d");
	const size_t pos_pc_pc = time_format.find("%%");
	if (nr_percent!=time_idx.size() || nr_percent!=nr_placeholders || pos_pc_pc!=std::string::npos)
		throw InvalidFormatException("Badly formatted time specification '"+time_format+"': argument appearing twice or using '%%'", AT);
}

Date CsvParameters::parseDate(const std::string& date_str, const std::string& time_str) const
{
	int args[6] = {0, 0, 0, 0, 0 ,0};

	bool status = false;
	switch( datetime_idx.size() ) {
		case 6:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], &args[ datetime_idx[5] ])==6);
			break;
		case 5:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ]) ==5);
			break;
		case 4:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ])==4);
			break;
		case 3:
			status = (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ])==3);
			break;
	}
	if (!status) return Date(); //we MUST have read successfuly at least the date part

	if (!time_idx.empty()) {
		//there is a +3 offset because the first 3 positions are used by the date part
		switch( time_idx.size() ) {
			case 3:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ], &args[ time_idx[1]+3 ], &args[ time_idx[2]+3 ])==3);
				break;
			case 2:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ], &args[ time_idx[1]+3 ])==2);
				break;
			case 1:
				status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0]+3 ])==1);
				break;
		}
	}

	if (!status) return Date();
	return Date(args[0], args[1], args[2], args[3], args[4], args[5], csv_tz);
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
      : cfg(configfile), indexer(), csvparam(), vecStations(),
        coordin(), coordinparam() { parseInputOutputSection(); }

CsvIO::CsvIO(const Config& cfgreader)
      : cfg(cfgreader), indexer(), csvparam(), vecStations(),
        coordin(), coordinparam() { parseInputOutputSection(); }

void CsvIO::parseInputOutputSection()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	
	const double in_TZ = cfg.get("TIME_ZONE", "Input");
	const std::string meteopath = cfg.get("METEOPATH", "Input");
	const std::vector< std::pair<std::string, std::string> > vecFilenames( cfg.getValues("STATION", "INPUT") );
	
	for (size_t ii=0; ii<vecFilenames.size(); ii++) {
		const std::string idx( vecFilenames[ii].first.substr(string("STATION").length()) );
		static const std::string dflt("CSV_"); //the prefix for a key for ALL stations
		const std::string pre( "CSV"+idx+"_" ); //the prefix for the current station only
		
		CsvParameters tmp_csv(in_TZ);
		std::string coords_specs;
		if (cfg.keyExists("POSITION"+idx, "INPUT")) cfg.getValue("POSITION"+idx, "INPUT", coords_specs);
		else cfg.getValue("POSITION", "INPUT", coords_specs);
		const Coords loc(coordin, coordinparam, coords_specs);
		const std::string name( FileUtils::removeExtension(vecFilenames[ii].second) );
		tmp_csv.setLocation(loc, name, "ID"+idx);
		
		
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
		
		tmp_csv.setFile(meteopath + "/" + vecFilenames[ii].second, vecMetaSpec);
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
	streampos fpointer = indexer.getIndex(dateStart);

	if (fpointer!=static_cast<streampos>(-1))
		fin.seekg(fpointer); //a previous pointer was found, jump to it
	else {
		//skip the headers (they have been read already, so we know this works)
		while (!fin.eof() && linenr<params.header_lines){
			line.clear();
			getline(fin, line, params.eoln);
			linenr++;
		}
	}
	
	//and now, read the data and fill the vector vecMeteo
	std::vector<MeteoData> vecMeteo;
	std::vector<std::string> tmp_vec;
	while (!fin.eof()){
		const streampos current_fpointer = fin.tellg();
		getline(fin, line, params.eoln);
		linenr++;
		if (line.empty()) continue; //Pure comment lines and empty lines are ignored
		
		const size_t nr_curr_data_fields = IOUtils::readLineToVec(line, tmp_vec, params.csv_delim);
		if (nr_of_data_fields==0) nr_of_data_fields = nr_curr_data_fields;
		if (nr_curr_data_fields!=nr_of_data_fields) {
			std::ostringstream ss;
			ss << "File \'" << filename << "\' declares (either as first data line or columns headers or units offset/multiplier) " << nr_of_data_fields << " columns ";
			ss << "but this does not match the following line:\n" << line << "\n";
			throw InvalidFormatException(ss.str(), AT);
		}
		
		const Date dt( params.parseDate(tmp_vec[params.date_col], tmp_vec[params.time_col]) );
		if (dt.isUndef()) {
			const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
			throw InvalidFormatException("Date or time could not be read in file \'"+filename+"' at line "+linenr_str, AT);
		}

		if ( (linenr % streampos_every_n_lines)==0 && (current_fpointer != static_cast<streampos>(-1)) )
			indexer.setIndex(dt, current_fpointer);
		if (dt<dateStart) continue;
		if (dt>dateEnd) break;
		
		MeteoData md(template_md);
		md.setDate(dt);
		for (size_t ii=0; ii<tmp_vec.size(); ii++){
			if (ii==params.date_col || ii==params.time_col) continue;
			if (tmp_vec[ii].empty()) { //treat empty value as nodata
				md( params.csv_fields[ii] ) = IOUtils::nodata;
				continue;
			}
			double tmp;
			if (!IOUtils::convertString(tmp, tmp_vec[ii])) {
				const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
				throw InvalidFormatException("Could not parse field '"+tmp_vec[ii]+"' in file \'"+filename+"' at line "+linenr_str, AT);
			}
			if (use_multiplier) tmp *= params.units_multiplier[ii];
			if (use_offset) tmp += params.units_offset[ii];
			md( params.csv_fields[ii] ) = tmp;
		}
		vecMeteo.push_back( md );
	}
	
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
