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
 * *Put here the informations about the standard format that is implemented*
 *
 * @section csvio_units Units
 *
 *
 * @section csvio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

//parse the user provided special fields specification
std::map< size_t, std::pair<size_t, std::string> > CsvParameters::parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec) const
{
	std::map< size_t, std::pair<size_t, std::string> > meta_spec;
	for (size_t ii=0; ii<vecMetaSpec.size(); ii++) {
		bool error = false;
		const size_t count_delims = std::count(vecMetaSpec[ii].begin(), vecMetaSpec[ii].end(), ':');
		if (count_delims==2) {
			const size_t delim1 = vecMetaSpec[ii].find(":");
			const size_t delim2 = vecMetaSpec[ii].find(":", delim1+1);
			if (delim2!=delim1 && delim1>0) {
				const int linenr = atoi( vecMetaSpec[ii].substr(delim1+1, delim2-delim1-1).c_str() );
				const int colnr = atoi( vecMetaSpec[ii].substr(delim2+1).c_str() );
				if (linenr>0 && colnr>0)
					meta_spec[ linenr ] = make_pair( colnr, vecMetaSpec[ii].substr(0, delim1));
				else
					error = true;
			} else {
				error = true;
			}
		} else {
			error = true;
		}
		
		if (error)
			throw InvalidFormatException("Wrong format for Metadata specification '"+vecMetaSpec[ii]+"'", AT);
	}
	
	return meta_spec;
}

inline bool isQuote(const char& c) { if (c=='"' || c=='\'') return true; return false;}

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
}

//read and parse the file's headers in order to extract all possible information
void CsvParameters::setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec)
{
	//HACK: support multiple fields per line, ie use multimap
	std::map< size_t, std::pair<size_t, std::string> > meta_spec( parseHeadersSpecs(vecMetaSpec) );
	
	//read and parse the file's headers
	file_and_path = i_file_and_path;
	if (!FileUtils::fileExists(file_and_path)) throw AccessException("File '"+file_and_path+"' does not exists", AT); //prevent invalid filenames
	errno = 0;
	std::ifstream fin(file_and_path.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << file_and_path << "\" for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	eoln = FileUtils::getEoln(fin);
	
	size_t linenr=0;
	std::string line;
	std::vector<std::string> vecStr;
	for (size_t ii=0; ii<header_lines; ii++) {
		getline(fin, line, eoln); //read complete signature line
		linenr++;
		if (meta_spec.count(linenr)>0) {
			IOUtils::readLineToVec(line, vecStr, csv_delim);
			if (meta_spec[linenr].first>vecStr.size())
				throw InvalidArgumentException("Metadata specification for '"+meta_spec[linenr].second+"' refers to a non-existent field for file '"+file_and_path+"'", AT);
			const std::string field_type( IOUtils::strToUpper(meta_spec[linenr].second) );
			//for the user, the first field is 1, not 0 -> "-1"
			if (field_type=="NAME") name=vecStr[meta_spec[linenr].first-1]; 
			if (field_type=="ID") id=vecStr[meta_spec[linenr].first-1];
		}
		if (linenr==columns_headers) {
			IOUtils::readLineToVec(line, csv_fields, csv_delim);
		}
	}
	fin.close(); //HACK try/catch to close the file even after a throw
	
	if (!csv_fields.empty()) { //cleanup potential '\r' char at the end of the line
		std::string &tmp = csv_fields.back();
		if (*tmp.rbegin()=='\r') tmp.erase(tmp.end()-1); //BUG this should be \r\n ?! so \r second from last...
	}
	
	parseFields(csv_fields, date_col, time_col);
}

struct sort_pred {
	bool operator()(const std::pair<size_t,size_t> &left, const std::pair<size_t,size_t> &right) {
		return left.first < right.first;
	}
};
	
//from a SPEC string such as "DD.MM.YYYY HH24:MIN:SS", build the format string for scanf as well as the parameters indices
//the indices are based on ISO timestamp, so year=0, month=1, etc
void CsvParameters::setDateTimeSpec(const std::string& datetime_spec, const double& tz_in)
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
	IOUtils::replace_all(datetime_format, "DD", "%u");
	IOUtils::replace_all(datetime_format, "MM", "%u");
	IOUtils::replace_all(datetime_format, "YYYY", "%u");
	IOUtils::replace_all(datetime_format, "HH24", "%u");
	IOUtils::replace_all(datetime_format, "MI", "%u");
	IOUtils::replace_all(datetime_format, "SS", "%u");
	
	csv_tz = tz_in;
}

Date CsvParameters::parseDate(const std::string& date_str, const std::string& /*time_str*/) const
{
	const size_t nrArgs = datetime_idx.size();
	unsigned int args[6];
	
	Date t;
	if (nrArgs==6) {
		if (sscanf(date_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], &args[ datetime_idx[5] ]) == 6) {
			t.setDate(args[0], args[1], args[2], args[3], args[4], args[5], csv_tz);
		}
	}
	
	return t;
}


///////////////////////////////////////////////////// Now the real CsvIO class starts //////////////////////////////////////////
CsvIO::CsvIO(const std::string& configfile) 
      : cfg(configfile), csvparam(), vecStations(),
        coordin(), coordinparam(), coordout(), coordoutparam()
{
	parseInputOutputSection();
}

CsvIO::CsvIO(const Config& cfgreader)
      : cfg(cfgreader), csvparam(), vecStations(),
        coordin(), coordinparam(), coordout(), coordoutparam()
{
	parseInputOutputSection();
}

void CsvIO::parseInputOutputSection()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	const double in_TZ = cfg.get("TIME_ZONE", "Input");
	std::string meteopath;
	cfg.getValue("METEOPATH", "Input", meteopath);
	std::vector<std::string> vecFilenames;
	cfg.getValues("STATION", "INPUT", vecFilenames);
	const std::vector< std::pair<std::string, std::string> > vecCoords( cfg.getValues("META", "INPUT") );
	if (vecFilenames.size()!=vecCoords.size())
		throw InvalidArgumentException("The declared stations and metadata must match!", AT);
	
	csvparam.resize(vecFilenames.size());
	for (size_t ii=0; ii<vecFilenames.size(); ii++) {
		csvparam[ii].location.setLatLon(vecCoords[ii].second, IOUtils::nodata);
		
		//HACK several of these should be per station!
		cfg.getValue("CSV_DELIMITER", "Input", csvparam[ii].csv_delim, IOUtils::nothrow);
		cfg.getValue("CSV_HEADER_LINES", "Input", csvparam[ii].header_lines, IOUtils::nothrow);
		cfg.getValue("CSV_COLUMNS_HEADERS", "Input", csvparam[ii].columns_headers, IOUtils::nothrow);
		cfg.getValue("CSV_FIELDS", "Input", csvparam[ii].csv_fields, IOUtils::nothrow);
		cfg.getValue("CSV_UNITS_OFFSET", "Input", csvparam[ii].units_offset, IOUtils::nothrow);
		cfg.getValue("CSV_UNITS_MULTIPLIER", "Input", csvparam[ii].units_multiplier, IOUtils::nothrow);
		
		std::string datetime_spec;
		cfg.getValue("CSV_DATETIME_SPEC", "Input", datetime_spec, IOUtils::nothrow);
		csvparam[ii].setDateTimeSpec(datetime_spec, in_TZ);
		
		std::vector<std::string> vecMetaSpec;
		cfg.getValue("CSV_SPECIAL_HEADERS", "Input", vecMetaSpec, IOUtils::nothrow);
		csvparam[ii].setFile(meteopath + "/" + vecFilenames[ii], vecMetaSpec);
	}
	
}

void CsvIO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	for (size_t ii=0; ii<csvparam.size(); ii++)
		vecStation.push_back( csvparam[ii].getStation() );
}

std::vector<MeteoData> CsvIO::readCSVFile(CsvParameters& params, const Date& /*dateStart*/, const Date& /*dateEnd*/)
{
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
	
	//skip the headers
	while (!fin.eof() && linenr<params.header_lines){
		line.clear();
		getline(fin, line, params.eoln);
		linenr++;
	}
	
	std::vector<MeteoData> vecMeteo;
	std::vector<std::string> tmp_vec;
	if (params.csv_fields.empty())
		throw InvalidArgumentException("No columns names could be retrieve, please provide them through the configuration file", AT);
	
	size_t nr_of_data_fields = params.csv_fields.size();
	const bool use_offset = !params.units_offset.empty();
	const bool use_multiplier = !params.units_multiplier.empty();
	if ((use_offset && params.units_offset.size()!=nr_of_data_fields) || (use_multiplier && params.units_multiplier.size()!=nr_of_data_fields)) {
		throw InvalidFormatException("The declared units_offset / units_multiplier must match the number of columns in the file!", AT);
	}
	
	StationData sd(params.location, params.id, params.name);
	while (!fin.eof()){
		line.clear();
		getline(fin, line, params.eoln);
		linenr++;
		
		if (line.empty()) continue; //Pure comment lines and empty lines are ignored
		
		const size_t nr_curr_data_fields = IOUtils::readLineToVec(line, tmp_vec, params.csv_delim);
		if (nr_of_data_fields==0) nr_of_data_fields=nr_curr_data_fields;
		if (nr_curr_data_fields!=nr_of_data_fields) {
			std::ostringstream ss;
			ss << "File \'" << filename << "\' declares (either as first data line or columns headers or units offset/multiplier) " << nr_of_data_fields << " columns ";
			ss << "but this does not match the following line:\n" << line << "\n";
			throw InvalidFormatException(ss.str(), AT);
		}
		
		const Date dt( params.parseDate(tmp_vec[params.date_col], tmp_vec[params.time_col]) );
		if (dt.isUndef()) {
			const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
			throw InvalidFormatException("Date could not be read in file \'"+filename+"' at line "+linenr_str, AT);
		}
		
		//MeteoData md( dt, vecStations[stat_idx]);
		MeteoData md( dt, sd);
		for (size_t ii=0; ii<tmp_vec.size(); ii++){
			if (ii==params.date_col) continue;
			double tmp;
			if (!IOUtils::convertString(tmp, tmp_vec[ii])) {
				const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
				throw InvalidFormatException("Could not parse field '"+tmp_vec[ii]+"' in file \'"+filename+"' at line "+linenr_str, AT);
			}
			if (use_multiplier) tmp *= params.units_multiplier[ii];
			if (use_offset) tmp += params.units_offset[ii];
			md.addParameter( params.csv_fields[ii] );
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
