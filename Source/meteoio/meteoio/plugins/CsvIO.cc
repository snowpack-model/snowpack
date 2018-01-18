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

bool isQuote(const char& c)
{
	if (c=='"' || c=='\'') return true;
	return false;
}

CsvIO::CsvIO(const std::string& configfile) 
      : cfg(configfile), indexer(), vecStations(), vecFilenames(), csv_fields(), units_offset(), units_multiplier(),
        datetime_idx(), coordin(), coordinparam(), coordout(), coordoutparam(), 
        meteopath(), datetime_format(), csv_tz(), header_lines(1), columns_headers(1), csv_delim(',')
{
	parseInputOutputSection();
}

CsvIO::CsvIO(const Config& cfgreader)
      : cfg(cfgreader), indexer(), vecStations(), vecFilenames(), csv_fields(), units_offset(), units_multiplier(),
        datetime_idx(), coordin(), coordinparam(), coordout(), coordoutparam(), 
        meteopath(), datetime_format(), csv_tz(), header_lines(1), columns_headers(1), csv_delim(',')
{
	parseInputOutputSection();
}

void CsvIO::parseInputOutputSection()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	//HACK several of these should be per station!
	cfg.getValue("CSV_DELIMITER", "Input", csv_delim, IOUtils::nothrow);
	cfg.getValue("CSV_HEADER_LINES", "Input", header_lines, IOUtils::nothrow);
	cfg.getValue("CSV_COLUMNS_HEADERS", "Input", columns_headers, IOUtils::nothrow);
	cfg.getValue("CSV_FIELDS", "Input", csv_fields, IOUtils::nothrow);
	std::vector<std::string> vecMetaSpec;
	cfg.getValue("CSV_SPECIAL_HEADERS", "Input", vecMetaSpec, IOUtils::nothrow);
	cfg.getValue("CSV_UNITS_OFFSET", "Input", units_offset, IOUtils::nothrow);
	cfg.getValue("CSV_UNITS_MULTIPLIER", "Input", units_multiplier, IOUtils::nothrow);
	
	cfg.getValue("TIME_ZONE", "Input", csv_tz);
	std::string datetime_spec;
	cfg.getValue("CSV_DATETIME_SPEC", "Input", datetime_spec, IOUtils::nothrow);
	datetime_format = setDateParsing(datetime_spec);
	
	cfg.getValue("METEOPATH", "Input", meteopath);
	cfg.getValues("STATION", "INPUT", vecFilenames);
	const std::vector< std::pair<std::string, std::string> > vecMeta( cfg.getValues("META", "INPUT") );
	if (vecFilenames.size()!=vecMeta.size())
		throw InvalidArgumentException("The declared stations and metadata must match!", AT);
	vecStations = initStations(vecFilenames, vecMeta, vecMetaSpec);
}

struct sort_pred {
	bool operator()(const std::pair<size_t,size_t> &left, const std::pair<size_t,size_t> &right) {
		return left.first < right.first;
	}
};
	
//from a SPEC string such as "DD.MM.YYYY HH24:MIN:SS", build the format string for scanf as well as the parameters indices
//the indices are based on ISO timestamp, so year=0, month=1, etc
std::string CsvIO::setDateParsing(const std::string& datetime_spec)
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
	
	std::string format( datetime_spec );
	IOUtils::replace_all(format, "DD", "%u");
	IOUtils::replace_all(format, "MM", "%u");
	IOUtils::replace_all(format, "YYYY", "%u");
	IOUtils::replace_all(format, "HH24", "%u");
	IOUtils::replace_all(format, "MI", "%u");
	IOUtils::replace_all(format, "SS", "%u");
	
	return format;
}

std::vector<StationData> CsvIO::initStations(const std::vector<std::string>& i_vecFilenames, const std::vector< std::pair<std::string, std::string> >& vecMeta, const std::vector<std::string>& vecMetaSpec) const
{
	//HACK add special coding on how to extract meta information out of the header: like id:1:2 name:1:6 (as line:col)
	//HACK: support multiple fields per line, ie use multimap
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
	
	if (!meta_spec.empty()) {
		std::vector<std::string> vecStr;
		for (size_t ii=0; ii<i_vecFilenames.size(); ii++) {
			const std::string filename( meteopath+"/"+i_vecFilenames[ii] );
			if (!FileUtils::fileExists(filename)) throw AccessException("File '"+filename+"' does not exists", AT); //prevent invalid filenames
			errno = 0;
			std::ifstream fin(filename.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
			if (fin.fail()) {
				ostringstream ss;
				ss << "Error opening file \"" << filename << "\" for reading, possible reason: " << strerror(errno);
				ss << " Please check file existence and permissions!";
				throw AccessException(ss.str(), AT);
			}
			const char eoln = FileUtils::getEoln(fin);
			size_t linenr=0;
			std::string line;
	
			for (size_t jj=0; jj<header_lines; jj++) {
				getline(fin, line, eoln); //read complete signature line
				linenr++;
				if (meta_spec.count(linenr)>0) {
					IOUtils::readLineToVec(line, vecStr, csv_delim);
					if (meta_spec[linenr].first>vecStr.size())
						throw InvalidArgumentException("Metadata specification for '"+meta_spec[linenr].second+"' refers to a non-existent field for file '"+filename+"'", AT);
					const std::string field_type( IOUtils::strToUpper(meta_spec[linenr].second) );
					if (field_type=="NAME") std::cout << "found name='" << vecStr[meta_spec[linenr].first-1] << "\n";
				}
			}
		}
	}
		
		
	std::vector<StationData> vecStats;
	for (size_t ii=0; ii<vecMeta.size(); ii++) {
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		const Coords curr_point(coordin, coordinparam, vecMeta[ii].second);
		const std::string id_num( vecMeta[ii].first.substr(string("Meta").length()) );
		StationData sd(curr_point, "ST"+id_num, "Station_"+id_num);
		
		vecStats.push_back( sd );
	}
	
	return vecStats;
}

//return the columns' names, identify the timestamp or separate date and time columns indices
std::vector<std::string> CsvIO::readHeaders(std::ifstream& fin, const char& eoln, size_t &date_col, size_t &time_col) const
{
	date_col = time_col = 0; //default: first column for datetime
	size_t linenr=0;
	std::string line;
	std::vector<std::string> vecNames;
	
	for (size_t ii=0; ii<header_lines; ii++) {
		getline(fin, line, eoln); //read complete signature line
		linenr++;
		if (linenr==columns_headers) {
			IOUtils::readLineToVec(line, vecNames, csv_delim);
		}
	}
	
	for (size_t ii=0; ii<vecNames.size(); ii++) {
		std::string &tmp = vecNames[ii];
		IOUtils::toUpper( tmp );
		tmp.erase(std::remove_if(tmp.begin(), tmp.end(), &isQuote), tmp.end());
		if (tmp.empty()) continue;
		
		//HACK let the user provide the columns indices if it is not in the header.
		if (tmp.compare("TIMESTAMP")==0) {
			date_col = time_col = ii;
		} else if (tmp.compare("DATE")==0) {
			date_col = ii;
		} else if (tmp.compare("TIME")==0) {
			time_col = ii;
		}
	}
	
	if (!vecNames.empty()) { //cleanup potential '\r' char at the end of the line
		std::string &tmp = vecNames.back();
		if (*tmp.rbegin()=='\r') tmp.erase(tmp.end()-1);
	}
	
	return vecNames;
}

void CsvIO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation = vecStations;
}

Date CsvIO::parseDate(const std::string& date_str, const std::string& /*time_str*/) const
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

std::vector<MeteoData> CsvIO::readCSVFile(const std::string& filename, const size_t& stat_idx, const Date& /*dateStart*/, const Date& /*dateEnd*/)
{
	if (!FileUtils::fileExists(filename)) throw AccessException("File '"+filename+"' does not exists", AT); //prevent invalid filenames
	errno = 0;
	std::ifstream fin(filename.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << filename << "\" for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	const char eoln = FileUtils::getEoln(fin);
	
	std::string line;
	std::vector<MeteoData> vecMeteo;
	std::vector<std::string> tmp_vec;
	size_t date_col, time_col;
	std::vector<std::string> vecColNames( readHeaders(fin, eoln, date_col, time_col) );
	if (vecColNames.empty()) vecColNames = csv_fields;
	if (vecColNames.empty())
		throw InvalidArgumentException("No columns names could be retrieve, please provide them through the configuration file", AT);
	
	size_t nr_of_data_fields = vecColNames.size();
	const bool use_offset = !units_offset.empty();
	const bool use_multiplier = !units_multiplier.empty();
	if ((use_offset && units_offset.size()!=nr_of_data_fields) || (use_multiplier && units_multiplier.size()!=nr_of_data_fields)) {
		throw InvalidFormatException("The declared units_offset / units_multiplier must match the number of columns in the file!", AT);
	}
	
	size_t linenr=header_lines;
	while (!fin.eof()){
		line.clear();
		getline(fin, line, eoln);
		linenr++;
		
		if (line.empty()) continue; //Pure comment lines and empty lines are ignored
		
		const size_t nr_curr_data_fields = IOUtils::readLineToVec(line, tmp_vec, csv_delim);
		if (nr_of_data_fields==0) nr_of_data_fields=nr_curr_data_fields;
		if (nr_curr_data_fields!=nr_of_data_fields) {
			std::ostringstream ss;
			ss << "File \'" << filename << "\' declares (either as first data line or columns headers or units offset/multiplier) " << nr_of_data_fields << " columns ";
			ss << "but this does not match the following line:\n" << line << "\n";
			throw InvalidFormatException(ss.str(), AT);
		}
		
		const Date dt( parseDate(tmp_vec[date_col], tmp_vec[time_col]) );
		if (dt.isUndef()) {
			const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
			throw InvalidFormatException("Date could not be read in file \'"+filename+"' at line "+linenr_str, AT);
		}
		
		MeteoData md( dt, vecStations[stat_idx]);
		for (size_t ii=0; ii<tmp_vec.size(); ii++){
			if (ii==date_col) continue;
			double tmp;
			if (!IOUtils::convertString(tmp, tmp_vec[ii])) {
				const std::string linenr_str( static_cast<ostringstream*>( &(ostringstream() << linenr) )->str() );
				throw InvalidFormatException("Could not parse field '"+tmp_vec[ii]+"' in file \'"+filename+"' at line "+linenr_str, AT);
			}
			if (use_multiplier) tmp *= units_multiplier[ii];
			if (use_offset) tmp += units_offset[ii];
			md.addParameter( vecColNames[ii] );
			md( vecColNames[ii] ) = tmp;
		}
		vecMeteo.push_back( md );
	}
	
	return vecMeteo;
}

void CsvIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	for (size_t ii=0; ii<vecFilenames.size(); ii++) {
		vecMeteo.push_back( readCSVFile(meteopath+"/"+vecFilenames[ii], ii, dateStart, dateEnd) );
	}
}

} //namespace
