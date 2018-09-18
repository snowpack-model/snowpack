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
#ifndef CSVIO_H
#define CSVIO_H

#include <meteoio/IOInterface.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/Coords.h>

#include <string>
#include <vector>

namespace mio {

class CsvParameters {
	public:
		CsvParameters(const double& tz_in)
		: csv_fields(), units_offset(), units_multiplier(), nodata("NAN"), date_col(0), time_col(0), header_lines(1), columns_headers(IOUtils::npos), units_headers(IOUtils::npos), csv_delim(','), eoln('\n'), location(), datetime_idx(), time_idx(), file_and_path(), datetime_format(), time_format(), single_field(), name(), id(), slope(IOUtils::nodata), azi(IOUtils::nodata), csv_tz(static_cast<float>(tz_in)), has_tz(false) {}
		
		void setDateTimeSpec(const std::string& datetime_spec);
		void setTimeSpec(const std::string& time_spec);
		void setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec, const std::string& filename_spec, const std::string& station_idx="");
		void setLocation(const Coords i_location, const std::string& i_name, const std::string& i_id) {location=i_location; name=i_name; id=i_id;}
		Date parseDate(const std::string& datetime_str, const std::string& time_str) const;
		std::string getFilename() const {return file_and_path;}
		StationData getStation() const;
		
		std::vector<std::string> csv_fields;		///< the user provided list of field names
		std::vector<double> units_offset, units_multiplier;		///< offsets and multipliers to convert the data to SI
		
		std::string nodata;
		size_t date_col, time_col;
		size_t header_lines, columns_headers, units_headers;
		char csv_delim;
		char eoln;
	private:
		void assignMetadataVariable(const std::string& field_type, const std::string& field_val, double &lat, double &lon);
		void parseFileName(std::string filename, const std::string& filename_spec, double &lat, double &lon);
		void parseFields(const std::vector<std::string>& headerFields, std::vector<std::string>& fieldNames, size_t &dt_col, size_t &tm_col);
		static std::multimap< size_t, std::pair<size_t, std::string> > parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec);
		void parseSpecialHeaders(const std::string& line, const size_t& linenr, const std::multimap< size_t, std::pair<size_t, std::string> >& meta_spec, double &lat, double &lon);
		void parseUnits(const std::string& line);
		static void checkSpecString(const std::string& spec_string, const size_t& nr_params);
		
		Coords location;
		std::vector<size_t> datetime_idx, time_idx;		///< order of the datetime fields for use in parseDate
		std::string file_and_path, datetime_format, time_format, single_field; 		///< the scanf() format string for use in parseDate, the parameter in case of a single value contained in the Csv file
		std::string name, id;
		double slope, azi;
		double csv_tz;		///< timezone to apply to parsed dates
		bool has_tz;		///< does the user-provided date/time format contains a TZ?
};

/**
 * @class CsvIO
 * @brief Reads meteo data from a comma separated file.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2018-01
 */
class CsvIO : public IOInterface {
	public:
		CsvIO(const std::string& configfile);
		CsvIO(const CsvIO&);
		CsvIO(const Config& cfgreader);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void parseInputOutputSection();
		void cleanup() throw();
		std::string setDateParsing(const std::string& datetime_spec);
		std::vector<std::string> readHeaders(std::ifstream& fin, CsvParameters& params) const;
		Date parseDate(const std::string& date_str, const std::string& time_str) const;
		std::vector<MeteoData> readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd);
		
		const Config cfg;
		std::map<std::string, mio::FileUtils::FileIndexer> indexer_map;
		std::vector<CsvParameters> csvparam;
		std::vector<StationData> vecStations;
		std::string coordin, coordinparam; //projection parameters
		static const size_t streampos_every_n_lines; //save current stream pos every n lines of data
		bool silent_errors; ///< when reading a file, should errors throw or just be ignored?
		bool errors_to_nodata;    //unparseable values are treated as nodata, but the dataset is kept
};

} //namespace
#endif
