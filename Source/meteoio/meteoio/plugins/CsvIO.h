// SPDX-License-Identifier: LGPL-3.0-or-later
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

///class to contain date and time parsing information
class CsvDateTime {
	public:
		CsvDateTime() : decimal_date(IOUtils::npos), date_str(IOUtils::npos), time_str(IOUtils::npos), year(IOUtils::npos), jdn(IOUtils::npos), month(IOUtils::npos), day(IOUtils::npos), time(IOUtils::npos), hours(IOUtils::npos), minutes(IOUtils::npos), seconds(IOUtils::npos), max_dt_col(0), year_cst(IOUtils::inodata), decimal_date_type(JULIAN), auto_wrap(true) {}
		
		//this matches the formats that are supported in Date
		typedef enum DECIMAL_DATE_FORMATS {
			EXCEL, ///< Excel date
			JULIAN, ///< standard julian date
			MJULIAN, ///< Modified julian date
			MATLAB, ///< Matlab  date
			RFC868, ///< RFC 868 date
			UNIX ///< Unix date
		} decimal_date_formats;
		
		void updateMaxCol();
		int getFixedYear(const double& i_jdn);
		int getFixedYear(const int& i_month);
		bool isSet() const;
		std::string toString() const;
		
		//time is a field that contains numerical time, for example 0920
		size_t decimal_date, date_str, time_str, year, jdn, month, day, time, hours, minutes, seconds;
		size_t max_dt_col;
		int year_cst;
		decimal_date_formats decimal_date_type; ///< in case of decimal date, which representation is associated with it
		bool auto_wrap; ///< if true, dates >= October will be assumed to belong to (year_cst-1) until a date < October is encountered
	};

class CsvParameters {
	public:
		CsvParameters(const double& tz_in)
		: csv_fields(), units_offset(), units_multiplier(), skip_fields(), nodata("NAN"), header_repeat_mk(), filter_ID(), ID_col(IOUtils::npos), header_lines(1), columns_headers(IOUtils::npos), units_headers(IOUtils::npos), single_param_idx(IOUtils::npos), csv_delim(','), header_delim(','), eoln('\n'), comments_mk('\n'), header_repeat_at_start(false), asc_order(true), purgeQuotes(false),  location(), datetime_idx(), time_idx(), file_and_path(), datetime_format(), time_format(), single_field(), name(), id(), date_cols(), slope(IOUtils::nodata), azi(IOUtils::nodata), csv_tz(tz_in), has_tz(false), dt_as_components(false), dt_as_year_and_jdn(false), dt_as_decimal(false) {}
		
		void setPurgeQuotes(const bool& i_purgeQuotes) {purgeQuotes=i_purgeQuotes;}
		void setHeaderRepeatMk(const std::string& marker) {header_repeat_mk=marker;}
		void setDelimiter(const std::string& delim);
		void setHeaderDelimiter(const std::string& delim);
		void setSkipFields(const std::vector<size_t>& vecSkipFields);
		void setUnits(const std::string& csv_units,  const char& delim);
		void setDateTimeSpec(const std::string& datetime_spec);
		void setTimeSpec(const std::string& time_spec);
		void setDecimalDateType(std::string decimaldate_type);
		void setFixedYear(const int& i_year, const bool& auto_wrap);
		void setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec, const std::string& filename_spec, const std::string& station_idx="");
		void setLocation(const Coords i_location, const std::string& i_name, const std::string& i_id) {location=i_location; name=i_name; id=i_id;}
		void setSlope(const double& i_slope, const double& i_azimuth) {slope=i_slope; azi=i_azimuth;}
		Date parseDate(const std::vector<std::string>& vecFields);
		std::string getFilename() const {return file_and_path;}
		StationData getStation() const;
		
		std::vector<std::string> csv_fields;		///< the user provided list of field names
		std::vector<double> units_offset, units_multiplier;		///< offsets and multipliers to convert the data to SI
		std::map<size_t, bool> skip_fields;		///< Fields that should not be read
		
		std::string nodata, header_repeat_mk, filter_ID;
		size_t ID_col;
		size_t header_lines, columns_headers, units_headers;
		size_t single_param_idx;
		char csv_delim, header_delim;
		char eoln, comments_mk;
		bool header_repeat_at_start, asc_order, purgeQuotes;
	private:
		static std::string identifyField(const std::string& fieldname);
		void assignMetadataVariable(const std::string& field_type, const std::string& field_val, double &lat, double &lon, double &easting, double &northing);
		void parseFileName(std::string filename, const std::string& filename_spec, double &lat, double &lon, double &easting, double &northing);
		void parseFields(const std::vector<std::string>& headerFields, std::vector<std::string>& fieldNames);
		static std::multimap< size_t, std::pair<size_t, std::string> > parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec);
		void parseSpecialHeaders(const std::string& line, const size_t& linenr, const std::multimap< size_t, std::pair<size_t, std::string> >& meta_spec, double &lat, double &lon, double &easting, double &northing);
		static Date createDate(const float args[6], const double i_tz);
		static bool parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, int& value);
		static bool parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, double& value);
		Date parseJdnDate(const std::vector<std::string>& vecFields);
		Date parseDate(const std::string& date_str, const std::string& time_str) const;
		Date parseDate(const std::string& value_str, const CsvDateTime::decimal_date_formats& format) const;
		static void checkSpecString(const std::string& spec_string, const size_t& nr_params);
		
		Coords location;
		std::vector<size_t> datetime_idx;		///< order of the datetime fields for use in parseDate: Year Month Day Hour Minutes Seconds
		std::vector<size_t> time_idx;		///< order of the time fields for use in parseDate for split date / time
		std::string file_and_path, datetime_format, time_format, single_field; 		///< the scanf() format string for use in parseDate, the parameter in case of a single value contained in the Csv file
		std::string name, id;
		CsvDateTime date_cols;		///< index of each column containing the a date/time component
		double slope, azi;
		double csv_tz;		///< timezone to apply to parsed dates
		bool has_tz;		///< does the user-provided date/time format contains a TZ?
		bool dt_as_components; 	///< is date/time provided as components each in its own column (ie an hour column, a day column, etc)?
		bool dt_as_year_and_jdn;	///< is date provided as year + julian day?
		bool dt_as_decimal;	///< is date provided as a single decimal number?
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
		void cleanup() noexcept;
		std::string setDateParsing(const std::string& datetime_spec);
		std::vector<std::string> readHeaders(std::ifstream& fin, CsvParameters& params) const;
		static MeteoData createTemplate(const CsvParameters& params);
		static Date getDate(CsvParameters& params, const std::vector<std::string>& vecFields, const bool& silent_errors, const std::string& filename, const size_t& linenr);
		std::vector<MeteoData> readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd);
		
		const Config cfg;
		std::map<std::string, FileUtils::FileIndexer> indexer_map;
		std::vector<CsvParameters> csvparam;
		std::vector<StationData> vecStations;
		std::string coordin, coordinparam; //projection parameters
		static const size_t streampos_every_n_lines; //save current stream pos every n lines of data
		bool silent_errors; ///< when reading a file, should errors throw or just be ignored?
		bool errors_to_nodata;    //unparseable values are treated as nodata, but the dataset is kept
};

} //namespace
#endif
