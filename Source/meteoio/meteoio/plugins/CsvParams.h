// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2023 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef CSVPARAMS_H
#define CSVPARAMS_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOInterface.h> //for LinesRange

#include <string>
#include <vector>

namespace mio {

///class to contain date and time parsing information
class CsvDateTime {
	public:
		CsvDateTime(const double& tz_in);
		
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
		int getFixedHour();
		bool isSet() const;
		
		void setDateTimeSpec(const std::string& datetime_spec);
		void setDateSpec(const std::string& date_spec);
		void setTimeSpec(const std::string& time_spec);
		void setDecimalDateType(std::string i_decimaldate_type);
		void setFixedYear(const int& i_year, const bool& i_auto_wrap);
		void setFixedHour(const int& i_hour);
		bool parseField(const std::string& fieldname, const size_t &ii);
		Date parseDate(const std::vector<std::string>& vecFields);
		std::string toString() const;
		
		size_t max_dt_col; ///< Maximum index of a date/time field (for optimized parsing)
		bool auto_wrap; ///< if true, dates >= October will be assumed to belong to (year_cst-1) until a date < October is encountered
		
	private:
		static bool parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, int& value);
		static bool parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, double& value);
		Date parseDate(const std::string& date_time_str) const;
		Date parseDate(const std::string& value_str, const CsvDateTime::decimal_date_formats& format) const;
		bool parseDate(const std::string& date_str, float args[3]) const;
		bool parseTime(const std::string& time_str, float args[3], double& tz) const;
		double parseTime(const std::string& time_str, double& tz) const;
		
		static int castToInt(const float &val);
		static void checkSpecString(const std::string& spec_string, const size_t& nr_params);
		
		std::vector<size_t> datetime_idx;		///< order of the datetime fields for use in parseDate: Year Month Day Hour Minutes Seconds
		std::vector<size_t> date_idx;		///< order of the datetime fields for use in parseDate: Year Month Day
		std::vector<size_t> time_idx;		///< order of the time fields for use in parseDate for split date / time
		std::string datetime_format, date_format, time_format; 		///< the scanf() format string for use in parseDate, the parameter in case of a single value contained in the Csv file
		decimal_date_formats decimal_date_type; ///< in case of decimal date, which representation is associated with it
		double csv_tz;		///< timezone to apply to parsed dates
		//time is a field that contains numerical time, for example 0920
		size_t idx_decimal_date, idx_date_time_str, idx_date_str, idx_time_str, idx_year, idx_jdn, idx_month, idx_day, idx_ntime, idx_hours, idx_minutes, idx_seconds;
		static const int cutoff_year = 40;
		int year_cst, hour_cst; ///< When the year/hour is not provided as such but set from the folder name as "fallback year"
		bool has_tz;		///< does the user-provided date/time format contains a TZ?
		bool dt_as_decimal;	///< is date provided as a single decimal number?
		bool dt_2digits_year;	///< is the year only provided with 2 digits?
	};

class CsvParameters {
	public:
		CsvParameters(const double& tz_in);
		
		void setHeaderRepeatMk(const std::string& marker) {header_repeat_mk=marker;}
		void setDelimiter(const std::string& delim);
		void setHeaderDelimiter(const std::string& delim);
		void setSkipFields(const std::string& skipFieldSpecs, const bool& negate);
		void setUnits(const std::string& csv_units,  const char& delim=' ');
		void setLinesExclusions(const std::vector< LinesRange >& linesSpecs) {linesExclusions=linesSpecs;}
		void setNodata(const std::string& nodata_markers);
		void setPurgeChars(const std::string& chars_to_purge);
		void setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec, const std::string& filename_spec, const std::string& station_idx="");
		void setLocation(const Coords i_location, const std::string& i_name, const std::string& i_id) {location=i_location; name=i_name; id=i_id;}
		void setSlope(const double& i_slope, const double& i_azimuth) {slope=i_slope; azi=i_azimuth;}
		void setDateTimeSpecs(const std::string &datetime_spec, const std::string &date_spec, const std::string &time_spec, const std::string &decimaldate_type);
		void setFixedYear(const int& i_year, const bool& i_auto_wrap) {date_cols.setFixedYear(i_year, i_auto_wrap);}
		void setFixedHour(const int& i_hour) {date_cols.setFixedHour(i_hour);}
		
		std::string toString() const;
		std::string getFilename() const {return file_and_path;}
		StationData getStation() const;
		Date getDate(const std::vector<std::string>& vecFields) {return date_cols.parseDate(vecFields);}
		bool excludeLine(const size_t& linenr, bool& hasExclusions);
		bool skipField(const size_t& fieldnr) const;
		bool hasPurgeChars() const {return !purgeCharsSet.empty();}
		void purgeChars(std::string &line) {IOUtils::removeChars(line, purgeCharsSet);}
		bool isNodata(const std::string& value) const;
		
		std::vector<std::string> csv_fields;		///< the user provided list of field names
		std::vector<double> units_offset, units_multiplier;		///< offsets and multipliers to convert the data to SI
		std::vector<double> field_offset, field_multiplier;		///< offsets and multipliers to apply to each field
		
		std::string header_repeat_mk, filter_ID, fields_postfix;
		size_t ID_col;
		size_t header_lines, columns_headers, units_headers;
		char csv_delim, header_delim;
		char eoln, comments_mk;
		bool header_repeat_at_start, asc_order;
		bool number_fields;	///< include a column number in the field names as well as an optional field_postfix (this helps when debugging invalid/changing column order)
		
	private:
		static std::string identifyField(const std::string& fieldname);
		void assignMetadataVariable(const std::string& field_type, const std::string& field_val, double &lat, double &lon, double &easting, double &northing);
		void parseFileName(std::string filename, const std::string& filename_spec, double &lat, double &lon, double &easting, double &northing);
		void parseFields(const std::vector<std::string>& headerFields, std::vector<std::string>& fieldNames);
		static std::multimap< size_t, std::pair<size_t, std::string> > parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec);
		void parseSpecialHeaders(const std::string& line, const size_t& linenr, const std::multimap< size_t, std::pair<size_t, std::string> >& meta_spec, double &lat, double &lon, double &easting, double &northing);
		
		CsvDateTime date_cols;		///< index of each column containing the a date/time component
		Coords location;
		std::set<std::string> nodata;		///< representations of nodata with their variants (with quotes, with double quotes, etc)
		std::set<size_t> skip_fields;		///< Fields that should not be read
		std::set<char> purgeCharsSet;			///< characters to purge from each line (such as quotes, double quotes, etc)
		std::vector< LinesRange > linesExclusions;	///< lines to exclude from reading
		std::string file_and_path, single_field; 		///< the scanf() format string for use in parseDate, the parameter in case of a single value contained in the Csv file
		std::string name, id;
		double slope, azi;
		size_t exclusion_idx;		///< pointer to the latest exclusion period that has been found, if using lines exclusion
		size_t exclusion_last_linenr; ///< pointer to the last line number that has been checked for exclusions
		size_t last_allowed_field;	///< index of the last allowed field (as set by the user with setSkipFields(negate=true)
};

} //namespace
#endif
