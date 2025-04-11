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
#include <meteoio/plugins/CsvParams.h>
#include <meteoio/FileUtils.h>

#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <regex>
#include <iomanip>

using namespace std;

namespace mio {

//helper function to sort the static keys used for specifying the date/time formats
inline bool sort_dateKeys(const std::pair<size_t,size_t> &left, const std::pair<size_t,size_t> &right) { return left.first < right.first;}

CsvDateTime::CsvDateTime(const double& tz_in) 
           : max_dt_col(0), csv_tz(tz_in), auto_wrap(false), datetime_idx(), date_idx(), time_idx(), datetime_format(), date_format(), time_format(), decimal_date_type(JULIAN),
           idx_decimal_date(IOUtils::npos), idx_date_time_str(IOUtils::npos), idx_date_str(IOUtils::npos), idx_time_str(IOUtils::npos), idx_year(IOUtils::npos), idx_jdn(IOUtils::npos), idx_month(IOUtils::npos), idx_day(IOUtils::npos), idx_ntime(IOUtils::npos), idx_hours(IOUtils::npos), idx_minutes(IOUtils::npos), idx_seconds(IOUtils::npos), year_cst(IOUtils::inodata), hour_cst(IOUtils::inodata),
           has_tz(false), dt_as_decimal(false), dt_2digits_year(false)
{}

void CsvDateTime::updateMaxCol() 
{
	if (idx_decimal_date!=IOUtils::npos && idx_decimal_date>max_dt_col) max_dt_col=idx_decimal_date;
	if (idx_date_time_str!=IOUtils::npos && idx_date_time_str>max_dt_col) max_dt_col=idx_date_time_str;
	if (idx_date_str!=IOUtils::npos && idx_date_str>max_dt_col) max_dt_col=idx_date_str;
	if (idx_time_str!=IOUtils::npos && idx_time_str>max_dt_col) max_dt_col=idx_time_str;
	if (idx_year!=IOUtils::npos && idx_year>max_dt_col) max_dt_col=idx_year;
	if (idx_jdn!=IOUtils::npos && idx_jdn>max_dt_col) max_dt_col=idx_jdn;
	if (idx_month!=IOUtils::npos && idx_month>max_dt_col) max_dt_col=idx_month;
	if (idx_day!=IOUtils::npos && idx_day>max_dt_col) max_dt_col=idx_day;
	if (idx_ntime!=IOUtils::npos && idx_ntime>max_dt_col) max_dt_col=idx_ntime;
	if (idx_hours!=IOUtils::npos && idx_hours>max_dt_col) max_dt_col=idx_hours;
	if (idx_minutes!=IOUtils::npos && idx_minutes>max_dt_col) max_dt_col=idx_minutes;
	if (idx_seconds!=IOUtils::npos && idx_seconds>max_dt_col) max_dt_col=idx_seconds;
}

int CsvDateTime::getFixedYear(const double& i_jdn)
{
	if (i_jdn<273.) auto_wrap = false;
	if (auto_wrap) return year_cst - 1;
	return year_cst;
}

int CsvDateTime::getFixedYear(const int& i_month)
{
	if (i_month<10) auto_wrap = false;
	if (auto_wrap) return year_cst - 1;
	return year_cst;
}

int CsvDateTime::getFixedHour()
{
	return hour_cst;
}

bool CsvDateTime::isSet() const 
{
	//date and time strings
	if ((idx_date_str!=IOUtils::npos && idx_time_str!=IOUtils::npos) || idx_date_time_str!=IOUtils::npos) return true;
	
	//purely decimal timestamps
	if (idx_decimal_date!=IOUtils::npos) return true;
	
	//as components, possibly mixed with some string representations
	const bool has_component_date = (idx_year!=IOUtils::npos || year_cst!=IOUtils::inodata) && ((idx_month!=IOUtils::npos && idx_day!=IOUtils::npos) || (idx_jdn!=IOUtils::npos));
	const bool has_date = has_component_date || idx_date_str!=IOUtils::npos;
	const bool has_component_time = (idx_hours!=IOUtils::npos) || (idx_ntime!=IOUtils::npos);
	const bool has_time = has_component_time || idx_time_str!=IOUtils::npos || hour_cst!=IOUtils::inodata ;
	if (has_date && has_time) return true;

	return false;
}

//check that the format is usable (and prevent parameters injection / buffer overflows)
void CsvDateTime::checkSpecString(const std::string& spec_string, const size_t& nr_params)
{
	const size_t nr_percent = (unsigned)std::count(spec_string.begin(), spec_string.end(), '%');
	const size_t nr_placeholders0 = IOUtils::count(spec_string, "%f");
	const size_t nr_placeholders2 = IOUtils::count(spec_string, "%2f");
	const size_t nr_placeholders4 = IOUtils::count(spec_string, "%4f");
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
//the indices are based on ISO timestamp, so year=0, month=1, ..., ss=5 while tz is handled separately
void CsvDateTime::setDateTimeSpec(const std::string& datetime_spec)
{
	//support for 2 digits year is a little hacky: it requires special processing in order to avoid machting both
	//YYYY and YY as well as having a boolean to know if we should add an offset to the year or not
	static const std::vector< std::pair<std::string, unsigned short> > keys( {{"YYYY", 0}, {"YY", 0}, {"MM", 1}, {"DD", 2}, {"HH24", 3}, {"MI", 4}, {"SS", 5}} );
	dt_2digits_year = (datetime_spec.find("YYYY")==std::string::npos && datetime_spec.find("YY")!=std::string::npos);
	
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (const auto& spec_key : keys) {
		if (!dt_2digits_year && spec_key.first=="YY") continue;	//skip looking for 2 digit years if not applicable
		const size_t key_pos = datetime_spec.find( spec_key.first );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, spec_key.second) );
	}
	
	//fill datetime_idx as a vector of [0-5] indices (for ISO fields) in the order they appear in the user-provided format string
	std::sort(sorting_vector.begin(), sorting_vector.end(), &sort_dateKeys);
	for (const auto& spec_key : sorting_vector)
		datetime_idx.push_back( spec_key.second ); //push the matching position in the static const "keys"
	
	datetime_format = datetime_spec;
	const size_t tz_pos = datetime_format.find("TZ");
	if (tz_pos!=std::string::npos) {
		if (tz_pos!=(datetime_format.length()-2))
			throw InvalidFormatException("When providing TZ in a date/time format, it must be at the very end of the string", AT);
		has_tz = true;
		datetime_format.replace(tz_pos, 2, "%32s");
	}
	//the order of the replacements is important in order to match YYYY first and then YY
	IOUtils::replace_all(datetime_format, "DD",   "%2f");
	IOUtils::replace_all(datetime_format, "MM",   "%2f");
	IOUtils::replace_all(datetime_format, "YYYY", "%4f");
	IOUtils::replace_all(datetime_format, "YY",   "%2f");
	IOUtils::replace_all(datetime_format, "HH24", "%2f");
	IOUtils::replace_all(datetime_format, "MI",   "%2f");
	IOUtils::replace_all(datetime_format, "SS",   "%f");
	
	const size_t nr_params_check = (has_tz)? datetime_idx.size()+1 : datetime_idx.size();
	checkSpecString(datetime_format, nr_params_check);
}

//from a SPEC string such as "DD.MM.YYYY", build the format string for scanf as well as the parameters indices
//the indices are based on ISO timestamp, so year=0, month=1, day=2
void CsvDateTime::setDateSpec(const std::string& date_spec)
{
	//support for 2 digits year is a little hacky: it requires special processing in order to avoid machting both
	//YYYY and YY as well as having a boolean to know if we should add an offset to the year or not
	static const std::vector< std::pair<std::string, unsigned short> > keys( {{"YYYY", 0}, {"YY", 0}, {"MM", 1}, {"DD", 2}, {"HH24", 3}, {"MI", 4}, {"SS", 5}} );
	dt_2digits_year = (date_spec.find("YYYY")==std::string::npos && date_spec.find("YY")!=std::string::npos);
	
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (const auto& spec_key : keys) {
		if (!dt_2digits_year && spec_key.first=="YY") continue;	//skip looking for 2 digit years if not applicable
		const size_t key_pos = date_spec.find( spec_key.first );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, spec_key.second) );
	}
	
	//fill date_idx as a vector of [0-2] indices (for ISO fields) in the order they appear in the user-provided format string
	std::sort(sorting_vector.begin(), sorting_vector.end(), &sort_dateKeys);
	for (const auto& spec_key : sorting_vector)
		date_idx.push_back( spec_key.second ); //push the matching position in the static const "keys"
	
	date_format = date_spec;	
	//the order of the replacements is important in order to match YYYY first and then YY
	IOUtils::replace_all(date_format, "DD",   "%2f");
	IOUtils::replace_all(date_format, "MM",   "%2f");
	IOUtils::replace_all(date_format, "YYYY", "%4f");
	IOUtils::replace_all(date_format, "YY",   "%2f");
	
	const size_t nr_params_check = date_idx.size();
	checkSpecString(date_format, nr_params_check);
}

void CsvDateTime::setTimeSpec(const std::string& time_spec)
{
	if (time_spec.empty()) return;
	static const std::vector< std::pair<std::string, unsigned short> > keys( {{"HH24", 0}, {"MI", 1}, {"SS", 2}} );
	std::vector< std::pair<size_t, size_t> > sorting_vector;
	for (const auto& spec_key : keys) {
		const size_t key_pos = time_spec.find( spec_key.first );
		if (key_pos!=std::string::npos)
			sorting_vector.push_back( make_pair( key_pos, spec_key.second) );
	}

	//fill time_idx as a vector of [0-3] indices (for ISO fields) in the order they appear in the user-provided format string
	std::sort(sorting_vector.begin(), sorting_vector.end(), &sort_dateKeys);
	for (const auto& spec_key : sorting_vector)
		time_idx.push_back( spec_key.second ); //push the matching position in the static const "keys"

	time_format = time_spec;
	const size_t tz_pos = time_format.find("TZ");
	if (tz_pos!=std::string::npos) {
		if (tz_pos!=(time_format.length()-2))
			throw InvalidFormatException("When providing TZ in a date/time format, it must be at the very end of the string", AT);
		has_tz = true;
		time_format.replace(tz_pos, 2, "%32s");
	}
	IOUtils::replace_all(time_format, "HH24", "%2f");
	IOUtils::replace_all(time_format, "MI", "%2f");
	IOUtils::replace_all(time_format, "SS", "%f");

	const size_t nr_params_check = (has_tz)? time_idx.size()+1 : time_idx.size();
	checkSpecString(time_format, nr_params_check);
}

void CsvDateTime::setDecimalDateType(std::string i_decimaldate_type)
{
	IOUtils::toUpper( i_decimaldate_type );
	if (i_decimaldate_type=="EXCEL") {
		decimal_date_type = CsvDateTime::EXCEL;
	} else if (i_decimaldate_type=="JULIAN") {
		decimal_date_type = CsvDateTime::JULIAN;
	} else if (i_decimaldate_type=="MJULIAN") {
		decimal_date_type = CsvDateTime::MJULIAN;
	} else if (i_decimaldate_type=="MATLAB") {
		decimal_date_type = CsvDateTime::MATLAB;
	} else if (i_decimaldate_type=="RFC868") {
		decimal_date_type = CsvDateTime::RFC868;
	} else if (i_decimaldate_type=="UNIX") {
		decimal_date_type = CsvDateTime::UNIX;
	} else
		throw InvalidArgumentException("Unknown decimal date type '"+i_decimaldate_type+"'", AT);
	
	dt_as_decimal = true;
}

void CsvDateTime::setFixedYear(const int& i_year, const bool& i_auto_wrap)
{
	year_cst = i_year;
	auto_wrap = i_auto_wrap;
}

void CsvDateTime::setFixedHour(const int& i_hour)
{
	hour_cst = i_hour;
}

// return true if the field must be skiped (all special fields are marked as SKIP since they are read in a special way)
bool CsvDateTime::parseField(const std::string& fieldname, const size_t &ii)
{
	if (fieldname.compare("TIMESTAMP")==0 || fieldname.compare("TS")==0 || fieldname.compare("DATETIME")==0) {
		if (dt_as_decimal) 
			idx_decimal_date = ii;
		else 
			idx_date_time_str = ii;
		return true;
	} else if (fieldname.compare("DATE")==0 || fieldname.compare("GIORNO")==0 || fieldname.compare("FECHA")==0) {
		idx_date_str = ii;
		return true;
	} else if (fieldname.compare("TIME")==0 || fieldname.compare("ORA")==0 || fieldname.compare("HORA")==0) {
		idx_time_str = ii;
		return true;
	} else if (fieldname.compare("YEAR")==0) {
		idx_year = ii;
		return true;
	} else if (fieldname.compare("YEAR_2DIGITS")==0) {
		idx_year = ii;
		dt_2digits_year = true;
		return true;
	} else if (fieldname.compare("JDAY")==0 || fieldname.compare("JDN")==0 || fieldname.compare("YDAY")==0 || fieldname.compare("DAY_OF_YEAR")==0 || fieldname.compare("DOY")==0) {
		idx_jdn = ii;
		return true;
	} else if (fieldname.compare("MONTH")==0) {
		idx_month = ii;
		return true;
	} else if (fieldname.compare("DAY")==0) {
		idx_day = ii;
		return true;
	} else if (fieldname.compare("NTIME")==0) {
		idx_ntime = ii;
		return true;
	} else if (fieldname.compare("HOUR")==0 || fieldname.compare("HOURS")==0) {
		idx_hours = ii;
		return true;
	} else if (fieldname.compare("MINUTE")==0 || fieldname.compare("MINUTES")==0) {
		idx_minutes = ii;
		return true;
	} else if (fieldname.compare("SECOND")==0 || fieldname.compare("SECONDS")==0) {
		idx_seconds = ii;
		return true;
	}
	
	//check that the string parsing specs have been set if necessary
	if (idx_date_time_str!=IOUtils::npos && datetime_idx.empty()) throw InvalidArgumentException("Please define how to parse DATETIME strings with key DATETIME_SPEC", AT);
	if (idx_date_str!=IOUtils::npos && date_idx.empty()) throw InvalidArgumentException("Please define how to parse DATE strings with key DATE_SPEC", AT);
	if (idx_time_str!=IOUtils::npos && time_idx.empty()) throw InvalidArgumentException("Please define how to parse TIME strings with key TIME_SPEC", AT);
	
	return false;
}

int CsvDateTime::castToInt(const float &val)
{
	const int ival = (int)val;
	if ((float)ival!=val) return IOUtils::inodata;
	return ival;
}

bool CsvDateTime::parseDate(const std::string& date_str, float args[3]) const
{
	//parse the date information and return the ymd components (the date must be exactly read, if there is "rest" it is wrong!)
	char rest[32] = "";
	bool status = (sscanf(date_str.c_str(), date_format.c_str(), &args[ date_idx[0] ], &args[ date_idx[1] ], &args[ date_idx[2] ], rest)>=3);
	if (!status || rest[0]) return false;
	return true;
}

bool CsvDateTime::parseTime(const std::string& time_str, float args[3], double& tz) const
{
	//parse the time information and return the hms components
	char rest[32] = "";
	bool status = false;
	switch( time_idx.size() ) {
		case 3:
			status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0] ], &args[ time_idx[1] ], &args[ time_idx[2] ], rest)>=3);
			break;
		case 2:
			status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0] ], &args[ time_idx[1] ], rest)>=2);
			break;
		case 1:
			status = (sscanf(time_str.c_str(), time_format.c_str(), &args[ time_idx[0] ], rest)>=1);
			break;
		default: // do nothing;
			break;
	}
	if (!status) return false;
	
	tz = (has_tz)? Date::parseTimeZone(rest) : csv_tz;
	return true;
}

double CsvDateTime::parseTime(const std::string& time_str, double& tz) const
{
	static const double seconds_to_days = 1. / (24.*3600.);
	//parse the time information and return the fractional day
	float args_tm[3] = {0., 0., 0.};
	if (!parseTime(time_str, args_tm, tz)) return IOUtils::nodata;
	
	const double fractional_day = (static_cast<double>(args_tm[0])*3600. + static_cast<double>(args_tm[1])*60. + static_cast<double>(args_tm[2])) * seconds_to_days;
	return fractional_day;
}

Date CsvDateTime::parseDate(const std::string& date_time_str) const
{
	float args[6] = {0., 0., 0., 0., 0., 0.};
	char rest[32] = "";
	bool status = false;
	switch( datetime_idx.size() ) {
		case 6:
			status = (sscanf(date_time_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], &args[ datetime_idx[5] ], rest)>=6);
			break;
		case 5:
			status = (sscanf(date_time_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], &args[ datetime_idx[4] ], rest)>=5);
			break;
		case 4:
			status = (sscanf(date_time_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], &args[ datetime_idx[3] ], rest)>=4);
			break;
		case 3:
			status = (sscanf(date_time_str.c_str(), datetime_format.c_str(), &args[ datetime_idx[0] ], &args[ datetime_idx[1] ], &args[ datetime_idx[2] ], rest)>=3);
			break;
		default: // do nothing;
			break;
	}
	if (!status) return Date(); //we MUST have read successfuly at least the date part

	const double tz = (has_tz)? Date::parseTimeZone(rest) : csv_tz;
	
	//now build a Date object
	int year = castToInt(args[0]);
	if (dt_2digits_year) year = (year>cutoff_year)? 1900+year : 2000+year;
	
	return Date(year, castToInt(args[1]), castToInt(args[2]), castToInt(args[3]), castToInt(args[4]), static_cast<double>(args[5]), tz);
}

Date CsvDateTime::parseDate(const std::string& value_str, const CsvDateTime::decimal_date_formats& format) const
{
	Date dt;
	
	if (format==CsvDateTime::UNIX) {
		time_t value;
		if (!IOUtils::convertString(value, value_str)) return dt;
		dt.setUnixDate(value);
		return dt;
	} else {
		double value=0.;
		if (!IOUtils::convertString(value, value_str)) return dt;
		
		if (format==CsvDateTime::EXCEL) {
			dt.setExcelDate(value, csv_tz);
		} else if (format==CsvDateTime::JULIAN) {
			dt.setDate(value, csv_tz);
		} else if (format==CsvDateTime::MJULIAN) {
			dt.setModifiedJulianDate(value, csv_tz);
		} else if (format==CsvDateTime::MATLAB) {
			dt.setMatlabDate(value, csv_tz);
		} else if (format==CsvDateTime::RFC868) {
			dt.setRFC868Date(value, csv_tz);
		} 
	}
	
	return dt;
}

bool CsvDateTime::parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, int& value)
{
	if (idx==IOUtils::npos) {
		value=0;
		return true;
	}
	
	return IOUtils::convertString(value, vecFields[ idx ]);
}

bool CsvDateTime::parseDateComponent(const std::vector<std::string>& vecFields, const size_t& idx, double& value)
{
	if (idx==IOUtils::npos) {
		value=0.;
		return true;
	}
	
	return IOUtils::convertString(value, vecFields[ idx ]);
}

Date CsvDateTime::parseDate(const std::vector<std::string>& vecFields)
{
	//datetime as pure string
	if (idx_date_time_str!=IOUtils::npos) {
		return parseDate(vecFields[ idx_date_time_str ]);
	}

	//datetime as purely decimal
	if (dt_as_decimal) {
		return parseDate(vecFields[ idx_decimal_date ], decimal_date_type);
	}
	
	//date and time are provided as some combination of components
	int year=IOUtils::inodata, month=IOUtils::inodata, day=IOUtils::inodata, hour=IOUtils::inodata, minute=0;
	double seconds = 0., fractional_day = IOUtils::nodata, jdn = IOUtils::nodata, ntime = IOUtils::nodata, tz = csv_tz;

	//date as string
	if (idx_date_str != IOUtils::npos) {
		if (date_idx.size()!=3)
			throw InvalidFormatException("String date representation can only contain year, month and day when reading date/time as component", AT);

		float args[3] = {0., 0., 0.};
		if (!parseDate(vecFields[ idx_date_str ], args))
			throw InvalidFormatException("Could not parse date", AT);

		year = castToInt( args[0] );
		month = castToInt( args[1] );
		day = castToInt( args[2] );
	}

	//time as string
	if (idx_time_str != IOUtils::npos) {
		//parse the time information and compute the decimal jdn
		const std::string time_str( vecFields[ idx_time_str ] );
		fractional_day = parseTime(time_str, tz);
	}
	
	//attempt to read all possible components
	if (idx_jdn!=IOUtils::npos && !parseDateComponent(vecFields, idx_jdn, jdn)) return Date();
	if (idx_year!=IOUtils::npos && !parseDateComponent(vecFields, idx_year, year)) return Date();
	if (idx_month!=IOUtils::npos && !parseDateComponent(vecFields, idx_month, month)) return Date();
	if (idx_day!=IOUtils::npos && !parseDateComponent(vecFields, idx_day, day)) return Date();
	if (idx_hours!=IOUtils::npos && !parseDateComponent(vecFields, idx_hours, hour)) return Date();
	if (idx_minutes!=IOUtils::npos && !parseDateComponent(vecFields, idx_minutes, minute)) return Date();
	if (idx_seconds!=IOUtils::npos && !parseDateComponent(vecFields, idx_seconds, seconds)) return Date();
	if (idx_ntime!=IOUtils::npos && !parseDateComponent(vecFields, idx_ntime, ntime)) return Date();
	
	//special handling of year: fixed year provided by the user, 2 digits year
	if (year==IOUtils::inodata && year_cst!=IOUtils::inodata) {
		if (jdn!=IOUtils::nodata) year = getFixedYear( jdn );
		else if (month!=IOUtils::nodata) year = getFixedYear( month );
	}
	if (year!=IOUtils::inodata && dt_2digits_year) {
		if (year < cutoff_year) year += 2000;
		else year += 1900;
	}
	
	//specaial handling of time: fixed hour provided by the user
	if (hour==IOUtils::inodata && hour_cst!=IOUtils::inodata) {
		hour = getFixedHour();
	}

	//special handling of numerical time
	if (ntime!=IOUtils::nodata) {
		hour = (int)( ntime / 100. );
		minute = (int)( ntime - hour*100. );
	}

	//now build a Date object
	if (jdn!=IOUtils::nodata) {
		//jdn only represents the date, other components represent the time
		if (hour!=IOUtils::inodata) {
			static const double seconds_to_days = 1. / (24.*3600.);
			return Date(year, static_cast<double>(jdn)+(hour*3600.+minute*60.+seconds)*seconds_to_days, tz);
		}
		//jdn only represents the date, otherwise we have a fractional day
		if (fractional_day!=IOUtils::nodata) {
			Date tmp(year, month, day, 0, 0, tz);
			tmp += fractional_day;
			return tmp;
		}
		
		//jdn represents both date and time
		return Date(year, static_cast<double>(jdn), tz);
	}
	
	//we already have a fractional_day, so use it
	if (fractional_day!=IOUtils::nodata) {
		Date tmp(year, month, day, 0, 0, tz);
		tmp += fractional_day;
		return tmp;
	}
	
	return Date(year, month, day, hour, minute, seconds, tz);
}

std::string CsvDateTime::toString() const 
{
	std::ostringstream os;
	os << "[";
	if (idx_decimal_date!=IOUtils::npos) os << "idx_decimal_date→" << idx_decimal_date << " ";
	if (idx_date_time_str!=IOUtils::npos) os << "idx_date_time_str→" << idx_date_time_str << " ";
	if (idx_date_str!=IOUtils::npos) os << "idx_date_str→" << idx_date_str << " ";
	if (idx_time_str!=IOUtils::npos) os << "idx_time_str→" << idx_time_str << " ";
	if (idx_year!=IOUtils::npos) os << "idx_year→" << idx_year << " ";
	if (dt_2digits_year) os << "(cutoff_year=" << cutoff_year << ") ";
	if (year_cst!=IOUtils::nodata) os << "year_cst→" << year_cst << " ";
	if (hour_cst!=IOUtils::nodata) os << "hour_cst→" << hour_cst << " ";
	if (idx_jdn!=IOUtils::npos) os << "idx_jdn→" << idx_jdn << " ";
	if (idx_month!=IOUtils::npos) os << "idx_month→" << idx_month << " ";
	if (idx_day!=IOUtils::npos) os << "idx_day→" << idx_day << " ";
	if (idx_ntime!=IOUtils::npos) os << "idx_ntime→" << idx_ntime << " ";
	if (idx_hours!=IOUtils::npos) os << "idx_hours→" << idx_hours << " ";
	if (idx_minutes!=IOUtils::npos) os << "idx_minutes→" << idx_minutes << " ";
	if (idx_seconds!=IOUtils::npos) os << "idx_seconds→" << idx_seconds << " ";
	os << "tz→" << csv_tz << " ";
	if (has_tz) os << "tz_in_data ";
	if (auto_wrap) os << "auto_wrap";
	
	os << "]";
	return os.str();
}

///////////////////////////////////////////////////// Start of the CsvParameters class //////////////////////////////////////////

CsvParameters::CsvParameters(const double& tz_in) : csv_fields(), units_offset(), units_multiplier(), field_offset(), field_multiplier(), header_repeat_mk(), filter_ID(), fields_postfix(), ID_col(IOUtils::npos), header_lines(1), columns_headers(IOUtils::npos), units_headers(IOUtils::npos), csv_delim(','), header_delim(','), eoln('\n'), comments_mk('\n'), header_repeat_at_start(false), asc_order(true), number_fields(false), date_cols(tz_in), location(), nodata(), skip_fields(), purgeCharsSet(), linesExclusions(), file_and_path(), single_field(), name(), id(), start_hint(), end_hint(), slope(IOUtils::nodata), azi(IOUtils::nodata), exclusion_idx(0), exclusion_last_linenr(0), last_allowed_field(IOUtils::npos)
{
	//prepare default values for the nodata markers
	setNodata( "NAN NULL" );
}

std::string CsvParameters::toString() const
{
	std::ostringstream os;
	os << "<CsvParameters>\n";
	os << "\t" << file_and_path << " - " << name << " - " << id << " - " << location.toString(Coords::FULL) << " - " << date_cols.toString() << "\n";
	os << "</CsvParameters>\n";
	return os.str();
}

//parse the user provided special headers specification. It is stored as <line_nr, <column, field_type>> in a multimap
//(since there can be multiple keys on the same line)
std::multimap< size_t, std::pair<size_t, std::string> > CsvParameters::parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec)
{
	std::multimap< size_t, std::pair<size_t, std::string> > meta_spec;
	for (const auto& metaSpec : vecMetaSpec) {
		std::vector<std::string> vecArgs;
		if (IOUtils::readLineToVec(metaSpec, vecArgs, ':') !=3)
			throw InvalidFormatException("Wrong format for Metadata specification '"+metaSpec+"'", AT);
		const int linenr = atoi( vecArgs[1].c_str() );
		const int colnr = atoi( vecArgs[2].c_str() );
		if (linenr<=0 || colnr<=0)
			throw InvalidFormatException("Line numbers and column numbers must be >0 in Metadata specification '"+metaSpec+"'", AT);
		
		meta_spec.insert( make_pair( linenr, make_pair( colnr, vecArgs[0]) ) );
	}
	
	return meta_spec;
}

//Given a list of fields to skip, fill the skip_fields set
void CsvParameters::setSkipFields(const std::string& skipFieldSpecs, const bool& negate)
{
	const std::string where = (negate)? "INPUT::CSV#_ONLY_FIELDS" : "INPUT::CSV#_SKIP_FIELDS";
	//HACK temportarily look for old, space delimited syntax
	static const std::regex old_syntax_regex("[^;|#]*[0-9]+(\\s+)[0-9]+.*", std::regex::optimize); //space delimited list of ints
	if (std::regex_match(skipFieldSpecs, old_syntax_regex))
		throw InvalidArgumentException("Using old, space delimited list for " + where + ". It should now be a comma delimited list (ranges are also supported)", AT);

	const std::vector< LinesRange > fieldRange( IOInterface::initLinesRestrictions(skipFieldSpecs, where, negate) );
	
	//keep single fields as such, enumerate ranges so "1, 12-15" will generate "1 12 13 14 15"
	for (const auto& skipField : fieldRange) {
		if (skipField.end==static_cast<size_t>(-1)) { //convert an open-ended skip to last_allowed_field 
			last_allowed_field = skipField.start - 2;	//last valid position is -1, real idx start at 0, so -1 again
			continue;
		}

		for (size_t ii=skipField.start; ii<=skipField.end; ii++)
			skip_fields.insert( ii-1 );
	}
}

void CsvParameters::setDelimiter(const std::string& delim)
{
	if (delim.size()==1) {
		csv_delim = delim[0];
	} else {
		if (delim.compare("SPACE")==0 || delim.compare("TAB")==0) 
			csv_delim=' ';
		else 
			throw InvalidArgumentException("The CSV delimiter must be a single character or SPACE or TAB", AT);
	}
}

void CsvParameters::setHeaderDelimiter(const std::string& delim)
{
	if (delim.size()==1) {
		header_delim = delim[0];
	} else {
		if (delim.compare("SPACE")==0 || delim.compare("TAB")==0)
			header_delim=' ';
		else
			throw InvalidArgumentException("The CSV header delimiter must be a single character or SPACE or TAB", AT);
	}
}

std::string CsvParameters::identifyField(const std::string& fieldname)
{
	if (fieldname.compare(0, 15, "TEMPERATURE_AIR")==0 || fieldname.compare(0, 7, "AIRTEMP")==0 || fieldname.compare(0, 16, "TEMPERATURA_ARIA")==0) return "TA";
	else if (fieldname.compare(0, 16, "SOIL_TEMPERATURE")==0 || fieldname.compare(0, 8, "SOILTEMP")==0) return "TSG";
	else if (fieldname.compare(0, 13, "PRECIPITATION")==0 || fieldname.compare(0, 4, "PREC")==0 || fieldname.compare(0, 14, "PRECIPITAZIONE")==0) return "PSUM";
	else if (fieldname.compare(0, 19, "REFLECTED_RADIATION")==0 || fieldname.compare(0, 26, "RADIAZIONE_SOLARE_RIFLESSA")==0) return "RSWR";
	else if (fieldname.compare(0, 18, "INCOMING_RADIATION")==0 || fieldname.compare(0, 26, "INCOMINGSHORTWAVERADIATION")==0 || fieldname.compare(0, 27, "RADIAZIONE_SOLARE_INCIDENTE")==0) return "RSWR";
	else if (fieldname.compare(0, 14, "WIND_DIRECTION")==0 || fieldname.compare(0, 2, "WD")==0 || fieldname.compare(0, 15, "DIREZIONE_VENTO")==0) return "DW";
	else if (fieldname.compare(0, 17, "RELATIVE_HUMIDITY")==0 || fieldname.compare(0, 16, "RELATIVEHUMIDITY")==0 || fieldname.compare(0, 15, "UMIDIT_RELATIVA")==0) return "RH";
	else if (fieldname.compare(0, 13, "WIND_VELOCITY")==0 || fieldname.compare(0, 2, "WS")==0 || fieldname.compare(0, 13, "VELOCIT_VENTO")==0) return "VW";
	else if (fieldname.compare(0, 8, "PRESSURE")==0 || fieldname.compare(0, 15, "STATIONPRESSURE")==0) return "P";
	else if (fieldname.compare(0, 17, "INCOMING_LONGWAVE")==0 || fieldname.compare(0, 25, "INCOMINGLONGWAVERADIATION")==0) return "ILWR";
	else if (fieldname.compare(0, 22, "SNOWSURFACETEMPERATURE")==0 ) return "TSS";
	else if (fieldname.compare(0, 6, "WS_MAX")==0) return "VW_MAX";
	
	return fieldname;
}

//Given a provided field_type, attribute the value to the proper metadata variable.
void CsvParameters::assignMetadataVariable(const std::string& field_type, const std::string& field_val, double &lat, double &lon, double &easting, double &northing)
{
	if (field_type=="ID") {
		if (id.empty()) id = field_val;
	} else if (field_type=="NAME") {
		if (name.empty()) name = field_val;
	} else if (field_type=="NODATA") {
		setNodata( field_val );
	} else if (field_type=="SKIP") {
		return;
	} else if (field_type=="PARAM") {
		std::string param( IOUtils::strToUpper( field_val ) );
		if (MeteoData::getStaticParameterIndex( param )!=IOUtils::npos) {
			single_field = param;
			return;
		}
		
		IOUtils::replaceInvalidChars(param); //remove accentuated characters, etc
		param = identifyField( param ); //try to map non-standard names to mio's names
		
		single_field = param;
	} else {
		if (field_type=="ALT" || field_type=="LON" || field_type=="LAT" || field_type=="SLOPE" || field_type=="AZI" || field_type=="EASTING" || field_type=="NORTHING") {
			double tmp;
			if (!IOUtils::convertString(tmp, field_val))
				throw InvalidArgumentException("Could not extract metadata '"+field_type+"' for "+file_and_path, AT);
			
			if (field_type=="ALT") location.setAltitude( tmp, false);
			if (field_type=="LON") lon = tmp;
			if (field_type=="LAT") lat = tmp;
			if (field_type=="SLOPE") slope = tmp;
			if (field_type=="AZI") azi = tmp;
			if (field_type=="EASTING") easting = tmp;
			if (field_type=="NORTHING") northing = tmp;
		} else 
			throw InvalidFormatException("Unknown parsing key '"+field_type+"' when extracting metadata", AT);
	}
}

//Using the special headers parsed specification (done in parseHeadersSpecs), some metadata is extracted from the headers
void CsvParameters::parseSpecialHeaders(const std::string& line, const size_t& linenr, const std::multimap< size_t, std::pair<size_t, std::string> >& meta_spec, double &lat, double &lon, double &easting, double &northing)
{
	std::vector<std::string> vecStr;
	IOUtils::readLineToVec(line, vecStr, header_delim);
	
	const bool readID = (id.empty()); //if the user defined CSV_ID, it has priority
	const bool readName = (name.empty()); //if the user defined CSV_NAME, it has priority
	std::string prev_ID, prev_NAME;
	std::multimap<size_t, std::pair<size_t, std::string> >::const_iterator it;
	for (it=meta_spec.equal_range(linenr).first; it!=meta_spec.equal_range(linenr).second; ++it) {
		const size_t colnr = (*it).second.first;
		const std::string field_type( IOUtils::strToUpper( (*it).second.second) );
		if (colnr>vecStr.size() || colnr==0)
			throw InvalidArgumentException("Metadata specification for '"+field_type+"' refers to a non-existent field for file (either 0 or too large) "+file_and_path, AT);
		
		//remove the quotes from the field
		std::string field_val( vecStr[colnr-1] );
		IOUtils::removeQuotes(field_val);
		
		//we handle ID and NAME differently in order to support appending
		if (field_type=="ID" && readID) {
			id = prev_ID+field_val;
			prev_ID = id+"-";
		} else if (field_type=="NAME" && readName) {
			name = prev_NAME+field_val;
			prev_NAME = name+"-";
		} else {
			assignMetadataVariable(field_type, field_val, lat, lon, easting, northing);
		}
	}
}

//Extract metadata from the filename, according to a user-provided specification.
//filename_spec in the form of {ID}_{NAME}-{PARAM}_-_{SKIP} where {ID} is a variable and '_-_' is a constant pattern
void CsvParameters::parseFileName(std::string filename, const std::string& filename_spec, double &lat, double &lon, double &easting, double &northing)
{
	filename = FileUtils::removeExtension( FileUtils::getFilename(filename) );
	size_t pos_fn = 0, pos_mt = 0; //current position in the filename and in the filename_spec
	if (filename_spec[0]!='{') { //there is a constant pattern at the beginning, getting rid of it
		const size_t start_var = filename_spec.find('{');
		if (start_var==std::string::npos) throw InvalidFormatException("No variables defined for filename parsing", AT);
		const std::string pattern( filename_spec.substr(0, start_var) );
		if (filename.substr(0, start_var)!=pattern) throw InvalidFormatException("The filename pattern '"+filename_spec+"' does not match with the given filename ('"+filename+"') for metadata extraction", AT);
		pos_mt = start_var;
		pos_fn = start_var;
	}

	const bool readID = (id.empty()); //if the user defined CSV_ID, it has priority
	const bool readName = (name.empty()); //if the user defined CSV_NAME, it has priority
	std::string prev_ID, prev_NAME;
	//we now assume that we start with a variable
	do { //HACK redo with regex
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
			if (pos_pattern_fn==std::string::npos) throw InvalidFormatException("The filename pattern '"+filename_spec+"' does not match with the given filename ('"+filename+"') for metadata extraction", AT);
			len_var = pos_pattern_fn - pos_fn;
		}

		//read the variable type and value
		const std::string field_type( IOUtils::strToUpper(filename_spec.substr(pos_mt+1, start_pattern-pos_mt-1)) ); //skip { and }
		const std::string value( filename.substr(pos_fn, len_var) );
		//we handle ID and NAME differently in order to support appending
		if (field_type=="ID" && readID) {
			id = prev_ID+value;
			prev_ID = id+"-";
		} else if (field_type=="NAME" && readName) {
			name = prev_NAME+value;
			prev_NAME = name+"-";
		} else {
			assignMetadataVariable(field_type, value, lat, lon, easting, northing);
		}

		if (end_pattern==std::string::npos) break; //nothing more to parse
		pos_mt = end_pattern;
		pos_fn = pos_fn + len_var + pattern_len;
	} while (true);

}

//very basic units parsing: a few hard-coded units are recognized and provide the necessary
//offset and multiplier to convert the values back to SI
void CsvParameters::setUnits(const std::string& csv_units, const char& delim)
{
	static const std::set<std::string> noConvUnits = {"TS", "S", "RN", "W/M2", "M/S", "K", "M", "N", "PA", "V", "VOLT", "DEG", "°", "KG/M2", "KG/M3", "DB", "DBM"};

	std::vector<std::string> units;
	if (delim!=' ') 
		IOUtils::readLineToVec(csv_units, units, delim);
	else 
		IOUtils::readLineToVec(csv_units, units);
	units_offset.resize(units.size(), 0.);
	units_multiplier.resize(units.size(), 1.);
	
	for (size_t ii=0; ii<units.size(); ii++) {
		std::string tmp( units[ii] );
		IOUtils::toUpper( tmp );
		IOUtils::removeQuotes(tmp);
		if (tmp.empty() || tmp=="1" || tmp=="-" || tmp=="0 OR 1" || tmp=="0/1" || tmp=="??") continue; //empty unit
		if (noConvUnits.count(tmp)>0) continue; //this unit does not need conversion
		
		if (tmp=="%" || tmp=="PC" || tmp=="CM") units_multiplier[ii] = 0.01;
		else if (tmp=="‰") units_multiplier[ii] = 0.001;
		else if (tmp=="C" || tmp=="DEGC" || tmp=="GRAD C" || tmp=="°C") units_offset[ii] = Cst::t_water_freezing_pt;
		else if (tmp=="HPA") units_multiplier[ii] = 1e2;
		else if (tmp=="MM" || tmp=="MV" || tmp=="MA" || tmp=="MS") units_multiplier[ii] = 1e-3;
		else if (tmp=="MIN") units_multiplier[ii] = 60.;
		else if (tmp=="IN") units_multiplier[ii] = 0.0254;
		else if (tmp=="FT") units_multiplier[ii] = 0.3048;
		else if (tmp=="F") { units_multiplier[ii] = 5./9.; units_offset[ii] = -32.*5./9.;}
		else if (tmp=="KM/H") units_multiplier[ii] = 1./3.6;
		else if (tmp=="MPH") units_multiplier[ii] = 1.60934 / 3.6;
		else if (tmp=="KT") units_multiplier[ii] = 1.852 / 3.6;
		else {
			std::cerr << "CsvIO: Can not parse unit '" << tmp << "', please inform the MeteoIO developers\n";
		}
	}
}

void CsvParameters::setNodata(const std::string& nodata_markers)
{
	//by NOT calling clear(), we append the provided markers to potential previous ones
	//(such as default values)
	
	std::vector<std::string> vecNodata;
	IOUtils::readLineToVec(nodata_markers, vecNodata);
	
	for (const auto& marker : vecNodata) {
		nodata.insert( marker );
		nodata.insert( "\""+marker+"\"" );
		nodata.insert( "'"+marker+"'" );
	}
}

void CsvParameters::setPurgeChars(const std::string& chars_to_purge)
{
	std::vector<std::string> vecPurge;
	IOUtils::readLineToVec(chars_to_purge, vecPurge);
	char rest[32] = "";
	
	for (const auto& purgeChar : vecPurge) {
		unsigned int c;
		const char* c_str( purgeChar.c_str() );
		
		if (purgeChar.length()==1) { //the character has been directly given
			purgeCharsSet.insert( purgeChar[0] );
		} else if (sscanf(c_str, "0x%2x%31s", &c, rest) == 1) { //hexadecimal
			purgeCharsSet.insert( static_cast<char>(c) );
		} else if (sscanf(c_str, "%3u%31s", &c, rest) == 1) { //hexadecimal
			purgeCharsSet.insert( static_cast<char>(c) );
		} else {
			throw InvalidArgumentException("Invalid purge chars specification '"+std::string(c_str)+"' for the CSV plugin, please use either single chars or decimal notation", AT);
		}
	}
}

/** @brief parse a simplified ISO date and build a Date object
 * @param[in] Date_str ths string to parse for an ISO date
 * @param[in] early_interpretation interpret missing date components toward the earliest possible date (if set to true), otherwise toward the latest possible date
 * @return Date object
 */
Date CsvParameters::parseDateHint(const std::string& Date_str, const double tz, const bool& early_interpretation)
{
	static const std::regex ISO_date_regex("^([0-9]{4})(?:-([0-9]{2}))?(?:-([0-9]{2}))?$", std::regex::optimize); //either only year, or year and month or year, month, day
	std::smatch match;

	if (!std::regex_match(Date_str, match, ISO_date_regex))
		throw InvalidArgumentException("Invalid temporal coverage hint '"+Date_str+"' for the CSV plugin, please use proper ISO notation", AT);

	int year=IOUtils::inodata, month=IOUtils::inodata, day=IOUtils::inodata;
	if (!IOUtils::convertString(year, match.str(1)))
		throw InvalidArgumentException("Could not parse temporal coverage hint '"+Date_str+"' for the CSV plugin", AT);

	//since we can not easily count the number of groups that have effectively been captured, we need to check their content
	if (!match.str(2).empty()) {
		if (!IOUtils::convertString(month, match.str(2)))
			throw InvalidArgumentException("Could not parse temporal coverage hint '"+Date_str+"' for the CSV plugin", AT);
		if (!match.str(3).empty()) {
			if (!IOUtils::convertString(day, match.str(3)))
				throw InvalidArgumentException("Could not parse temporal coverage hint '"+Date_str+"' for the CSV plugin", AT);
		}
	}

	if (early_interpretation) {
		if (month==IOUtils::inodata) month = 1;
		if (day==IOUtils::inodata) day = 1;
		return Date(year, month, day, 0, 0, 0., tz);
	} else {
		if (month==IOUtils::inodata) month = 12;
		if (day==IOUtils::inodata) day = Date::getNumberOfDaysInMonth(year, month);
		return Date(year, month, day, 23, 59, 59.999, tz);
	}
}

/** @brief Set the date range from a provided date range string
 * @details Such a string can either contain only one date (at least the year must be provided, everything else is optional)
 * or two dates separated by a dash. Examples of valid ranges:
 *  * 2024 - 2025 (this means from 2024-01-01 until 2025-12-31 at midnight)
 *  * 2024-10 (this means the whole month of October 2024)
 *  * 2024-10 - 2025 (this means from 2024-10-01 until 2025-12-31)
 * 
 * @param[in] range_spec the range string
 */
void CsvParameters::setCoverageHint(const std::string& range_spec)
{
	const double tz_in = date_cols.csv_tz;
	static const std::regex ISO_range_regex("^([0-9]{4}(?:\\-[0-1][0-9](?:\\-[0-3][0-9])?)?)(?:\\s+\\-\\s+([0-9]{4}(?:\\-[0-1][0-9](?:\\-[0-3][0-9])?)?)?)?$", std::regex::optimize); //2 capturing group: we either have a single date or a range
	std::smatch match;

	//do we have two dates or only one?
	if (!std::regex_match(range_spec, match, ISO_range_regex))
	 	throw InvalidArgumentException("Invalid temporal coverage hint format '"+range_spec+"' for the CSV plugin, please use proper ISO notation with daily resolution", AT);

	//since we can not easily count the number of groups that have effectively been captured, we need to check their content
	if (!match.str(1).empty()) start_hint = parseDateHint(match.str(1), tz_in, true);
	if (!match.str(2).empty()) 
		end_hint = parseDateHint(match.str(2), tz_in, false);
	else
		end_hint = parseDateHint(match.str(1), tz_in, false);

	if (start_hint.isUndef() || end_hint.isUndef()) return;
	if (start_hint>=end_hint) return;
}

bool CsvParameters::hasDates(const Date& start, const Date& end) const 
{
	if (start_hint.isUndef()) return true; 
	const bool overlap = (start_hint<=end && end_hint>=start);
	return overlap;
}

bool CsvParameters::excludeLine(const size_t& linenr, bool& hasExclusions)
{
	//As an optimimzation, we reuse the exclusion periods index over calls.
	//as long as line numbers are always increasing, this is OK.
	//So we check that this is the case with exclusion_last_linenr
	if (linenr<exclusion_last_linenr) exclusion_idx=0;
	exclusion_last_linenr = linenr;

	if (linesExclusions.empty() || linenr>linesExclusions.back().end) {
		hasExclusions = false;
		return false; //no more exclusions to handle
	}
	
	for (; exclusion_idx<linesExclusions.size(); exclusion_idx++) {
		if (linesExclusions[ exclusion_idx ].in( linenr ))
			return true;
		if (linenr<linesExclusions[ exclusion_idx ].start)
			break;
	}
	
	return false;
}

bool CsvParameters::skipField(const size_t& fieldnr) const
{
	if (skip_fields.count(fieldnr)>0) return true;
	if (last_allowed_field!=IOUtils::npos && fieldnr>last_allowed_field) return true;
	
	return false;
}

bool CsvParameters::isNodata(const std::string& value) const
{
	if (value.empty()) return true;
	if (nodata.find( value ) != nodata.end()) return true;
	
	return false;
}

void CsvParameters::setDateTimeSpecs(const std::string &datetime_spec, const std::string &date_spec, const std::string &time_spec, const std::string &decimaldate_type)
{
	if (!decimaldate_type.empty() && (!datetime_spec.empty() || !date_spec.empty() || !time_spec.empty() ))
		throw InvalidArgumentException("It is not possible to define both decimaldate_type and other date / time specifications", AT);
	if (!datetime_spec.empty() && (!date_spec.empty() || !time_spec.empty()) )
		throw InvalidArgumentException("It is not possible to define both datetime_spec and date_spec or time_spec", AT);
	
	if (!datetime_spec.empty()) date_cols.setDateTimeSpec( datetime_spec );
	if (!date_spec.empty()) date_cols.setDateSpec( date_spec );
	if (!time_spec.empty()) date_cols.setTimeSpec( time_spec );
	if (!decimaldate_type.empty()) date_cols.setDecimalDateType( decimaldate_type );
}

//user provided field names are in fieldNames, header field names are in headerFields
//and user provided fields have priority (headerFields only used if fieldNames is empty).
void CsvParameters::parseFields(const std::vector<std::string>& headerFields, std::vector<std::string>& fieldNames)
{
	const bool user_provided_field_names = (!fieldNames.empty());
	if (headerFields.empty() && !user_provided_field_names) {
		if (single_field.empty())
			throw InvalidArgumentException("No columns names could be found. Please either provide CSV_COLUMNS_HEADERS or CSV_FIELDS (or a PARAM metadata)", AT);
	}
	if (!user_provided_field_names) fieldNames = headerFields;
	
	std::map<std::string, size_t> data_fields;	// this is only used to prevent multiple use of the same name and handle single_field if necessary
	bool single_field_found = false;
	const size_t nFields = fieldNames.size();
	const int fwidth = static_cast<int>( ceil( log10( nFields ) ) ); //required if number_fields==true
	
	for (size_t ii=0; ii<nFields; ii++) {
		//if this has been given in the list of indices to skip by the user, don't even try to read the field name,
		//so there won't be an error if the same name is used multiple times but skipped
		if (skipField(ii)) continue;

		std::string &tmp = fieldNames[ii];
		IOUtils::trim( tmp ); //there could still be leading/trailing whitespaces in the individual field name
		IOUtils::toUpper( tmp );
		IOUtils::removeDuplicateWhitespaces(tmp); //replace internal spaces by '_'
		IOUtils::replaceWhitespaces(tmp, '_');
		if (tmp.empty()) continue;
		
		if (tmp.compare("PARAM")==0) {
			tmp = single_field;
			if (single_field_found)
				throw InvalidArgumentException("It is not possible to have more than one PARAM field!", AT);
			single_field_found = true;
		} else if (date_cols.parseField(tmp, ii)) { //this is a date/time component
			skip_fields.insert( ii );
		} else if (tmp.compare("SKIP")==0 || tmp.compare("-")==0) {
			skip_fields.insert( ii );
		} else if (tmp.compare("ID")==0 || tmp.compare("STATIONID")==0) {
			ID_col = ii;
			skip_fields.insert( ii );
		} else {
			//number the field names if so requested by the user
			if (number_fields) {
				std::ostringstream os;
				os << std::setw(fwidth) << std::setfill('0') << ii+1 << "_" << tmp << fields_postfix;
				tmp = os.str();
			}
			
			if (data_fields.count( tmp ) > 0)
				throw InvalidArgumentException("Multiple definitions of the same field name ('"+tmp+"') either in column headers or user-provided CSV_FIELDS", AT);
			data_fields[ tmp ] = ii;
		}
		
		//tmp = identifyField( tmp ); //try to identify known fields
	}
	date_cols.updateMaxCol();
	
	//check for time handling consistency
	if (!date_cols.isSet()) throw UnknownValueException("Please fully define how to parse the date and time information. Check that all date/time data is available (including CSV_FALLBACK_YEAR if not year information is provided in the file)!", AT);

	//the user wants to keep only one column, find the one he wants...
	//if there is a parameter name from the filename or header it has priority:
	if (!single_field.empty() && data_fields.size()==1 && !single_field_found)
		fieldNames[ data_fields.begin()->second ] = data_fields.begin()->first;
}

//read and parse the file's headers in order to extract all possible information (including how to interpret the date/time information)
void CsvParameters::setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec, const std::string& filename_spec, const std::string& station_idx)
{
	file_and_path = i_file_and_path;
	const std::multimap< size_t, std::pair<size_t, std::string> > meta_spec( parseHeadersSpecs(vecMetaSpec) );
	double lat = IOUtils::nodata;
	double lon = IOUtils::nodata;
	double easting = IOUtils::nodata, northing = IOUtils::nodata;
	
	//read and parse the file's headers
	if (!FileUtils::fileExists(file_and_path)) throw AccessException("File "+file_and_path+" does not exists", AT); //prevent invalid filenames
	errno = 0;
	if (!filename_spec.empty()) parseFileName( file_and_path, filename_spec, lat, lon, easting, northing);
	std::ifstream fin(file_and_path.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		const std::string msg( "Error opening file " + file_and_path + " for reading, possible reason: " + std::strerror(errno) + " Please check file existence and permissions!" );
		throw AccessException(msg, AT);
	}
	
	const bool user_auto_wrap = date_cols.auto_wrap; //we might trigger it, so we must reset it after the pre-reading
	const bool read_units = (units_headers!=IOUtils::npos && units_offset.empty() && units_multiplier.empty());
	size_t linenr=0, repeat_markers=0;
	std::string line;
	std::vector<std::string> headerFields; //this contains the column headers from the file itself
	std::vector<std::string> tmp_vec; //to read a few lines of data
	Date prev_dt;
	size_t count_asc=0, count_dsc=0; //count how many ascending/descending timestamps are present
	static const size_t min_valid_lines = 10; //we want to correctly parse at least that many lines before quitting our sneak peek into the file
	const bool delimIsNoWS = (csv_delim!=' ');
	const bool hasHeaderRepeatMk = (!header_repeat_mk.empty());
	bool fields_ready = false;

	try {
		eoln = FileUtils::getEoln(fin);
		for (size_t ii=0; ii<(header_lines+1000); ii++) {
			getline(fin, line, eoln); //read complete line
			linenr++;
			
			IOUtils::trim(line);
			if (fin.eof()) {
				if ((linenr-repeat_markers)>header_lines) break; //eof while reading the data section
				std::ostringstream ss;
				ss << "Declaring " << header_lines << " header line(s) for file " << file_and_path << ", but it only contains " << linenr << " lines";
				throw InvalidArgumentException(ss.str(), AT);
			}
			
			if (hasHeaderRepeatMk && !header_repeat_at_start && line.find(header_repeat_mk)!=std::string::npos) {
				header_repeat_at_start = true; //so we won't match another header_repeat_mk marker
				repeat_markers++;
				continue;
			}
			
			if (comments_mk!='\n') IOUtils::stripComments(line, comments_mk);
			if (line.empty()) continue;
			if (*line.rbegin()=='\r') line.erase(line.end()-1); //getline() skipped \n, so \r comes in last position
			
			if (meta_spec.count(linenr-repeat_markers)>0) 
				parseSpecialHeaders(line, linenr-repeat_markers, meta_spec, lat, lon, easting, northing); //do not count repeat_markers if any
			if ((linenr-repeat_markers)==columns_headers) { //so user provided csv_fields have priority. If columns_headers==npos, this will also never be true
				if (delimIsNoWS) { //even if header_delim is set, we expect the fields to be separated by csv_delim
					IOUtils::cleanFieldName(line, false); //we'll handle whitespaces when parsing
					IOUtils::readLineToVec(line, headerFields, csv_delim);
				} else {
					IOUtils::cleanFieldName(line, false); //don't touch whitespaces
					IOUtils::readLineToVec(line, headerFields);
				}
			}
			if (read_units && (linenr-repeat_markers)==units_headers)
				setUnits(line, csv_delim);

			if ((linenr-repeat_markers)<=header_lines) continue; //we are still parsing the header
			if (!fields_ready) { //we should now have all the information from the headers, so build what we need for data parsing
				parseFields(headerFields, csv_fields);
				fields_ready=true;
				continue;
			}
			
			const size_t nr_curr_data_fields = (delimIsNoWS)? IOUtils::readLineToVec(line, tmp_vec, csv_delim) : IOUtils::readLineToVec(line, tmp_vec);
			if (nr_curr_data_fields>date_cols.max_dt_col) {
				const Date dt( date_cols.parseDate(tmp_vec) );
				if (dt.isUndef()) continue;
				if (!prev_dt.isUndef()) {
					if (dt>prev_dt) count_asc++;
					else count_dsc++;
				}
				prev_dt = dt;
			}
			if (count_asc+count_dsc >= min_valid_lines) break; //we've add enough valid lines to understand the file, quitting
		}
	} catch (...) {
		fin.close();
		throw;
	}
	fin.close();
	date_cols.auto_wrap = user_auto_wrap; //resetting it since we might have triggered it
	if (count_dsc>count_asc) asc_order=false;
	
	if (lat!=IOUtils::nodata || lon!=IOUtils::nodata) {
		const double alt = location.getAltitude(); //so we don't change previously set altitude
		location.setLatLon(lat, lon, alt, false); //we let Coords handle possible missing data / wrong values, etc
	}
	if (easting!=IOUtils::nodata || northing!=IOUtils::nodata) {
		const double alt = location.getAltitude();
		location.setXY(easting, northing, alt, false); //coord system was set on keyword parsing
	}
	//location is either coming from POSITIONxx ini keys or from file name parsing or from header parsing
	if (location.isNodata()) 
		throw NoDataException("Missing geographic coordinates for '"+i_file_and_path+"', please consider providing the POSITION ini key", AT);
	location.check("Inconsistent geographic coordinates in file \"" + file_and_path + "\": ");
	
	if (name.empty()) name = FileUtils::removeExtension( FileUtils::getFilename(i_file_and_path) ); //fallback if nothing else could be find
	if (id.empty()) {
		if (station_idx.empty()) 
			id = name; //really nothing, copy "name"
		else
			id = "ID"+station_idx; //automatic numbering of default IDs
	}
}

StationData CsvParameters::getStation() const 
{
	StationData sd(location, id, name);
	if (slope==0. || (slope!=IOUtils::nodata && azi!=IOUtils::nodata))
		sd.setSlope(slope, azi);
	return sd;
}

} //namespace
