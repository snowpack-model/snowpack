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
#include <meteoio/plugins/CsvIO.h>
#include <meteoio/plugins/plugin_utils.h>

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <utility>
#include <cerrno>

using namespace std;

namespace mio {
	using namespace PLUGIN;
/**
 * @page csvio CsvIO
 * @section csvio_format Format
 * This plugins offers a flexible way to read Comma Separated Values (<A HREF="https://en.wikipedia.org/wiki/Comma-separated_values">CSV</A>) files. 
 * It is however assumed that:
 *     - each line contains a data record (or is an empty line)
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
 * **The final units MUST be <a href="https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf">coherent derived SI units</a>** 
 * (section 2.3.4 in the SI-Brochure). If not, the conversion offsets/factors must be provided to convert the data back to SI (see required keywords below)
 * or the units declared (in the headers) and supported by this plugin.
 *
 * @section csvio_keywords Keywords
 * This plugin uses the keywords described below, in the [Input] section. First, there are some general plugin options:
 * - COORDSYS: coordinate system (see Coords);
 * - COORDPARAM: extra coordinates parameters (see Coords);
 * - TIME_ZONE: the timezone that should be used to interpret the dates/times (default: 0);
 * - METEOPATH: the directory where the data files are available (mandatory);
 * - METEOPATH_RECURSIVE: if set to true, the scanning of METEOPATH is performed recursively (default: false);
 * - CSV_FILE_EXTENSION: When scanning the whole directory, look for these files (default: .csv). Note that this matching isn't restricted to the end of the file name so if you had files stat1_jan.csv, stat1_feb.csv and stat2_jan.csv you could select January's data by putting "_jan" here;
 * - CSV_SILENT_ERRORS: if set to true, lines that can not be read will be silently ignored (default: false, has priority over CSV_ERRORS_TO_NODATA);
 * - CSV_ERRORS_TO_NODATA: if true, unparseable fields (like text fields) are set to nodata, but the rest of the line is kept (default: false).
 * 
 * You can now describe the specific format for all files (prefixing the following keys by \em "CSV_") or for each particular file (prefixing the following 
 * keys by \em "CSV#_" where \em "#" represents the station index). Of course, you can mix keys that are defined for all files with some keys only defined for a 
 * few specific files (the keys defined for a particular station have priority over the global version).
 * - CSV\#_DELIMITER: field delimiter to use (default: ','), use SPACE or TAB for whitespaces (in this case, multiple whitespaces directly following each other are considered to be only one whitespace);
 * - CSV\#_NODATA: a space delimited list of strings (of course, this also contains numbers such as -6999) that should be interpreted as \em nodata (default: NAN NULL);
 * - <b>Headers handling</b>
 *    - CSV\#_NR_HEADERS: how many lines should be treated as headers? (default: 1);
 *    - CSV\#_HEADER_DELIMITER: different field delimiter to use in header lines; optional
 *    - CSV\#_HEADER_REPEAT_MK: a string that is used to signal another copy of the headers mixed with the data in the file (the matching is done anywhere in the line) (default: empty);
 *    - CSV\#_UNITS_HEADER: header line providing the measurements units (the subset of recognized units is small, please inform us if one is missing for you); optional
 *    - CSV\#_UNITS_OFFSET: offset to add to each value in order to convert it to SI (if also providing CSV\#_UNITS, would be applied afterwards); optional
 *    - CSV\#_UNITS_MULTIPLIER: factor to multiply each value by, in order to convert it to SI (if also providing CSV\#_UNITS, would be applied afterwards); optional
 * - <b>Data parsing restrictions</b>
 *    - CSV\#_COMMENTS_MK: a single character to use as comments delimiter, everything after this char until the end of the line will be skipped (default: no comments);
 *    - CSV\#_PURGE_CHARS: space delimited list of ascii characters to purge from the input, either directly given or as decimal representation or as hexadecimal representation (prefixed with <i>0x</i>). Example: 0x40 13 \" ;
 *    - CSV\#_EXCLUDE_LINES: a comma delimited list of line ranges (numbers separated by a dash) or line numbers to exclude from parsing (ie the lines within these ranges will be read and discarded immediately). Example:  <i>18 - 36, 52, 55, 167 - 189</i>. Please note that it is not possible to mix CSV\#_EXCLUDE_LINES and CSV\#_ONLY_LINES and that additional spaces (for more clarity in the input, as in the provided example) can be used although they are not mandatory.
 *    - CSV\#_ONLY_LINES: a comma delimited list of line ranges (numbers separated by a dash enclosed in spaces) or line numbers to restrict the parsing to (ie the lines outside of these ranges will be read and discarded immediately). Example:  <i>18 - 36, 52, 55, 167 - 189</i>. Please note that it is not possible to mix CSV\#_EXCLUDE_LINES and CSV\#_ONLY_LINES.
 *    - CSV\#_COVERAGE_HINT: a simplified date range (at most daily resolution) fully encompassing the data contained in the file. This is a useful optimization when a dataset is made of many CSV files covering a very large temporal range and MeteoIO is used to generate subsets of the dataset on-demand. In this case, MeteoIO can skip some files if their temporal coverage does not overlap with the requested temporal period. Examples: <i>2021</i> (this means, all of 2021) or <i>2020-10 - 2023</i> (this means from 2020-10-01T00:00 until the last second of 2023).
 * - <b>Fields parsing</b>
 *    - CSV\#_COLUMNS_HEADERS: header line to interpret as columns headers (default: 1, see also \ref csvio_special_fields "special field names");
 *    - CSV\#_FIELDS: one line providing the columns headers (if they don't exist in the file or to overwrite them). If a field is declared as "ID" then only the lines that have the proper ID for the current station will be kept; if a field is declared as "SKIP" it will be skipped; otherwise date/time parsing fields are supported according to <i>Date/Time parsing</i> below (see also the \ref csvio_special_fields "special field names" for more); optional
 *    - CSV\#_FILTER_ID: if the data contains an "ID" column, which ID should be kept (all others will be rejected); default: station ID
 *    - CSV\#_UNITS: one line providing space delimited units for each column (including the timestamp), no units is represented as "-". This is an alternative to relying on a units line in the file itself or relying on units_offset / units_multiplier. Please keep in mind that the choice of recognized units is very limited... (C, degC, cm, in, ft, F, deg, pc, % and a few others). If CSV\#UNITS_OFFSET / MULTIPLIER were also provided, CSV\#_UNITS would be applied first.
 *    - CSV\#_SKIP_FIELDS: a comma-delimited list of fields to skip (first field is numbered 1, ranges such as <i>12 - 17</i> are supported as well). Keep in mind that when using parameters such as UNITS_OFFSET, the skipped field MUST be taken into consideration (since even if a field is skipped, it is still present in the file!); optional
 *    - CSV\#_ONLY_FIELDS: a comma-delimited list of fields to keep, all others will be skipped (so this is the opposite of CSV\#_SKIP_FIELDS. Please note that if using both CSV\#_SKIP_FIELDS and CSV\#_ONLY_FIELDS, the most restrictive interpretation will be applied: only fields that are included in the ONLY list and NOT in the SKIP list will be kept); optional
 *    - CSV\#_NUMBER_FIELDS: prefix every given field name by its column index in the original file (this is useful to "debug" CSV files when the columns' content don't match what was expected. Please note that special fields related to date/time, station ID or SKIP are left unchanged); optional
 *    - CSV\#_FIELDS_POSTFIX: postfix every given field name by the provided string (this is also to "debug" CSV files); optional
 * - <b>Date/Time parsing</b>. There are several possibilities: the date/time is provided as one or two strings; as a purely decimal number following a given representation; as each component as a separate column.
 *    - Date/Time as string(s):
 *       - CSV\#_DATETIME_SPEC: mixed date and time format specification (default is ISO_8601: YYYY-MM-DDTHH24:MI:SS);
 *       - CSV\#_DATE_SPEC: date format specification (default: YYYY_MM_DD);
 *       - CSV\#_TIME_SPEC: time format specification (default: HH24:MI:SS);
 *    - Date/Time decimal representation:
 *       - CSV\#_DECIMALDATE_TYPE: the numerical representation that is used, one of EXCEL, JULIAN, MJULIAN, MATLAB, RFC868 or UNIX (see \ref decimal_date_representation "decimal date representations");
 *    - Date/Time as separate components: 
 *       - the fields must be named (either from the headers or through the CSV\#_FIELDS key) as YEAR, YEAR_2DIGITS (only the last 2 digits of the year, numbers before 40 will be converted to years after 2000), JDAY (number of days since the beginning of the year), MONTH, DAY, NTIME (numerical representation of time, for example 952 for 09:52), HOURS, MINUTES, SECONDS (if minutes or seconds are missing, they will be assumed to be zero). See \ref csvio_special_fields "special field names" for accepted synonyms;
 *       - if/when no year component is provided, it is possible to define a fallback year with the CSV\#_FALLBACK_YEAR key;
 *       - when using CSV\#_FALLBACK_YEAR, it will by default assume that all data for times greater than 1st October that appear 
 * before data belonging to times before 1st of October are actually data from the year before. Please set CSV\#_FALLBACK_AUTO_WRAP to false if this is not desired.
 * - <b>Metadata</b>
 *    - CSV\#_NAME: a descriptive station name to use (if provided, has priority over the special headers);
 *    - CSV\#_ID: the (short) station id to use (if provided, has priority over the special headers);
 *    - CSV\#_SLOPE: the slope angle in degrees at the station (if providing a slope, also provide an azimuth);
 *    - CSV\#_AZIMUTH: the slope azimuth in degrees from North at the station ;
 *    - CSV\#_SPECIAL_HEADERS: description of how to extract more metadata out of the headers; optional
 *    - CSV\#_FILENAME_SPEC: pattern to parse the filename and extract metadata out of it; optional
 *    - The following two keys provide mandatory data for each station, therefore there is no "global" version and they must be defined:
 *       - METEOFILE\#: input filename (in METEOPATH). As many meteofiles as needed may be specified. If nothing is specified, the METEOPATH directory will be scanned for files with the extension specified in CSV_FILE_EXTENSION;
 *       - POSITION\#: coordinates of the station (default: reading key "POSITION", see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax). This key can only be omitted if lat/lon/altitude are provided in the file name or in the headers (see CSV\#_FILENAME_SPEC and CSV\#_SPECIAL_HEADERS);
 *
 * If no ID has been provided, an automatic station ID will be generated as "ID{n}" where *n* is the current station's index. Regarding the units handling, 
 * it is only performed through either the CSV_UNITS_OFFSET key or the CSV_UNITS_OFFSET / CSV_UNITS_MULTIPLIER keys. These keys expect a value for each
 * column of the file, including the date and time.
 * 
 * @subsection csvio_special_fields Special field names
 * When reading the field names, either from a file header or as provided in the configuration file with the CSV\#_FIELDS key, the fields will be attributed to 
 * variables bearing the same names. But some field names will be recognized and automatically interpreted as either known internal 
 * parameter names (see \ref meteoparam "this list") or date/time parameters or stationID the data belongs to. Besides MeteoIO's internal parameter names, the 
 * following field names are also automatically recognized (synonyms are separated by ',' while different parameters are separated by ';'):
 *    - TIMESTAMP, TS, DATETIME; DATE; TIME; YEAR; JDAY, JDN, YDAY, DAY_OF_YEAR, DOY; MONTH; DAY; NTIME; HOUR, HOURS; MINUTE, MINUTES; SECOND, SECONDS; ID, STATIONID;
 *    - TEMPERATURE_AIR, AIRTEMP; SOIL_TEMPERATURE, SOILTEMP; PRECIPITATION, PREC; REFLECTED_RADIATION; INCOMING_RADIATION, INCOMINGSHORTWAVERADIATION; WIND_DIRE
CTION, WD; RELATIVE_HUMIDITY, RELATIVEHUMIDITY; WIND_VELOCITY, WS; PRESSURE, STATIONPRESSURE; INCOMING_LONGWAVE, INCOMINGLONGWAVERADIATION; SNOWSURFACETEMPERATURE; WS_MAX;
 * 
 * @note Since most parameter won't have names that are recognized by MeteoIO, it is advised to map them to \ref meteoparam "MeteoIO's internal names". 
 * This is done either by using the CSV_FIELDS key or using the EditingMove feature of the 
 * \ref data_editing "Input Data Editing" stage.
 * 
 * @section csvio_date_specs Date and time specification
 * In order to be able to read any date and time format, the format has to be provided in the configuration file. This is provided as a string containing
 * the following special markers:
 * - YYYY: the 4 digits year;
 * - YY: the 2 digits year (using 40 as cutoff year: date greater than 40 get converted to 1900+year, otherwise to 2000+year);
 * - MM: the two digits month;
 * - DD: the two digits day;
 * - HH24: the two digits hour of the day (0-24);
 * - MI: the two digits minutes (0-59);
 * - SS: the number of seconds (0-59.98), that can be decimal;
 * - TZ: the numerical timezone as offset to GMT (see note below).
 *
 * Any other character is interpreted as itself, present in the string. It is possible to either provide a combined datetime field (so date and time are combined into
 * one single field) or date and time as two different fields. For example:
 * - YYYY-MM-DDTHH24:MI:SS described an <A HREF="https://en.wikipedia.org/wiki/ISO_8601">ISO 8601</A> datetime field;
 * - MM/DD/YYYY described an anglo-saxon date;
 * - DD.MM.YYYY HH24:MI:SS is for a Swiss formatted datetime.
 * 
 * @note When providing a timezone field, it \em must appear at the end of the string. it can either be numerical (such as "+1.") or an abbreviation
 * such as "Z" or "CET" (see https://en.wikipedia.org/wiki/List_of_time_zone_abbreviations).
 * 
 * When this plugin identifies the fields by their column headers, it will look for TIMESTAMP or DATETIME for a combined date and time field, or DATE or TIME for (respectively) a
 * date and time field. Usually, other labels will not be recognized.
 * 
 * @section csvio_metadata_extraction Metadata extraction
 * Since there is no unified way of providing metadata (such as the location, station name, etc) in CSV files, this information has to
 * be either provided in the configuration file (see \ref csvio_keywords "Configuration keywords") or extracted out of either the file
 * name or the file headers. A specific syntax allows to describe where to find which metadata field type.
 * 
 * @subsection csvio_metadata_field_types Metadata fields types
 * The following field types are supported:
 * - NAME;
 * - ID (this will be used as a handle for the station);
 * - ALT (for the altitude);
 * - LON (for the longitude);
 * - LAT (for the latitude);
 * - EASTING (as per your input coordinate system);
 * - NORTHING (if LON/LAT is not used);
 * - SLOPE (in degrees);
 * - AZI (for the slope azimuth, in degree as read from a compass);
 * - NODATA (string to interpret as nodata);
 * - PARAM (the extracted metadata will replace the PARAM field either as found in the file's headers or in the CSV_FIELDS user configuration key);
 * - SKIP or - (skip this field).
 * 
 * If ID or NAME appear more than once in one specification string, their multiple values will be appended.
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
 * For example, to parse the filename "H0118_Generoso-Calmasino_-_Precipitation.csv" use (please note that the extension is NOT provided):
 * @code
 * CSV_FILENAME_SPEC = {ID}_{NAME}-{SKIP}_-_{PARAM}
 * @endcode
 * 
 * If the CSV_FIELDS key is also present, it will have priority. Therefore, it is possible to define one CSV_FILENAME_SPEC for several files and 
 * only define CSV\#_FIELDS for the files that would require a different handling (for example because their parameter would not be recognized).
 * Moreover, it is possible to set \em "AUTOMERGE" to "true" in the [InputEditing] section, so all files leading to the same station ID will be merged together
 * into one single station.
 * 
 * @note Obviously, the {PARAM} metadata field type can only be used for files that contain the time information (either as datetime or separate date and time) and one
 * meteorological parameter. If there are multiple (potentially unimportant) parameters in your file you have to set CSV_SINGLE_PARAM_INDEX to the column number
 * matching your parameter.
 * 
 * @section csvio_examples Examples
 * This section contains some exemplary CSV files together with the INI configuration to read them.
 * 
 * In order to read a bulletin file downloaded from IDAWEB, you need the following configuration:
 * 
 * CSV
 * @code
 * "COLUMN TO SKIP" "TIMESTAMP" "COL 1" "COL 2" "COL 3" "COL TO SKIP" "COL TO SKIP" "COL 4" "COL TO SKIP" "COL TO SKIP" "COL 5"
 * 000 2023082201 1 2 3 000 000 4 000 000 5
 * 000 2023082202 1 2 3 000 000 4 000 000 5
 * 000 2023082203 1 2 3 000 000 4 000 000 5
 * @endcode
 * INI
 * @code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_DELIMITER = SPACE
 * CSV_NR_HEADERS = 1
 * CSV_COLUMNS_HEADERS = 1
 * CSV_DATETIME_SPEC = YYYYMMDDHH24
 * CSV_NODATA = -
 *
 * METEOFILE1 = IDA_station1.csv
 * POSITION1 = latlon (46.80284, 9.77726, 2418)
 * CSV1_NAME = TEST
 * CSV1_ID = myID
 * CSV1_FIELDS = SKIP TIMESTAMP HS RSWR TA SKIP SKIP RH SKIP SKIP ILWR
 * CSV1_UNITS_OFFSET = 0 0 0 0 273.15 0 0 0 0 0 0
 * CSV1_UNITS_MULTIPLIER = 1 1 0.01 1 1 1 1 0.01 1 1 1
 * @endcode
 *
 * In order to read a CSV file produced by a Campbell data logger with Swiss-formatted timestamps, you need the following configuration:
 * 
 * CSV
 * @code
 * An arbitrary header line
 * "Some","Information","About","the","data","arbitrary","no of columns","not used for parsing"
 * "TIMESTAMP","COLUMN_NAME 1","COLUMN_NAME 2","COLUMN_NAME 3","COLUMN_NAME 4","COLUMN_NAME 5","COLUMN_NAME 6"
 * "TS","unit 1","unit 2","unit 3","unit 4","unit 5","unit 6"
 * 04.05.2023 13:00:00,0,300,600,900,1200,1500
 * 04.05.2023 14:00:00,0,300,600,900,1200,1500
 * 04.05.2023 15:00:00,0,300,600,900,1200,1500
 * 04.05.2023 16:00:00,0,300,600,900,1200,1500
 * @endcode
 * INI
 * @code
 * METEO = CSV
 * METEOPATH = ./input/meteo
 * CSV_NR_HEADERS = 4
 * CSV_COLUMNS_HEADERS = 2
 * CSV_UNITS_HEADERS = 3
 * CSV_DATETIME_SPEC = DD.MM.YYYY HH24:MI:SS
 * CSV_SPECIAL_HEADERS = name:1:2 id:1:4
 *
 * METEOFILE1 = DisMa_DisEx.csv
 * POSITION1 = latlon 46.810325 9.806657 2060
 * CSV1_ID = DIS4
 * @endcode
 *
 * For a logger produced csv file with repeating headers and quoted timestamps, you need the following configuration:
 * 
 * CSV
 * @code
 * # a information header line that is skipped
 * "Some","Information","About","the","data","arbitrary","no of columns","not used for parsing"
 * "TIMESTAMP","COLUMN_NAME 1","COLUMN_NAME 2","COLUMN_NAME 3","COLUMN_NAME 4","COLUMN_NAME 5","COLUMN_NAME 6"
 * "TS","unit 1","unit 2","unit 3","unit 4","unit 5","unit 6"
 * "2023-01-11 13:30:00",0,300,600,900,1200,1500
 * "2023-01-11 13:40:00",0,300,600,900,1200,1500
 * "2023-01-11 13:50:00",0,300,600,900,1200,1500
 * "2023-01-11 14:00:00",0,300,600,900,1200,1500
 * # a information header line that is skipped
 * "Some","Information","About","the","data","arbitrary","no of columns","not used for parsing"
 * "TIMESTAMP","COLUMN_NAME 1","COLUMN_NAME 2","COLUMN_NAME 3","COLUMN_NAME 4","COLUMN_NAME 5","COLUMN_NAME 6"
 * "TS","unit 1","unit 2","unit 3","unit 4","unit 5","unit 6"
 * "2023-01-11 15:40:00",0,300,600,900,1200,1500
 * "2023-01-11 15:50:00",0,300,600,900,1200,1500
 * "2023-01-11 16:00:00",0,300,600,900,1200,1500
 * @endcode
 * INI
 * @code
 * METEO = CSV
 * METEOPATH = /path/to/input
 * CSV_DELIMITER = ,
 * CSV_NR_HEADERS = 3
 * CSV_HEADER_REPEAT_MK = #
 * CSV_UNITS_SOURCE = FROM HEADERS
 * CSV_UNITS_HEADERS = 3
 * CSV_COLUMNS_HEADERS = 2
 * CSV_TIMESTAMP = COMBINED
 * CSV_DATETIME_SPEC = "YYYY-MM-DD HH24:MI:SS"
 * POSITION = xy(198754, 723458,2200)
 * @endcode
 * 
 * In order to read a set of files each containing only one parameter and merge them together (see \ref data_editing "input data editing" for more
 * on the merge feature), extracting the station ID, name and meteorological parameter from the filename:
 *@code
 * [Input]
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
 * METEOFILE1 = H0118_Generoso_-_Calmasino_precipitation.csv
 *
 * METEOFILE2 = H0118_Generoso_-_Calmasino_temperature.csv    #the parameter name is ambiguous, it will not be recognized
 * CSV2_FIELDS = DATE TIME TA                               #so we define the parameter manually
 * CSV2_UNITS_OFFSET = 0 0 273.15
 *
 * METEOFILE3 = H0118_Generoso_-_Calmasino_reflected_solar_radiation.csv
 * METEOFILE4 = H0118_Generoso_-_Calmasino_relative_humidity.csv
 * METEOFILE5 = H0118_Generoso_-_Calmasino_wind_velocity.csv
 *
 * [InputEditing]
 * AUTOMERGE = true
 * @endcode
 * Please note that here the file's headers will look like 'DATE;TIME;PARAM' which will allow PARAM to be replaced by
 * the PARAM special value extracted from the file name.
 * 
 * In order to read a set of files and merge them together (see \ref data_editing "input data editing" for more
 * on the merge feature):
 *@code
 * [Input]
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
 * METEOFILE1 = H0118_lg23456.csv
 *
 * CSV2_FIELDS = DATE TIME TA
 * METEOFILE2 = H0118_lg7850.csv
 * CSV2_UNITS_OFFSET = 0 0 273.15
 *
 * CSV3_FIELDS = DATE TIME RSWR ISWR
 * METEOFILE3 = H0118_lg64520.csv
 *
 * CSV4_FIELDS = DATE TIME RH
 * METEOFILE4 = H0118_lg45302.csv
 * CSV4_UNITS_MULTIPLIER = 1 1 0.01
 *
 * CSV5_FIELDS = DATE TIME VW
 * METEOFILE5 = H0118_wind_velocity.csv
 *
 * [InputEditing]
 * ID1::MERGE = ID2 ID3 ID4 ID5
 * @endcode
 * 
 * When reading a file containing the data from multiple stations, each line containing the station ID it applies to, it is necessary
 * to provide the requested station ID for each new station as well as declare which is the ID field (if the headers declare an "ID" column or
 * a "STATIONID" column, this will also work)
 * @code
 * METEO = CSV
 * METEOPATH = ./input
 * CSV_DELIMITER = ,   
 * CSV_NR_HEADERS = 1
 * CSV_COLUMNS_HEADERS = 1
 * CSV_DATETIME_SPEC = YYYY-MM-DD HH24:MI:SS
 * CSV_FIELDS = ID TIMESTAMP SKIP TA RH DW VW SKIP HS SKIP SKIP P SKIP SKIP SKIP ISWR SKIP SKIP SKIP TSG SKIP SKIP ISWR ILWR TSS
 * CSV_UNITS_OFFSET = 0 0 0 273.15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 273.15
 * CSV_UNITS_MULTIPLIER = 1 1 1 1 0.01 1 1 1 0.01 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 * 
 * METEOFILE1 = Extracted_data.csv
 * POSITION1 = latlon (46.8, 9.81, 1511.826)
 * CSV1_ID = 109
 * CSV1_NAME = Station109
 * 
 * METEOFILE2 = Extracted_data.csv
 * POSITION2 = xy (45.8018, 9.82, 111.826)
 * CSV2_ID = 105
 * CSV2_NAME = Station105
 * @endcode
 *
 * @section csvio_debugging Debugging CSV files reading
 * Unfortunately, there are some cases when the data in one or multiple files does not match the expected field content. This can happen with
 * the data being spread over multiple files and the column order changing between files or even with the column order changing within
 * a given file. Depending on the amount of data and the number of files, this can be quite cumbersome to correct, if even possible at all.
 * Here is some rough procedure to help correcting such issues:
 *     1. Split the data at the time of changes. If the data is provided for example in yearly files and such changes appear
 *        between files, created one station ID per file (such as {base_id}_{year}). If such changes appear within one file,
 *        split the file by reading it multiple times (with multiple METEOFILE# statements) with different lines exclusions
 *        commands and attribute a different station ID for each.
 *     2. Name each column with its column index, using the CSV\#_NUMBER_FIELDS key. You can combine it with CSV\#_FIELDS_POSTFIX
 *        in order to also have an indication of the origin of the data (such as the year or the line ranges).
 *     3. Plot all the fields for each station ID (like a contact print or as individual plots) with the field name clearly
 *        visible. It is often not that hard to distinguish between different meteorological parameters (such as between snow height
 *        and wind speed), so it will be possible to quickly spot obvious mis-labeled fields.
 *     4. Load successive data periods together on the same plot, select one parameter from one period and cycle through the
 *        parameters of the other period until you find that there is a smooth transition between the two periods. The column index
 *        in the field names (as advised to do in (2)) allow to know which field is at which index for which period.
 *
 */


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

	//handle the deprecated STATION# syntax
	std::string meteofiles_key( "METEOFILE" );
	const std::vector< std::string > vecDeprecated( cfg.getKeys("STATION", "INPUT") );
	if (!vecDeprecated.empty()) {
		meteofiles_key = "STATION";
		std::cerr << "[W] The STATION# syntax has been deprecated for the CSV input plugin, please rename these keys as METEOFILE#!\n";
		//throw InvalidArgumentException("The STATION# syntax has been deprecated for the CSV plugin, please rename these keys as METEOFILE#!", AT);
	}
	
	const double in_TZ = cfg.get("TIME_ZONE", "Input");
	const std::string meteopath = cfg.get("METEOPATH", "Input");
	std::vector< std::pair<std::string, std::string> > vecFilenames( cfg.getValues(meteofiles_key, "INPUT") );

	if (vecFilenames.empty()) { //scan all of the data path for a given file extension if no stations are specified
		std::string csvext(".csv");
		cfg.getValue("CSV_FILE_EXTENSION", "Input", csvext, IOUtils::nothrow); 
		
		std::vector<std::string> tmpFilenames;
		scanMeteoPath(cfg, meteopath, tmpFilenames, csvext);

		size_t hit = 0; //human readable iterator
		for (const std::string& filename : tmpFilenames)	{
			hit++;
			std::stringstream ss;
			ss << meteofiles_key << hit; //assign alphabetically ordered ID to station
			
			const std::pair<std::string, std::string> stat_id_and_name(ss.str(), filename);
			vecFilenames.push_back(stat_id_and_name);
		}
	} 
	
	for (size_t ii=0; ii<vecFilenames.size(); ii++) {
		const std::string idx( vecFilenames[ii].first.substr(meteofiles_key.length()) );
		static const std::string dflt("CSV_"); //the prefix for a key for ALL stations
		const std::string pre( "CSV"+idx+"_" ); //the prefix for the current station only
		
		CsvParameters tmp_csv(in_TZ);
		std::string coords_specs;
		if (cfg.keyExists("POSITION"+idx, "INPUT")) cfg.getValue("POSITION"+idx, "INPUT", coords_specs);
		else cfg.getValue("POSITION", "INPUT", coords_specs, IOUtils::nothrow);
		
		std::string name;
		if (cfg.keyExists(pre+"NAME", "Input")) cfg.getValue(pre+"NAME", "Input", name);
		else cfg.getValue(dflt+"NAME", "Input", name, IOUtils::nothrow);
		
		std::string id;
		if (cfg.keyExists(pre+"ID", "Input")) cfg.getValue(pre+"ID", "Input", id);
		else cfg.getValue(dflt+"ID", "Input", id, IOUtils::nothrow);
		if (!coords_specs.empty()) {
			const Coords loc(coordin, coordinparam, coords_specs);
			tmp_csv.setLocation(loc, name, id);
		} else if (!coordin.empty()) { //coord system parameters in [INPUT], coordinates in csv files
			tmp_csv.setLocation(Coords(coordin, coordinparam), name, id);
		} else {
			tmp_csv.setLocation(Coords(), name, id);
		}
		
		double slope=IOUtils::nodata;
		if (cfg.keyExists(pre+"SLOPE", "INPUT")) cfg.getValue(pre+"SLOPE", "INPUT", slope);
		else cfg.getValue(dflt+"SLOPE", "INPUT", slope, IOUtils::nothrow);
		double azimuth=IOUtils::nodata;
		if (cfg.keyExists(pre+"AZIMUTH", "INPUT")) cfg.getValue(pre+"AZIMUTH", "INPUT", azimuth);
		else cfg.getValue(dflt+"AZIMUTH", "INPUT", azimuth, IOUtils::nothrow);
		tmp_csv.setSlope(slope, azimuth);
		
		std::string csv_nodata;
		if (cfg.keyExists(pre+"NODATA", "Input")) cfg.getValue(pre+"NODATA", "Input", csv_nodata);
		else cfg.getValue(dflt+"NODATA", "Input", csv_nodata, IOUtils::nothrow);
		tmp_csv.setNodata( csv_nodata );
		
		std::string delim_spec(","); //default delimiter
		if (cfg.keyExists(pre+"DELIMITER", "Input")) cfg.getValue(pre+"DELIMITER", "Input", delim_spec);
		else cfg.getValue(dflt+"DELIMITER", "Input", delim_spec, IOUtils::nothrow);
		tmp_csv.setDelimiter(delim_spec);
		
		std::string PurgeCharsSpecs;
		if (cfg.keyExists(pre+"PURGE_CHARS", "INPUT")) cfg.getValue(pre+"PURGE_CHARS", "INPUT", PurgeCharsSpecs);
		else cfg.getValue(dflt+"PURGE_CHARS", "INPUT", PurgeCharsSpecs, IOUtils::nothrow);
		if (!PurgeCharsSpecs.empty()) {
			tmp_csv.setPurgeChars( PurgeCharsSpecs );
		}
		
		char comments_mk='\n';
		if (cfg.keyExists(pre+"COMMENTS_MK", "Input")) cfg.getValue(pre+"COMMENTS_MK", "Input", comments_mk);
		else cfg.getValue(dflt+"COMMENTS_MK", "Input", comments_mk, IOUtils::nothrow);
		if (comments_mk!='\n') tmp_csv.comments_mk = comments_mk;

		std::string header_delim_spec;
		if (cfg.keyExists(pre+"HEADER_DELIMITER", "Input"))
			cfg.getValue(pre+"HEADER_DELIMITER", "Input", header_delim_spec);
		else if (cfg.keyExists(dflt+"HEADER_DELIMITER", "Input"))
			cfg.getValue(dflt+"HEADER_DELIMITER", "Input", header_delim_spec);
		else
			header_delim_spec = delim_spec;
		tmp_csv.setHeaderDelimiter(header_delim_spec);

		std::string hdr_repeat_mk;
		if (cfg.keyExists(pre+"HEADER_REPEAT_MK", "Input")) cfg.getValue(pre+"HEADER_REPEAT_MK", "Input", hdr_repeat_mk);
		else cfg.getValue(dflt+"HEADER_REPEAT_MK", "Input", hdr_repeat_mk, IOUtils::nothrow);
		tmp_csv.setHeaderRepeatMk(hdr_repeat_mk);
		
		if (cfg.keyExists(pre+"NR_HEADERS", "Input")) cfg.getValue(pre+"NR_HEADERS", "Input", tmp_csv.header_lines);
		else cfg.getValue(dflt+"NR_HEADERS", "Input", tmp_csv.header_lines, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"COLUMNS_HEADERS", "Input")) cfg.getValue(pre+"COLUMNS_HEADERS", "Input", tmp_csv.columns_headers);
		else cfg.getValue(dflt+"COLUMNS_HEADERS", "Input", tmp_csv.columns_headers, IOUtils::nothrow);
		if (tmp_csv.columns_headers>tmp_csv.header_lines) tmp_csv.columns_headers = IOUtils::npos;
		
		if (cfg.keyExists(pre+"FIELDS", "Input")) cfg.getValue(pre+"FIELDS", "Input", tmp_csv.csv_fields);
		else cfg.getValue(dflt+"FIELDS", "Input", tmp_csv.csv_fields, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"FIELDS_POSTFIX", "Input")) cfg.getValue(pre+"FIELDS_POSTFIX", "Input", tmp_csv.fields_postfix);
		else cfg.getValue(dflt+"FIELDS_POSTFIX", "Input", tmp_csv.fields_postfix, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"NUMBER_FIELDS", "Input")) cfg.getValue(pre+"NUMBER_FIELDS", "Input", tmp_csv.number_fields);
		else cfg.getValue(dflt+"NUMBER_FIELDS", "Input", tmp_csv.number_fields, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"FILTER_ID", "Input")) cfg.getValue(pre+"FILTER_ID", "Input", tmp_csv.filter_ID);
		else cfg.getValue(dflt+"FILTER_ID", "Input", tmp_csv.filter_ID, IOUtils::nothrow);
		
		std::string skipFieldSpecs;
		if (cfg.keyExists(pre+"SKIP_FIELDS", "Input")) cfg.getValue(pre+"SKIP_FIELDS", "Input", skipFieldSpecs);
		else cfg.getValue(dflt+"SKIP_FIELDS", "Input", skipFieldSpecs, IOUtils::nothrow);
		if (!skipFieldSpecs.empty()) tmp_csv.setSkipFields( skipFieldSpecs, false );
		
		std::string onlyFieldSpecs;
		if (cfg.keyExists(pre+"ONLY_FIELDS", "Input")) cfg.getValue(pre+"ONLY_FIELDS", "Input", onlyFieldSpecs);
		else cfg.getValue(dflt+"ONLY_FIELDS", "Input", onlyFieldSpecs, IOUtils::nothrow);
		if (!onlyFieldSpecs.empty()) tmp_csv.setSkipFields( onlyFieldSpecs, true );
		
		if (cfg.keyExists(pre+"UNITS_HEADERS", "Input")) cfg.getValue(pre+"UNITS_HEADERS", "Input", tmp_csv.units_headers);
		else cfg.getValue(dflt+"UNITS_HEADERS", "Input", tmp_csv.units_headers, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"UNITS_OFFSET", "Input")) cfg.getValue(pre+"UNITS_OFFSET", "Input", tmp_csv.field_offset);
		else cfg.getValue(dflt+"UNITS_OFFSET", "Input", tmp_csv.field_offset, IOUtils::nothrow);
		
		if (cfg.keyExists(pre+"UNITS_MULTIPLIER", "Input")) cfg.getValue(pre+"UNITS_MULTIPLIER", "Input", tmp_csv.field_multiplier);
		else cfg.getValue(dflt+"UNITS_MULTIPLIER", "Input", tmp_csv.field_multiplier, IOUtils::nothrow);
		
		std::string csv_units;
		if (cfg.keyExists(pre+"UNITS", "Input")) cfg.getValue(pre+"UNITS", "Input", csv_units);
		else cfg.getValue(dflt+"UNITS", "Input", csv_units, IOUtils::nothrow);
		if (!csv_units.empty()) tmp_csv.setUnits( csv_units );
		
		//Date and time formats. The defaults will be set when parsing the column names (so they are appropriate for the available columns)
		std::string datetime_spec, date_spec, time_spec, decimaldate_type;
		if (cfg.keyExists(pre+"DECIMALDATE_TYPE", "Input")) cfg.getValue(pre+"DECIMALDATE_TYPE", "Input", decimaldate_type);
		if (cfg.keyExists(pre+"DATETIME_SPEC", "Input")) cfg.getValue(pre+"DATETIME_SPEC", "Input", datetime_spec);
		if (cfg.keyExists(pre+"DATE_SPEC", "Input")) cfg.getValue(pre+"DATE_SPEC", "Input", date_spec);
		if (cfg.keyExists(pre+"TIME_SPEC", "Input")) cfg.getValue(pre+"TIME_SPEC", "Input", time_spec);
		if (decimaldate_type.empty() && datetime_spec.empty() && date_spec.empty() && time_spec.empty() && datetime_spec.empty()) {
			cfg.getValue(dflt+"DECIMALDATE_TYPE", "Input", decimaldate_type, IOUtils::nothrow);
			cfg.getValue(dflt+"DATETIME_SPEC", "Input", datetime_spec, IOUtils::nothrow);
			cfg.getValue(dflt+"DATE_SPEC", "Input", date_spec, IOUtils::nothrow);
			cfg.getValue(dflt+"TIME_SPEC", "Input", time_spec, IOUtils::nothrow);
		}
		tmp_csv.setDateTimeSpecs(datetime_spec, date_spec, time_spec, decimaldate_type);
		
		int fixed_year=IOUtils::inodata;
		int fixed_hour=IOUtils::inodata;
		if (cfg.keyExists(pre+"FALLBACK_YEAR", "Input")) cfg.getValue(pre+"FALLBACK_YEAR", "Input", fixed_year);
		else cfg.getValue(dflt+"FALLBACK_YEAR", "Input", fixed_year, IOUtils::nothrow);
		bool auto_wrap_year=true;
		if (cfg.keyExists(pre+"FALLBACK_AUTO_WRAP", "Input")) cfg.getValue(pre+"FALLBACK_AUTO_WRAP", "Input", auto_wrap_year);
		else cfg.getValue(dflt+"FALLBACK_AUTO_WRAP", "Input", auto_wrap_year, IOUtils::nothrow);
		if (fixed_year!=IOUtils::inodata) tmp_csv.setFixedYear( fixed_year, auto_wrap_year );

		if (cfg.keyExists(pre+"FALLBACK_TIME", "Input")) cfg.getValue(pre+"FALLBACK_TIME", "Input", fixed_hour);
		else cfg.getValue(dflt+"FALLBACK_TIME", "Input", fixed_hour, IOUtils::nothrow);
		if (fixed_hour!=IOUtils::inodata) tmp_csv.setFixedHour( fixed_hour );

		std::string dateRangeHint;
		if (cfg.keyExists(pre+"COVERAGE_HINT", "INPUT")) cfg.getValue(pre+"COVERAGE_HINT", "INPUT", dateRangeHint);
		else cfg.getValue(dflt+"COVERAGE_HINT", "INPUT", dateRangeHint, IOUtils::nothrow);
		if (!dateRangeHint.empty()) tmp_csv.setCoverageHint( dateRangeHint );

		std::vector<std::string> vecMetaSpec;
		if (cfg.keyExists(pre+"SPECIAL_HEADERS", "Input")) cfg.getValue(pre+"SPECIAL_HEADERS", "Input", vecMetaSpec);
		else cfg.getValue(dflt+"SPECIAL_HEADERS", "Input", vecMetaSpec, IOUtils::nothrow);
		
		//handling of lines restrictions (either as ONLY or EXCLUDE statements)
		std::string linesExclusionsSpecs;
		if (cfg.keyExists(pre+"EXCLUDE_LINES", "INPUT")) cfg.getValue(pre+"EXCLUDE_LINES", "INPUT", linesExclusionsSpecs);
		else cfg.getValue(dflt+"EXCLUDE_LINES", "INPUT", linesExclusionsSpecs, IOUtils::nothrow);
		
		std::string linesRestrictionsSpecs;
		if (cfg.keyExists(pre+"ONLY_LINES", "INPUT")) cfg.getValue(pre+"ONLY_LINES", "INPUT", linesRestrictionsSpecs);
		else cfg.getValue(dflt+"ONLY_LINES", "INPUT", linesRestrictionsSpecs, IOUtils::nothrow);
		
		if (!linesExclusionsSpecs.empty() || !linesRestrictionsSpecs.empty()) {
			if (!linesExclusionsSpecs.empty() && !linesRestrictionsSpecs.empty()) 
				throw InvalidArgumentException("It is not possible to provide both CSV_EXCLUDE_LINES and CSV_ONLY_LINES", AT);
			if (!linesExclusionsSpecs.empty()) {
				const std::vector< LinesRange > lrRange( initLinesRestrictions(linesExclusionsSpecs, "INPUT::CSV_EXCLUDE_LINES", false) );
				tmp_csv.setLinesExclusions( lrRange );
			} else {
				const std::vector< LinesRange > lrRange( initLinesRestrictions(linesRestrictionsSpecs, "INPUT::CSV_ONLY_LINES", true) );
				tmp_csv.setLinesExclusions( lrRange );
			}
		}
		
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
	for (const CsvParameters &params : csvparam)
		vecStation.push_back( params.getStation() );
}

//build MeteoData template, based on parameters available in the csv file
MeteoData CsvIO::createTemplate(const CsvParameters& params)
{
	const size_t nr_of_data_fields = params.csv_fields.size(); //this has been checked by CsvParameters
	
	//build MeteoData template
	MeteoData template_md( Date(0., 0.), params.getStation() );
	for (size_t ii=0; ii<nr_of_data_fields; ii++) {
		if (params.skipField( ii)) continue;
		template_md.addParameter( params.csv_fields[ii] );
	}

	return template_md;
}

Date CsvIO::getDate(CsvParameters& params, const std::vector<std::string>& vecFields, const bool& silent_errors, const std::string& filename, const size_t& linenr)
{
	try {
		const Date dt( params.getDate(vecFields) );
		if (dt.isUndef()) {
			const std::string err_msg( "Date or time could not be read in file \'"+filename+"' at line "+IOUtils::toString(linenr) );
			throw InvalidFormatException(err_msg, AT);
		}

		return dt;
	} catch (...) {
		const std::string err_msg( "Date or time could not be read in file \'"+filename+"' at line "+IOUtils::toString(linenr) );
		std::cerr << err_msg << "\n";
		if (!silent_errors)
			throw;
	}
	return Date();
}

std::vector<MeteoData> CsvIO::readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd)
{
	const std::string filename( params.getFilename() );
	size_t nr_of_data_fields = params.csv_fields.size(); //this has been checked by CsvParameters
	const bool use_unit_offset = !params.units_offset.empty();
	const bool use_unit_multiplier = !params.units_multiplier.empty();
	if ((use_unit_offset && params.units_offset.size()!=nr_of_data_fields) || (use_unit_multiplier && params.units_multiplier.size()!=nr_of_data_fields)) {
		const std::string msg( "in file '"+filename+"', the declared units_offset ("+IOUtils::toString(params.units_offset.size())+") / units_multiplier ("+IOUtils::toString(params.units_multiplier.size())+") must match the number of columns ("+IOUtils::toString(nr_of_data_fields)+") in the file!" );
		throw InvalidFormatException(msg, AT);
	}
	
	const bool use_field_offset = !params.field_offset.empty();
	const bool use_field_multiplier = !params.field_multiplier.empty();
	if ((use_field_offset && params.field_offset.size()!=nr_of_data_fields) || (use_field_multiplier && params.field_multiplier.size()!=nr_of_data_fields)) {
		const std::string msg( "in file '"+filename+"', the declared field_offset ("+IOUtils::toString(params.field_offset.size())+") / field_multiplier ("+IOUtils::toString(params.field_multiplier.size())+") must match the number of columns ("+IOUtils::toString(nr_of_data_fields)+") in the file!" );
		throw InvalidFormatException(msg, AT);
	}
	
	const MeteoData template_md( createTemplate(params) );

	//now open the file
	if (!FileUtils::fileExists(filename)) throw AccessException("File '"+filename+"' does not exists", AT); //prevent invalid filenames
	errno = 0;
	std::ifstream fin(filename.c_str(), ios::in|ios::binary); //ascii does end of line translation, which messes up the pointer code
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Error opening file \"" << filename << "\" for reading, possible reason: " << std::strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}
	
	std::string line;
	size_t linenr=0;
	streampos fpointer = indexer_map[filename].getIndex(dateStart, linenr);
	if (fpointer!=static_cast<streampos>(-1) && params.asc_order) {
		fin.seekg(fpointer); //a previous pointer was found, jump to it
	} else {
		//skip the headers (they have been read already, so we know this works)
		const size_t skip_count = (params.header_repeat_at_start)? params.header_lines+1 : params.header_lines;
		FileUtils::skipLines(fin, skip_count);
		linenr += skip_count;
	}
	
	//and now, read the data and fill the vector vecMeteo
	std::vector<MeteoData> vecMeteo;
	std::vector<std::string> tmp_vec;
	const std::string filterID = (params.filter_ID.empty())? template_md.getStationID() : params.filter_ID; //necessary if filtering on stationID field
	const char comments_mk = params.comments_mk;
	const bool delimIsNoWS = (params.csv_delim!=' ');
	const bool hasHeaderRepeat = (!params.header_repeat_mk.empty());
	const bool purgeChars = params.hasPurgeChars();
	bool has_exclusions = true; //this wil be set at the first call to params.excludeLine
	Date prev_dt;
	while (!fin.eof()){
		getline(fin, line, params.eoln);
		linenr++;
		
		if (has_exclusions && params.excludeLine( linenr, has_exclusions )) continue;
		
		if (comments_mk!='\n') IOUtils::stripComments(line, comments_mk);
		if (purgeChars) params.purgeChars(line);
		IOUtils::trim( line );
		if (line.empty()) continue; //Pure comment lines and empty lines are ignored
		if (hasHeaderRepeat && line.find(params.header_repeat_mk)!=std::string::npos) {
			FileUtils::skipLines(fin, params.header_lines);
			linenr += params.header_lines;
			continue;
		}
		
		const size_t nr_curr_data_fields = (delimIsNoWS)? IOUtils::readLineToVec(line, tmp_vec, params.csv_delim) : IOUtils::readLineToVec(line, tmp_vec);
		if (nr_of_data_fields==0) nr_of_data_fields = nr_curr_data_fields;
		
		//filter on ID
		if (params.ID_col!=IOUtils::npos) {
			if (tmp_vec.size()<=params.ID_col) { //we can not filter on the ID although it has been requested so we have to stop!
				std::ostringstream ss;
				ss << "File \'" << filename << "\' declares station ID in column " << params.ID_col << " but only has " << tmp_vec.size() << " columns at line ";
				ss << linenr << " :\n'" << line << "'\n";
				throw InvalidFormatException(ss.str(), AT);
			}
			
			if (tmp_vec[params.ID_col]!=filterID) continue;
		}
		
		//check that we have the expected number of fields
		if (nr_curr_data_fields!=nr_of_data_fields) {
			std::ostringstream ss;
			ss << "File \'" << filename << "\' declares (either as first data line or columns headers or units offset/multiplier) " << nr_of_data_fields << " columns ";
			ss << "but this does not match line " << linenr << " with " << nr_curr_data_fields << " fields :\n'" << line << "'\n";
			if (silent_errors) {
				std::cerr << ss.str();
				continue;
			} else throw InvalidFormatException(ss.str(), AT);
		}

		const Date dt( getDate(params, tmp_vec, silent_errors, filename, linenr) );
		if (dt.isUndef() && silent_errors) continue;
		if (!prev_dt.isUndef()) {
			if (dt==prev_dt)
				std::cerr << "File \'" << filename << "\' has duplicated timestamps for " << dt.toString(Date::ISO) << " at line " << linenr << "\n";
			if (dt<prev_dt)
				std::cerr << "File \'" << filename << "\' has out of order timestamps for " << dt.toString(Date::ISO) << " at line " << linenr << "\n";
		}
		prev_dt = dt;

		if (linenr % streampos_every_n_lines == 0) {
			fpointer = fin.tellg();
			if (fpointer != static_cast<streampos>(-1)) indexer_map[filename].setIndex(dt, fpointer, linenr);
		}
		if (params.asc_order) {
			if (dt<dateStart) continue;
			if (dt>dateEnd) break;
		} else {
			if (dt<dateStart) break;
			if (dt>dateEnd) continue;
		}
		
		MeteoData md(template_md);
		md.setDate(dt);
		bool no_errors = true;
		for (size_t ii=0; ii<tmp_vec.size(); ii++) {
			if (params.skipField( ii )) continue; //the user has requested this field to be skipped or this is a special field
			if (params.isNodata( tmp_vec[ii] )) continue; //recognize nodata
			
			double tmp;
			if (!IOUtils::convertString(tmp, tmp_vec[ii])) {
				const std::string err_msg( "Could not parse field '"+tmp_vec[ii]+"' in file \'"+filename+"' at line "+IOUtils::toString(linenr) );
				if (silent_errors) {
					std::cerr << err_msg << "\n";
					no_errors = false;
					continue;
				} else if (errors_to_nodata) {
					tmp = IOUtils::nodata;
				} else throw InvalidFormatException(err_msg, AT);
			}
			
			if (tmp!=IOUtils::nodata) {
				if (use_unit_multiplier) tmp *= params.units_multiplier[ii];
				if (use_unit_offset) tmp += params.units_offset[ii];
				if (use_field_multiplier) tmp *= params.field_multiplier[ii];
				if (use_field_offset) tmp += params.field_offset[ii];
			}
			
			md( params.csv_fields[ii] ) = tmp;
		}
		if (no_errors) vecMeteo.push_back( md );
	}

	if (!params.asc_order) std::reverse(vecMeteo.begin(), vecMeteo.end());
	
	return vecMeteo;
}

void CsvIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	vecMeteo.reserve( csvparam.size() ); //max possible size

	for (size_t ii=0; ii<csvparam.size(); ii++) { //loop over CSV files
		const std::vector<MeteoData> tmpVec( readCSVFile(csvparam[ii], dateStart, dateEnd) );
		if (tmpVec.empty()) continue;
		const std::string tmpID( tmpVec.front().getStationID() );

		//does this stationID already exist?
		size_t append_index = IOUtils::npos;
		for (size_t jj=0; jj<std::min(ii, vecMeteo.size()); jj++) {
			if (vecMeteo[jj].back().getStationID()==tmpID) {
				append_index = jj;
				break;
			}
		}

		if (append_index==IOUtils::npos) { //this is a new stationID, add to vecMeteo
			vecMeteo.push_back( tmpVec );
		} else { //this stationID already exists, appending
			vecMeteo[append_index].insert( vecMeteo[append_index].end(), tmpVec.begin(), tmpVec.end() );
		}
	}

	vecMeteo.shrink_to_fit();	//free non-used storage
}

} //namespace
