// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <fstream>
#include <iostream>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/plugins/iCSVHelper.h>
#include <regex>

// clang-format off
static const std::regex POINT("^POINT\\s*\\(\\s*-?\\d+(\\.\\d+)?\\s+-?\\d+(\\.\\d+)?\\s*\\)$");
static const std::regex FALSEPOINT("^POINT\\s*\\(\\s*-?\\d+(\\.\\d+)?\\s+-?\\d+(\\.\\d+)?\\s+-?\\d+(\\.\\d+)?\\s*\\)$");
static const std::regex POINTZ("^POINTZ\\s*\\(\\s*-?\\d+(\\.\\d+)?\\s+-?\\d+(\\.\\d+)?\\s+-?\\d+(\\.\\d+)?\\s*\\)$");
// clang-format on

/**
 * @file iCSVHelper.cc
 * @brief Contains helper functions and classes for processing iCSV data.
 */
namespace mio {
namespace iCSV {
    
// ----------------- iCSV Helpers -----------------
static void custom_assert(const std::string& value, const std::string& asssert_val) {
    if (value.empty() && !asssert_val.empty()) {
        std::cerr << "Assertion failed: No input value" << std::endl;
        std::abort();
    }
    bool condition = value == asssert_val;
    if (!condition) {
        std::cerr << "Assertion failed: " << value << " != " << asssert_val << std::endl;
        std::abort();
    }
}

/**
* Checks if the first line of iCSV data is conforming to the standard.
*
* @param firstline The first line of iCSV data.
* @return True if the first line is valid, false otherwise.
*/
static bool isValidFirstLine(const std::string &firstline) {
    std::istringstream iss(firstline);
    std::string hash, icsv, version, encoding;
    iss >> hash >> icsv >> version >> encoding;
    custom_assert(icsv, "iCSV");
    custom_assert(version, "1.0");
    custom_assert(encoding, "UTF-8");
    return true;
}

/**
* @brief Extracts values out of geometry strings.
*/
static std::vector<double> extractMatches(const std::string &match) {
    static const std::regex coord_regex("-?\\d+(\\.\\d+)?");
    auto words_begin = std::sregex_iterator(match.begin(), match.end(), coord_regex);
    auto words_end = std::sregex_iterator();

    std::vector<double> temp_coordinates;
    for (std::sregex_iterator ii = words_begin; ii != words_end; ++ii) {
        const std::smatch smatch = *ii;
        const std::string match_str( smatch.str() );
        temp_coordinates.push_back(std::stod(match_str));
    }

    return temp_coordinates;
}

/**
* @brief process geometry strings to return a geoLocation object.
*
* For Now only returs the last point in the geometry string.
*/
static geoLocation processGeometryRegexes(const std::string &geometry, const std::regex &reg, const bool& has_Z) {
    size_t stride = has_Z ? 3 : 2;
    std::smatch matches;
    geoLocation location;
    if (std::regex_search(geometry, matches, reg)) {
        const std::vector<double> temp_coordinates( extractMatches(matches.str()) );
        for (size_t ii = 0; ii < temp_coordinates.size(); ii += stride) {
            location = geoLocation(temp_coordinates[ii], temp_coordinates[ii + 1], stride == 3 ? temp_coordinates[ii + 2] : default_nodata);
        }
    }
    return location;
}

/**
* @brief Extracts the EPSG code from a string in the format "EPSG:XXXX".
*
* @param epsgStr The string containing the EPSG code.
* @return The extracted EPSG code as a string.
*/
static std::string extractEpsgCode(const std::string &epsgStr) {
    static const std::regex pattern("^EPSG:(\\d{4,5})$");
    std::smatch matches;

    if (std::regex_search(epsgStr, matches, pattern) && matches.size() > 1) {
        return matches[1].str(); // matches[1] contains the first capture group
    } else {
        return "";
    }
}

static bool isValidEpsgCode(const std::string &epsgStr) { return !extractEpsgCode(epsgStr).empty(); }

/**
* Converts a vector of strings to a vector of doubles.
*
* @param vec The vector of strings to be converted.
* @return The converted vector of doubles.
*/
std::vector<double> convertVector(const std::vector<std::string> &vec) {
    std::vector<double> result;
    for (const auto &str : vec) {
        double val;
        IOUtils::convertString(val, str);
        result.push_back(val);
    }
    return result;
}

std::vector<Coords> convertVector(const std::vector<geoLocation> &vec, const int& epsg) {
    std::vector<Coords> result;
    for (auto loc : vec) {
        loc.standardizeNodata();
        Coords coordinate;
        coordinate.setPoint(loc.x, loc.y, loc.z, epsg);
        result.push_back(coordinate);
    }
    return result;
}

/**
* Converts the given Coords object to a geoLocation object.
* @param[in] loc The Coords object to convert.
* @param[in] epsg <a href="https://epsg.org/home.html">EPSG code</a> of the provided coordinates
* @return The converted geoLocation object.
*/
geoLocation toiCSVLocation(Coords loc, const int& epsg) {
    geoLocation location;
    if (Coords::latlon_epsgs.find(epsg) != Coords::latlon_epsgs.end()) {
        location = geoLocation(loc.getLat(), loc.getLon(), loc.getAltitude());
    } else {
        loc.setEPSG(epsg);
        location = geoLocation(loc.getEasting(), loc.getNorthing(), loc.getAltitude());
    }
    return location;
}

/**
* @brief extracts a geoLocation object out of a <a href="https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry">WKT geometry string</a>.
* @param[in] geometry WKT geometry string
* @return the matching geoLocation object
*/
geoLocation extractCoordinates(const std::string &geometry) {
    geoLocation coordinates;

    // In your main function, replace the two for-loops with calls to processRegexes:
    coordinates = processGeometryRegexes(geometry, POINT, false);

    if (coordinates.isEmpty()) {
        coordinates = processGeometryRegexes(geometry, POINTZ, true);
    }

    if (coordinates.isEmpty()) {
        coordinates = processGeometryRegexes(geometry, FALSEPOINT, true);
    }
    coordinates.standardizeNodata();
    return coordinates;
}


// ----------------- iCSVFile -----------------
iCSVFile::iCSVFile()
    : skip_lines_to_data(0), filename(""), firstline(""), station_location(), METADATA(), FIELDS(), location_in_header(true),
        timezone_in_data(false), timestamp_present(false), julian_present(false), time_id(IOUtils::npos), location_id(IOUtils::npos),
        dates_in_file(), row_data(), locations_in_data() 
{}

iCSVFile::iCSVFile(const iCSVFile &other)
    : skip_lines_to_data(other.skip_lines_to_data), filename(other.filename), firstline(other.firstline),
        station_location(other.station_location), METADATA(other.METADATA), FIELDS(other.FIELDS),
        location_in_header(other.location_in_header), timezone_in_data(other.timezone_in_data),
        timestamp_present(other.timestamp_present), julian_present(other.julian_present), time_id(other.time_id),
        location_id(other.location_id), dates_in_file(other.dates_in_file), row_data(other.row_data),
        locations_in_data(other.locations_in_data) 
{}

/**
* iCSVFile constructor
*
* @brief Reads a iCSV file, and checks the format
*
* The data is only read if read_sequential is set to false. Otherwise, only the HEADER is processed
*
* @param[in] infile The path to the input file.
* @param[in] read_sequential If true, only HEADER will be processed.
*/
iCSVFile::iCSVFile(const std::string &infile, const bool& read_sequential)
    : skip_lines_to_data(0), filename(infile), firstline(""), station_location(), METADATA(), FIELDS(), location_in_header(true),
        timezone_in_data(false), timestamp_present(false), julian_present(false), time_id(IOUtils::npos), location_id(IOUtils::npos),
        dates_in_file(), row_data(), locations_in_data() 
{
    readFile(infile, read_sequential);
    parseGeometry();
    checkFormatValidity();
    checkMeteoIOCompatibility();
}

/**
* Reads the contents of a iCSV file.
*
* @param[in] infile The path of the input file to be read.
* @param[in] sequentially Flag indicating whether to read the data
* @throws IOException If there is an error opening or reading the file.
* @throws InvalidFormatException If the iCSV file format is invalid.
*/
void iCSVFile::readFile(const std::string &infile, const bool& sequentially) {
    size_t hash_count = 0;
    std::string line;

    std::ifstream file(infile);

    if (!file.is_open()) {
        throw IOException("Unable to open file " + infile, AT);
    }

    std::getline(file, firstline);
    if (!isValidFirstLine(firstline)) {
        throw IOException("Invalid first line: " + firstline, AT);
    }
    hash_count++;

    std::string section( "meta" );

    // In your main function or where the while loop was
    while (std::getline(file, line)) {
        IOUtils::trim(line);
        hash_count++;
        if (!processLine(line, section, sequentially))
            break;
        ;
    }
    file.close();
    skip_lines_to_data = hash_count;
}

bool iCSVFile::processLine(const std::string &line, std::string &section, const bool &sequentially) {
    if (line == "#")
        return true;
    if (line.find("[METADATA]") != std::string::npos) {
        section = "meta";
        return true;
    }
    if (line.find("[DATA]") != std::string::npos) {
        section = "data";
        if (METADATA.field_delimiter == '0' || FIELDS.fields.empty()) {
            throw InvalidFormatException("DATA section has to appear after METADATA and FIELDS sections!", AT);
        }
        findTime();
        findLocation();
        if (sequentially)
            return true;
        return true;
    }
    if (line.find("[FIELDS]") != std::string::npos) {
        section = "fields";
        if (METADATA.field_delimiter == '0') {
            throw InvalidFormatException("FIELDS section has to appear after METADATA section!", AT);
        }
        return true;
    }
    if (line.empty())
        return true;
    if (section == "data" && sequentially)
        return false;
    if (section != "data" && line[0] != '#')
        throw InvalidFormatException("Invalid line: " + line + section, AT);

    const std::string content = section == "data" ? line : line.substr(1);
    processContent(content, section);
    return true;
}

void iCSVFile::processContent(const std::string &content, const std::string &section) {
    std::string key, value;
    if (section != "data") {
        const std::vector<std::string> key_value( IOUtils::split(content, '=') );

        if (key_value.size() != 2) {
            throw IOException("Invalid line: " + content, AT);
        }

        key = key_value[0];
        value = key_value[1];
    }

    if (section == "meta") {
        populateMetaData(key, value);
    } else if (section == "fields") {
        populateFields(key, value);
    } else if (section == "data") {
        processData(content);
    } else {
        throw InvalidFormatException("Invalid iCSV format", AT);
    }
}

void iCSVFile::processData(const std::string &content) {
    std::vector<std::string> row( IOUtils::split(content, METADATA.field_delimiter) );
    Date tmp_date;
    // TODO: need to handle geometry = multiple columns
    if (!location_in_header) {
        locations_in_data.push_back(extractCoordinates(row[location_id]));
        if (locations_in_data.back().isEmpty()) {
            throw IOException("Invalid location: " + row[location_id], AT);
        }
        row[location_id] = std::to_string(getNoData());
    }
    if (!IOUtils::convertString(tmp_date, row[time_id], METADATA.timezone)) {
        std::string err_msg = getTimeZone() == IOUtils::nodata ? ". No timezone provided" : ". Invalid date";
        throw IOException("Cannot parse date: " + row[time_id] + err_msg, AT);
    }
    row[time_id] = std::to_string(getNoData());
    const std::vector<double> row_vals( convertVector(row) );
    row_data.push_back(row_vals);
    dates_in_file.push_back(tmp_date);
}

void iCSVFile::populateMetaData(const std::string &key, const std::string &value) {
    if (key == "field_delimiter") {
        if (value.size() != 1) {
            throw IOException("Invalid field delimiter: " + value, AT);
        }
        METADATA.field_delimiter = value[0];
    } else if (key == "geometry")
        METADATA.geometry = value;
    else if (key == "srid")
        METADATA.srid = value;
    else if (key == "station_id")
        METADATA.station_id = value;
    else if (key == "timestamp_meaning")
        METADATA.timestamp_meaning = value;
    else if (key == "nodata")
        IOUtils::convertString(METADATA.nodata, value);
    else if (key == "timezone")
        IOUtils::convertString(METADATA.timezone, value);
    else if (key == "doi" || key == "reference")
        METADATA.doi = value;
    else if (key == "slope_angle")
        IOUtils::convertString(station_location.slope_angle, value);
    else if (key == "slope_azi")
        IOUtils::convertString(station_location.slope_azi, value);
    else
        METADATA.optional_metadata[key] = value;
}

void iCSVFile::populateFields(const std::string &key, const std::string &value) {
    if (key == "fields")
        FIELDS.fields = IOUtils::split(value, METADATA.field_delimiter);
    else if (key == "units_multipliers") {
        const std::vector<std::string> str_units_multipliers( IOUtils::split(value, METADATA.field_delimiter) );
        FIELDS.units_multipliers = convertVector(str_units_multipliers);
    } else if (key == "units_offsets") {
        const std::vector<std::string> str_units_offsets( IOUtils::split(value, METADATA.field_delimiter) );
        FIELDS.units_offsets = convertVector(str_units_offsets);
    } else if (key == "units")
        FIELDS.units = IOUtils::split(value, METADATA.field_delimiter);
    else if (key == "long_name")
        FIELDS.long_name = IOUtils::split(value, METADATA.field_delimiter);
    else if (key == "standard_name")
        FIELDS.standard_name = IOUtils::split(value, METADATA.field_delimiter);
    else
        FIELDS.other_fields[key] = IOUtils::split(value, METADATA.field_delimiter);
}

/**
* Reads the data for a given date and fieldname from the iCSV file.
*
* @param[in] r_date The date for which to read the data.
* @param[in] fieldname The name of the field for which to read the data.
* @return The data value for the given date and fieldname.
* @throws IOException If there is no data available or if the date is out of range.
*/
double iCSVFile::readData(const Date &r_date, const std::string &fieldname) {
    if (row_data.empty()) {
        throw IOException("No data available", AT);
    }
    if (r_date < dates_in_file.front() || r_date > dates_in_file.back()) {
        throw IOException("Date out of range", AT);
    }
    for (size_t ii = 0; ii < dates_in_file.size(); ii++) {
        if (dates_in_file[ii] == r_date) {
            return row_data[ii][FIELDS.getFieldIndex(fieldname)];
        }
    }
    return IOUtils::nodata;
}

/**
* Aggregates data from the given vector of MeteoData objects.
*
* @param[in] vecMeteo The vector of MeteoData objects to aggregate.
* @throws IOException if the vector is empty.
*/
void iCSVFile::aggregateData(const std::vector<MeteoData> &vecMeteo) {
    if (vecMeteo.empty()) {
        throw IOException("No data available", AT);
    }
    const std::set<std::string> available_params( MeteoData::listAvailableParameters(vecMeteo) );
    for (auto md : vecMeteo) {
        md.date.setTimeZone(METADATA.timezone);
        dates_in_file.push_back(md.date);
        if (!location_in_header) {
            const Coords &loc = md.meta.position;
            locations_in_data.push_back(toiCSVLocation(loc, METADATA.epsg));
        }
        std::vector<double> row;
        for (const auto &field : FIELDS.fields) {
            if (field == "timestamp" || field == "julian")
                continue;
            if (field == METADATA.geometry) 
                continue;
            if (field == "OSWR" && available_params.count("RSWR") == 1) {
                row.push_back(md("RSWR"));
            } else if (available_params.count(field) == 0) {
                row.push_back(IOUtils::nodata);
            } else {
                row.push_back(md(field));
            }
        }
        row_data.push_back(row);
    }
}

// ----------------- HEADER Format -----------------
/**
* @brief Checks the validity of the format for the iCSVFile.
*
* @return true if the format is valid, false otherwise.
*/
bool iCSVFile::checkFormatValidity() {
    static const std::string format_type = "iCSV format: ";
    if (!isValidFirstLine(firstline))
        throw InvalidFormatException(format_type + "Invalid first line: " + firstline, AT);
    if (METADATA.field_delimiter == '0')
        throw InvalidFormatException(format_type + "Key: field_delimiter is required", AT);
    else if (std::find(METADATA.valid_delimiters.begin(), METADATA.valid_delimiters.end(), METADATA.field_delimiter) ==
                METADATA.valid_delimiters.end()) {
        throw InvalidFormatException(format_type + "Invalid field delimiter: " + METADATA.field_delimiter, AT);
    }
    if (METADATA.geometry.empty())
        throw InvalidFormatException(format_type + "Key: geometry is required", AT);
    if (METADATA.srid.empty())
        throw InvalidFormatException(format_type + "Key: srid is required", AT);
    if (!isValidEpsgCode(METADATA.srid))
        throw InvalidFormatException(format_type + "Invalid EPSG code: " + METADATA.srid, AT);
    // get the actual EPSG code if it is valid
    IOUtils::convertString(METADATA.epsg, extractEpsgCode(METADATA.srid));

    if (FIELDS.fields.empty())
        throw InvalidFormatException(format_type + "Key: fields is required", AT);

    if (METADATA.epsg == -1)
        throw InvalidFormatException(
            format_type + "Invalid EPSG conversion: " + METADATA.srid + "; " + std::to_string(METADATA.epsg), AT);

    // all fields must have the same number of entries if they are not empty
    if (!FIELDS.units_multipliers.empty() && FIELDS.units_multipliers.size() != FIELDS.fields.size()) {
        throw InvalidFormatException(format_type + "Invalid units_multipliers: " + std::to_string(FIELDS.units_multipliers.size()) +
                                            " != " + std::to_string(FIELDS.fields.size()), AT);
    }
    if (!FIELDS.units_offsets.empty() && FIELDS.units_offsets.size() != FIELDS.fields.size()) {
        throw InvalidFormatException(format_type + "Invalid units_offsets: " + std::to_string(FIELDS.units_offsets.size()) +
                                            " != " + std::to_string(FIELDS.fields.size()), AT);
    }
    if (!FIELDS.units.empty() && FIELDS.units.size() != FIELDS.fields.size()) {
        throw InvalidFormatException(format_type + "Invalid units: " + std::to_string(FIELDS.units.size()) +
                                            " != " + std::to_string(FIELDS.fields.size()), AT);
    }
    if (!FIELDS.long_name.empty() && FIELDS.long_name.size() != FIELDS.fields.size()) {
        throw InvalidFormatException(format_type + "Invalid long_name: " + std::to_string(FIELDS.long_name.size()) +
                                            " != " + std::to_string(FIELDS.fields.size()), AT);
    }
    if (!FIELDS.standard_name.empty() && FIELDS.standard_name.size() != FIELDS.fields.size()) {
        throw InvalidFormatException(format_type + "Invalid standard_name: " + std::to_string(FIELDS.standard_name.size()) +
                                            " != " + std::to_string(FIELDS.fields.size()), AT);
    }
    // iterate through the other_fields and check if they have the same number of entries as fields
    for (const auto &field : FIELDS.other_fields) {
        if (!field.second.empty() && field.second.size() != FIELDS.fields.size()) {
            throw InvalidFormatException(format_type + "Invalid " + field.first + ": " + std::to_string(field.second.size()) +
                                                " != " + std::to_string(FIELDS.fields.size()), AT);
        }
    }

    return true;
}

/**
* Checks the compatibility of the iCSV file with MeteoIO.
*
* @return true if the iCSV file is compatible with MeteoIO, false otherwise.
*/
bool iCSVFile::checkMeteoIOCompatibility() const {
    static const std::string format_type( "MeteoIO format: " );
    if (station_location.isEmpty() && location_in_header ) {
        throw InvalidFormatException(format_type + "Please provide a location for a MeteoStation", AT);
    }
    if (std::find(FIELDS.fields.begin(), FIELDS.fields.end(), "timestamp") == FIELDS.fields.end() &&
        std::find(FIELDS.fields.begin(), FIELDS.fields.end(), "julian") == FIELDS.fields.end()) {
        throw InvalidFormatException(format_type + "Please provide a ISO86 timestamps or julian dates", AT);
    }
    return true;
}

// ----------------- iCSVFile Helperss -----------------
bool iCSVFile::isValidLocation(const geoLocation& location) { return location.x != getNoData() && location.y != getNoData(); }

bool iCSVFile::isColumnName(const std::string &geometry) { // TODO: if column name can be multiple columns need to handle that
    return std::find(FIELDS.fields.begin(), FIELDS.fields.end(), geometry) != FIELDS.fields.end();
}

/**
* Finds parameters, that are used in MeteoIO but not given in a iCSV file.
*
* @param[in] vecMeteo The vector of MeteoData objects to compare with.
* @return A vector of strings containing the parameters to append.
* @throws IOException if the vector is empty.
*/
std::vector<std::string> iCSVFile::columnsToAppend(const std::vector<MeteoData> &vecMeteo) const {
    if (vecMeteo.empty()) {
        throw IOException("No data available", AT);
    }
    const std::set<std::string> available_params( MeteoData::listAvailableParameters(vecMeteo) );

    auto contains = [&](const std::string &s) {
        return std::find(FIELDS.fields.begin(), FIELDS.fields.end(), s) != FIELDS.fields.end();
    };

    std::vector<std::string> columns_to_append;

    // check if all available params are in the fields
    for (const auto &param : available_params) {
        // if the param is "RSWR" both "RSWR" and "OSWR" can be in the fields
        if (param == "RSWR" && !contains("RSWR") && !contains("OSWR")) {
            columns_to_append.push_back("RSWR");
        } else if (!contains(param)) {
            columns_to_append.push_back(param);
        }
    }
    return columns_to_append;
}

/**
* @brief Parses the geometry information of the iCSV file.
*
* This function extracts the coordinates from the metadata and checks if they are valid.
* If the coordinates are empty or not a known WKTS station location, an exception is thrown.
* If the z-coordinate is set to the default nodata value, it is updated to the actual nodata value.
* Finally, the station location is updated with the parsed coordinates.
*
* @throws IOException if the geometry is neither a column name nor a known WKTS station location,
* or if the location is invalid.
*/
void iCSVFile::parseGeometry() {
    geoLocation coordinates( extractCoordinates(METADATA.geometry) );
    if (coordinates.isEmpty()) {
        if (!isColumnName(METADATA.geometry)) {
            throw IOException("Geometry is neither a column name, nor a known WKTS station location (POINT/POINTZ)", AT);
        }
        location_in_header = false;
        station_location = coordinates;
    }

    if (coordinates.z == default_nodata) {
        coordinates.z = getNoData();
    }

    if (!isValidLocation(coordinates)) {
        throw IOException("Invalid location: " + coordinates.toString(), AT);
    }

    station_location = coordinates;
}

void iCSVFile::findTime() {
    if (std::find(FIELDS.fields.begin(), FIELDS.fields.end(), "timestamp") == FIELDS.fields.end()) {
        timestamp_present = false;
    } else {
        timestamp_present = true;
    }
    if (std::find(FIELDS.fields.begin(), FIELDS.fields.end(), "julian") == FIELDS.fields.end()) {
        julian_present = false;
    } else {
        julian_present = true;
    }
    if (!timestamp_present && !julian_present) {
        throw InvalidFormatException("MeteoIO format: No timestamp or julian field found", AT);
    }
    time_id = timestamp_present ? FIELDS.getFieldIndex("timestamp") : FIELDS.getFieldIndex("julian");
}

void iCSVFile::findLocation() {
    if (location_in_header)
        return;
    // TODO: if geometry can be multiple columns, needs parsing
    if (std::find(FIELDS.fields.begin(), FIELDS.fields.end(), METADATA.geometry) == FIELDS.fields.end()) {
        throw InvalidFormatException("MeteoIO format: No geometry field found", AT);
    }
    location_id = FIELDS.getFieldIndex(METADATA.geometry);
}

// ----------------- OPERATORS -----------------
bool operator==(const geoLocation &lhs, const geoLocation &rhs) { return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z; }
bool operator!=(const geoLocation &lhs, const geoLocation &rhs) { return !(lhs == rhs); }

geoLocation &geoLocation::operator=(const geoLocation &rhs) {
    this->x = rhs.x;
    this->y = rhs.y;
    this->z = rhs.z;
    this->slope_angle = rhs.slope_angle;
    this->slope_azi = rhs.slope_azi;
    return *this;
}

geoLocation::geoLocation(const geoLocation &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->slope_angle = other.slope_angle;
    this->slope_azi = other.slope_azi;
}

geoLocation::geoLocation(const double& in_x, const double& in_y, const double& in_z, const double& in_slope_angle, const double& in_slope_azi)
    : x(in_x), y(in_y), z(in_z), slope_angle(in_slope_angle), slope_azi(in_slope_azi) {}

bool operator==(const MetaDataSection &lhs, const MetaDataSection &rhs) {
    if (lhs.field_delimiter != rhs.field_delimiter || lhs.geometry != rhs.geometry || lhs.srid != rhs.srid ||
        lhs.epsg != rhs.epsg) {
        return false;
    }

    if (lhs.station_id != rhs.station_id && !lhs.station_id.empty() && !rhs.station_id.empty()) {
        return false;
    }

    if (lhs.timestamp_meaning != rhs.timestamp_meaning && !lhs.timestamp_meaning.empty() && !rhs.timestamp_meaning.empty()) {
        return false;
    }

    if (lhs.nodata != rhs.nodata && lhs.nodata != default_nodata && rhs.nodata != default_nodata) {
        return false;
    }

    if (lhs.timezone != rhs.timezone && lhs.timezone != default_nodata && rhs.timezone != default_nodata) {
        return false;
    }

    if (lhs.doi != rhs.doi && !lhs.doi.empty() && !rhs.doi.empty()) {
        return false;
    }

    for (const auto &pair : lhs.optional_metadata) {
        auto found = rhs.optional_metadata.find(pair.first);
        if (found->second != pair.second) {
            return false;
        }
    }
    return true;
}

bool roughlyEqual(const MetaDataSection &lhs, const MetaDataSection &rhs) {
    if (lhs.field_delimiter != rhs.field_delimiter && lhs.field_delimiter != '0' && rhs.field_delimiter != '0') {
    }
    if (lhs.geometry != rhs.geometry && !lhs.geometry.empty() && !rhs.geometry.empty()) {
    }
    if (lhs.srid != rhs.srid && !lhs.srid.empty() && !rhs.srid.empty()) {
    }
    if (lhs.epsg != rhs.epsg && lhs.epsg != -1 && rhs.epsg != -1) {
        return false;
    }

    if (lhs.station_id != rhs.station_id && !lhs.station_id.empty() && !rhs.station_id.empty()) {
        return false;
    }

    if (lhs.timestamp_meaning != rhs.timestamp_meaning && !lhs.timestamp_meaning.empty() && !rhs.timestamp_meaning.empty()) {
        return false;
    }

    if (lhs.nodata != rhs.nodata && lhs.nodata != default_nodata && rhs.nodata != default_nodata) {
        return false;
    }

    if (lhs.timezone != rhs.timezone && lhs.timezone != default_nodata && rhs.timezone != default_nodata) {
        return false;
    }

    if (lhs.doi != rhs.doi && !lhs.doi.empty() && !rhs.doi.empty()) {
        return false;
    }

    for (const auto &pair : lhs.optional_metadata) {
        const auto& found = rhs.optional_metadata.find(pair.first);
        if (found->second != pair.second) {
            return false;
        }
    }
    return true;
}

MetaDataSection &operator+=(MetaDataSection &lhs, const MetaDataSection &rhs) {
    if (roughlyEqual(lhs, rhs)) {
        if (lhs.field_delimiter == '0' && rhs.field_delimiter != '0') {
            lhs.field_delimiter = rhs.field_delimiter;
        }
        if (lhs.geometry.empty() && !rhs.geometry.empty()) {
            lhs.geometry = rhs.geometry;
        }
        if (lhs.srid.empty() && !rhs.srid.empty()) {
            lhs.srid = rhs.srid;
        }
        if (lhs.epsg == -1 && rhs.epsg != -1) {
            lhs.epsg = rhs.epsg;
        }
        if (lhs.station_id.empty() && !rhs.station_id.empty()) {
            lhs.station_id = rhs.station_id;
        }
        if (lhs.timestamp_meaning.empty() && !rhs.timestamp_meaning.empty()) {
            lhs.timestamp_meaning = rhs.timestamp_meaning;
        }
        if (lhs.nodata == default_nodata && rhs.nodata != default_nodata) {
            lhs.nodata = rhs.nodata;
        }
        if (lhs.timezone == default_nodata && rhs.timezone != default_nodata) {
            lhs.timezone = rhs.timezone;
        }
        if (lhs.doi.empty() && !rhs.doi.empty()) {
            lhs.doi = rhs.doi;
        }

        for (const auto &pair : rhs.optional_metadata) {
            lhs.optional_metadata[pair.first] = pair.second;
        }
    } else {
        throw IOException("Cannot merge incompatible MetaDataSections in iCSVFiles");
    }

    return lhs;
}

bool operator==(const fieldsSection &lhs, const fieldsSection &rhs) {
    if (lhs.fields != rhs.fields) {
        return false;
    }

    if (lhs.units_multipliers != rhs.units_multipliers && !lhs.units_multipliers.empty() && !rhs.units_multipliers.empty()) {
        return false;
    }

    if (lhs.units_offsets != rhs.units_offsets && !lhs.units_offsets.empty() && !rhs.units_offsets.empty()) {
        return false;
    }
    if (lhs.units != rhs.units && !lhs.units.empty() && !rhs.units.empty()) {
        return false;
    }
    if (lhs.long_name != rhs.long_name && !lhs.long_name.empty() && !rhs.long_name.empty()) {
        return false;
    }
    if (lhs.standard_name != rhs.standard_name && !lhs.standard_name.empty() && !rhs.standard_name.empty()) {
        return false;
    }

    for (const auto &pair : lhs.other_fields) {
        const auto& found = rhs.other_fields.find(pair.first);
        if (found->second != pair.second) {
            return false;
        }
    }
    return true;
}

bool roughlyEqual(const fieldsSection &lhs, const fieldsSection &rhs) {
    if (lhs.fields != rhs.fields && !lhs.fields.empty() && !rhs.fields.empty()) {
        return false;
    }
    if (lhs.units_multipliers != rhs.units_multipliers && !lhs.units_multipliers.empty() && !rhs.units_multipliers.empty()) {
        return false;
    }

    if (lhs.units_offsets != rhs.units_offsets && !lhs.units_offsets.empty() && !rhs.units_offsets.empty()) {
        return false;
    }
    if (lhs.units != rhs.units && !lhs.units.empty() && !rhs.units.empty()) {
        return false;
    }
    if (lhs.long_name != rhs.long_name && !lhs.long_name.empty() && !rhs.long_name.empty()) {
        return false;
    }
    if (lhs.standard_name != rhs.standard_name && !lhs.standard_name.empty() && !rhs.standard_name.empty()) {
        return false;
    }

    for (const auto &pair : lhs.other_fields) {
        const auto& found = rhs.other_fields.find(pair.first);
        if (found->second != pair.second) {
            return false;
        }
    }
    return true;
}

fieldsSection &operator+=(fieldsSection &lhs, const fieldsSection &rhs) {
    if (roughlyEqual(lhs, rhs)) {
        if (lhs.fields.empty() && !rhs.fields.empty()) {
            lhs.fields = rhs.fields;
        }
        if (lhs.units_multipliers.empty() && !rhs.units_multipliers.empty()) {
            lhs.units_multipliers = rhs.units_multipliers;
        }
        if (lhs.units_offsets.empty() && !rhs.units_offsets.empty()) {
            lhs.units_offsets = rhs.units_offsets;
        }
        if (lhs.units.empty() && !rhs.units.empty()) {
            lhs.units = rhs.units;
        }
        if (lhs.long_name.empty() && !rhs.long_name.empty()) {
            lhs.long_name = rhs.long_name;
        }
        if (lhs.standard_name.empty() && !rhs.standard_name.empty()) {
            lhs.standard_name = rhs.standard_name;
        }

        for (const auto &pair : rhs.other_fields) {
            lhs.other_fields[pair.first] = pair.second;
        }
    } else {
        throw IOException("Cannot merge incompatible fieldsSections in iCSVFiles");
    }

    return lhs;
}
} // namespace iCSV
} // namespace mio
