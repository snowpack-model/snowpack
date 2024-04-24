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
#ifndef iCSVHELP_H
#define iCSVHELP_H

#include <meteoio/IOUtils.h>

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

static const double default_nodata = M_PI;
namespace mio {
namespace iCSV {


// ----------------- Datat structures -----------------
/**
    * @brief Represents a geographic location extracted from iCSV files. (POINT(X Y)...)
    */
struct geoLocation {
    double x = default_nodata; /**< The x-coordinate of the location. */
    double y = default_nodata; /**< The y-coordinate of the location. */
    double z = default_nodata; /**< The z-coordinate of the location. */
    double slope_angle = default_nodata; /**< The slope angle at the location. */
    double slope_azi = default_nodata; /**< The slope azimuth at the location. */

    /**
        * @brief Checks if the geoLocation is empty.
        * @return True if the geoLocation is empty, false otherwise.
        */
    bool isEmpty() const { return x == default_nodata && y == default_nodata && z == default_nodata; }

    /**
        * @brief Converts the geoLocation to a string representation.
        * @return The string representation of the geoLocation.
        */
    std::string toString() const {
        std::stringstream ss;
        ss << "x: " << x << " y: " << y << " z: " << z
            << (slope_angle != default_nodata ? " slope_angle: " + std::to_string(slope_angle) : "")
            << (slope_azi != default_nodata ? " slope_azi: " + std::to_string(slope_azi) : "");
        return ss.str();
    }

    void standardizeNodata() {
        if (x == default_nodata) {
            x = IOUtils::nodata;
        }
        if (y == default_nodata) {
            y = IOUtils::nodata;
        }
        if (z == default_nodata) {
            z = IOUtils::nodata;
        }
    }

    Coords toCoords(const int& epsg) const {
        Coords loc;
        loc.setPoint(x, y, z, epsg);
        return loc;
    }

    geoLocation& operator=(const geoLocation& rhs);
    geoLocation(const geoLocation& rhs);
    geoLocation() = default;
    geoLocation(const double& in_x, const double& in_y, const double& in_z, const double& in_slope_angle=default_nodata, const double& in_slope_azi=default_nodata);

};

bool operator==(const geoLocation& lhs, const geoLocation& rhs);
bool operator!=(const geoLocation& lhs, const geoLocation& rhs);

/**
    * @brief Stores information on the METADATA section of a iCSV file.
    */
struct MetaDataSection {
    char field_delimiter = '0'; /**< The field delimiter character. */
    std::string geometry = ""; /**< The geometry information. */
    std::string srid = ""; /**< The spatial reference system identifier. */
    int epsg = -1; /**< The EPSG code. */
    std::string station_id = ""; /**< The station identifier. */
    std::string timestamp_meaning = ""; /**< The meaning of the timestamp. */
    double nodata = default_nodata; /**< The nodata value. */
    double timezone = default_nodata; /**< The timezone value. */
    std::string doi = ""; /**< The Digital Object Identifier. */

    std::map<std::string, std::string> optional_metadata = {}; /**< The optional metadata. */

    const std::vector<char> valid_delimiters = {',', '/', '\\', '|', ':', ';'}; /**< The valid delimiters. */

    /**
        * @brief Set the EPSG code and update the SRID.
        * @param in_epsg The EPSG code.
        */
    void setEPSG(const int &in_epsg) {
        epsg = in_epsg;
        srid = "EPSG:" + std::to_string(in_epsg);
    }

    /**
        * @brief Get the value of an optional metadata field.
        *
        * @param key The key of the optional metadata field.
        * @return The value of the optional metadata field.
        */
    std::string getOptionalMetaData(const std::string &key) const {
        const auto it = optional_metadata.find(key);
        if (it != optional_metadata.end()) {
            return it->second;
        }
        return "";
    }

    /**
        * @brief Convert the metadata section to a map.
        *
        * @return The metadata section as a map.
        */
    std::map<std::string, std::string> toMetaMap() const {
        std::map<std::string, std::string> dict;
        if (!doi.empty()) {
            dict["doi"] = doi;
        }

        // Append optional_metadata to dict
        dict.insert(optional_metadata.begin(), optional_metadata.end());

        return dict;
    }

    /**
        * @brief Convert the metadata section to an output map.
        *
        * @return The metadata section as an output map.
        */
    std::map<std::string, std::string> toOutputMap() const {
        std::map<std::string, std::string> dict;
        dict["field_delimiter"] = field_delimiter;
        dict["geometry"] = geometry;
        dict["srid"] = srid;
        if (!station_id.empty()) {
            dict["station_id"] = station_id;
        }
        if (!timestamp_meaning.empty()) {
            dict["timestamp_meaning"] = timestamp_meaning;
        }
        if (nodata != default_nodata) {
            dict["nodata"] = std::to_string(nodata);
        }
        if (timezone != default_nodata) {
            dict["timezone"] = std::to_string(timezone);
        }
        if (!doi.empty()) {
            dict["doi"] = doi;
        }

        dict.insert(optional_metadata.begin(), optional_metadata.end());
        return dict;
    }
};

bool operator==(const MetaDataSection &lhs, const MetaDataSection &rhs);

/**
    * @brief Checks if two MetaDataSection objects are roughly equal.
    *
    * @param lhs The first MetaDataSection object.
    * @param rhs The second MetaDataSection object.
    * @return true if the objects are roughly equal, false otherwise.
    */
bool roughlyEqual(const MetaDataSection &lhs, const MetaDataSection &rhs);

MetaDataSection& operator+=(MetaDataSection& lhs, const MetaDataSection& rhs);

/**
    * @brief Stores Information on the FIELDS section of a iCSV file.
    *
    */
struct fieldsSection {
    std::vector<std::string> fields = {};
    std::vector<double> units_multipliers = {};
    std::vector<double> units_offsets = {};
    std::vector<std::string> units = {};
    std::vector<std::string> long_name = {};
    std::vector<std::string> standard_name = {};

    std::map<std::string, std::vector<std::string>> other_fields = {};

    size_t getFieldIndex(const std::string &fieldname) const {
        for (size_t id = 0; id < fields.size(); id++) {
            if (fieldname == fields[id])
                return id;
        }
        return IOUtils::npos;
    }

    std::vector<std::string> getOtherFields(const std::string &key) {
        auto it = other_fields.find(key);
        if (it != other_fields.end()) {
            return it->second;
        }
        return {};
    }
    std::map<std::string, std::vector<std::string>> toMap() const {
        std::map<std::string, std::vector<std::string>> dict;
        dict["fields"] = fields;

        if (!units.empty()) {
            dict["units"] = units;
        }
        if (!long_name.empty()) {
            dict["long_name"] = long_name;
        }
        if (!standard_name.empty()) {
            dict["standard_name"] = standard_name;
        }

        // Convert vector<double> to vector<string>
        if (!units_multipliers.empty()) {
            std::vector<std::string> str_units_multipliers(units_multipliers.size());
            std::transform(units_multipliers.begin(), units_multipliers.end(), str_units_multipliers.begin(),
                            [](double d) { return std::to_string(d); });
            dict["units_multipliers"] = str_units_multipliers;
        }

        if (!units_offsets.empty()) {
            std::vector<std::string> str_units_offsets(units_offsets.size());
            std::transform(units_offsets.begin(), units_offsets.end(), str_units_offsets.begin(),
                            [](double d) { return std::to_string(d); });
            dict["units_offsets"] = str_units_offsets;
        }

        // Append other_fields to dict
        if (!other_fields.empty()) {
            dict.insert(other_fields.begin(), other_fields.end());
        }
        return dict;
    }
};

bool operator==(const fieldsSection &lhs, const fieldsSection &rhs);

/**
    * @brief Checks if two fieldsSection objects are roughly equal.
    *
    * @param lhs The first MetaDataSection object.
    * @param rhs The second MetaDataSection object.
    * @return true if the objects are roughly equal, false otherwise.
    */
bool roughlyEqual(const fieldsSection &lhs, const fieldsSection &rhs);
fieldsSection& operator+=(fieldsSection& lhs, const fieldsSection& rhs);

// ----------------- Helper functions -----------------
std::vector<double> convertVector(const std::vector<std::string> &vec);
std::vector<Coords> convertVector(const std::vector<geoLocation> &vec, const int& epsg);
geoLocation extractCoordinates(const std::string &geometry);
geoLocation toiCSVLocation(Coords coords, const int& epsg);



// ----------------- iCSVFile class -----------------

/**
* @class iCSVFile
*
* @brief This class is responsible for handling and storing iCSV files.
*
* General information about the file is stored as public members, the station_location is an additional member.
* The MetaData is stored in the METADATA member, and the fields are stored in the FIELDS member.
*
* It is possible to store all the data as well with read_sequential= false
*
* Format checking is done on initialization.
*
* Getters are provided to access the data, Metadata is directly accessed.
*
*/
class iCSVFile {
    public:
        // Essential information
        size_t skip_lines_to_data;
        std::string filename;
        std::string firstline;
        geoLocation station_location;

        MetaDataSection METADATA;
        fieldsSection FIELDS;

        // helper flags
        bool location_in_header = true;
        bool timezone_in_data = false;
        bool timestamp_present = false;
        bool julian_present = false;
        size_t time_id = IOUtils::npos;
        size_t location_id = IOUtils::npos;

    private:
        // File data
        std::vector<Date> dates_in_file;
        std::vector<std::vector<double>> row_data;
        std::vector<geoLocation> locations_in_data;


    public:
        // Constructors
        iCSVFile();
        iCSVFile(const iCSVFile &);
        iCSVFile(const std::string &infile, const bool& read_sequential);

        // Main methods
        void readFile(const std::string &infile, const bool& sequentially);
        bool checkFormatValidity();
        bool checkMeteoIOCompatibility() const;
        void parseGeometry();

        // Getters
        double getNoData() const {
            double nodat = METADATA.nodata == default_nodata ? IOUtils::nodata : METADATA.nodata;
            return nodat;
        }
        double getTimeZone() const { return METADATA.timezone == default_nodata ? IOUtils::nodata : METADATA.timezone; }
        std::vector<std::vector<double>> getRowData() const { return row_data; }
        std::vector<Date> getAllDatesInFile() const { return dates_in_file; }
        std::vector<geoLocation> getAllLocationsInData() const { return locations_in_data; }
        geoLocation getLocationAt(size_t index) const {
            if (index >= locations_in_data.size()) {
                throw std::out_of_range("Index out of range for locations");
            }
            return locations_in_data[index];
        }
        std::vector<Date> getDatesInFile(const Date &start_date, const Date &end_date) const {
            std::vector<Date> dates;
            for (const Date &date : dates_in_file) {
                if (date >= start_date && date <= end_date) {
                    dates.push_back(date);
                }
            }
            return dates;
        }
        std::vector<geoLocation> getLocationsInData(const Date &start_date, const Date &end_date) const {
            if (location_in_header)
                return std::vector<geoLocation>();

            std::vector<geoLocation> locations;
            for (size_t i = 0; i < dates_in_file.size(); i++) {
                if (dates_in_file[i] >= start_date && dates_in_file[i] <= end_date) {
                    locations.push_back(locations_in_data[i]);
                }
            }
            return locations;
        }
        // reading helpers
        bool processLine(const std::string& line, std::string& section, const bool& sequentially);
        void processContent(const std::string& content, const std::string& section);
        void processData(const std::string& content);

        // Format helpers
        bool isValidLocation(const geoLocation& location);
        bool isColumnName(const std::string &geometry);
        void findTime();
        void findLocation();
        std::vector<std::string> columnsToAppend(const std::vector<MeteoData> &vecMeteo) const;

        // Data methods
        double readData(const Date &r_date, const std::string &fieldname);
        void aggregateData(const std::vector<MeteoData> &vecMeteo);
        void populateMetaData(const std::string &key, const std::string &value);
        void populateFields(const std::string &key, const std::string &value);
};

} // namespace iCSV
} // namespace mio
#endif // iCSVHELP_H
