// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

/* https://codes.ecmwf.int/grib/param-db/?ordering=id-- parameter database

*/

/*
NOTES:
- needs documentatino
*/
#include "meteoio/IOUtils.h"
#include <meteoio/FStream.h>
#include <meteoio/FileUtils.h>
#include <meteoio/plugins/libcodes.h>

#include <cstring>
#include <string>

namespace mio {
    /**
     * @namespace codes
     * @brief This namespace handles all the low level manipulation of GRIB and BUFR files with ecCodes
     *
     * @section codes_intro Introduction
     * This namespace provides a set of functions to handle GRIB and BUFR files using the ecCodes C library by the ECMWF (https://confluence.ecmwf.int/display/ECC/ecCodes+installation) at least in
     * version 2.27.0. Either a file is indexed (contained in CodesIndexPtr) and messages (contained as CodesHandlePtr) are read from this index by using key/value pairs. Or all messages are read from
     * a file directly, without any additional information. The index does not use much memory, the handles do??. Therefore, the handles should be deleted as soon as they are not needed anymore.
     *
     *
     * @ingroup plugins
     * @author Patrick Leibersperge
     * @date   2024-04-17
     *
     */
    namespace codes {
        const int WMO_BUFR_TABLE_NO = 41;

        // ------------------- BUFR CONSTANTS -------------------

        const std::string BUFR_HEIGHT_KEY = "height"; // height // verticalDistanceOfSensor is a descriptor that will be coming soon
        // mapping of BUFR parameters to MeteoIO parameters
        const std::map<std::string, std::string> BUFR_PARAMETER{
            {"P", "pressure"},
            {"TA", "airTemperature"}, // Can also be found as airTemperature under heightOfSensor(somthing like this) = 2
            {"RH", "relativeHumidity"},
            {"TSG", "groundTemperature"},
            {"TSOIL", "soilTemperature"},
            {"TSNOW", "snowTemperature"},
            {"HS", "totalSnowDepth"}, // totalSnowDepth // snowDepth will be available soon
            {"VW", "windSpeed"},
            {"DW", "windDirection"},
            {"VW_MAX", "maximumWindGustSpeed"},
            {"RSWR", ""}, // TODO: Reflected shortwave radiation when its available
            {"ISWR", "instantaneousShortWaveRadiation"},
            {"ILWR", "instantaneousLongWaveRadiation"},
            {"TAU_CLD", "cloudCoverTotal"},
            {"PSUM", "totalPrecipitationOrTotalWaterEquivalent"},
            {"PSUM_PH", "precipitationType"}, // should the type be mapped to the phase?
            // Cryo specific
            {"SURFACEQUALIFIER","surfaceQualifierForTemperatureData"},
            {"TSS","skinTemperature"},
            {"ICE_THICKNESS","iceThickness"},
            {"GROUNDSTATE","stateOfGround"},
            {"SENSORTYPE","temperatureSensorType"},
            {"TICE","iceSurfaceTemperature"}, // iceSurfaceTemperature // iceTemperature will be available soon
            {"TWATER","waterTemperature"},
            {"T","airTemperature"} // To allow for general temperatures
        };

        // an alternative mapping of BUFR parameters to MeteoIO parameters
        const std::map<std::string, std::string> BUFR_PARAMETER_ALT{{"TA", "airTemperatureAt2m"}};

        // flags for the possible reference systems are 0-4
        const std::vector<int> FLAG_TO_EPSG = {4326, 4258, 4269, 4314};

        // Surface qualifier for snow
        const long SNOW_SURFACE_QUALIFIER = 6;

        // TODO: exchange the cryos Parameters when available
        // the descriptors should only be available in this file
        static const std::map<std::string, long> BUFR_DESCRIPTORS{
            {"Date", 301011},
            {"Time", 301013},
            {"StationId", 1018},
            {"StationNumber", 1002},
            {"StationName", 1015},
            {"LatLon", 301021},
            {"Alt", 7030},
            {"TA2m", 12004},
            {"TSG", 12120},
            {"TSOIL", 12130},
            {"HS", 13013},// 13013 // 13119 will be available soon but is not yet in the wmo table
            {"VW", 11012},
            {"DW", 11011},
            {"VW_MAX", 11041},
            {"ILWR", 14017},
            {"ISWR", 14018},
            {"TAU_CLD", 20010},
            {"PSUM", 13011},
            {"PSUM_PH", 20021},
            {"TSNOW", 12131},
            {"TSS", 12161}, // we need to set this with surface qualifier and skin temperature, as there is no snow surface temperature kw
            {"P", 7004},
            {"TA", 12101},
            {"RH", 13003},
            {"HSensor", 7007}, // 7007 // 7034 will be available soon
            // Cryos specific
            {"WIGOS_ID", 301150},
            {"LongStationName", 1019},
            {"StationType", 2001},
            {"SurfType", 8029},
            {"IceThickness", 13115},
            {"GroundState", 20062},
            {"SnoDepthMethod", 2177},
            {"SurfQualifier", 8010},
            {"SensorType", 2096},
            {"ChangeDataWidth", 201131},
            {"ChangeScale", 202129},
            {"TIce", 12132}, // 12132// 12133 will be available soon
            {"TWater", 13082}};

        static long getDescriptor(const std::string &key) {
            if (BUFR_DESCRIPTORS.find(key) == BUFR_DESCRIPTORS.end()) {
                throw InvalidArgumentException("Descriptor for " + key + " not found", AT);
            }
            return BUFR_DESCRIPTORS.at(key);
        }

        static std::vector<long> getDescriptors(std::vector<std::string> keys) {
            std::vector<long> descriptors;
            for (const std::string &key : keys) {
                descriptors.push_back(getDescriptor(key));
            }
            return descriptors;
        }

        // ------------------- GRIB CONSTANTS-------------------
        // TODO: Find out how to incorporate RSWR, PSUM_PH
        const long NOPARAMID = -1;
        const std::map<std::string, long> GRIB_DEFAULT_PARAM_TABLE{{"DEM", 129}, {"P", 134},       {"TA", 167},  {"TSS", 228},        {"TSG", 235},          {"RH", 157},        {"PSUM", 228},
                                                                   {"HS", 3066}, {"RSNO", 33},     {"SWE", 141}, {"ISWR", 169},       {"ILWR", 175},         {"VW_MAX", 201187}, {"DW", 3031},
                                                                   {"VW", 10},   {"TAU_CLD", 164}, {"ROT", 205}, {"ALB", 243},        {"SLOP", 163},         {"AZI", 162},       {"W", 135},
                                                                   {"V", 166},   {"U", 165},       {"TD", 3017}, {"RSWR", NOPARAMID}, {"PSUM_PH", NOPARAMID}};

        const std::map<std::string, std::string> GRIB_DEFAULT_LEVELTYPE_TABLE{{"DEM", "surface"},
                                                                              {"P", "surface"},
                                                                              {"TA", "surface"},
                                                                              {"TSS", "surface"},
                                                                              {"TSG", "surface"},
                                                                              {"RH", "heightAboveGround"},
                                                                              {"PSUM", "surface"},
                                                                              {"HS", "surface"},
                                                                              {"RSNO", "surface"},
                                                                              {"SWE", "surface"},
                                                                              {"ISWR", "surface"},
                                                                              {"ILWR", "surface"},
                                                                              {"VW_MAX", "heightAboveGround"},
                                                                              {"DW", "heightAboveGround"},
                                                                              {"VW", "heightAboveGround"},
                                                                              {"TAU_CLD", "surface"},
                                                                              {"ROT", "surface"},
                                                                              {"ALB", "surface"},
                                                                              {"SLOP", "surface"},
                                                                              {"AZI", "surface"},
                                                                              {"W", "hybrid"},
                                                                              {"V", "surface"},
                                                                              {"U", "surface"},
                                                                              {"TD", "heightAboveGround"},
                                                                              {"RSWR", ""},
                                                                              {"PSUM_PH", ""}};

        const std::map<std::string, long> GRIB_DEFAULT_LEVELNO_TABLE{{"DEM", 0},  {"P", 0},    {"TA", 0},   {"TSS", 0},     {"TSG", 0}, {"RH", 2},  {"PSUM", 0},    {"HS", 3066},  {"RSNO", 0},
                                                                     {"SWE", 0},  {"ISWR", 0}, {"ILWR", 0}, {"VW_MAX", 10}, {"DW", 10}, {"VW", 10}, {"TAU_CLD", 0}, {"ROT", 0},    {"ALB", 0},
                                                                     {"SLOP", 0}, {"AZI", 0},  {"W", 10},   {"V", 0},       {"U", 0},   {"TD", 2},  {"RSWR", 0},    {"PSUM_PH", 0}};

        // ------------------- POINTER HANDLING -------------------
        CodesHandlePtr makeUnique(codes_handle *h) {
            CodesHandlePtr ptr(h);
            h = nullptr;
            return ptr;
        }
        CodesIndexPtr makeUnique(codes_index *i) {
            CodesIndexPtr ptr(i);
            i = nullptr;
            return ptr;
        }

        // ------------------- FILE HANDLING -------------------
        // Index a file by a given set of keys, return the index
        CodesIndexPtr indexFile(const std::string &filename, const std::vector<std::string> &index_keys, bool verbose) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            // the keys have to be concatenated to a single string as "key1,key2,key3"
            std::string indexing_string;
            for (const std::string &key : index_keys) {
                if (!indexing_string.empty()) {
                    indexing_string += ",";
                }
                indexing_string += key;
            }

            int ret;
            // create the index
            codes_index *index = codes_index_new_from_file(0, filename.c_str(), indexing_string.c_str(), &ret);

            CODES_CHECK(ret, 0);

            if (verbose) {
                std::cerr << "Indexing " << filename << " with keys " << indexing_string << "\n";
                for (const std::string &key : index_keys) {
                    size_t size;
                    CODES_CHECK(codes_index_get_size(index, key.c_str(), &size), 0);
                    std::cerr << "Found " << size << " " << key << " in " << filename << "\n";

                    std::vector<char *> values(size);
                    CODES_CHECK(codes_index_get_string(index, key.c_str(), values.data(), &size), 0);

                    std::cerr << "Values:\n";
                    for (char *cstr : values) {
                        if (cstr != nullptr) {
                            std::cerr << cstr << "\n";
                            delete[] cstr;
                        }
                    }
                }
            }
            return makeUnique(index);
        }

        // ------------------- GETTERS -------------------
        // multiple overloads for different parameter types
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, double &parameterValue, const IOUtils::ThrowOptions &throwError) {
            int err = codes_get_double(h.get(), parameterName.c_str(), &parameterValue);
            if (err != 0) {
                if (throwError == IOUtils::dothrow) {
                    std::cout << "Error reading parameter " << parameterName << ": Errno " << err << "\n";
                    CODES_CHECK(err, 0);
                }
                return false;
            }
            if (parameterValue == CODES_MISSING_DOUBLE) {
                parameterValue = IOUtils::nodata;
            }
            return true;
        }
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, long &parameterValue, const IOUtils::ThrowOptions &throwError) {
            int err = codes_get_long(h.get(), parameterName.c_str(), &parameterValue);
            if (err != 0) {
                if (throwError == IOUtils::dothrow) {
                    std::cout << "Error reading parameter " << parameterName << ": Errno " << err << "\n";
                    CODES_CHECK(err, 0);
                }
                return false;
            }
            if (parameterValue == CODES_MISSING_LONG) {
                parameterValue = IOUtils::nodata;
            }
            return true;
        }
        // casts long to int
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, int &parameterValue, const IOUtils::ThrowOptions &throwError) {
            long tmp;
            int err = codes_get_long(h.get(), parameterName.c_str(), &tmp);
            if (err != 0) {
                if (throwError == IOUtils::dothrow) {
                    std::cout << "Error reading parameter " << parameterName << ": Errno " << err << "\n";
                    CODES_CHECK(err, 0);
                }
                return false;
            }
            if (tmp == CODES_MISSING_LONG) {
                std::cout << "Parameter " << parameterName << " is missing\n";
                parameterValue = IOUtils::nodata;
            }
            parameterValue = static_cast<int>(tmp);
            return true;
        }

        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, std::string &parameterValue, const IOUtils::ThrowOptions &throwError) {
            size_t len = 500;
            char name[500] = {'\0'};
            int err = codes_get_string(h.get(), parameterName.c_str(), name, &len);
            if (err != 0) {
                if (throwError == IOUtils::dothrow) {
                    std::cout << "Error reading parameter " << parameterName << ": Errno " << err << "\n";
                    CODES_CHECK(err, 0);
                }
                return false;
            }
            parameterValue = std::string(name);
            return true;
        }

        // ------------------- MESSAGE HANDLING -------------------

        // wrapper that opens a file, reads all messages and returns them as a vector of handles
        std::vector<CodesHandlePtr> getMessages(const std::string &filename, ProductKind product) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            FILE *fp = fopen(filename.c_str(), "r");
            return getMessages(fp, product);
        }

        // wrapper that reads all messages from a file and returns them as a vector of handles, closes on crash
        std::vector<CodesHandlePtr> getMessages(FILE *in_file, ProductKind product) {
            codes_handle *h = nullptr;
            int err = 0;
            std::vector<CodesHandlePtr> handles;
            while ((h = codes_handle_new_from_file(0, in_file, product, &err)) != nullptr) {
                if (!h) {
                    fclose(in_file);
                    throw IOException("Unable to create handle from file", AT);
                }
                if (err != 0) {
                    fclose(in_file);
                    CODES_CHECK(err, 0);
                    codes_handle_delete(h);
                    throw IOException("Error reading message: Errno " + std::to_string(err), AT);
                }
                handles.push_back(makeUnique(h));
            }
            fclose(in_file);
            return handles;
        }

        // ************************
        // GRIB
        // ************************

        // Return the timepoint a message is valid for
        Date getMessageDateGrib(CodesHandlePtr &h, const double &tz_in) {
            Date base;
            long validityDate, validityTime;
            getParameter(h, "validityDate", validityDate);
            getParameter(h, "validityTime", validityTime);

            const int year = static_cast<int>(validityDate / 10000), month = static_cast<int>(validityDate / 100 - year * 100), day = static_cast<int>(validityDate - month * 100 - year * 10000);
            const int hour = static_cast<int>(validityTime / 100), minutes = static_cast<int>(validityTime - hour * 100); // HACK: handle seconds!
            base.setDate(year, month, day, hour, minutes, tz_in);

            return base;
        }

        std::map<std::string, double> getGridParameters(CodesHandlePtr &h_unique) {
            // getting transformation parameters
            long Ni, Nj;
            getParameter(h_unique, "Ni", Ni);
            getParameter(h_unique, "Nj", Nj);

            double angleOfRotationInDegrees, latitudeOfSouthernPole, longitudeOfSouthernPole, latitudeOfNorthernPole, longitudeOfNorthernPole;
            bool success = getParameter(h_unique, "angleOfRotationInDegrees", angleOfRotationInDegrees, IOUtils::nothrow);
            if (!success)
                angleOfRotationInDegrees = 0.; // default value for when it is not rotated

            success = getParameter(h_unique, "latitudeOfSouthernPoleInDegrees", latitudeOfSouthernPole, IOUtils::nothrow);
            if (!success)
                latitudeOfSouthernPole = -90.;
            success = getParameter(h_unique, "longitudeOfSouthernPoleInDegrees", longitudeOfSouthernPole, IOUtils::nothrow);
            if (!success)
                longitudeOfSouthernPole = 0.;

            latitudeOfNorthernPole = -latitudeOfSouthernPole;
            longitudeOfNorthernPole = longitudeOfSouthernPole + 180.;

            // determining llcorner, urcorner and center coordinates
            double ll_latitude, ll_longitude, ur_latitude, ur_longitude;
            getParameter(h_unique, "latitudeOfFirstGridPointInDegrees", ll_latitude);
            getParameter(h_unique, "longitudeOfFirstGridPointInDegrees", ll_longitude);

            getParameter(h_unique, "latitudeOfLastGridPointInDegrees", ur_latitude);
            getParameter(h_unique, "longitudeOfLastGridPointInDegrees", ur_longitude);

            std::map<std::string, double> gridParams = {{"Ni", static_cast<double>(Ni)},
                                                        {"Nj", static_cast<double>(Nj)},
                                                        {"Nj", Nj},
                                                        {"ll_latitude", ll_latitude},
                                                        {"ll_longitude", ll_longitude},
                                                        {"ur_latitude", ur_latitude},
                                                        {"ur_longitude", ur_longitude},
                                                        {"angleOfRotationInDegrees", angleOfRotationInDegrees},
                                                        {"latitudeOfSouthernPole", latitudeOfSouthernPole},
                                                        {"longitudeOfSouthernPole", longitudeOfSouthernPole},
                                                        {"latitudeOfNorthernPole", latitudeOfNorthernPole},
                                                        {"longitudeOfNorthernPole", longitudeOfNorthernPole}};
            return gridParams;
        }

        std::map<std::string, double> getGriddedValues(CodesHandlePtr &h, std::vector<double> &values) {
            std::map<std::string, double> gridParams = getGridParameters(h);
            getGriddedValues(h, values, gridParams);
            return gridParams;
        }

        void getGriddedValues(CodesHandlePtr &h, std::vector<double> &values, std::map<std::string, double> &gridParams) {
            if (gridParams.empty()) {
                gridParams = getGridParameters(h);
            }

            size_t values_len;
            CODES_CHECK(codes_get_size(h.get(), "values", &values_len), 0);
            double Ni = gridParams.at("Ni"), Nj = gridParams.at("Nj");
            if (values_len != (unsigned)(Ni * Nj)) {
                std::ostringstream ss;
                ss << "Declaring grid of size " << Ni << "x" << Nj << "=" << Ni * Nj << " ";
                ss << "but containing " << values_len << " values. This is inconsistent!";
                throw InvalidArgumentException(ss.str(), AT);
            }

            values.resize(values_len);
            GRIB_CHECK(codes_get_double_array(h.get(), "values", values.data(), &values_len), 0);
        }

        void getNearestValues_grib(CodesHandlePtr &h, const std::vector<double> &in_lats, const std::vector<double> &in_lons, std::vector<double> &out_lats, std::vector<double> &out_lons,
                                   std::vector<double> &distances, std::vector<double> &values, std::vector<int> &indexes) {
            size_t npoints = in_lats.size();
            CODES_CHECK(codes_grib_nearest_find_multiple(h.get(), 0, in_lats.data(), in_lons.data(), static_cast<long>(npoints), out_lats.data(), out_lons.data(), values.data(), distances.data(),
                                                         indexes.data()),
                        0);
        }

        // ************************
        // BUFR
        // ************************

        void unpackMessage(CodesHandlePtr &m) {
            /* We need to instruct ecCodes to expand the descriptors
            i.e. unpack the data values */
            CODES_CHECK(codes_set_long(m.get(), "unpack", 1), 0);
        }

        /**
         * Returns either an empty string or a prefix to index the subset in BUFR messages.
         * As /subsetNumber=id/ where a key can follow. If no subset is present, an empty string is returned.
         *
         * @param subsetNumber The subset number.
         * @return The subset prefix as a string.
         */
        std::string getSubsetPrefix(const size_t &subsetNumber) {
            if (subsetNumber > 0) {
                return "/subsetNumber=" + std::to_string(subsetNumber) + "/";
            }
            return "";
        }

        // assumes UTC 0, need to convert to local time after
        Date getMessageDateBUFR(CodesHandlePtr &h, const size_t &subsetNumber, const double &tz_in) {
            std::string subset_prefix = getSubsetPrefix(subsetNumber);
            std::vector<std::string> parameters = {"year", "month", "day", "hour", "minute"}; // second is optional
            std::vector<int> values(6, -1);                                                   // year, month, day, hour, minute, second

            for (size_t i = 0; i < parameters.size(); ++i) {
                getParameter(h, subset_prefix + parameters[i], values[i]);
            }

            int second;
            bool success = getParameter(h, subset_prefix + "second", second, IOUtils::nothrow);
            if (!success) {
                second = -1; // second is optional, and if not present is set to weird number
            }
            values[5] = second;

            if (values[0] == -1 || values[1] == -1 || values[2] == -1 || values[3] == -1 || values[4] == -1)
                return Date();
            if (values[5] == -1)
                return Date(values[0], values[1], values[2], values[3], values[4], tz_in);

            Date base(values[0], values[1], values[2], values[3], values[4], values[5], tz_in);
            return base;
        };

        // ------------------- SETTERS -------------------
        void setMissingValue(CodesHandlePtr &message, double missingValue) { CODES_CHECK(codes_set_double(message.get(), "missingValue", missingValue), 0); }

        // multiple overloads for different parameter types
        bool selectParameter(codes_index *raw, const std::string &param_key, const std::string &paramId) { return codes_index_select_string(raw, param_key.c_str(), paramId.c_str()) == 0; };
        bool selectParameter(codes_index *raw, const std::string &param_key, const double &paramId) { return codes_index_select_double(raw, param_key.c_str(), paramId) == 0; };
        bool selectParameter(codes_index *raw, const std::string &param_key, const long &paramId) { return codes_index_select_long(raw, param_key.c_str(), paramId) == 0; };

        // ------------------- WRITE -------------------
        void writeToFile(CodesHandlePtr &h, const std::string &filename) {
            if (!FileUtils::fileExists(filename)) {
                ofilestream ofs(filename);
                ofs.close();
            }
            codes_write_message(h.get(), filename.c_str(), "a");
        }

        static void setHeader(codes_handle *ibufr, long num_subsets) {
            CODES_CHECK(codes_set_long(ibufr, "edition", 4), 0);
            CODES_CHECK(codes_set_long(ibufr, "masterTableNumber", 0), 0);
            CODES_CHECK(codes_set_long(ibufr, "masterTablesVersionNumber", WMO_BUFR_TABLE_NO), 0);
            CODES_CHECK(codes_set_long(ibufr, "dataCategory", 0), 0);
            CODES_CHECK(codes_set_long(ibufr, "internationalDataSubCategory", 2), 0);
            CODES_CHECK(codes_set_long(ibufr, "dataSubCategory", 2), 0);
            CODES_CHECK(codes_set_long(ibufr, "numberOfSubsets", num_subsets), 0);
        }

        // ------------------- HELPERS -------------------
        // adds the descriptors to the list of descriptors
        static void addDescriptors(std::vector<long> &descriptors, const std::vector<long> &to_add) { descriptors.insert(descriptors.end(), to_add.begin(), to_add.end()); }

        // adds the number of repitions to the repeated descriptors (not using delayed descriptors, as we know the number of repitions already and it shouldnt be too much)
        static void addRepeatedDescriptors(std::vector<long> &descriptors, std::vector<long> &replication_factors_in_subset, const std::vector<long> &repeated_descriptors, long num) {
            addDescriptors(descriptors, repeated_descriptors);
            replication_factors_in_subset.push_back(num);
        }

        // adds the standard descriptor, if it found in the available params and not repeated
        static void addStandardDescriptors(std::vector<long> &descriptors, long num, const std::set<std::string> &available_params, const std::string &param, long descriptor) {
            if (num == 0 || (num > 0 && available_params.find(param) != available_params.end())) {
                descriptors.push_back(descriptor);
            }
        }

        // ------------------- BUFR MESSAGE CREATION -------------------
        static void setCryosDescriptors(std::vector<long> &descriptors, std::vector<long> &replication_factors_in_subset, long num_heights) {
            // even though it seems similar we need to have to functions setting the descriptors, to later allow for the cryos sequence
            std::vector<long> info_descriptors(getDescriptors({"WIGOS_ID", "LongStationName", "StationType", "Date", "Time", "LatLon", "Alt"}));
            std::vector<long> single_occurence_descriptors(getDescriptors({"SurfType", "IceThickness", "GroundState", "HS", "SnoDepthMethod", "SurfQualifier", "TSS"}));
            std::vector<long> repeated_descriptors = {109000,
                                                      31001,
                                                      getDescriptor("HSensor"),
                                                      getDescriptor("SensorType"),
                                                      getDescriptor("TA"),
                                                      getDescriptor("TSNOW"),
                                                      getDescriptor("TIce"),
                                                      getDescriptor("TWater"),
                                                      getDescriptor("ChangeDataWidth"),
                                                      getDescriptor("ChangeScale"),
                                                      getDescriptor("TSOIL")};

            addDescriptors(descriptors, info_descriptors);
            addDescriptors(descriptors, single_occurence_descriptors);

            addRepeatedDescriptors(descriptors, replication_factors_in_subset, repeated_descriptors, num_heights);
            // TODO: When the descriptor is available, we can use this:
            // descriptors.push_back(307104);
            // replication_factors_in_subset.push_back(num_heights);
        }

        static void setMeteoIODesrciptors(std::vector<long> &descriptors, std::vector<long> &replication_factors_in_subset, const std::map<MeteoParam, size_t> &multi_param_occurences,
                                          const std::set<std::string> &available_params, const std::vector<MeteoParam> &POSSIBLE_MULTIPLE_PARAMETERS) {
            // these are all the descriptors needed to include the meteoio parameters
            std::vector<long> info_descriptors(getDescriptors({"Date", "Time", "StationId", "StationNumber", "StationName", "LatLon", "Alt"}));
            std::vector<long> single_occurence_descriptors(getDescriptors({"TA2m", "TSG", "SurfQualifier", "TSS", "HS", "VW", "DW", "VW_MAX", "ILWR", "ISWR", "TAU_CLD", "PSUM", "PSUM_PH"}));
            std::vector<long> repeated_descriptors_p = {102000, 31001, getDescriptor("HSensor"), getDescriptor("P")};
            std::vector<long> repeated_descriptors_ta = {102000, 31001, getDescriptor("HSensor"), getDescriptor("TA")};
            std::vector<long> repeated_descriptors_rh = {102000, 31001, getDescriptor("HSensor"), getDescriptor("RH")};
            std::vector<long> repeated_descriptors_tsoil = {102000, 31001, getDescriptor("HSensor"), getDescriptor("TSOIL")};
            std::vector<long> repeated_descriptors_tsnow = {102000, 31001, getDescriptor("HSensor"), getDescriptor("TSNOW")};

            // this is to avoid wrong ordering of the descriptors
            // the descriptors need to be hardcoded, so to avoid any confusion with setting new parameters or similar i do this wrapper
            std::map<MeteoParam, std::vector<long> *> repeated_descriptors_mapping = {{MeteoParam::P, &repeated_descriptors_p},
                                                                                      {MeteoParam::TA, &repeated_descriptors_ta},
                                                                                      {MeteoParam::RH, &repeated_descriptors_rh},
                                                                                      {MeteoParam::TSOIL, &repeated_descriptors_tsoil},
                                                                                      {MeteoParam::TSNOW, &repeated_descriptors_tsnow}};

            addDescriptors(descriptors, info_descriptors);
            addDescriptors(descriptors, single_occurence_descriptors);

            // need to add those before the replications, cause otherwise it looks horrible
            addStandardDescriptors(descriptors, multi_param_occurences.at(MeteoParam::P), available_params, "P", 7004);
            addStandardDescriptors(descriptors, multi_param_occurences.at(MeteoParam::RH), available_params, "RH", 13003);

            // we use delayed descriptors, so that the number of repitions can easily be read
            for (const MeteoParam &param : POSSIBLE_MULTIPLE_PARAMETERS) {
                if (repeated_descriptors_mapping.find(param) == repeated_descriptors_mapping.end()) {
                    throw InvalidArgumentException("The parameter " + std::to_string(static_cast<int>(param)) + " does not have a descriptor mapping", AT);
                }
                addRepeatedDescriptors(descriptors, replication_factors_in_subset, *repeated_descriptors_mapping.at(param), multi_param_occurences.at(param));
            }
        }

        // repeats the replication factors for all subsets
        static void setReplicationFactors(codes_handle *ibufr, const std::vector<long> &replication_factors_in_subset, long num_subsets) {
            if (!replication_factors_in_subset.empty()) {
                std::vector<long> all_replication_factors;
                for (int i = 0; i < num_subsets; i++) {
                    all_replication_factors.insert(all_replication_factors.end(), replication_factors_in_subset.begin(), replication_factors_in_subset.end());
                }
                codes_set_long_array(ibufr, "inputDelayedDescriptorReplicationFactor", all_replication_factors.data(), all_replication_factors.size());
            }
        }

        CodesHandlePtr createBUFRMessageFromSample(long num_subsets, const std::map<MeteoParam, size_t> &multi_param_occurences, const std::set<std::string> &available_params,
                                                   const std::vector<MeteoParam> &POSSIBLE_MULTIPLE_PARAMETERS, const bool &write_cryos_station, const long &num_cryo_heights) {

            if (write_cryos_station && num_cryo_heights == 0) {
                throw InvalidArgumentException("The number of cryo heights needs to be greater than 0, to write a Cryo Station File", AT);
            }

            std::cout << "Creating BUFR message with " << num_subsets << " subsets" << std::endl;
            codes_handle *ibufr = codes_handle_new_from_samples(NULL, "BUFR4");
            if (!ibufr) {
                throw IOException("Unable to create handle from sample", AT);
            }
            setHeader(ibufr, num_subsets);

            std::vector<long> descriptors;
            std::vector<long> replication_factors_in_subset;


            if (write_cryos_station)
                setCryosDescriptors(descriptors, replication_factors_in_subset, num_cryo_heights);
            else
                setMeteoIODesrciptors(descriptors, replication_factors_in_subset, multi_param_occurences, available_params, POSSIBLE_MULTIPLE_PARAMETERS);

            // the snow profile replication factors need to be accounted here.

            setReplicationFactors(ibufr, replication_factors_in_subset, num_subsets);

            CODES_CHECK(codes_set_long_array(ibufr, "unexpandedDescriptors", descriptors.data(), descriptors.size()), 0);
            return makeUnique(ibufr);
        }

        void setTime(CodesHandlePtr &ibufr, const Date &date, const std::string &subset_prefix) {
            int year, month, day, hour, minute, second;
            date.getDate(year, month, day, hour, minute, second);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "year").c_str(), year), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "month").c_str(), month), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "day").c_str(), day), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "hour").c_str(), hour), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "minute").c_str(), minute), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), (subset_prefix + "second").c_str(), second), 0);
        }

        void setTypicalTime(CodesHandlePtr &ibufr, const Date &date) {
            int year, month, day, hour, minute, second;
            date.getDate(year, month, day, hour, minute, second);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalYear", year), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalMonth", month), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalDay", day), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalHour", hour), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalMinute", minute), 0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalSecond", second), 0);
        }

        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const double &parameterValue) {
            if (parameterValue == IOUtils::nodata) return true;
            int err = codes_set_double(ibufr.get(), parameterName.c_str(), parameterValue);
            return err == 0;
        }
        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const long &parameterValue) {
            if (parameterValue == static_cast<long>(IOUtils::nodata)) return true;
            int err = codes_set_long(ibufr.get(), parameterName.c_str(), parameterValue);
            return err == 0;
        }

        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const std::vector<long> &parameterValues) {
            std::vector<long> values = parameterValues;
            double missingValue;
            CODES_CHECK(codes_get_double(ibufr.get(), "missingValue", &missingValue), 0);
            std::replace(values.begin(), values.end(), static_cast<long>(IOUtils::nodata), static_cast<long>(missingValue));
            int err = codes_set_long_array(ibufr.get(), parameterName.c_str(), values.data(), values.size());
            return err == 0;
        }

        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const std::string &parameterValue) {
            size_t len = parameterValue.size();
            int err = codes_set_string(ibufr.get(), parameterName.c_str(), parameterValue.c_str(), &len);
            return err == 0;
        }

        void packMessage(CodesHandlePtr &m) {
            /* We need to instruct ecCodes to pack the descriptors
            i.e. pack the data values */
            CODES_CHECK(codes_set_long(m.get(), "pack", 1), 0);
        }

    } // namespace codes

} // namespace mio
