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
     * This namespace provides a set of functions to handle GRIB and BUFR files using the ecCodes C library by the ECMWF (https://confluence.ecmwf.int/display/ECC/ecCodes+installation) at least in version 2.27.0.
     * Either a file is indexed (contained in CodesIndexPtr) and messages (contained as CodesHandlePtr) are read from this index by using key/value pairs.
     * Or all messages are read from a file directly, without any additional information.
     * The index does not use much memory, the handles do??. Therefore, the handles should be deleted as soon as they are not needed anymore.
     *
     *
     * @ingroup plugins
     * @author Patrick Leibersperge
     * @date   2024-04-17
     *
     */
    namespace codes {

        // ------------------- BUFR -------------------

        // mapping of BUFR parameters to MeteoIO parameters
        const std::map<std::string, std::string> BUFR_PARAMETER{
            {"P", "pressure"},
            {"TA", "airTemperatureAt2M"}, // Can also be found as airTemperature under heightOfSensor(somthing like this) = 2
            {"RH", "relativeHumidity"},
            {"TSG", "groundTemperature"},
            {"TSS", "snowTemperature"},
            {"HS", "totalSnowDepth"},
            {"VW", "windSpeed"},
            {"DW", "windDirection"},
            {"VW_MAX", "maximumWindGustSpeed"},
            {"RSWR", ""},                           // TODO: there is no parameter for it, should we save albedo? --> also need to find descriptor
            {"ISWR", "instantaneousShortWaveRadiation"},
            {"ILWR", "instantaneousLongWaveRadiation"},
            {"TAU_CLD", "cloudCoverTotal	"},
            {"PSUM", "totalPrecipitationOrTotalWaterEquivalent	"},
            {"PSUM_PH", "precipitationType"} // should the type be mapped to the phase?
        };

        // an alternative mapping of BUFR parameters to MeteoIO parameters
        const std::map<std::string, std::string> BUFR_PARAMETER_ALT { 
            {"TA", "airTemperature"}
        };

        // flags for the possible reference systems are 0-4
        const std::vector<int> FLAG_TO_EPSG = {4326, 4258, 4269, 4314};

        // ------------------- GRIB -------------------
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

        void unpackMessage(CodesHandlePtr &m) {
            /* We need to instruct ecCodes to expand the descriptors
            i.e. unpack the data values */
            CODES_CHECK(codes_set_long(m.get(), "unpack", 1), 0);
        }

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

        CodesHandlePtr createBUFRMessageFromSample() {
            codes_handle *ibufr = codes_handle_new_from_samples(NULL, "BUFR4");
            if (!ibufr) {
                throw IOException("Unable to create handle from sample", AT);
            }
            CODES_CHECK(codes_set_long(ibufr, "edition", 4),0);
            CODES_CHECK(codes_set_long(ibufr, "masterTableNumber", 0),0);
            CODES_CHECK(codes_set_long(ibufr, "masterTablesVersionNumber", 40),0);
            CODES_CHECK(codes_set_long(ibufr, "dataCategory", 0),0);
            CODES_CHECK(codes_set_long(ibufr, "internationalDataSubCategory", 2), 0);
            CODES_CHECK(codes_set_long(ibufr, "dataSubCategory", 2),0);
            // these are all the descriptors needed to include the meteoio parameters
            std::vector<long> descriptors = {301011, 301013, 1018,1002, 1015, 301021, 7030, 7004,12004,13003,12120,12131,13013,11012,11011,11041,14017,14018,20010,13011,20021};
            CODES_CHECK(codes_set_long_array(ibufr, "unexpandedDescriptors", descriptors.data(), descriptors.size()),0);
            return makeUnique(ibufr);
        }

        void setTime(CodesHandlePtr &ibufr, const Date &date) {
            int year, month, day, hour, minute, second;
            date.getDate(year, month, day, hour, minute, second);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalYear", year),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "year", year),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalMonth", month),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "month", month),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalDay", day),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "day", day),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalHour", hour),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "hour", hour),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalMinute", minute),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "minute", minute),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "typicalSecond", second),0);
            CODES_CHECK(codes_set_long(ibufr.get(), "second", second),0);
        }

        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const double &parameterValue) { 
            int err = codes_set_double(ibufr.get(), parameterName.c_str(), parameterValue); 
            return err == 0;
        }
        bool setParameter(CodesHandlePtr &ibufr, const std::string &parameterName, const long &parameterValue) {
            int err =  codes_set_long(ibufr.get(), parameterName.c_str(), parameterValue); 
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
