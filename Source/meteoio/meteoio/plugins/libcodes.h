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
#ifndef LIBCODES_H
#define LIBCODES_H

#include <meteoio/IOInterface.h>

#include <eccodes.h>
#include <vector>

#include <iostream>

namespace mio {
    /**
     * @namespace codes
     * @brief This namespace handles all the low level manipulation of GRIB and BUFR files with ecCodes
     *
     * @ingroup plugins
     * @author Patrick Leibersperge
     * @date   2024-04-17
     * 
     */
    namespace codes {

        // ------------------------- POINTER HANDLING -------------------------
        struct HandleDeleter {
            void operator()(codes_handle *h) const { codes_handle_delete(h); }
        };

        struct IndexDeleter {
            void operator()(codes_index *i) const { codes_index_delete(i); }
        };

        using CodesHandlePtr = std::unique_ptr<codes_handle, HandleDeleter>;
        using CodesIndexPtr = std::unique_ptr<codes_index, IndexDeleter>;

        CodesHandlePtr makeUnique(codes_handle *h);
        CodesIndexPtr makeUnique(codes_index *i);

        // ------------------------- FILE AND MESSAGE HANDLING -------------------------
        CodesIndexPtr indexFile(const std::string &filename, const std::vector<std::string> &index_keys, bool verbose);

        // ------------------------- GETTERS -------------------------
        template <typename T> std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::string &param_key, const T &paramID, const std::string &level_key, const std::string &levelType);
        std::vector<CodesHandlePtr> getMessages(const std::string &filename, ProductKind product = PRODUCT_GRIB);
        std::vector<CodesHandlePtr> getMessages(FILE* in_file, ProductKind product = PRODUCT_GRIB);
        void unpackMessage(CodesHandlePtr &m);

        Date getMessageDateGrib(CodesHandlePtr &h, const double &tz_in);
        std::string getSubsetPrefix(const size_t &subsetNumber);
        Date getMessageDateBUFR(CodesHandlePtr &h, const size_t& subsetNumber, const double &tz_in = 0);

        std::map<std::string, double> getGridParameters(CodesHandlePtr &h_unique);
        void getGriddedValues(CodesHandlePtr &h, std::vector<double> &values, std::map<std::string, double> &gridParams);
        std::map<std::string, double> getGriddedValues(CodesHandlePtr &h, std::vector<double> &values);

        void getNearestValues_grib(CodesHandlePtr &h, const std::vector<double> &in_lats, const std::vector<double> &in_lons, std::vector<double> &out_lats, std::vector<double> &out_lons,
                                   std::vector<double> &distances, std::vector<double> &values, std::vector<int> &indexes);

        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, double &param_value, const IOUtils::ThrowOptions& throwError = IOUtils::dothrow);
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, long &param_value, const IOUtils::ThrowOptions& throwError = IOUtils::dothrow);
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, int &param_value, const IOUtils::ThrowOptions& throwError = IOUtils::dothrow);
        bool getParameter(CodesHandlePtr &h, const std::string &parameterName, std::string &param_value, const IOUtils::ThrowOptions& throwError = IOUtils::dothrow);
        template <typename T> bool getParameter(CodesHandlePtr &h, const std::vector<std::string> &paramNames, T &param_value, const size_t& subset_number = -1);

        // ------------------------- SETTERS -------------------------
        void setMissingValue(CodesHandlePtr &message, double missingValue);

        bool selectParameter(codes_index *raw, const std::string &param_key, const std::string &paramId);
        bool selectParameter(codes_index *raw, const std::string &param_key, const double &paramId);
        bool selectParameter(codes_index *raw, const std::string &param_key, const long &paramId);

        // ------------------------- WRITE -------------------------
        void writeToFile(CodesHandlePtr &h, const std::string &filename);
        CodesHandlePtr createBUFRMessageFromSample();
        void setTime(CodesHandlePtr& ibufr, const Date& date);
        bool setParameter(CodesHandlePtr& ibufr, const std::string& parameterName, const double& parameterValue);
        bool setParameter(CodesHandlePtr& ibufr, const std::string& parameterName, const long& parameterValue);
        bool setParameter(CodesHandlePtr& ibufr, const std::string& parameterName, const std::string& parameterValue);
        void packMessage(CodesHandlePtr &m);

        // ------------------------- CONSTANTS -------------------------
        extern const std::map<std::string, std::string> BUFR_PARAMETER;
        extern const std::map<std::string, std::string> BUFR_PARAMETER_ALT;
        extern const std::vector<int> FLAG_TO_EPSG;
        extern const std::map<std::string, long> GRIB_DEFAULT_PARAM_TABLE;
        extern const std::map<std::string, std::string> GRIB_DEFAULT_LEVELTYPE_TABLE;
        extern const std::map<std::string, long> GRIB_DEFAULT_LEVELNO_TABLE;

        // ------------------------- TEMPLATE FUNCTIONS -------------------------
        // definition of the template functions
        // -1 to indicate that the subset number is not used
        template <typename T> bool getParameter(CodesHandlePtr &h, const std::vector<std::string> &paramNames, T &param_value, const size_t &subset_number) {
            std::string subset_prefix = getSubsetPrefix(subset_number);
            T tmp = param_value;
            for (const auto &paramName : paramNames) {
                getParameter(h, subset_prefix + paramName, param_value, IOUtils::nothrow);
                if (param_value != tmp) {
                    return true;
                }
            }
            return false;
        }

        template <typename T>
        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::string &param_key, const T &paramID, const std::string &level_key, const std::string &levelType) {
            codes_index *raw = index.get();
            if (!selectParameter(raw, param_key, paramID)) {
                return {};
            };

            CODES_CHECK(codes_index_select_string(raw, level_key.c_str(), levelType.c_str()), 0);

            codes_handle *h = nullptr;
            int ret;
            std::vector<CodesHandlePtr> handles;
            while ((h = codes_handle_new_from_index(raw, &ret)) != nullptr) {
                if (!h)
                    throw IOException("Unable to create grib handle from index", AT);
                if (ret != 0 && ret != CODES_END_OF_INDEX) {
                    throw IOException("Error reading message: Errno " + std::to_string(ret), AT);
                }
#ifdef DEBUG
                size_t len = 500;
                char name[len] = {'\0'};
                std::cerr << "Found message ";
                CODES_CHECK(codes_get_string(h, "name", name, &len), 0);
                std::cerr << "With name " << name << "\n";
                CODES_CHECK(codes_get_string(h, "shortName", name, &len), 0);
                std::cerr << "With shortName " << name << "\n";
                if (levelType != 1) {
                    long level;
                    CODES_CHECK(codes_get_long(h, "level", &level), 0);
                    std::cerr << "With level " << level << "\n";
                }
#endif
                handles.push_back(makeUnique(h));
            }
            return handles;
        }
    };
} // namespace mio
#endif // LIBCODES_H