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

/**
 * @note For reading profiles, it will be possible to read all iterations in a subset (only way to allow subsets)
 **/

#include <meteoio/plugins/BUFRFile.h>

namespace mio {

    using namespace codes;
    // initialize start_date very low, and end_date very high
    BUFRFile::BUFRFile(const std::string &in_filename, const std::string &ref_coords, bool verbose_in, double tz_in)
        : filename(in_filename), meta_data(), station_timezones(), start_date(), end_date(), default_timezone(tz_in), isMeteoTimeseries(false), isProfile(false), subset_numbers(), station_ids_in_file(), verbose(verbose_in)  {

        readMetaData(ref_coords);
        if (default_timezone == IOUtils::nodata) {
            throw IOException("Default timezone cannot be nodata", AT);
        }
    };


    // ------------------ STATIC HELPERS ------------------
    //--------- read the metadata
    static void getStationIdfromMessage(CodesHandlePtr &h, std::string &station_id, std::string &station_name, const size_t& subsetNumber) {
        std::vector<std::string> id_keys = {"stationNumber", "wigosIdentifierSeries", "shortStationName"};
        std::vector<std::string> name_keys = {"stationOrSiteName", "longStationName"};
        std::string stationID, stationName;
        getParameter(h, id_keys, stationID, subsetNumber);
        getParameter(h, name_keys, stationName, subsetNumber);

        if (stationID.empty()) {
            if (stationName.empty()) {
                throw AccessException("No station ID or name found in BUFR file, need to identify stations", AT);
            } else {
                stationID = stationName;
            }
        }
        station_id = stationID;
        station_name = stationName;
        return;
    };

    static void getPosition(CodesHandlePtr &h, double &latitude, double &longitude, double &altitude, const size_t &subsetNumber) {
        std::string subset_prefix = getSubsetPrefix(subsetNumber);
        getParameter(h, subset_prefix + "latitude", latitude);
        getParameter(h, subset_prefix + "longitude", longitude);

        std::vector<std::string> height_keys = {"heightOfStation", "height", "elevation", "heightOfStationGroundAboveMeanSeaLevel"};
        bool success = getParameter(h, height_keys, altitude, subsetNumber);
        if (!success) {
            std::cerr << "No altitude for station found in BUFR file" << std::endl;
            altitude = IOUtils::nodata;
        }
    }

    static void handleReferenceCoordinates(CodesHandlePtr &h, Coords &position, const double& latitude, const double& longitude, const double& altitude,const std::string &subset_prefix, const std::string &ref_coords, std::string &info_msg) {
        long ref_flag;
        bool found_ref = getParameter(h, subset_prefix + "coordinateReferenceSystem", ref_flag, IOUtils::nothrow);
        if (!found_ref) { // weird behaviour where value is changed to random values if not found
            ref_flag = 999;
        }

        // multiple possible ways to get the reference coordinates
        if (ref_flag == 999 || ref_flag == 65535) {
            if (ref_coords.empty()) {
                std::string error = (ref_flag == 999) ? "No reference coordinates found in BUFR file and none provided in the configuration"
                                                      : "Missing reference coordinates in BUFR file, and none provided in configuration";
                throw InvalidFormatException(error, AT);
            } else {
                Coords ref_coords_obj(ref_coords);
                ref_coords_obj.setLatLon(latitude, longitude, altitude);
                position = ref_coords_obj;
                info_msg = (ref_flag == 999) ? "Using reference coordinates from configuration. Could not find any in BUFR file."
                                            : "Using reference coordinates from configuration. Missing coordinates reference in file";
            }
        } else if (ref_flag > 3) {
            std::ostringstream ss;
            ss << "Unsuppported reference flag " << ref_flag << " in BUFR file";
            throw InvalidFormatException(ss.str(), AT);
        } else if (ref_flag <= 3) {
            info_msg = "Using coordinates reference from BUFR file";
            position.setEPSG(FLAG_TO_EPSG[ref_flag]);
            position.setPoint(latitude, longitude, altitude);
        } 
    }

    static void getTimezone(CodesHandlePtr &h, double &timezone, const std::string &subset_prefix) {
        // get the timezone if available otherwise default to UTC
        bool hasTz = getParameter(h, subset_prefix + "timeDifferenceUtcLmt", timezone, IOUtils::nothrow);
        if (!hasTz) {
            timezone = IOUtils::nodata; // default to UTC
        } else {
            timezone /= 60.0; // convert to hours
        }
    }

    // !!!Read all the necessary metadata from a BUFR file, returns timezone if possible
    static StationData getStationDataBUFR(CodesHandlePtr &h, double& timezone,const std::string &ref_coords, const size_t &subsetNumber, std::string &info_msg) {
        std::string error;
        std::string subset_prefix = getSubsetPrefix(subsetNumber);
        StationData stationData;

        // get the station ID and station name
        std::string stationID, stationName;
        getStationIdfromMessage(h, stationID, stationName, subsetNumber);


        // get the position 
        double latitude, longitude, altitude;
        getPosition(h, latitude, longitude, altitude, subsetNumber);


        Coords position;

        handleReferenceCoordinates(h, position, latitude, longitude, altitude, subset_prefix, ref_coords, info_msg);

        stationData.setStationData(position, stationID, stationName);

        getTimezone(h, timezone, subset_prefix);
        return stationData;
    };


    // -----------Read the date from a BUFR message
    static void fillAdditionalParams(MeteoData &md, CodesHandlePtr &message, const std::vector<std::string> &additional_params, const std::string &subset_prefix) {
        for (const auto &params : additional_params) {
            double value = IOUtils::nodata;
            getParameter(message, subset_prefix + params, value); // should throw an error, as user defined parameters should be in the file
            
            size_t param_id = md.addParameter(params);
            md(param_id) = value;
        }
    };

    static void handleSpecialCases(CodesHandlePtr &/* message */, MeteoData &md, const std::string &/* subset_prefix */, const std::string &param_name, bool& /* success */, const size_t &id) {
        // handle special cases, where the parameter name is not the same as in the MeteoData object
        // e.g. airTemperatureAt2M is not defined, but maybe airTemperature is and conver TAU_CLD to decimal
        if (param_name == "TAU_CLD" && md(id) != IOUtils::nodata)
            md(id) = md(id) / 100.0;
    };

    // !!!Read the date from a BUFR message, by looping through all meteodata parameters
    static std::vector<std::string> fillFromMessage(MeteoData &md, CodesHandlePtr &message, const std::vector<std::string> &additional_params, const size_t &subsetNumber) {
        std::string subset_prefix = getSubsetPrefix(subsetNumber);
        // fill the MeteoData object with data from the message and additional_params
        std::vector<std::string> param_names_found;
        for (size_t id = 0; id < md.getNrOfParameters(); id++) {
            std::string param_name = md.getParameterName(id);

            bool success = getParameter(message, subset_prefix + BUFR_PARAMETER.at(param_name), md(id), IOUtils::nothrow);
            if (!success)
                success = getParameter(message, subset_prefix + BUFR_PARAMETER_ALT.at(param_name), md(id), IOUtils::nothrow);
            handleSpecialCases(message, md, subset_prefix, param_name, success, id);

            if (!success) 
                md(id) = IOUtils::nodata;
            else {
                param_names_found.push_back(param_name);
            }

        }
        
        fillAdditionalParams(md, message, additional_params, subset_prefix);

        return param_names_found;
    };

    static bool areMultipleSubsetsInMessage(CodesHandlePtr &h, double &num_subsets) {
        getParameter(h, "numberOfSubsets", num_subsets);
        return num_subsets > 0;
    };

    // ------------------ GET ALL STATION DATA ------------------
    void BUFRFile::readMetaData(const std::string &ref_coords) {
        // a bufr file can contain multiple messages, each message can contain multiple subsets
        // the order of stations and parameters is arbitrary.
        // So we go through all messages and its subsets to find all stations, and then check that there are no duplicate dates for any station 
        // Also we check that for each station we have the necessary metadata, and that the metadata is consistent
        std::vector<CodesHandlePtr> messages = getMessages(filename, PRODUCT_BUFR);
        const double NO_SUBSETS_FOUND = -1;
        size_t message_counter = 1;
        std::map<std::string, std::set<Date>> station_dates;
        for (auto &message : messages) {
            double num_subsets = getNumSubsets(message, NO_SUBSETS_FOUND);
            subset_numbers.push_back(static_cast<size_t>(num_subsets));

            unpackMessage(message);

            processSubsets(message, ref_coords, subset_numbers.back(), station_dates);
            message_counter++;
        }
        if (meta_data.empty()) {
            throw IOException("No valid stations found in file " + filename, AT);
        }

        isMeteoTimeseries = true;
    }

    // ------------------ STATION HELPERS ------------------
    double BUFRFile::getNumSubsets(CodesHandlePtr &message, double default_value) {
        bool hasMultipleSubsets = areMultipleSubsetsInMessage(message, default_value);
        return hasMultipleSubsets ? default_value : 0;
    }

    void BUFRFile::processSubsets(CodesHandlePtr &message, const std::string &ref_coords, const size_t& num_subsets, std::map<std::string, std::set<Date>> &station_dates) {
        std::string info_msg;
        bool warned_tz = false;
        for (size_t sb_id = 1; sb_id <= num_subsets; sb_id++) {
            std::string new_info;
            double tz_for_station;
            StationData new_meta = getStationDataBUFR(message, tz_for_station,ref_coords, sb_id, new_info);
            if (!info_msg.empty() && verbose && new_info != info_msg){
                std::cout << new_info << std::endl;
                info_msg = new_info;
            }

            if (tz_for_station == IOUtils::nodata) {
                if (!warned_tz) {
                    std::cerr << "No timezone found for station in file " << filename << ". Using User provided Timezone: " << default_timezone <<" (default: GMT+0)" << std::endl;
                    warned_tz = true;
                }
                tz_for_station = default_timezone;
            }
            
            Date new_date = getMessageDateBUFR(message, sb_id, tz_for_station);
            std::string station_id = new_meta.getStationID();

            if (isNewStation(station_id)) {
                processNewStation(station_id, new_meta, new_date, station_dates);
            } else {
                processExistingStation(station_id, new_meta, new_date, station_dates);
            }
            updateDateRange(new_date);
        }
    }

    void BUFRFile::processNewStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates) {
        if (!new_meta.isValid())
            throw IOException("No valid station metadata for station " + station_id + " file " + filename, AT);
        station_ids_in_file.push_back(station_id);
        meta_data[station_id] = new_meta;

        // create a set to be sure we done have duplicate timepoints
        std::set<Date> dates;
        dates.insert(new_date);
        station_dates[station_id] = dates;
        station_timezones[station_id] = new_date.getTimeZone(); 
    }

    void BUFRFile::processExistingStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates) {
        // check if the date is already in the set
        auto success = station_dates[station_id].insert(new_date);
        if (!success.second) {
            throw IOException("Duplicate date for station " + station_id + " in file " + filename);
        }

        // check consistency between station data
        if (!new_meta.isEmpty() && meta_data[station_id] != new_meta) {
            throw IOException("Inconsistent station data for station " + station_id + " in file " + filename);
        }
    }

    void BUFRFile::updateDateRange(const Date &new_date) {
        start_date.isUndef() ? start_date = new_date : start_date = std::min(start_date, new_date);
        end_date.isUndef() ? end_date = new_date : end_date = std::max(end_date, new_date);
    }

    // ------------------ READ DATA ------------------
    void BUFRFile::readData(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping, const std::vector<std::string> &additional_params) {
        // a bufr file can contain multiple messages, each message can contain multiple subsets
        // the order of where stations are stored at what date is quite arbitrary, so we need to map the station ids to the correct index in vecMeteo
        // as well as keeping track between calls, what station ids we already have
        auto messages = getMessages(filename, PRODUCT_BUFR);
        updateStationIdsMapping(vecMeteo, station_ids_mapping);

        for (size_t msg_id = 0; msg_id < messages.size(); msg_id++) {
            unpackMessage(messages[msg_id]);
            processMessage(vecMeteo, station_ids_mapping, additional_params, messages[msg_id], msg_id);
        }
    }

    // ------------------ HELPERS ------------------
    void BUFRFile::updateStationIdsMapping(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping) {
        for (const auto& station : station_ids_in_file) {
            if (station_ids_mapping.find(station) == station_ids_mapping.end()) {
                METEO_SET empty_set;
                station_ids_mapping[station] = vecMeteo.size();
                vecMeteo.push_back(empty_set);
            }
        }
    }

    void BUFRFile::processMessage(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping, const std::vector<std::string> &additional_params, CodesHandlePtr &message, size_t msg_id) {
        // TODO: if we need speedup cache the date and meta data (station id) for each message and subset
        std::set<std::string> param_names_in_file;
        for (size_t sub_id = 1; sub_id <= subset_numbers[msg_id]; sub_id++) {
            std::string station_id, stationName;
            getStationIdfromMessage(message, station_id, stationName, sub_id);
            size_t index_in_vecMeteo = station_ids_mapping[station_id];
            
            Date date = getMessageDateBUFR(message, sub_id,station_timezones[station_id]);
            MeteoData md;
            md.setDate(date);
            md.meta = meta_data[station_id];
            std::vector<std::string> found_params = fillFromMessage(md, message, additional_params, sub_id);
            if (verbose) {
                param_names_in_file.insert(found_params.begin(), found_params.end());
            }
            vecMeteo[index_in_vecMeteo].push_back(md);
        }
        if (verbose) {
            std::cout << "Found parameters in file " << filename << ": ";
            for (const auto &param : param_names_in_file) {
                std::cout << param << " ";
            }
            std::cout << std::endl;
        }
    }
} // namespace mio