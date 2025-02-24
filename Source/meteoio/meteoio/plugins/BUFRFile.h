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

#ifndef BUFRFILE_H
#define BUFRFILE_H

#include <meteoio/IOInterface.h>
#include <meteoio/plugins/libcodes.h>

#include <string>

namespace mio {

using namespace codes;

class BUFRFile {
    public:
        BUFRFile(const std::string &filename, const std::string& ref_coords, bool verbose = false, double default_timezone = IOUtils::nodata);
        
        void readMetaData(const std::string& ref_coords);
        // takes a map that contatins the bufr keys as keys, and names as values. Should map to itself, if no renaming is needed
        // returns a vector of meteo data timeseries, one for each station
        void readData(std::vector<METEO_SET> &vecMeteo, std::map<std::string,size_t>& station_ids, const std::vector<std::string> &additional_params);

        const std::map<std::string, StationData>& getMetadata() const { return meta_data; };

        bool hasMultipleSubets(const size_t &msg_index) const { return subset_numbers[msg_index] > 1; };

    private:
        std::string filename;
        std::map<std::string, StationData> meta_data;
        std::map<std::string, double> station_timezones;


        Date start_date;
        Date end_date;
        double default_timezone;

        bool isMeteoTimeseries;
        bool isProfile;
        // for each message see if we need subset indexing
        std::vector<size_t> subset_numbers;

        std::vector<std::string> station_ids_in_file;

        bool verbose;

        // HELPERS
        bool isNewStation(const std::string &station_id) { return station_ids_in_file.end() == std::find(station_ids_in_file.begin(), station_ids_in_file.end(), station_id); };
        double getNumSubsets(CodesHandlePtr &message, double default_value);
        void processSubsets(CodesHandlePtr &message, const std::string &ref_coords, const size_t& num_subsets, std::map<std::string, std::set<Date>> &station_dates);
        void processNewStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates);
        void processExistingStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates);
        void updateDateRange(const Date &new_date);
        void updateStationIdsMapping(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping);
        void processMessage(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping, const std::vector<std::string> &additional_params, CodesHandlePtr &message, size_t msg_id);



};

} //namespace


#endif // BUFRFILE_H
