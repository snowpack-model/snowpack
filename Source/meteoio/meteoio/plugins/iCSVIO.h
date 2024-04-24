// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2024 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef iCSV_H
#define iCSV_H

#include <meteoio/FileUtils.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/plugins/iCSVHelper.h>
#include <meteoio/plugins/libacdd.h>

#include <vector>

namespace mio {

    using namespace iCSV;
    /**
     * @class iCSVIO
     * @brief A class to read and write iCSV files
     *
     * @section TODO (when available)
     * - if meteoparam metadata is available use it to pass along long_name ... and also write it out
     *
     *
     * @ingroup plugins
     * @author Patrick Leibersperger
     * @date   2024-02-015
     */
    class iCSVIO : public IOInterface {
    public:
        iCSVIO(const std::string &configfile);
        iCSVIO(const iCSVIO &);
        iCSVIO(const Config &cfgreader);

        virtual void readStationData(const Date &date, std::vector<StationData> &vecStation);
        virtual void readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecMeteo);

        virtual void writeMeteoData(const std::vector<std::vector<MeteoData>> &vecMeteo, const std::string &name = "");

    private:
        // input section
        const Config cfg;
        std::string coordin, coordinparam, coordout, coordoutparam; // projection parameters
        bool snowpack_slopes;
        bool read_sequential;

        // file information
        std::vector<iCSVFile> stations_files;

        // output section
        ACDD acdd_metadata;
        double TZ_out;
        std::string outpath;
        bool allow_overwrite;
        bool allow_append;
        char out_delimiter;
        std::string file_extension_out;


        // constants
        const std::string iCSV_version = "1.0";
        const std::string iCSV_firstline = "# iCSV " + iCSV_version + " UTF-8";
        static const double snVirtualSlopeAngle;

        // constructor helpers
        void parseInputSection();
        void parseOutputSection();

        // read helpers
        void identify_fields(const std::vector<std::string> &fields, std::vector<size_t> &indexes, MeteoData &md,
                                     const std::string &geometry_field);
        double getSnowpackSlope(const std::string &id);
        void read_meta_data(const iCSVFile &current_file, StationData &meta);
        void setMetaDataPosition(const iCSVFile &current_file, StationData &meta, const double &nodata_value);
        void setMetaDataSlope(const iCSVFile &current_file, StationData &meta, const double &nodata_value);

        void readDataSequential(iCSVFile &current_file);
        std::vector<MeteoData> createMeteoDataVector(iCSVFile &current_file, std::vector<Date> &date_vec,
                                                             std::vector<geoLocation> &location_vec);
        void setMeteoDataLocation(MeteoData &tmp_md, geoLocation &loc, iCSVFile &current_file, double nodata);
        void setMeteoDataFields(MeteoData &tmp_md, iCSVFile &current_file, Date &date, std::vector<size_t> &indexes, double nodata);

        // write helpers
        void prepareOutfile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleNewFile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleFileAppend(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo);

        std::string getGeometry(const geoLocation &loc);
        bool checkLocationConsistency(const std::vector<MeteoData> &vecMeteo);
        bool createFilename(iCSVFile &outfile, const StationData &station, size_t station_num);
        void createMetaDataSection(iCSVFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void createFieldsSection(iCSVFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void writeToFile(const iCSVFile &outfile);
    };

} // namespace
#endif
