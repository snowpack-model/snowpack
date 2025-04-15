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

		virtual void readStationData(const Date &date, std::vector<StationData> &vecStation) override;
		virtual void readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecMeteo) override;

		virtual void writeMeteoData(const std::vector<std::vector<MeteoData>> &vecMeteo, const std::string &name = "") override;

	private:
		// constructor helpers
		void parseInputSection();
		void parseOutputSection();

		// read helpers
		static void identify_fields(const std::vector<std::string> &fields, std::vector<size_t> &indexes, MeteoData &md, const std::string &geometry_field);
		static double getSnowpackSlope(const std::string &id);
		void read_meta_data(const iCSVFile &current_file, StationData &meta) const;
		void setMetaDataPosition(const iCSVFile &current_file, StationData &meta, const double &nodata_value) const;
		void setMetaDataSlope(const iCSVFile &current_file, StationData &meta, const double &nodata_value) const;

		void readDataSequential(iCSVFile &current_file) const;
		std::vector<MeteoData> createMeteoDataVector(iCSVFile &current_file, std::vector<Date> &date_vec, std::vector<geoLocation> &location_vec) const;
		static void setMeteoDataLocation(MeteoData &tmp_md, geoLocation &loc, iCSVFile &current_file, double nodata);
		static void setMeteoDataFields(MeteoData &tmp_md, iCSVFile &current_file, const Date &date, std::vector<size_t> &indexes, double nodata);

		// write helpers
		void prepareOutfile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists) const;
		void handleNewFile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists) const;
		static void handleFileAppend(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo);

		static std::string getGeometry(const geoLocation &loc);
		static bool checkLocationConsistency(const std::vector<MeteoData> &vecMeteo);
		bool createFilename(iCSVFile &outfile, const std::vector<MeteoData> &vecvecMeteo, size_t station_num) const;
		void createMetaDataSection(iCSVFile &current_file, const std::vector<MeteoData> &vecMeteo) const;
		static void createFieldsSection(iCSVFile &current_file, const std::vector<MeteoData> &vecMeteo);
		void writeToFile(const iCSVFile &outfile) const;


		// constants
		const std::string iCSV_version = "1.0";
		const std::string iCSV_firstline = "# iCSV " + iCSV_version + " UTF-8";

		std::vector<iCSVFile> stations_files;
		const Config cfg;
		ACDD acdd_metadata;
		std::string coordin, coordinparam, coordout, coordoutparam; // projection parameters
		std::string file_extension_out, outpath, versioning_str;
		double TZ_out;
		VersioningType outputVersioning; //this is usefull when generating multiple versions of the same dataset, for example with forecast data
		char out_delimiter;
		bool snowpack_slopes, read_sequential, allow_overwrite, allow_append;
};

} // namespace
#endif
