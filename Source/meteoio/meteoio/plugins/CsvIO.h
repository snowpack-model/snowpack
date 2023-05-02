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
#ifndef CSVIO_H
#define CSVIO_H

#include <meteoio/IOInterface.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/plugins/CsvParams.h>

#include <string>
#include <vector>

namespace mio {
/**
 * @class CsvIO
 * @brief Reads meteo data from a comma separated file.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2018-01
 */
class CsvIO : public IOInterface {
	public:
		CsvIO(const std::string& configfile);
		CsvIO(const CsvIO&);
		CsvIO(const Config& cfgreader);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void parseInputOutputSection();
		void cleanup() noexcept;
		std::string setDateParsing(const std::string& datetime_spec);
		std::vector<std::string> readHeaders(std::ifstream& fin, CsvParameters& params) const;
		static MeteoData createTemplate(const CsvParameters& params);
		static Date getDate(CsvParameters& params, const std::vector<std::string>& vecFields, const bool& silent_errors, const std::string& filename, const size_t& linenr);
		std::vector<MeteoData> readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd);
		
		const Config cfg;
		std::map<std::string, FileUtils::FileIndexer> indexer_map;
		std::vector<CsvParameters> csvparam;
		std::vector<StationData> vecStations;
		std::string coordin, coordinparam; //projection parameters
		static const size_t streampos_every_n_lines; //save current stream pos every n lines of data
		bool silent_errors; ///< when reading a file, should errors throw or just be ignored?
		bool errors_to_nodata;    //unparseable values are treated as nodata, but the dataset is kept
};

} //namespace
#endif
