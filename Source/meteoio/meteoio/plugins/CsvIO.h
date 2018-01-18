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

#include <string>
#include <vector>

namespace mio {

/**
 * @class CsvIO
 * @brief This (empty) class is to be used as a template for developing new plugins
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2010-06-14
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
		void cleanup() throw();
		std::string setDateParsing(const std::string& datetime_spec);
		std::vector<StationData> initStations(const std::vector<std::string>& i_vecFilenames, const std::vector< std::pair<std::string, std::string> >& vecMeta, const std::vector<std::string>& vecMetaSpec) const;
		std::vector<std::string> readHeaders(std::ifstream& fin, const char& eoln, size_t &date_col, size_t &time_col) const;
		Date parseDate(const std::string& date_str, const std::string& time_str) const;
		std::vector<MeteoData> readCSVFile(const std::string& filename, const size_t& stat_idx, const Date& dateStart, const Date& dateEnd);
		
		const Config cfg;
		mio::FileUtils::FileIndexer indexer; //in order to save file pointers
		std::vector<StationData> vecStations;
		std::vector<std::string> vecFilenames, csv_fields;
		std::vector<double> units_offset, units_multiplier;
		std::vector<size_t> datetime_idx;
		
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::string meteopath, datetime_format;
		double csv_tz;
		static const size_t streampos_every_n_lines; //save current stream pos every n lines of data
		size_t header_lines, columns_headers;
		char csv_delim;
};

} //namespace
#endif
