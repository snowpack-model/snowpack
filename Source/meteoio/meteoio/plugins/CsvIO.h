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

#include <string>
#include <vector>

namespace mio {

class CsvParameters {
	public:
		CsvParameters() : csv_fields(), units_offset(), units_multiplier(), location(), name(), id(), header_lines(1), columns_headers(1), date_col(0), time_col(0), csv_delim(','), eoln('\n'), datetime_idx(), file_and_path(), datetime_format(), csv_tz(0.) {}
		
		void setDateTimeSpec(const std::string& datetime_spec, const double& tz_in);
		void setFile(const std::string& i_file_and_path, const std::vector<std::string>& vecMetaSpec);
		Date parseDate(const std::string& datetime_str, const std::string& /*time_str*/) const;
		std::string getFilename() const {return file_and_path;}
		StationData getStation() const {return StationData(location, id, name);}
		
		//mio::FileUtils::FileIndexer indexer;
		std::vector<std::string> csv_fields;		///< the user provided list of field names
		std::vector<double> units_offset, units_multiplier;		///< offsets and multipliers to convert the data to SI
		
		Coords location;
		std::string name, id;
		size_t header_lines, columns_headers;
		size_t date_col, time_col;
		char csv_delim, eoln;
	private:
		static void parseFields(std::vector<std::string>& fieldNames, size_t &dt_col, size_t &tm_col);
		std::map< size_t, std::pair<size_t, std::string> > parseHeadersSpecs(const std::vector<std::string>& vecMetaSpec) const;
		
		std::vector<size_t> datetime_idx;		///< order of the datetime fields for use in parseDate
		std::string file_and_path, datetime_format; 		///< the scanf() format string for use in parseDate
		double csv_tz;		///< timezone to apply to parsed dates
		
};

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
		std::vector<std::string> readHeaders(std::ifstream& fin, CsvParameters& params) const;
		Date parseDate(const std::string& date_str, const std::string& time_str) const;
		std::vector<MeteoData> readCSVFile(CsvParameters& params, const Date& dateStart, const Date& dateEnd);
		
		const Config cfg;
		std::vector<CsvParameters> csvparam;
		std::vector<StationData> vecStations;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		static const size_t streampos_every_n_lines; //save current stream pos every n lines of data
};

} //namespace
#endif
