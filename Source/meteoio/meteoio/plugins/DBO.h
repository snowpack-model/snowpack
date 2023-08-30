// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2017 SLF                                                                                                                                */
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
#ifndef DBO_H
#define DBO_H

#include <meteoio/IOInterface.h>

#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

class JsonWrapper; //forward declaration, it is defined in an include called in the .cc file

/**
 * @class DBO
 * @brief This class enables the access to the DBO RESTful web service
 *
 * @ingroup plugins
 * @date   2017-01-26
 */

class DBO : public IOInterface {
	public:
		DBO(const std::string& configfile);
		DBO(const Config&);
		DBO(const DBO&);
		~DBO();

		DBO& operator=(const mio::DBO&);
		
		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);
		
	private:
		typedef struct ts_Meta {
			ts_Meta() : since(), until(), param_dbo(), parname(), agg_type(), id(0), units_factor(1.), units_offset(0), interval(IOUtils::unodata), ts_offset(IOUtils::unodata), sequence(IOUtils::unodata) {} //required by std::map
			ts_Meta(const std::string& i_param_str, const Date& i_since, const Date& i_until, const std::string& i_agg_type, const size_t i_id, const unsigned int& i_interval, const unsigned int& i_ts_offset, const unsigned int& i_sequence)
			                : since(i_since), until(i_until), param_dbo(i_param_str), parname(), agg_type(i_agg_type), id(i_id), units_factor(1.), units_offset(0), interval(i_interval), ts_offset(i_ts_offset), sequence(i_sequence) {}

			std::string toString() const {
				std::ostringstream os;
				os << "   " << id << " " << param_dbo << " [";
				os << ((since.isUndef())? "-∞" : since.toString(Date::ISO)) << " - ";
				os << ((until.isUndef())? "∞" : until.toString(Date::ISO)) << "] ";
				os << agg_type << " - " << interval << " s - #" << sequence;
				//os << " - x" << units_factor << " +" << units_offset;
				return os.str();
			}

			Date since, until;
			std::string param_dbo, parname, agg_type;
			size_t id;
			double units_factor, units_offset;
			unsigned int interval, ts_offset, sequence;
		} tsMeta;

		typedef struct ts_Data {
			ts_Data() : date(), val(IOUtils::nodata) {}
			ts_Data(const Date& i_date, const double& i_val) : date(i_date), val(i_val){}

			Date date;
			double val;
		} tsData;
		
		void fillStationMeta();
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		std::vector<DBO::tsMeta> getTsProperties() const;
		std::vector<DBO::tsData> getTimeSerie(const size_t& tsID, const double& factor, const double& offset) const;
		static void setUnitsConversion(DBO::tsMeta& ts);
		static std::string getParameter(const std::string& param_str, const std::string& agg_type);
		void mergeTimeSeries(const MeteoData& md_pattern, const size_t& param, const std::vector<DBO::tsData>& vecData, std::vector<MeteoData>& vecMeteo) const;
		void initDBOConnection(const Config& cfg);
		
		std::vector<std::string> vecStationName;
		std::vector<StationData> vecMeta;
		std::vector< std::vector<DBO::tsMeta> > vecTsMeta; ///< for every station, a map that contains for each timeseries the relevant timeseries properties
		std::string coordin, coordinparam; ///< projection parameters
		std::string endpoint; ///< Variables for endpoint configuration
		JsonWrapper *json;

		bool dbo_debug;

		static const int http_timeout_dflt;
		static const std::string endpoint_default, metadata_api, data_api;
};

} //end namespace mio

#endif
