// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 SLF                                                                                                                                */
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
#ifndef METEOBLUE_H
#define METEOBLUE_H

#include <meteoio/IOInterface.h>
#include <meteoio/thirdParty/picojson.h>

#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

/**
 * @class MeteoBlue
 * @brief This class enables the access to the MeteoBlue RESTful web service
 *
 * @ingroup plugins
 * @date   2021-11-24
 */

class MeteoBlue : public IOInterface {
	public:
		MeteoBlue(const std::string& configfile);
		MeteoBlue(const MeteoBlue&);
		MeteoBlue(const Config&);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation) override;
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo) override;

		typedef struct METEOPARAM {
			METEOPARAM() : param(), units_multiplier(1.), units_offset(0.) {}
			METEOPARAM(const MeteoGrids::Parameters& i_param)
			                : param(i_param), units_multiplier(1.), units_offset(0.) {}
			METEOPARAM(const MeteoGrids::Parameters& i_param, const double& multiplier, const double& offset)
			                : param(i_param), units_multiplier(multiplier), units_offset(offset) {}

			std::string getParameterName() const {return MeteoGrids::getParameterName( param );}
			double convertValue(const double& val) const {return (units_multiplier*val + units_offset);}
			
			std::string toString() const {
				std::ostringstream os;
				os << "[ " << MeteoGrids::getParameterName( param ) << " *" << units_multiplier << " +" << units_offset << "]";
				return os.str();
			}

			MeteoGrids::Parameters param;
			double units_multiplier, units_offset; ///< first we apply the multiplier, THEN the offset
		} meteoParam;
		
	private:
		void init();
		void readData(const StationData& sd, std::vector<MeteoData> &vecMeteo);
		void readTime(const picojson::value &v, const StationData& sd, std::vector<MeteoData> &vecMeteo) const;
		void readParameter(const picojson::value &v, const std::string& paramname, std::vector<MeteoData> &vecMeteo) const;
		static picojson::value goToJSONPath(const std::string& path, const picojson::value& v);
		static size_t data_write(void* buf, const size_t size, const size_t nmemb, void* userp);
		bool curl_read(const std::string& url, std::ostream& os) const;

		const Config cfg;
		std::vector<StationData> vecMeta;
		std::string coordin, coordinparam; ///< projection parameters
		std::string endpoint, apikey, packages;
		int http_timeout; //time out for http connections
		bool debug;

		static std::map< std::string, meteoParam > params_map; ///< parameters to extract from the files
		static const std::string dflt_endpoint;
		static const int http_timeout_dflt;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

} //end namespace mio

#endif
