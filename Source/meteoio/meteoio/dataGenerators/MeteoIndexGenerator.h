// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef METEOINDEXGENERATOR_H
#define METEOINDEXGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>
#include <meteoio/meteoLaws/Sun.h>

namespace mio {

/**
 * @class MeteoIndex
 * @ingroup parametrizations
 * @brief Mteorological indices data generator.
 * @details
 * This generator offers several common meteorological indices that aggregate several 
 * parameters into one convenient index. It is possible to choose which index to compute
 * with the TYPE argument:
 *    - WINDCHILL: human-percived feeling of air temperature on exposed skin due to wind;
 *    - HEATINDEX: human-perceived air temperature due to humidity;
 *    - WET_BULB: lowest temperature that could be reached by water evaporation.
 * 
 * @code
 * *::edit1 = CREATE
 * *::arg1::algorithm = METEOINDEX
 * *::arg1::type = WINDCHILL
 * *::arg1::param = CHILL
 * @endcode
 */
class MeteoIndex : public GeneratorAlgorithm {
	public:
		MeteoIndex(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ)
			: GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), sun(), model(WINDCHILL) { parse_args(vecArgs); }
		bool generate(const size_t& param, MeteoData& md);
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo);
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		static bool windChill(const size_t& param, MeteoData& md);
		static bool heatIndex(const size_t& param, MeteoData& md);
		static bool wetBulbTemperature(const size_t& param, MeteoData& md);
		bool WBGT_index(const size_t& param, MeteoData& md);
		
		typedef enum PARAMETRIZATION {
			WINDCHILL,
			HEATINDEX,
			WET_BULB,
			WBGT_INDEX
		} parametrization;
		
		SunObject sun;
		static const double soil_albedo, snow_albedo, snow_thresh;
		parametrization model;
};

} //end namespace mio

#endif
