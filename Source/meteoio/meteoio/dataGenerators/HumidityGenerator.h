// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef HUMIDITYGENERATOR_H
#define HUMIDITYGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class HumidityGenerator
 * @ingroup parametrizations
 * @brief (Relative/Specifc) Humidity or Dew Point temperature generator.
 * @details
 * Since Absolute Humidity (AH), Relative Humidity (RH), Specific Humidity (QI) or Dew Point Temperature (TD)
 * are all measures of the atmospheric humidity, this generates any one of these from any other one (depending on
 * what is available).The parameter that should be generated is provided as argument (default: RH).
 * @code
 * [Generators]
 * RH::generator1 = HUMIDITY
 * RH::arg1::type = RH
 * @endcode
 */
class HumidityGenerator : public GeneratorAlgorithm {
	public:
		HumidityGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ)
			: GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), type(GEN_RH) { parse_args(vecArgs); }
		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo) override;
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo) override;

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs) override;

		static bool generateAH(double& value, MeteoData& md);
		static bool generateRH(double& value, MeteoData& md);
		static bool generateTD(double& value, MeteoData& md);
		static bool generateQI(double& value, MeteoData& md);

		typedef enum GEN_TYPE {
			GEN_AH,
			GEN_RH,
			GEN_TD,
			GEN_QI
		} gen_type;
		gen_type type;
};

} //end namespace mio

#endif
