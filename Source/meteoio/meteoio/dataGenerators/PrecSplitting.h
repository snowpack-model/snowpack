// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013-2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
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
#ifndef PRECSPLITTING_H
#define PRECSPLITTING_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class PrecSplitting
 * @ingroup parametrizations
 * @brief Generate precipitation phase or splitting according to the selected method
 * @details
 * In order to handle solid and liquid precipitation, there are two possibilities: either as total amount and phase (respectively "PSUM" and "PSUM_PH")
 * or as solid and liquid amounts (respectively "PSUM_S" and "PSUM_L"). This generator can convert from one representation to the other one or 
 * generate the missing components (for example, the precipitation phase associated with a given precipitation amout based on a splitting model).
 * 
 * The component that will be generated depends on the parameter name, so only the following parameters are supported: PSUM, PSUM_PH, PSUM_L, PSUM_S.
 * Whenever a component is missing, a splitting model is called to compute the splitting. It is therefore recommended to configure a splitting model.
 * 
 * It takes the following arguments:
 *  - TYPE: the splitting method to use, any of the following:
 *     - NONE: don't do any splitting (default). This can be usefull to overwrite an inherited generator (from an imported file)
 *     - THRESH: a provided fixed air temperature threshold splits precipitation as either fully solid or fully liquid
 *     - RANGE: two air temperature thresholds provide the lower and upper range for fully solid / fully liquid precipitation.
 *                 Within the provided range, a linear transition is assumed.
 *  - SNOW: when using a fixed air temperature threshold, this gives the snow/rain threshold (in K). When using two air temperatures
 * thresholds, this provides the temperature below which only solid precipitation is found (in K);
 *  - RAIN: when using two air temperatures thresholds, this provides the temperature above which only liquid precipitation is found (in K);
 *
 * For example, to generate the precipitation phase when only the precipitation amout is available, using a simple temperature
 * threshold for doing the splitting:
 * @code
 * [Generators]
 * PSUM_PH::generator1 = PRECSPLITTING
 * PSUM_PH::arg1::type = THRESH
 * PSUM_PH::arg1::snow = 274.35
 * @endcode
 * 
 * To generate the liquid and solid amounts from the precipitation sum and phase, relying on a simple temperature threshold for the splitting
 * (if the precipitation phase is available, it will be used instead of calling the splitting model):
 * @code
 * [InputEditing]
 * ;using the '*' station ID, this will apply to all stations
 * *::edit1           = CREATE
 * *::arg1::algorithm = PRECSPLITTING
 * *::arg1::param     = PSUM_L
 * *::arg1::type      = THRESH
 * *::arg1::snow      = 274.35
 * 
 * *::edit2           = CREATE
 * *::arg2::algorithm = PRECSPLITTING
 * *::arg2::param     = PSUM_S
 * *::arg2::type      = THRESH
 * *::arg2::snow      = 274.35
 * @endcode
 * 
 * To generate the precipitation sum and phase from the liquid and solid amounts without any fallback algorithm for the splitting:
 * @code
 * [Generators]
 * PSUM::generator1 = PRECSPLITTING
 * PSUM_PH::generator1 = PRECSPLITTING
 * @endcode
 * 
 * @note When generating PSUM_L / PSUM_S, you most probably also need to set their resampling to "accumulate", like for PSUM...
 */
class PrecSplitting : public GeneratorAlgorithm {
	public:
		PrecSplitting(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ)
			: GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), model(NONE), where( section+"::"+algo ), fixed_thresh(IOUtils::nodata),
			range_start(IOUtils::nodata), range_norm(IOUtils::nodata) { parse_args(vecArgs); }

		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo) override;
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo) override;

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs) override;
		bool generatePSUM_S(double& value, MeteoData& md) const;
		bool generatePSUM_L(double& value, MeteoData& md) const;
		bool generatePSUM_PH(double& value, MeteoData& md) const;
		bool generatePSUM(double& value, MeteoData& md) const;
		bool runModel(double &value, MeteoData& md) const;

		typedef enum PARAMETRIZATION {
			THRESH,
			RANGE,
			NONE
		} parametrization;
		parametrization model;
		const std::string where;
		double fixed_thresh, range_start, range_norm;
};

} //end namespace mio

#endif
