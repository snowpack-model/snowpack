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
#ifndef CLEARSKYLWGENERATOR_H
#define CLEARSKYLWGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class ClearSkyLWGenerator
 * @ingroup parametrizations
 * @brief ILWR clear sky parametrization
 * @details
 * Using air temperature (TA) and relative humidity (RH), this offers the choice of several clear sky parametrizations with the
 * following arguments:
 *  - TYPE: specify the parametrization that should be used, any of the following:
 *      - BRUTSAERT -- from Brutsaert, <i>"On a Derivable Formula for Long-Wave Radiation From Clear Skies"</i>,
 * Journal of Water Resources Research, <b>11</b>, No. 5, October 1975, pp 742-744.
 *      - DILLEY -- from Dilley and O'Brien, <i>"Estimating downward clear sky
 * long-wave irradiance at the surface from screen temperature and precipitable water"</i>, Q. J. R. Meteorolo. Soc., <b>124</b>, 1998, pp 1391-1401.
 *      - PRATA -- from Prata, <i>"A new long-wave formula for estimating downward clear-sky radiation at the surface"</i>, Q. J. R. Meteorolo. Soc., <b>122</b>, 1996, pp 1127-1151.
 *      - CLARK -- from Clark & Allen, <i>"The estimation of atmospheric radiation for clear and
 * cloudy skies"</i>, Proceedings of the second national passive solar conference, <b>2</b>, 1978, p 676.
 *      - TANG -- from Tang et al., <i>"Estimates of clear night sky emissivity in the
 * Negev Highlands, Israel"</i>, Energy Conversion and Management, <b>45.11</b>, 2004, pp 1831-1843.
 *      - IDSO -- from Idso, <i>"A set of equations for full spectrum and 8 to 14 um and
 * 10.5 to 12.5 um thermal radiation from cloudless skies"</i>, Water Resources Research, <b>17</b>, 1981, pp 295-304.
 *
 * Please keep in mind that for energy balance modeling, this significantly underestimate the ILWR input.
 * @code
 * [Generators]
 * ILWR::generator1 = clearsky_LW
 * ILWR::arg1::type = Dilley
 * @endcode
 * 
 * The graph below shows the comparison between measured and modeled ILWR depending on the chosen parametrization. The measured data (ISWR, TA, RH and the reference ILWR)
 * comes from the Weissfluhjoch *WFJ AWS (2691m, Davos, Switzerland) for the 2010-08-01 -- 2019-08-01 period with half-hourly resolution. The data has been binned every 5 W/m²,
 * the black dots represent the average of the bin, the greay area contains every data point (ie it shows the minimum and maximum data) while the brown area is defined as average±σ.
 * \image html clear_sky_ilwr_cmp.svg "Comparison between measured and parametrized ILWR at the Weissfluhjoch *WFJ station (2691m, Davos, Switzerland) for the 2010-08-01 – 2019-08-01 period"
 * \image latex clear_sky_ilwr_cmp.eps "Comparison between measured and parametrized ILWR at the Weissfluhjoch *WFJ station (2691m, Davos, Switzerland) for the 2010-08-01 – 2019-08-01 period" width=0.9\textwidth
 *
 */
class ClearSkyLWGenerator : public GeneratorAlgorithm {
	public:
		ClearSkyLWGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ)
			: GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), model(BRUTSAERT) { parse_args(vecArgs); }
		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo) override;
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo) override;
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs) override;
		typedef enum PARAMETRIZATION {
			BRUTSAERT,
			DILLEY,
			PRATA,
			CLARK,
			TANG,
			IDSO
		} parametrization;
		parametrization model;
};

} //end namespace mio

#endif
