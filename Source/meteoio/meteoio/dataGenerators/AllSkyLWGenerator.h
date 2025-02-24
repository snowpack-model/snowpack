// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013-2021 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
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
#ifndef ALLSKYGENERATOR_H
#define ALLSKYGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>
#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/meteoLaws/Sun.h>

#include <map>

namespace mio {

/**
 * @class AllSkyLWGenerator
 * @ingroup parametrizations
 * @brief ILWR all sky parametrization
 * @details
 * Using air temperature (TA) and relative humidity (RH) and optionnally cloud transmissivity (TAU_CLD),
 * this offers the choice of several all-sky parametrizations, with the following arguments:
 *  - TYPE: specify which parametrization should be used, from the following:
 *      - LHOMME -- from Lhomme et al. -- <i>"Estimating downward long-wave 
 * radiation on the Andean Altiplano"</i>, Agric. For. Meteorol., <b>145</b>, 2007, 
 * pp 139–148, doi:10.1016/j.agrformet.2007.04.007
 *      - CARMONA -- from Carmona et al., <i>"Estimation of daytime downward 
* longwave radiation under clear and cloudy skies conditions over a sub-humid region."</i> Theoretical and applied climatology <b>115.1-2</b> (2014): 281-295.
 *      - CRAWFORD -- from Crawford and Duchon, <i>"An Improved Parametrization for Estimating Effective Atmospheric Emissivity for Use in Calculating Daytime
 * Downwelling Longwave Radiation"</i>, Journal of Applied Meteorology, <b>38</b>, 1999, pp 474-480
 *      - OMSTEDT -- from Omstedt, <i>"A coupled one-dimensional sea ice-ocean model applied to a semi-enclosed basin"</i>,
 * Tellus, <b>42 A</b>, 568-582, 1990, DOI:10.1034/j.1600-0870.1990.t01-3-00007.
 *      - KONZELMANN -- from Konzelmann et al., <i>"Parameterization of global and longwave incoming radiation
 * for the Greenland Ice Sheet."</i> Global and Planetary change <b>9.1</b> (1994): 143-164.
 *      - UNSWORTH -- from Unsworth and Monteith, <i>"Long-wave radiation at the ground"</i>,
 * Q. J. R. Meteorolo. Soc., Vol. 101, 1975, pp 13-24 coupled with a clear sky emissivity following (Dilley, 1998).
 *  - it also takes all the arguments of TauCLDGenerator, including the option to provide an Horizon file (see the ProcShade filter for the format).
 *
 * If no cloud transmissivity is provided in the data, it is calculated from the solar index (ratio of measured iswr to potential iswr, therefore using
 * the current location (lat, lon, altitude) and ISWR to parametrize the cloud cover). This relies on (Kasten and Czeplak, 1980) by default 
 * except for Crawford and Lhomme that provide their own parametrizations (it can be forced through the TauCLDGenerator options).
 * The last evaluation of cloud transmissivity is used all along during the times when no ISWR is available (as it is night, or the station stands in 
 * the shade or the ISWR measurement is missing) and the last valid value is not too old (ie. no more than 1 day old).  The example below gives an
 * example use, providing a file with the horizons for the stations (see TauCLDGenerator) but not relying on a DEM. Of course, it is possible to
 * use the generator without any DEM or horizon file!
 * 
 * @code
 * [Generators]
 * ILWR::generator1 = allsky_LW
 * ILWR::arg1::type = Omstedt
 * ILWR::arg1::infile = input/horizons.txt
 * ILWR::arg1::shade_from_dem = FALSE
 * ILWR::arg1::use_rswr = FALSE
 * @endcode
 *
 * Finally, it is recommended to also use a clear sky generator (declared after this one) for the case of no available short wave measurement
 * (by declaring the ClearSky generator \em after AllSky).
 *
 * The graph below shows the comparison between measured and modeled ILWR depending on the chosen parametrization. The measured data (ISWR, TA, RH and the reference ILWR)
 * comes from the Weissfluhjoch *WFJ AWS (2691m, Davos, Switzerland) for the 2010-08-01 -- 2019-08-01 period with half-hourly resolution. The data has been binned every 5 W/m²,
 * the black dots represent the average of the bin, the greay area contains every data point (ie it shows the minimum and maximum data) while the brown area is defined as average±σ.
 * \image html all_sky_ilwr_cmp.svg "Comparison between measured and parametrized ILWR at the Weissfluhjoch *WFJ station (2691m, Davos, Switzerland) for the 2010-08-01 – 2019-08-01 period"
 * \image latex all_sky_ilwr_cmp.eps "Comparison between measured and parametrized ILWR at the Weissfluhjoch *WFJ station (2691m, Davos, Switzerland) for the 2010-08-01 – 2019-08-01 period" width=0.9\textwidth
 *
 */
class AllSkyLWGenerator : public TauCLDGenerator {
	public:
		//TauCLDGenerator will do its own arguments parsing, then AllSkyLWGenerator
		//so make sure that we don't use here the same name as an argument to TauCLDGenerator!
		AllSkyLWGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ, const Config &i_cfg) : TauCLDGenerator(vecArgs, i_algo, i_section, TZ, i_cfg), model(OMSTEDT) {parse_args(vecArgs); }
		
		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo) override;
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo) override;
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs) override;

		typedef enum PARAMETRIZATION {
			LHOMME,
			CARMONA,
			CRAWFORD,
			KONZELMANN,
			OMSTEDT,
			UNSWORTH
		} parametrization;

		parametrization model;
};

} //end namespace mio

#endif
