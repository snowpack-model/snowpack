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
#ifndef TAUCLDGENERATOR_H
#define TAUCLDGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>
#include <meteoio/meteoLaws/Sun.h>
#include <meteoio/dataClasses/DEMObject.h>

#include <map>

namespace mio {

/**
 * @class TauCLDGenerator
 * @ingroup parametrizations
 * @brief Atmospheric transmissivity generator.
 * @details
 * Generate the atmospheric transmissivity (or clearness index, see \ref meteoparam) from other parameters. If a parameter
 * named "CLD" is available, it will be interpreted as cloud cover / cloudiness: in okta between
 * 0 (fully clear) and 8 (fully cloudy). For synop reports, it is possible to include a value of exactly 9 (sky obstructed
 * from view by fog, heavy precipitation...) that will be transparently reset to 8 (fully cloudy).
 *
 * If no such parameter is available, the atmospheric transmissivity is calculated from the solar index
 * (ratio of measured iswr to potential iswr, therefore using the current location (lat, lon, altitude) and ISWR
 * to parametrize the cloud cover). This relies on (Kasten and Czeplak, 1980).
 *
 * It takes the following (optional) argument:
 *    - CLD_TYPE: cloudiness model, either LHOMME, KASTEN or CRAWFORD (default: KASTEN, see AllSkyLWGenerator for the references of the papers. 
 * Please also note that CRAWFORD and LHOMME are exactly identical as the both simply consider that the cloudiness is <em>1-clearness_index</em>);
 *    - SHADE_FROM_DEM: if set to true, the DEM defined in the [Input] section will be used to compute the shading;
 *    - INFILE: a file containing the horizon for some or all stations (see the ProcShade filter for the format);
 *    - OUTFILE: a file to write the computed horizons to. If some horizons have been read from INFILE, they will also be written out in OUTFILE;
 *    - USE_RSWR. If set to TRUE, when no ISWR is available but RSWR and HS are available, a ground albedo is estimated
 * (either soil or snow albedo) and ISWR is then computed from RSWR. Unfortunatelly, this is not very precise... (thus default is false)
 *    - USE_RAD_THRESHOLD: when relying on measured ISWR to parametrize the cloudiness, there is a risk that the measuring station would
 * stand in a place where it is shaded by the surrounding terrain at some point during the day. This would lead to an overestimation 
 * of the cloudiness that is undesirable. In this case, it is possible to set USE_RAD_THRESHOLD to TRUE in order to interpolate the cloudiness
 * over all periods of low radiation measured ISWR. This is less performant that only considering the solar elevation but improves things
 * in this specific scenario.
 * 
 * Please not that it is possible to combine SHADE_FROM_DEM and INFILE: in this case, stations that don't have any horizon in the provided
 * INFILE will be computed from DEM. If none is available, a fixed 5 degrees threshold is used.
 *
 * @code
 * [Generators]
 * TAU_CLD::generator1     = TAU_CLD
 * TAU_CLD::arg1::use_rswr = false
 * @endcode
 */
class TauCLDGenerator : public GeneratorAlgorithm {
	public:
		typedef enum CLF_PARAMETRIZATION {
			DEFAULT,	//will be mapped to KASTEN
			CLF_LHOMME,
			KASTEN,
			CLF_CRAWFORD
		} clf_parametrization;

		TauCLDGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ, const Config &i_cfg);
		~TauCLDGenerator();
		
		bool generate(const size_t& param, MeteoData& md);
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo);
	protected:
		double getCloudiness(const MeteoData& md, SunObject& sun, bool &is_night);
		double getClearness(const double& cloudiness) const;
		static std::vector< std::pair<double,double> > computeMask(const DEMObject& i_dem, const StationData& sd);
		
		std::map< std::string, std::pair<double, double> > last_cloudiness; //as < station_hash, <julian_gmt, cloudiness> >
		std::map< std::string , std::vector< std::pair<double,double> > > masks;
		std::string horizons_outfile;
		const Config &cfg;
		DEMObject dem;
		clf_parametrization cloudiness_model;
		bool use_rswr, use_rad_threshold;
		bool write_mask_out, has_dem;
};

} //end namespace mio

#endif
