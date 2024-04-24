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
 * (either soil or snow albedo) and ISWR is then computed from RSWR. Unfortunately, this is not very precise... (thus default is false)
 *    - USE_RAD_THRESHOLD: when relying on measured ISWR to parametrize the cloudiness, there is a risk that the measuring station would
 * stand in a place where it is shaded by the surrounding terrain at some point during the day. This would lead to an overestimation 
 * of the cloudiness that is undesirable. In this case, it is possible to set USE_RAD_THRESHOLD to TRUE in order to interpolate the cloudiness
 * over all periods of low radiation measured ISWR. This is less performant that only considering the solar elevation (with shading 
 * computed from DEM) but improves things in this specific scenario, when no shading is available.
 * 
 * \note
 * Please note that it is possible to combine SHADE_FROM_DEM and INFILE: in this case, stations that don't have any horizon in the provided
 * INFILE will be computed from DEM. It is also possible to define wildcard station ID in the horizon file. If SHADE_FROM_DEM has been set
 * to false and no INFILE has been provided, a fixed 5 degrees threshold is used.
 *
 * @code
 * [Generators]
 * TAU_CLD::generator1     = TAU_CLD
 * TAU_CLD::arg1::use_rswr = false
 * @endcode
 *
 * The horizon file contains on each line a station ID followed by an azimuth (in degrees, starting from North) and an elevation above the
 * horizontal (also in degrees). It is possible to define a wildcard station ID '*' to be used as fallback. The elevation for any given azimuth
 * will be linearly interpolated between the provided horizons before and after. If only one azimuth is given for a station ID, its horizon
 * elevation is assumed to be constant. See below an example horizon file, defining two station IDs (SLF2 and FLU2):
 * @code
 * SLF2 0 5
 * SLF2 45 25
 * SLF2 180 30
 * SLF2 245 10
 * FLU2 0 7
 * @endcode
 * 
 * \note 
 * In order to pre-compute the horizons for multiple stations, write an ini file reading all stations of interest, declare a dataGenerator 
 * that relies on TauCLDGenerator and run meteoio_timeseries on it for at least one time step. This will be enough to create the horizon file
 * for all stations...
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
		
		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo);
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo);
	protected:

		typedef struct CLOUDCACHE {
			CLOUDCACHE() : last_valid(std::make_pair(IOUtils::nodata, IOUtils::nodata)) {}
			CLOUDCACHE(const double& julian_gmt, const double& cloudiness) : last_valid(std::make_pair(IOUtils::nodata, IOUtils::nodata)) { addCloudiness(julian_gmt, cloudiness); }

			void addCloudiness(const double& julian_gmt, const double& cloudiness);

			std::pair<double, double> last_valid;
		} cloudCache;

		double interpolateCloudiness(const std::string& station_hash, const double& julian_gmt) const;
		double getCloudiness(const MeteoData& md);
		double computeCloudiness(const MeteoData& md, bool &is_night);
		double getClearness(const double& cloudiness) const;
		static std::vector< std::pair<double,double> > computeMask(const DEMObject& i_dem, const StationData& sd);
		double getHorizon(const MeteoData& md, const double& sun_azi);
		
		std::map< std::string, cloudCache > last_cloudiness; //as < station_hash, <julian_gmt, cloudiness> >
		std::map< std::string , std::vector< std::pair<double,double> > > masks;
		std::string horizons_outfile;
		const Config &cfg;
		DEMObject dem;
		SunObject sun;
		clf_parametrization cloudiness_model;
		bool use_rswr, use_rad_threshold;
		bool write_mask_out, use_horizons, from_dem;
};

} //end namespace mio

#endif
