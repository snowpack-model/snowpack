/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SNOWPACKINTERFACEWORKER_H
#define SNOWPACKINTERFACEWORKER_H

#include <iostream>
#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

class SnGrids {
	public:
		/// \anchor SnGrids this enum provides names for possible Snowpack grids
		enum Parameters {firstparam=0,
                                          TA=firstparam, ///< Air temperature
                                          RH, ///< Relative humidity
                                          VW, ///< Wind velocity
                                          ISWR, ///< Incoming short wave radiation
                                          ISWR_DIFF, ///< Incoming short wave, diffuse
                                          ISWR_DIR, ///< Incoming short wave, direct
                                          ILWR, ///< Incoming long wave radiation
                                          HS, ///< Height of snow
                                          PSUM, ///< Water equivalent of precipitations, either solid or liquid
                                          PSUM_PH, ///<  Precipitation phase, between 0 (fully solid) and 1 (fully liquid)
                                          TSG, ///< Temperature ground surface
                                          TSS, ///< Temperature snow surface
                                          TS0, ///< Temperature soil surface
                                          TSNOW, ///< Snow temperature at depth xxx m
                                          TSNOW_AVG, ///< Average snow temperature in the top xxx m
                                          RHOSNOW_AVG, ///< Average snow density in the top xxx m
                                          TSOIL, ///< Temperature within the soil, at a given depth
                                          SWE, ///< Snow Water Equivalent
                                          RSNO, ///< Snow mean density
                                          TOP_ALB, ///< Albedo from the top (ie above canopy)
                                          SURF_ALB, ///< Albedo of the surface (ie below canopy)
                                          SP, ///< sphericity
                                          RB, ///< bond radius
                                          RG, ///< grain radius
                                          N3, ///< grain Coordination number
                                          MS_SNOWPACK_RUNOFF, ///< runoff on the surface of the soil (vitual lysimeter)
                                          MS_SOIL_RUNOFF, ///< runoff at the bottom of the snow/soil column
                                          SFC_SUBL, ///< The mass loss or gain of the top element due to snow (ice) sublimating
                                          STORE, ///< internal usage (precipitation events that are delayed because they are too small)
                                          GLACIER, ///< mask showing the glaciated pixels
                                          GLACIER_EXPOSED, ///< mask showing the exposed glaciated pixels (ie not snow covered)
                                          lastparam=GLACIER_EXPOSED};


		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const size_t& parindex);
		static size_t getParameterIndex(const std::string& parname);

	private:
		static std::vector<std::string> paramname;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

#include <alpine3d/DataAssimilation.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/runoff/Runoff.h>

class SnowpackInterfaceWorker
{
	public:
		SnowpackInterfaceWorker(const mio::Config& io_cfg,
		                        const mio::DEMObject& dem_in,
		                        const mio::Grid2DObject& landuse_in,
		                        const std::vector< std::pair<size_t,size_t> >& pts_in,
		                        const std::vector<SnowStation*>& snow_stations,
		                        const size_t offset_in);

		~SnowpackInterfaceWorker();

		void setUseDrift(const bool useDrift_in) {useDrift = useDrift_in;}
		void setUseEBalance(const bool useEBalance_in) {useEBalance = useEBalance_in;}

		void getOutputSNO(std::vector<SnowStation*>& snow_station) const;
		void getOutputSpecialPoints(std::vector<SnowStation*>& ptr_snow_pixel, std::vector<CurrentMeteo*>& ptr_meteo_pixel, 
		                                                std::vector<SurfaceFluxes*>& ptr_surface_flux);
		void clearSpecialPointsData();

		mio::Grid2DObject getGrid(const SnGrids::Parameters& param) const;
		double& getGridPoint(const SnGrids::Parameters& param, const size_t& ii, const size_t& jj);

		void runModel(const mio::Date& julian,
		              const mio::Grid2DObject &psum,
		              const mio::Grid2DObject &psum_ph,
		              const mio::Grid2DObject &rh,
		              const mio::Grid2DObject &ta,
		              const mio::Grid2DObject &vw,
		              const mio::Grid2DObject &mns,
		              const mio::Grid2DObject &shortwave,
		              const mio::Grid2DObject &diffuse,
		              const mio::Grid2DObject &longwave,
		              const double solarElevation);

		static int round_landuse(const double& landuse_dbl);
		static bool skipThisCell(const double& landuse_val, const double& dem_val);
		static bool is_special(const std::vector< std::pair<size_t,size_t> >& pts_in, const size_t& ix, const size_t& iy);

	private:
		void initGrids(std::vector<std::string>& params);
		void gatherSpecialPoints(const CurrentMeteo& meteoPixel, const SnowStation& snowPixel, const SurfaceFluxes& surfaceFlux);
		void fillGrids(const size_t& ii, const size_t& jj, const CurrentMeteo& meteoPixel, const SnowStation& snowPixel, const SurfaceFluxes& surfaceFlux);

	private:
		SnowpackConfig sn_cfg; // created on element
		Snowpack sn;
		Meteo meteo;
		Stability stability;

		const mio::DEMObject dem;
		const size_t dimx, dimy, offset;
		std::vector<SnowStation*> SnowStations; // Save different Pixel values
		std::vector<bool> isSpecialPoint;

		const mio::Grid2DObject landuse;
		mio::Grid2DObject store;
		std::map< SnGrids::Parameters, mio::Grid2DObject > grids;

		// cache special point data for output on master process:
		std::vector<SnowStation> snow_pixel;
		std::vector<CurrentMeteo> meteo_pixel;
		std::vector<SurfaceFluxes> surface_flux;
		
		double calculation_step_length;
		double height_of_wind_value;
		double soil_temp_depth, snow_temp_depth, snow_avg_temp_depth, snow_avg_rho_depth;
		bool useDrift, useEBalance, useCanopy;
};

#endif
