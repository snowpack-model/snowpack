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

#include <alpine3d/DataAssimilation.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/runoff/Runoff.h>
#include <alpine3d/MeteoObj.h> //for the SnGrids

class SnowpackInterfaceWorker
{
	public:
		SnowpackInterfaceWorker(const mio::Config& io_cfg,
		                        const mio::DEMObject& dem_in,
		                        const mio::Grid2DObject& landuse_in,
		                        const std::vector< std::pair<size_t,size_t> >& pts_in,
		                        const std::vector<SnowStation*>& snow_stations,
		                        const std::vector<std::pair<size_t,size_t> >& snow_stations_coord,
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
		              const mio::Grid2DObject &psum_tech,
		              const mio::Grid2DObject &rh,
		              const mio::Grid2DObject &ta,
		              const mio::Grid2DObject &tsg,
		              const mio::Grid2DObject &vw,
		              const mio::Grid2DObject &vw_drift,
		              const mio::Grid2DObject &dw,
		              const mio::Grid2DObject &mns,
		              const mio::Grid2DObject &shortwave,
		              const mio::Grid2DObject &diffuse,
		              const mio::Grid2DObject &longwave,
		              const double solarElevation);

		void grooming(const mio::Date &current_date, const mio::Grid2DObject &grooming_map);

		static int round_landuse(const double& landuse_dbl);
		static bool skipThisCell(const double& landuse_val, const double& dem_val);
		static bool is_special(const std::vector< std::pair<size_t,size_t> >& pts_in, const size_t& ix, const size_t& iy);
		static void uniqueOutputGrids(std::vector<std::string>& output_grids);
		void getLateralFlow(std::vector<SnowStation*>& snow_station);
		void setLateralFlow(const std::vector<SnowStation*>& snow_station);

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
		std::vector<std::pair<size_t,size_t> > SnowStationsCoord;
		std::vector<bool> isSpecialPoint;

		const mio::Grid2DObject landuse;
		mio::Grid2DObject store;
		mio::Grid2DObject erodedmass;
		std::map< SnGrids::Parameters, mio::Grid2DObject > grids;

		// cache special point data for output on master process:
		std::vector<SnowStation> snow_pixel;
		std::vector<CurrentMeteo> meteo_pixel;
		std::vector<SurfaceFluxes> surface_flux;
		std::vector<double> soil_temp_depths;

		double calculation_step_length;
		double height_of_wind_value;
		double snow_temp_depth, snow_avg_temp_depth, snow_avg_rho_depth;
		bool enable_simple_snow_drift;
		bool useDrift, useEBalance, useCanopy;
};

#endif
