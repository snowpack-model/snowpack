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
#ifndef TERRAINRADIATION_H
#define TERRAINRADIATION_H


#include <meteoio/MeteoIO.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/ViewFactors.h>
#include <alpine3d/ebalance/VFSymetricMatrix.h>
#include <alpine3d/ebalance/ViewFactorsSectors.h>
#include <alpine3d/ebalance/ViewFactorsCluster.h>

#include <ctime>

/**
 * @class  TerrainRadiation
 * @brief Full computation of the terrain reflected radiation by ray-tracing.
 * Ray tracing is used between all cells to compute the terrain view factors. Once these are known,
 * a progressive refinement approach is used to compute the exchange of radiation between all cells.
 *
 * The progressive refinement approach is controlled by the following keys in the [EBalance] section:
 *   - itEps_SW
 *   - itEps1_SW
 *   - itEps_LW
 *   - sw_radius, the search radius for short wave
 *   - lw_radius, the search radius for long wave
 *
*/

class TerrainRadiation : public TerrainRadiationAlgorithm {
	public:
		TerrainRadiation(const mio::Config& i_cfg, const mio::DEMObject &dem_in, const int& i_nbworkers, const std::string& method);

		void getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain);
		void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta,
		              const mio::Array2D<double>& rh, const mio::Array2D<double>& ilwr);

	private:
		void Compute();
		int SWTerrainRadiationStep(const double threshold_itEps_SW, size_t& c, size_t& d,
		                           mio::Array2D<double>& sw_Q_out, double& Rmy);
		int LWTerrainRadiationStep(const double threshold_itEps_LW, const int itMax_LW, const int i_shoot,
		                           const int j_shoot, unsigned int n);
		int ComputeTerrainRadiation(const double& threshold_itEps_SW, const bool& day, size_t& c, size_t& d,
		                            mio::Array2D<double>& sw_Q_out, double& Rmy);
		void ComputeProgressiveRefinement();
		void InitializeTerrainSwSplitting(size_t& i_max_unshoot, size_t& j_max_unshoot, mio::Array2D<double>& sw_Q_out);
		void InitializeTerrainRadiation(const bool& day);
		void fillSWResultsGrids(const bool& day);
		void InitializeLW (const int i, const int j, double &lw_eps_stern);


		static inline void CalculateIndex(const int indice, const int distance_max, int dim, int * min, int * max);
		static inline void LWTerrainRadiationCore(const double bx2, const int j_shoot, const double z_shoot, const int j,
		                                          const double z, const double cellsize, const double t_snow_shoot,
		                                          const double t_snow_shoot_value, const double t_a, const double vf,
		                                          double * lwi, int * s);

		ViewFactors viewFactorsObj;
		mio::DEMObject dem;
		size_t dimx, dimy;
		double cellsize;

		mio::Array2D<double> albedo_grid, meteo2d_ilwr, meteo2d_ta, meteo2d_rh;
		mio::Array2D<double> total_diff, sw_t, glob_start, glob_h_isovf, glob_h, t_snowold;
		mio::Array2D<double> lw_t,lwi, lw_sky;
		std::vector<CellsList> lwt_byCell;

		mio::Array2D<double> tdir, tdiff;
		mio::Array2D<double> Mmy;
		mio::Array2D<double> Mt;

		double itEps1_SW;
		double mean_glob_start;      // mean reflectable direct and diffuse sky shortwave radiation
		double lw_eps_stern;			//stopping criterion
		double max_glob_start;       // maximum of reflectable direct and diffuse sky shortwave radiation
		double itEps_SW;
		double itEps_LW;
		double max_alb;              // max ground albedo
		double lw_start_l1;
		double sw_radius;

		double threshold_itEps_SW;

		size_t startx, nx;
		size_t starty, ny;

		int LW_distance_index;			//for LW: index of the maximum emitting distance

		const static int NB_UNROLL = 3; //Define the number of unrolled calculation in loop to LWTerrainRadiationStep
};

#endif
