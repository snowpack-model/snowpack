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
#ifndef TERRAINRADIATIONHELBIG_H
#define TERRAINRADIATIONHELBIG_H

#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/VFSymetricMatrix.h>
#include <alpine3d/ebalance/ViewFactorsHelbig.h>
#include <alpine3d/ebalance/ViewFactorsSectors.h>
#include <alpine3d/ebalance/ViewFactorsCluster.h>

#include <ctime>

//Optimisation #6 by GS : Function to sort CellList's array by radiation
inline int CellsRadComparator_Helbig(const void * cell1, const void * cell2)
{
	const double c1 = ((CellsList *)cell1)->radiation;
	const double c2 = ((CellsList *)cell2)->radiation;
	if (c1 > c2) return -1;
	if (c1 < c2) return 1;
	return 0;
}

/**
 * @class  TerrainRadiationHelbig
 * @brief Radiosity terrain radiation.
 * For each cell of the domain, a view factor to every other visible cell of the domain is computed. Then
 * the contributions from the terrain reflected radiations of every visible cell are added until a convergence 
 * criteria is reached. This is described in Helbig, Nora, Henning Löwe, and Michael Lehning, 
 * <i>"Radiosity approach for the shortwave surface radiation balance in complex terrain"</i>, 
 * Journal of the Atmospheric Sciences, <b>66.9</b> ,2009, pp 2900-2912.
 * 
 * The following keys in the [EBalance] section control the behavior of the model:
 *      - itEps_SW:  stopping tolerance/iteration error for shortwave radiation --(default: 0.4)--
 *            (make it larger to stop the iteration earlier, i.e. less terrain reflections are
 *            included --> with 10 % (0.1) a good accuracy is obtained)
 *      - itEps1_SW: stopping tolerance/iteration error for shortwave radiation --(default: 0.1)--
 *            (make it lower to stop the iteration earlier, i.e. less terrain reflections are
 *             included)
 *      - itEps_LW: stopping tolerance/iteration error for longwave radiation  --(default: 0.4)--
 *            note that two stopping critera are used
 *      - sw_radius: the distance radius (in m) around each grid cell until which terrain reflection
 *            is taken into account --(default: 3000.)--
 *      - lw_radius: the distance radius (in m) around each grid cell until which longwave emission
 *            is taken into account --(default: 1500. (Landl(2007)))--
 *      - vf_in_ram: Progressive Refinement iteration is used; but you can enable storing the view factors 
 *            in memory (in case of sufficient memory capacity --> faster terrain radiation computation)
 *      - sub_crit: substructuring threshold (in %) of the patches in the view factor computation --(default: 0.4)--
 *            note that mccluney(1994) proposes 0.1, i.e. for accurate view factors / radiation exchange
 *            computations a threshold of least 0.1 should be used
 * 
 * It is also possible to read the view factors from file or write them to a file. The following keys,
 * either in the [Input] or in the [Output] sections define the file names:
 *       - vf_file: file containing the sky view factors
 *       - tvfarea: file containing the terrain view factors x surface
 *
 */
class TerrainRadiationHelbig: public TerrainRadiationAlgorithm {

	public:
		TerrainRadiationHelbig(const mio::Config& i_cfg, const mio::DEMObject& dem_in, const int& i_nbworkers, const std::string& method);

		void getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain);
		void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta,
		              const mio::Array2D<double>& rh, const mio::Array2D<double>& ilwr);

	private:
		mio::DEMObject dem;
		//double iswr_ref;
		//double ea_ref;
		double sw_radius;

		int dimx, dimy;
		double cellsize;

		double itEps1_SW;
		double mean_glob_start;      // mean reflectable direct and diffuse sky shortwave radiation
		int LW_distance_index;			//for LW: index of the maximum emitting distance
		const static int NB_UNROLL = 3; //Define the number of unrolled calculation in loop to LWTerrainRadiationStep

		mio::Array2D<double> meteo2d_ta, meteo2d_rh;
		mio::Array2D<double> lw_t,lwi, lw_sky;

		double lw_eps_stern;			//stopping criterion
		double max_glob_start;       // maximum of reflectable direct and diffuse sky shortwave radiation

		double itEps_SW;
		double itEps_LW;
		double max_alb;              // max ground albedo

		mio::Array2D<double> total_diff, tdir, tdiff, sw_t, glob_start, glob_h_isovf, glob_h, t_snowold;
		//SnowpackInterface *snowpack;

		double lw_start_l1;

		//VFSymetricMatrix<float, double> vf; // view factor matrix with dynamic dimension

		ViewFactorsHelbig viewFactorsHelbigObj;
		ViewFactorsSectors viewSectorFactorsObj;
		ViewFactorsCluster viewFactorsClusterObj;
		//mio::Array2D<double> tdir;
		//mio::Array2D<double> tdiff;

		mio::Array2D<double> albedo_grid, meteo2d_ilwr;

		std::vector<CellsList> lwt_byCell;

		void Compute();
		int SWTerrainRadiationStep(const double threshold_itEps_SW, int *c, int *d, unsigned int n, const clock_t t0);
		int LWTerrainRadiationStep(const double threshold_itEps_LW, const int itMax_LW, const int i_shoot, const int j_shoot, unsigned int n, const clock_t t0);
		int ComputeTerrainRadiation(const bool& day, int c, int d);
		void ComputeRadiationBalance();
		void InitializeTerrainSwSplitting(const int& i, const int& j,
		                                  int *i_max_unshoot, int *j_max_unshoot, double *diffmax_sw);
		void InitializeTerrainRadiation(const bool& day, int *c, int *d);
		void fillSWResultsGrids(const bool& day);

		void InitializeLW (const int i, const int j);

		static inline void CalculateIndex(const int indice, const int distance_max, int dim, int * min, int * max);
		static inline void LWTerrainRadiationCore(const double bx2, const int j_shoot, const double z_shoot, const int j, const double z, const double cellsize, const double t_snow_shoot,const double t_snow_shoot_value ,const double t_a, const double vf, double * lwi, int * s);

};

#endif
