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
#include <alpine3d/ebalance/TerrainRadiationHelbig.h>

#include <cstdio>
#include <algorithm>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"

using namespace mio;
using namespace std;

TerrainRadiationHelbig::TerrainRadiationHelbig(const mio::Config &cfg, const mio::DEMObject &dem_in, const int &, const std::string &method)
	: TerrainRadiationAlgorithm(method), dem(dem_in), dimx(dem.getNx()), dimy(dem.getNy()), cellsize(dem.cellsize), viewFactorsHelbigObj(cfg, dem_in),
	  viewSectorFactorsObj(cfg, dem_in), viewFactorsClusterObj(cfg, dem_in), lwt_byCell(dimx * dimy)
{
	cfg.getValue("itEps_SW", "EBalance", itEps_SW);
	cfg.getValue("itEps_LW", "EBalance", itEps_LW);
	cfg.getValue("itEps1_SW", "EBalance", itEps1_SW);
	cfg.getValue("sw_radius", "EBalance", sw_radius);
	cfg.getValue("lw_radius", "EBalance", lw_radius);

	total_diff.resize(dimx, dimy);
	tdir.resize(dimx, dimy);
	tdiff.resize(dimx, dimy);
	sw_t.resize(dimx, dimy);
	glob_start.resize(dimx, dimy);
	glob_h_isovf.resize(dimx, dimy);
	glob_h.resize(dimx, dimy);
	t_snowold.resize(dimx, dimy);
	total_terrain.resize(dimx, dimy);

	lw_t.resize(dimx, dimy);
	lwi.resize(dimx, dimy);
	lw_sky.resize(dimx, dimy);

	LW_distance_index = (int)ceil(lw_radius / cellsize);

	bool write_sky_vf = false;
	cfg.getValue("WRITE_SKY_VIEW_FACTOR", "output", write_sky_vf, IOUtils::nothrow);

	if (MPIControl::instance().master() && write_sky_vf)
	{
		std::cout << "[i] Writing sky view factor grid" << std::endl;
		mio::Array2D<double> sky_vf(dimx, dimy);
		getSkyViewFactor(sky_vf);
		mio::IOManager io(cfg);
		io.write2DGrid(mio::Grid2DObject(dem_in.cellsize, dem_in.llcorner, sky_vf), "SKY_VIEW_FACTOR");
	}
}

void TerrainRadiationHelbig::getRadiation(mio::Array2D<double> &direct, mio::Array2D<double> &diffuse,
										  mio::Array2D<double> &terrain, const mio::Array2D<double> &direct_unshaded_horizontal, const mio::Array2D<double> &total_ilwr, mio::Array2D<double> &sky_ilwr,
										  mio::Array2D<double> &terrain_ilwr,
										  double solarAzimuth, double solarElevation)
{
	std::cout << "[i] Computing Helbig radiation" << std::endl;
	tdir = direct;
	tdiff = diffuse;
	tot_ilwr = total_ilwr;
	Compute();
	terrain = total_terrain;
	diffuse = tdiff;
}

void TerrainRadiationHelbig::getSkyViewFactor(mio::Array2D<double> &o_sky_vf)
{
	viewFactorsHelbigObj.getSkyViewFactor(o_sky_vf);
}

void TerrainRadiationHelbig::setMeteo(const mio::Array2D<double> &albedo, const mio::Array2D<double> &ta)
{
	meteo2d_ta = ta;
	albedo_grid = albedo;
	max_alb = albedo_grid.getMax();
}

void TerrainRadiationHelbig::Compute()
{
	ComputeRadiationBalance();
}

int TerrainRadiationHelbig::SWTerrainRadiationStep(const double threshold_itEps_SW, int &i_max_unshoot, int &j_max_unshoot, unsigned int n, const clock_t t0)
{
	// Computation of shortwave terrain radiation
	// At every iteration step, a reference grid cell reflects ('shoots') radiation to every other grid cell
	// The coordinates of the next cell with the most unshoot radiation is writen in (*c,*d)
	const int i_shoot = i_max_unshoot, j_shoot = j_max_unshoot; //variables for the grid cell with the most unshot radiation
	double diffmax_sw = 0.;										// reference product for detecting the grid cell with most unshot shortwave radiation
	const double diffmax_thres = 0.;							// threshold for when to stop looking for maximum cells
	double eps_stern = 0.;										// stopping criterion
	double be = 0.;												// matrix difference B^(k) - E resp. (L_sky + L_terrain) - L_sky
	unsigned int s = 0;											// counts gathering patches, i.e. those patches within limited distance radius

	//variable needed to optimize the computational speed
	const double sw_radius2 = sw_radius * sw_radius;
	const double z_shoot = dem.grid2D(i_shoot, j_shoot);
	const double sw_shoot = sw_t(i_shoot, j_shoot);

	// every grid cell ij receives radiation from the chosen one with coordinates ab
	for (int j = 0; j < dimy; j++)
	{
		for (int i = 0; i < dimx; i++)
		{
			const double bx = (i_shoot - i) * cellsize;	  // bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
			const double by = (j_shoot - j) * cellsize;	  // by = dy (m) going from (i_shoot,j_shoot) to (i,j)
			const double bz = z_shoot - dem.grid2D(i, j); // bz = dz (m) going from (i_shoot,j_shoot) to (i,j)
			// distance between the surfaces
			const double dist2 = ((bx * bx) + (by * by) + (bz * bz));
			// the shortwave radiation will only be emitted if the radius is lower than
			// a preassumed radius (sw_radius) -> changeable in: sw_t(a,b) * sw_radius
			//please note that it also excludes the current shooting cell itself!
			if (dist2 <= sw_radius2 && dist2 > 0.)
			{
				// This call either
				// computes the ViewFactor (vf_in_ram=false) or
				// take the pre-computed ViewFactor (vf_in_ram=true)
				// This differentian according to vf_in_ram is done in ViewFactorsHelbig.cc
				const double viewFactor = viewFactorsHelbigObj.GetViewfactor(i, j, i_shoot, j_shoot);
				// calculation of the reflected amount of radiation of ab that arrives at ij
				const double rad = viewFactor * albedo_grid(i, j) * sw_shoot;
				// the received amount is added to the reflectable amount of radiation (unshot)
				sw_t(i, j) += rad;
				// in addition the received amount is added to the total radiation at ij
				glob_h(i, j) += rad;
				total_terrain(i, j) += viewFactor * sw_shoot;

				//available radiation to reflect at the current cell
				const double reflected_i_j = viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) * sw_t(i, j);

				// the next ('shooting') grid cell with the most unshot radiation is selected
				// taking into account albedo, actual (slope) area size, reflected shortwave
				// radiation and the sum of total terrain view factor
				if (reflected_i_j > diffmax_sw)
				{
					if (i != i_shoot && j != j_shoot)
					{
						diffmax_sw = reflected_i_j;
						i_max_unshoot = i;
						j_max_unshoot = j;
					}
				}
				s++;
				// the first stopping criterion:
				// |Delta B^(k) * Sum_j(Fij) * A|_1
				eps_stern += fabs(reflected_i_j);
				// the second stopping criterion:
				// |B^(k) - E|_1
				be += fabs(glob_h(i, j) - glob_start(i, j));
			}
		}
	}

	// set the reflected 'shot' radiation at that cell to zero
	sw_t(i_shoot, j_shoot) = 0.;

	// check for stopping the iteration
	if (eps_stern <= threshold_itEps_SW || be >= (itEps1_SW *
												  (mean_glob_start * dimx * dimy)))
	// || diffmax_sw<=diffmax_thres)
	{
		printf("[i] EBALANCE: SW converged after n=%u steps eps_stern <= %2.10f be <= %2.10f with %u gathering patches\n", n, eps_stern, be, s);
		printf("    time for SW %f seconds\n", (double)((clock() - t0) / CLOCKS_PER_SEC));
		fflush(stdout);
		return 1;
	}

	return 0;
}

int TerrainRadiationHelbig::LWTerrainRadiationStep(const double threshold_itEps_LW, const int itMax_LW,
												   int &i_max_unshoot_lw, int &j_max_unshoot_lw, unsigned int n, const clock_t t0)
{

	// Computation of longwave terrain radiation
	// At every iteration step, one chosen reference grid cell (*e,*f) reflects ('shoots') radiation to every other grid cell
	// within a given distance (LW_distance_index in cells units). Every grid cell is only once the emitting cell (no multiple reflections)
	//counts gathering patches, i.e. those patches within limited distance radius
	int s = 0;
	// stopping criterion
	double eps_stern = 0;
	// reference product for detecting the grid cell with most unshot longwave radiation
	double diffmax_lw = 0.;

	// The coordinates of the next cell with the most unshoot radiation is writen in (*c,*d)
	//variables for the grid cell with the most unshot radiation
	const int i_shoot = i_max_unshoot_lw, j_shoot = j_max_unshoot_lw;

	//variable needed to optimize the computational speed
	const double lw_radius2 = lw_radius * lw_radius;
	const double z_shoot = dem.grid2D(i_shoot, j_shoot);

	// every grid cell ij receives radiation from the chosen one with coordinates ab
	for (int j = 0; j < dimy; j++)
	{
		for (int i = 0; i < dimx; i++)
		{
			const double bx = (i_shoot - i) * cellsize;	  // bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
			const double by = (j_shoot - j) * cellsize;	  // by = dy (m) going from (i_shoot,j_shoot) to (i,j)
			const double bz = z_shoot - dem.grid2D(i, j); // bz = dz (m) going from (i_shoot,j_shoot) to (i,j)
			// distance between the surfaces
			const double dist2 = ((bx * bx) + (by * by) + (bz * bz));
			// the longwave radiation will only be emitted if the radius is lower than
			// a preassumed radius (lw_radius) -> changeable in: lw_t(a,b) * lw_radius
			if (dist2 <= lw_radius2 && dist2 > 0.)
			{

				// This call either
				// computes the ViewFactor (vf_in_ram=false) or
				// take the pre-computed ViewFactor (vf_in_ram=true)
				// This differentian according to vf_in_ram is done in ViewFactorsHelbig.cc
				const double viewFactor = viewFactorsHelbigObj.GetViewfactor(i, j, i_shoot, j_shoot);

				// calculation of the reflected amount of radiation of ab that arrives at ij
				// formulation of the long-wave radiation coming from the terrain:
				// Model runs of Sebastian Hoch with MODTRAN4 to check out at what distances are terrain patches
				// seen still important; it is on one hand depending on the surface temperature of a specific patch
				// and on the other hand on the surface properties, respectively the albedo (?).
				// afterwards empirical formulation of Michi Lehning using R
				// the formulation gives LW in W / (m2 * sr)
				// as the view factor is for emitting in the hemisphere (divided by PI in the view factor formula)
				// we here have to multiply by PI

				double rad = (0.000009886 * (pow(log(sqrt(dist2)), 1.07)) * (meteo2d_ta(i, j) - t_snowold(i_shoot, j_shoot)) + 0.000003456 * meteo2d_ta(i, j) + 0.0001452 * t_snowold(i_shoot, j_shoot) - 0.0304) * viewFactor * 10000. * M_PI;

				// the received amount is added to the total radiation at ij
				lwi(i, j) += rad;

				//available radiation to emit at the current cell
				const double emitted_i_j = viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) * lw_t(i, j);

				// the next ('shooting') grid cell with the most unshot radiation is selected
				// taking into account an air column reduction factor,
				// the actual (slope) area size, emittable longwave radiation
				// and the sum of total terrain view factor
				if (i != i_shoot || j != j_shoot)
				{
					if (emitted_i_j > diffmax_lw)
					{
						diffmax_lw = emitted_i_j;
						i_max_unshoot_lw = i;
						j_max_unshoot_lw = j;
					}
				}
				s++;
				// stopping criterion : Delta B^(k) * Sum_j(Fij) * A
				eps_stern += fabs(emitted_i_j);
			}
		}
	}

	// set the emitted 'shot' radiation at that cell to zero, but in contrast to the
	// shortwave radiation exchange taking into account multiple reflections here
	// every grid cell emitts only once
	lw_t(i_shoot, j_shoot) = 0;

	// check for stopping the iteration
	if (eps_stern <= threshold_itEps_LW || (n >= itMax_LW))
	{
		printf("[i] EBALANCE: LW converged after n=%u steps eps_stern <= %2.10lf with %u gathering patches\n", n, eps_stern, s);
		printf("    time for LW %f seconds\n", (double)((clock() - t0) / CLOCKS_PER_SEC));
		fflush(stdout);
		return 1;
	}
	return 0;
}

void TerrainRadiationHelbig::ComputeTerrainRadiation(const bool &day, int i_max_unshoot, int j_max_unshoot, int i_max_unshoot_lw, int j_max_unshoot_lw)
{
	// computes long wave and short wave terrain radiation using
	// a Progressive Refinement Radiosity (e.g. Gortler et al. (1994)
	// view factors are either stored or computed on the fly (according to vf_in_ram)
	unsigned int n;			   // counts iteration steps
	double threshold_itEps_SW; // stopping threshold for SW PR iteration
	double threshold_itEps_LW; // stopping threshold for LW PR iteration
	int itMax_LW;			   // maximum number of iteration steps for longwave emission
	int converged = 0;		   // 0 until iterated solution reaches a certain accuracy
	clock_t t0;

	//computing SW terrain radiation
	if (day == true)
	{		   //We compute SWR terrain radiation only during the day
		n = 0; // iteration step counter
		t0 = clock();

		//limits redundant terrain radiation computation to a maximum
		//reflectable direct + diffuse sky of 40 W/m2 and albedo of 0.99
		if (max_glob_start < 40.)
			converged = 1;
		else
			converged = 0;
		// itEps_SW: SW radiosity stopping tolerance
		// Two options:
		// epsilon = itEps_SW * |unshot radiosity|_1 sum / faktor y
		//threshold_itEps_SW = itEps_SW * (mean_glob_start * dimx * dimy) * viewFactorsHelbigObj.min_area
		//          * viewFactorsHelbigObj.min_vterr * (1. - max_alb) / max_alb;
		threshold_itEps_SW = itEps_SW * (mean_glob_start * dimx * dimy) * viewFactorsHelbigObj.min_area * viewFactorsHelbigObj.min_vterr * (1 - max_alb * viewFactorsHelbigObj.max_vterr) /
							 (max_alb * viewFactorsHelbigObj.min_vterr);

		mio::Timer timer_sw;
		timer_sw.start();
		//HACK: do computation on square, as for LW
		std::cout << "converged: " << converged << std::endl;
		while (converged != 1)
		{
			n++;
			converged = SWTerrainRadiationStep(threshold_itEps_SW, i_max_unshoot, j_max_unshoot, n, t0);
		}
		timer_sw.stop();
		std::cout << "calc sw radiation: " << timer_sw.getElapsed() << std::endl;
	}
	//computing LW terrain radiation
	t0 = clock();
	// itEps_LW: LW radiosity stopping tolerance
	// Two options:
	//threshold_itEps_LW = itEps_LW * lw_start_l1 * viewFactorsHelbigObj.min_area * viewFactorsHelbigObj.min_vterr;
	threshold_itEps_LW = itEps_LW * lw_start_l1 * viewFactorsHelbigObj.min_area * viewFactorsHelbigObj.min_vterr * (1 - viewFactorsHelbigObj.max_vterr) / (viewFactorsHelbigObj.max_vterr);

	// maximum number of iteration steps for longwave emission is limited to the number of
	// grid cells per model domain
	itMax_LW = dimx * dimy;

	//---> Block to uncomment to put back LW
	// converged = 0;
	// n = 0; // iteration step counter = shooting cell counter
	// while (converged != 1)
	// {
	// 	n++;
	// 	converged = LWTerrainRadiationStep(threshold_itEps_LW, itMax_LW, i_max_unshoot_lw, j_max_unshoot_lw, n, t0);
	// }
}

void TerrainRadiationHelbig::ComputeRadiationBalance()
{
	// This routine computes the distributed radiation balance
	// correction parameter for the distributed theoretical direct radiation
	bool day;																			  // switch for daytime / nighttime
	int i_max_unshoot = 0, j_max_unshoot = 0, i_max_unshoot_lw = 0, j_max_unshoot_lw = 0; // variables for the grid cell with the most unshot shortwave radiation
	//double solarAzimuth;
	//double toa_h, direct_h, diffuse;

	//calculate atmosphere parameters at the station
	if (tdir.getMax() + tdiff.getMax() > 0)
		day = true;
	else
		day = false;

	// maximum of reflectable direct and diffuse sky shortwave radiation
	max_glob_start = 0.;
	// mean of reflectable direct and diffuse sky shortwave radiation
	mean_glob_start = 0.;
	// emittable longwave start distribution with vector norm l1
	lw_start_l1 = 0.;

	// calculation of the reflected radiation from the terrain
	InitializeTerrainRadiation(day, i_max_unshoot, j_max_unshoot, i_max_unshoot_lw, j_max_unshoot_lw);
	ComputeTerrainRadiation(day, i_max_unshoot, j_max_unshoot, i_max_unshoot_lw, j_max_unshoot_lw);

	fillSWResultsGrids(day);
}

void TerrainRadiationHelbig::InitializeTerrainSwSplitting(const int i, const int j,
														  int &i_max_unshoot, int &j_max_unshoot, double &diffmax_sw)
{

	// Correct diffuse for sky VF
	tdiff(i, j) = tdiff(i, j) * viewFactorsHelbigObj.getSkyViewFactor(i, j); //no reprojection for diffuse rad

	// local *reflected* direct and diffuse sky radiation (later glob_h will include terrain reflected radiation)
	glob_h(i, j) = albedo_grid(i, j) * (tdir(i, j) + tdiff(i, j)); //HACK: despite its name, it is on the slope, not horizontal!

	// storing the start shortwave radiation distribution
	glob_start(i, j) = glob_h(i, j);

	// computing the mean reflectable start shortwave radiation distribution
	mean_glob_start += fabs(glob_start(i, j)) / (dimx * dimy);

	// the possibly reflected shortwave radiation at the beginning
	sw_t(i, j) = glob_h(i, j);

	// the grid cell with the most unshot radiation is selected taking into account albedo,
	// actual (slope) area size, reflected shortwave radiation and sum of total terrain view factor
	if (sw_t(i, j) * viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) > diffmax_sw)
	{
		diffmax_sw = sw_t(i, j) * viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j);
		i_max_unshoot = i;
		j_max_unshoot = j;
	}

	// lower bound for iteration start:
	// no reflected SW terrain radiation needs to be accounted for if not a patch exists
	// with max_glob_start > 40 Wm2 and a view factor sum of > 0.1 (prevents to take mountain peaks)
	if (sw_t(i, j) > max_glob_start && viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) > 0.1)
	{
		max_glob_start = sw_t(i, j);
	}
}

void TerrainRadiationHelbig::InitializeTerrainRadiation(const bool &day, int &i_max_unshoot, int &j_max_unshoot, int &i_max_unshoot_lw, int &j_max_unshoot_lw)
{
	double diffmax_sw = 0.; // reference product for detecting the grid cell with most unshot shortwave radiation
	//int horizon_x, horizon_y; // corresponding horizon coordinates in solar azimuth direction
	double diffmax_lw = 0; // reference product for detecting the grid cell with most unshot longwave radiation

	total_terrain.resize(dimx, dimy, 0);

	lw_eps_stern = 0.; //Stop criterion for the LW terrain radiation

	for (int i = 0; i < dimx; i++)
	{
		for (int j = 0; j < dimy; j++)
		{
			total_terrain(i, j) = 0;
			if (dem.grid2D(i, j) == mio::IOUtils::nodata)
				continue;

			if (day == true)
			{ //then, compute short wave radiation
				//set initial values for each cell for SW diff (view caftor correction)
				InitializeTerrainSwSplitting(i, j, i_max_unshoot, j_max_unshoot, diffmax_sw);
			}
			//compute long wave radiation
			// saturation vapor pressure in Pa (routine in snowpack)
			//const double e_stern = mio::Atmosphere::vaporSaturationPressure( meteo2d_ta(i,j) );
			// wvp = water vapor pressure
			//const double wvp = meteo2d_rh(i,j) * e_stern;
			//Long Wave initialization
			InitializeLW(i, j, i_max_unshoot_lw, j_max_unshoot_lw, diffmax_lw);
		}
	}
	//Optimisation by GS #6 : Sort array by emitted LW radiation
	//if (TERRAIN_RADIATION!=0)
	//qsort(lwt_byCell, dimx*dimy, sizeof(CellsList), CellsRadComparator_Helbig );
	sort(lwt_byCell.begin(), lwt_byCell.end(), operator_greater);
}

//Function to calculate the min/max index around a given cell, knowing a given radius in cell units
void TerrainRadiationHelbig::CalculateIndex(const int indice, const int distance_max, int dim, int *min, int *max)
{
	int i_max = indice + distance_max;
	int i_min = indice - distance_max;

	if (i_max >= dim)
	{
		*max = dim - 1;
	}
	else
	{
		*max = i_max;
	}
	if (i_min < 0)
	{
		*min = 0;
	}
	else
	{
		*min = i_min;
	}
}

//Function to calculate radiation transmit by base cell to other cell and add them in the radiation received by the other cell
//Calculate also gathering patches
void TerrainRadiationHelbig::LWTerrainRadiationCore(const double bx2, const int j_shoot, const double z_shoot, const int j, const double z, const double cellsize, const double t_snow_shoot, const double t_snow_shoot_value, const double t_a, const double vf, double *lwi, int *s)
{
	// formulation of the long-wave radiation coming from the terrain:
	// Model runs of Sebastian Hoch with MODTRAN4 to check out at what distances are terrain patches
	// seen still important; it is on one hand depending on the surface temperature of a specific patch
	// and on the other hand on the surface properties, respectively the albedo (?).
	// afterwards empirical formulation of Michi Lehning using R
	// the formulation gives LW in W / (m2 * sr)
	// as the view factor is for emitting in the hemisphere (divided by PI in the view factor formula)
	// we here have to multiply by PI

	//Optimisation #4 by GS
	//Calculation of radiation only if vf != 0
	if (vf != 0.0)
	{
		const double by = (double)(j_shoot - j) * cellsize; // by = dy (m) going from (i_shoot,j_shoot) to (i,j)
		const double bz = z_shoot - z;						// bz = dz (m) going from (i_shoot,j_shoot) to (i,j)
		const double dist2 = (bx2 + (by * by) + (bz * bz));

		/*
		const double rad = ( 0.000009886 * (pow( 0.5 * log(dist2),1.07 )) * (t_a - t_snow_shoot)
		+ 0.000003456 * t_a + 0.0001452 * t_snow_shoot - 0.0304 )
		* vf * 10000. * PI;
		*/

		//Optimisation #3.3
		//Precalculation of some parts
		const double rad = vf * (4.7088896243311818e-06 * pow(log(dist2), 1.07) * (t_a - t_snow_shoot) + 0.000003456 * t_a + t_snow_shoot_value) * M_PI * 10000;

		// the received amount is added to the total radiation at ij
		(*lwi) += rad;
		(*s)++;
	}
}

void TerrainRadiationHelbig::fillSWResultsGrids(const bool &day)
{
	if (day == true)
	{
		//double sumrad_o = 0; // total outgoing radiation to compute effective albedos
		//double sumrad_i = 0; // total incoming radiation to compute effective albedos
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				if (dem.grid2D(i, j) == mio::IOUtils::nodata)
					continue;
				// Compute the effective albedo of the model domain
				//sumrad_i += glob_start(i,j) / meteo2D(i,j).alb;
				//sumrad_o += glob_h(i,j) * sky_vf(i,j);

				// Fill in the return values for SNOWPACK or for the ebalance OUTPUT
				// the INCIDENT global shortwave radiation:
				if (albedo_grid(i, j) > 0)
					glob_h(i, j) = glob_h(i, j) / albedo_grid(i, j); //HACK: this is BAD!
				// TOTAL incident shortwave diffuse radiation
				if (albedo_grid(i, j) > 0)
					total_diff(i, j) = tdiff(i, j) + glob_h(i, j) - (glob_start(i, j) / albedo_grid(i, j));
			}
		}
	}
	else
	{
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				if (dem.grid2D(i, j) == mio::IOUtils::nodata)
					continue;
				glob_h(i, j) = 0.;
				total_diff(i, j) = 0.;
			}
		}
	}
}

void TerrainRadiationHelbig::InitializeLW(const int i, const int j, int &i_max_unshoot_lw, int &j_max_unshoot_lw, double &diffmax_lw)
{
	// longwave sky radiation calculation
	lw_sky(i, j) = tot_ilwr(i, j) * viewFactorsHelbigObj.getSkyViewFactor(i, j);
	lwi(i, j) = lw_sky(i, j);

	// longwave radiation of ij emittable to the surrounding terrain without attenuation and view factor
	const double tss = t_snowold(i, j);
	lw_t(i, j) = Cst::stefan_boltzmann * (tss * tss * tss * tss);

	// computing the longwave emission start distribution with l1 vector norm
	lw_start_l1 += fabs(lw_t(i, j));

	// the grid cell with the most unshot radiation is selected taking into account an air column reduction factor,
	// the actual (slope) area size, emittable longwave radiation and sum of total terrain view factor
	const double rad_t = viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) * lw_t(i, j);

	// the grid cell with the most unshot radiation is selected taking into account an air column reduction factor,
	// the actual (slope) area size, emittable longwave radiation and sum of total terrain view factor
	if (lw_t(i, j) * viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j) > diffmax_lw)
	{
		diffmax_lw = lw_t(i, j) * viewFactorsHelbigObj.getSymetricTerrainViewFactor(i, j);
		i_max_unshoot_lw = i;
		j_max_unshoot_lw = j;
	}
	//Optmisation #6 by GS : put the current element in the table of LW TerrainRadiation
	//if (TERRAIN_RADIATION!=0) {
	CellsList &lwtrbc = lwt_byCell[i * dimy + j]; //temporary reference to the current element's location
	lwtrbc.radiation = rad_t;
	lwtrbc.x = i;
	lwtrbc.y = j;
}

#pragma GCC diagnostic pop
