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
#include <alpine3d/ebalance/TerrainRadiation.h>

#include <algorithm>

using namespace std;
using namespace mio;

TerrainRadiation::TerrainRadiation(const mio::Config& cfg, const mio::DEMObject& dem_in, const int&, const std::string& method)
     : TerrainRadiationAlgorithm(method), viewFactorsObj(cfg, dem_in), dem(dem_in), dimx(dem.getNx()), dimy(dem.getNy()),
       cellsize(dem.cellsize), albedo_grid(), meteo2d_ilwr(), meteo2d_ta(), meteo2d_rh(), total_diff(),
       sw_t(), glob_start(), glob_h_isovf(), glob_h(), t_snowold(), lw_t(), lwi(), lw_sky(), lwt_byCell(dimx*dimy),
       tdir(), tdiff(), Mmy(), Mt()
{
	double lw_radius;
	cfg.getValue("itEps_SW", "EBalance", itEps_SW );
	cfg.getValue("itEps_LW", "EBalance", itEps_LW );
	cfg.getValue("itEps1_SW", "EBalance", itEps1_SW);
	cfg.getValue("sw_radius", "EBalance", sw_radius);
	cfg.getValue("lw_radius", "EBalance", lw_radius);

	total_diff.resize( dimx, dimy, 0);
	tdir.resize( dimx, dimy, 0);
	tdiff.resize( dimx, dimy, 0);
	sw_t.resize( dimx, dimy,0);
	glob_start.resize( dimx, dimy, 0);
	glob_h_isovf.resize( dimx, dimy, 0);
	glob_h.resize( dimx, dimy, 0);
	t_snowold.resize( dimx, dimy, 0);

	lw_t.resize( dimx, dimy, 0);
	lwi.resize( dimx, dimy, 0 );
	lw_sky.resize( dimx, dimy, 0);

	LW_distance_index = (int)ceil(lw_radius / cellsize);

	startx = 0;
	nx = dimx;
	starty = 0;
	ny = dimy;

	if (MPIControl::instance().size() > 1) {
		MPIControl::instance().getArraySliceParams(dimx, startx, nx);

		const size_t endx = startx + nx;
		const size_t endy = starty + ny;
		viewFactorsObj.setBoundarys(startx, endx, starty, endy);
		Mmy.resize( dimx, dimy, 0);
		Mt.resize( dimx, dimy, 1);

		for (size_t ii=startx; ii<startx+nx; ii++) {
			for (size_t jj=starty; jj<starty+ny; jj++) {
				Mmy(ii,jj) = 1;
				Mt(ii,jj) = 0;
			}
		}
	} else {
		Mmy.resize( dimx, dimy, 1);
		Mt.resize( dimx, dimy, 0);
	}

	std::cout << "[TR] startx=" << startx << " nx=" << nx << std::endl;
	std::cout << "[TR] starty=" << starty << " ny=" << ny << std::endl;

	viewFactorsObj.generate();
}

void TerrainRadiation::getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain)
{
	if (MPIControl::instance().master())
		printf("[i] Calculating terrain radiation with standard method\n");

	tdir = Mmy*direct;
	tdiff = Mmy*diffuse;
	if (max_alb > 0) Compute();

	MPIControl::instance().allreduce_sum(total_diff);
	terrain = total_diff;
}

void TerrainRadiation::setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta,
                                const mio::Array2D<double>& rh, const mio::Array2D<double>& ilwr)
{
        meteo2d_ta = ta;
        meteo2d_rh = rh;
        meteo2d_ilwr = ilwr;
        albedo_grid = albedo;
        max_alb = albedo_grid.getMax();
}

void TerrainRadiation::Compute()
{	//we take sn_MdataT values as reference radiation values for all cells
	//iswr_ref = sn_MdataT.iswr;
	//ea_ref = sn_MdataT.ea;

	//distributed radiation calculation
	ComputeProgressiveRefinement();

	//sn_MdataT.elev = solarElevation; //needed for Canopy
	// Fill shortwave, longwave, total diffuse shortwave and solar elevation into
	// corresponding snowpack data structure
	//if ( snowpack != NULL ) {
		//snowpack->setRadiationComponents ( glob_h, lwi, total_diff, solarElevation, timestamp);
	//}
}

int TerrainRadiation::SWTerrainRadiationStep(const double i_threshold_itEps_SW, size_t& c, size_t& d,
                                             mio::Array2D<double>& sw_Q_out, double& Rmy)
{
	// Computation of shortwave terrain radiation
	// At every iteration step, a reference grid cell reflects ('shoots') radiation to every other grid cell
	// The coordinates of the next cell with the most unshoot radiation is written in (*c,*d)
	//std::cout << "c=" << *c << " d=" << *d << std::endl;
	const size_t i_shoot = c, j_shoot = d;	//variables for the grid cell with the most unshot radiation
	double diffmax_sw=0.;		// reference product for detecting the grid cell with most unshot shortwave radiation
	const double diffmax_thres=0.;	// threshold for when to stop looking for maximum cells
	double eps_stern=0.;		// stopping criterion
	double be=0.;			// matrix difference B^(k) - E resp. (L_sky + L_terrain) - L_sky
	unsigned int s=0;			// counts gathering patches, i.e. those patches within limited distance radius

	//variable needed to optimize the computational speed
	const double sw_radius2 = sw_radius * sw_radius;
	const double z_shoot=dem.grid2D(i_shoot,j_shoot);
	const double sw_shoot = sw_t(i_shoot,j_shoot);

	//if (!(i_shoot >= startx && i_shoot <startx+nx && j_shoot >= starty && j_shoot <starty+ny))
	//	return 1;

	// every grid cell ij receives radiation from the chosen one with coordinates ab
	for ( size_t j = 0; j < dimy; j++ ) {
		for ( size_t i = 0; i < dimx; i++ ) {
			if (albedo_grid(i,j) == mio::IOUtils::nodata) continue;
			if (dem.grid2D(i,j) == mio::IOUtils::nodata) continue;

			const double bx = (i_shoot - i) * cellsize;	// bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
			const double by = (j_shoot - j) * cellsize;	// by = dy (m) going from (i_shoot,j_shoot) to (i,j)
			const double bz = z_shoot - dem.grid2D(i,j);	// bz = dz (m) going from (i_shoot,j_shoot) to (i,j)

			// distance between the surfaces
			const double dist2 = ( (bx * bx) + (by * by) + (bz * bz) );


			if (i >= startx && i <startx+nx && j >= starty && j <starty+ny) {
				//available radiation to reflect at the current cell
				const double reflected_i_j = viewFactorsObj.getSymetricTerrainViewFactor(i,j) * sw_t(i,j);


				// the next ('shooting') grid cell with the most unshot radiation is selected
				// taking into account albedo, actual (slope) area size, reflected shortwave
				// radiation and the sum of total terrain view factor

				//if (i >= startx && i <startx+nx && j >= starty && j <starty+ny)
				if ( reflected_i_j  > diffmax_sw ) {
					if (!(i==i_shoot && j==j_shoot)) {
						diffmax_sw = reflected_i_j;
						c = i;
						d = j;
					}
				}

				// the first stopping criterion:
				// |Delta B^(k) * Sum_j(Fij) * A|_1
				eps_stern += fabs( reflected_i_j );
			}

			// the shortwave radiation will only be emitted if the radius is lower than
			// a preassumed radius (sw_radius) -> changeable in: sw_t(a,b) * sw_radius
			//please note that it also excludes the current shooting cell itself!


			if ( dist2 <= sw_radius2 && dist2>0.) {

				//const double viewFactor = viewFactorsObj.GetViewfactor(i, j, i_shoot, j_shoot);
				const double viewFactor = viewFactorsObj.GetViewfactor(i_shoot, j_shoot, i, j);

				// calculation of the reflected amount of radiation of ab that arrives at ij
				const double rad = viewFactor * albedo_grid(i,j) * sw_shoot;

				//std::cout << "rad=" << rad << " viewFactor=" << viewFactor << " albedo_grid(i,j)=" << albedo_grid(i,j)<< " sw_shoot=" << sw_shoot<< std::endl;

				// the received amount is added to the reflectable amount of radiation (unshot)
				sw_t(i,j) += rad;

				// in addition the received amount is added to the total radiation at ij
				sw_Q_out(i,j) += rad;

				//std::cout << "rad=" << rad << std::endl;

				s++;
			} // end of if only for distances lower than sw_radius

			// the second stopping criterion:
			// |B^(k) - E|_1
			be += fabs(sw_Q_out(i,j) - glob_start(i,j));

		} // end of for i loop
	} // end of for j loop

	// set the reflected 'shot' radiation at that cell to zero
	sw_t(i_shoot,j_shoot) = 0.;

	// check for stopping the iteration
	bool stop_criteria = false;
	if (MPIControl::instance().size() > 1) {
		//std::cout << "c=" << *c << " d=" << *d << std::endl;
		//std::cout << "i_threshold_itEps_SW=" << i_threshold_itEps_SW << " eps_stern=" << eps_stern << std::endl;
		//std::cout << sw_Q_out.toString() << std::endl;
		stop_criteria = (eps_stern <= i_threshold_itEps_SW);
		//const bool stop_criteria = (diffmax_sw < 0.0001);
		//const bool stop_criteria = true;
	} else {
		stop_criteria =((eps_stern <= i_threshold_itEps_SW) || be >= (itEps1_SW *
		               ( mean_glob_start * dimx * dimy )) || (diffmax_sw<=diffmax_thres));
	}

	if (stop_criteria) {
		Rmy = eps_stern;
		fflush( stdout );
		return 1;
	}

	return 0;
}


int TerrainRadiation::LWTerrainRadiationStep(const double threshold_itEps_LW, const int itMax_LW,
                                             const int i_shoot, const int j_shoot, unsigned int n)
{
// Computation of longwave terrain radiation. NOTE: This function takes up a lot of computational time globally!
// At every iteration step, one reference grid cell (*e,*f) reflects ('shoots') radiation to every other grid cell
// within a given distance (LW_distance_index in cells units). Every grid cell is only once the emitting cell (no multiple reflections)

//Optimizations by GS (see report for more details)
// inlines functions have been created to calculate radiation, distance and next shooting cell
// each shooting cell is only processed within a square around it (within LW_distance_index)
// #6.4 calculate the distance only when the view-factor is different from 0
// #6.5 use pointers to access Array2D elements
// We now use an array to find the next shooting cell

	const double diffmax_thresh=0.;		//when do we consider that the contributions are too small -> stop this iteration?
	//const unsigned int ab = i_shoot + j_shoot*dimx;	//1D addressing of the shooting cell for the view factors matrix
	int s = 0;				//counts gathering patches, i.e. those patches within limited distance radius

	// variables needed to optimize the computational speed
	const double z_shoot=dem.grid2D(i_shoot,j_shoot);
	const double t_snow_shoot=t_snowold(i_shoot,j_shoot);
	const double t_snow_shoot_value=0.0001452 * t_snow_shoot- 0.0304;

	//Optimisation #6 by GS : Calculation of max/min indice
	int distance_min_x, distance_max_x, distance_min_y, distance_max_y;
	CalculateIndex(i_shoot, LW_distance_index, dimx, &distance_min_x, &distance_max_x);
	CalculateIndex(j_shoot, LW_distance_index, dimy, &distance_min_y, &distance_max_y);

	//Optimisation #7 : Calculations to unroll loop
	const int rest_itr = (distance_max_y - distance_min_y) % NB_UNROLL; //leftover iterations after unrolling
	const int nb_itr = (distance_max_y - distance_min_y) - rest_itr; //iterations processed in the unrolling

	if ( viewFactorsObj.vf_in_ram ) {	//view factors are stored, we use them directly
		//Optimisation #3.1 by GS : precompute the invariant of the two loops
		int i = distance_min_x;
		for ( ; i <= distance_max_x; i++ ) {
			//Optimisation #3.1 by GS: precompute the invariant of the loop
			const double bx = (double)(i_shoot - i) * cellsize; // bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
			const double bx2 = bx*bx;

			//Optimization #6.5 by GS: Create pointers to table elements
			//WARNING: this assumes that the element are stored in consecutive locations
			double *pt_ta = (double *)(&meteo2d_ta(i, distance_min_y));
			double *pt_lwi = (double *)(&lwi(i,distance_min_y));
			double *pt_z = (double *)(&dem.grid2D(i,distance_min_y));

			//Optimisation #7.1 by GS: unrolling loops 8x
			int j = distance_min_y;
			for ( ; j < nb_itr; j+=3) {
				//const unsigned int ij_0 = j * dimx + i;
				//const unsigned int ij_1 = (j+1) * dimx + i;
				//const unsigned int ij_2 = (j+2) * dimx + i;

				const double sx_0 = dem.Nx(i,j);
				const double sx_1 = dem.Nx(i,j+1);
				const double sx_2 = dem.Nx(i,j+2);

				const double sy_0 = dem.Ny(i,j);
				const double sy_1 = dem.Ny(i,j+1);
				const double sy_2 = dem.Ny(i,j+2);

				double vf_0 = viewFactorsObj.GetViewfactor(i,j,i_shoot,j_shoot);
				double vf_1 = viewFactorsObj.GetViewfactor(i,j+1,i_shoot,j_shoot);
				double vf_2 = viewFactorsObj.GetViewfactor(i,j+1,i_shoot,j_shoot);

				//view factor is symetric part (from vf() ) divided by the "area" of the current cell
				vf_0 /= ( sqrt( sx_0*sx_0 + sy_0*sy_0 + 1. ) * cellsize * cellsize );
				vf_1 /= ( sqrt( sx_1*sx_1 + sy_1*sy_1 + 1. ) * cellsize * cellsize );
				vf_2 /= ( sqrt( sx_2*sx_2 + sy_2*sy_2 + 1. ) * cellsize * cellsize );

				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, vf_0, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j+1, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, vf_1, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j+2, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, vf_2, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;

			} // end of for j unrolled loop

			//Optimisation #7.1 by GS : Do the rest of unrolled loop
			for ( ; j <= distance_max_y; j++ ) {
				/*const int ij = j * dimx + i;*/

				//Optimisation #3.2 by GS
				//Use a variable to replace vf[0][0]
				//const double tmp_vf = vf(ij,ab);
				const double tmp_vf = viewFactorsObj.GetViewfactor(i, j, i_shoot, j_shoot);
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, tmp_vf, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;

			} // end of for j loop
		} // end of for i loop
	} else {		//view factors are not stored, they are calculated on the fly
		int i = distance_min_x;
		for ( ; i < distance_max_x; i++ ) {
			//Optimisation #3.1 by GS
			//Output the invariant of the loop
			const double bx = (double)(i_shoot - i) * cellsize;           // bx = dx (m) going from (i_shoot,j_shoot) to (i,j)
			const double bx2 = bx*bx;

			//Create pointers to O #6.5
			double * pt_ta = (double *)(&meteo2d_ta(i,distance_min_y));
			double * pt_lwi = (double *)(&lwi(i,distance_min_y));
			double * pt_z = (double *)(&dem.grid2D(i,distance_min_y));

			int j = distance_min_y;
			for ( ; j < nb_itr; j+=NB_UNROLL ) {
				double viewFactor;
				/*vf(0,0) = 0.; //our value will be stored in vf[0][0]
				ComputeHorizon (solar_elev, i, j, i_shoot, j_shoot);*/
				viewFactor = viewFactorsObj.GetViewfactor(i, j, i_shoot, j_shoot);
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, viewFactor,&(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;

				/*vf(0,0) = 0.;
				ComputeHorizon (solar_elev, i, j+1, i_shoot, j_shoot);*/
				viewFactor = viewFactorsObj.GetViewfactor(i, j+1, i_shoot, j_shoot);
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j+1, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, viewFactor, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;

				/*vf(0,0) = 0.;
				ComputeHorizon (solar_elev, i, j+2, i_shoot, j_shoot);*/
				viewFactor = viewFactorsObj.GetViewfactor(i, j+2, i_shoot, j_shoot);
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j+2, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, viewFactor, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;

			} // end of for j unrolled loop

			//Optimisation #7.1 by GS : Do the rest of unrolled loop
			for ( ; j <= distance_max_y; j++ ) {
				/*vf(0,0) = 0.;
				ComputeHorizon (solar_elev, i, j, i_shoot, j_shoot);*/
				const double viewFactor = viewFactorsObj.GetViewfactor(i, j, i_shoot, j_shoot);
				LWTerrainRadiationCore(bx2, j_shoot, z_shoot, j, *pt_z, cellsize, t_snow_shoot, t_snow_shoot_value, *pt_ta, viewFactor, &(*pt_lwi), &s);
				pt_z++;
				pt_ta++;
				pt_lwi++;
			} // end of for j loop
		} // end of for i loop
	} //end of stored/unstored view factors test

	//compute emitted radiation that will be used for the convergence criteria
	const double local_radiation = viewFactorsObj.getSymetricTerrainViewFactor(i_shoot, j_shoot) * lw_t(i_shoot, j_shoot);

	// set the emitted 'shot' radiation at that cell to zero, but in contrast to the
	// shortwave radiation exchange taking into account multiple reflections here
	// every grid cell emitts only once
	lw_t(i_shoot,j_shoot) = 0.;

	// check for stopping the iteration: physical threshold OR too many iterations OR no more cell that can emit left
	if ( lw_eps_stern <= threshold_itEps_LW || ((int)n >= itMax_LW) || local_radiation<=diffmax_thresh) {
		printf( "day with PR: LW converged after n %u steps error %2.10f with %u gathering patches\n",
				n, lw_eps_stern, s);
		//printf(" time for LW %f seconds\n", (double)((clock() - t0) / CLOCKS_PER_SEC) );
		fflush( stdout );
		return 1;
	}
	return 0;
}

int TerrainRadiation::ComputeTerrainRadiation(const double& i_threshold_itEps_SW, const bool& day, size_t& c,
                                              size_t& d, mio::Array2D<double>& sw_Q_out, double& Rmy)
{
	// computes long wave and short wave terrain radiation using
	// a Progressive Refinement Radiosity (e.g. Gortler et al. (1994)
	// view factors are either stored or computed on the fly (according to vf_in_ram)
	unsigned int n;				// counts iteration steps
	//double i_threshold_itEps_SW;	// stopping threshold for SW PR iteration
	//double threshold_itEps_LW;	// stopping threshold for LW PR iteration
	//int itMax_LW;			// maximum number of iteration steps for longwave emission
	int converged = 0;			// 0 until iterated solution reaches a certain accuracy

	//computing SW terrain radiation
	if (day==true && max_alb > 0) { //We compute SWR terrain radiation only during the day
		n = 0;	// iteration step counter

		//limits redundant terrain radiation computation to a maximum
		//reflectable direct + diffuse sky of 40 W/m2 and albedo of 0.99
		if ( max_glob_start < 40. )
			converged = 1;
		else
			converged = 0;

		mio::Timer timer_sw;
		timer_sw.start();
		//HACK: do computation on square, as for LW
		while (converged != 1) {
			n++;
			converged = SWTerrainRadiationStep(i_threshold_itEps_SW, c, d, sw_Q_out, Rmy);
		}

		timer_sw.stop();
	}

	//computing LW terrain radiation
	// itEps_LW: LW radiosity stopping tolerance
	/*threshold_itEps_LW = itEps_LW * lw_start_l1 * viewFactorsObj.min_area * viewFactorsObj.min_vterr;

	// maximum number of iteration steps for longwave emission is limited to the number of
	// grid cells per model domain
	itMax_LW = dimx * dimy;

	converged = 0;
	n = 0;	// iteration step counter = shooting cell counter
	while (converged != 1) {
		n++;
		converged = LWTerrainRadiationStep(threshold_itEps_LW, itMax_LW, lwt_byCell[n].x, lwt_byCell[n].y, n);
		lw_eps_stern -= fabs(lwt_byCell[n-1].radiation); //updating convergence criteria
	}*/

	return EXIT_SUCCESS;
} // end of ComputeTerrainRadiation


void TerrainRadiation::ComputeProgressiveRefinement()
{
	// This routine computes the distributed radiation balance
	//double corr_dir;		// correction parameter for the distributed theoretical direct radiation
	bool day;			     // switch for daytime / nighttime
	size_t c = 0, d = 0;	// variables for the grid cell with the most unshot shortwave radiation
	//double solarAzimuth;
	//double toa_h, direct_h, diffuse;

	//calculate atmosphere parameters at the station
	if ( tdir.getMax() + tdiff.getMax() > 0 )
		day=true;
	else
		day=false;

	// in case of usage of cloud cover as input the emissivity has to be converted in eighth
	//if (CLOUD==1) ea_ref = ea_ref/8.;

	// maximum of reflectable direct and diffuse sky shortwave radiation
	max_glob_start = 0.;
	// mean of reflectable direct and diffuse sky shortwave radiation
	mean_glob_start = 0.;
	// emittable longwave start distribution with vector norm l1
	lw_start_l1 = 0.;

	// calculation of the diffuse radiation from the terrain
	InitializeTerrainRadiation(day);

	//glob_h Qout

	/**
		Example Qout decomposition:

		Qout = Qt + Qmy

               /           |            \           /           |            \
               |     0     |            |           |    Q^     |            |
               |           |            |           |           |            |
          Qt = |----------              |      Qmy= |----------              |
               |                        |           |                        |
               |                Q*      |           |                  0     |
               \                        /           \                        /

		LOOP:

			Qmy input for shortwave calculation per node
			shortwave progressive refinement calculation on whole domain
			Qt reduced between nodes


	*/

	mio::Array2D<double> Qout=glob_h;
	//glob_h = 0;

	MPIControl::instance().allreduce_sum(Qout);

	double Rmy = 0;
	const int numNodes = MPIControl::instance().size();
	double Rt = (double)threshold_itEps_SW*(double)(numNodes-1)/(double)numNodes;

	bool converged = false;
	int i = 0;
	while (!converged) {
		mio::Array2D<double> Qmy = Qout * Mmy;

		if (day) InitializeTerrainSwSplitting(c, d, Qout);

		//thoems: std::cout << "[TR] starting local PR with e*=" << threshold_itEps_SW << " Rt=" << Rt << std::endl;

		mio::Array2D<double> Qout_s(dimx, dimy, 0);

		ComputeTerrainRadiation(threshold_itEps_SW-Rt,  day, c, d, Qout_s, Rmy);

		mio::Array2D<double> Qmy_s = Qout_s * Mmy;
		mio::Array2D<double> Qt_s = Qout_s * Mt;

		glob_h += Qout_s;

		MPIControl::instance().allreduce_sum(Qt_s);

		Qout = Qt_s + sw_t * Mmy;// not nice bro!

		double R = Rmy;

		MPIControl::instance().allreduce_sum(R);

		Rt = R - Rmy;

		for (unsigned int ii=startx; ii<startx+nx; ii++)
			for (unsigned int jj=starty; jj<starty+ny; jj++)
				R += abs(viewFactorsObj.getSymetricTerrainViewFactor(ii,jj) * Qt_s(ii,jj));

		//printf("[TR] GLOBAL ITERATION: e*=%d Max(Qout)=%d\n", R, Qout.getMax());
		converged = (R <= threshold_itEps_SW);
		//converged = (Qout.getMax() < 0.000001);

		if (MPIControl::instance().size() > 1) {
			double conv = 0;
			if (!converged) conv++;
			MPIControl::instance().allreduce_sum(conv);
			converged = !(conv > 0);
		}
		i++;
	}

	fillSWResultsGrids(day);
}

void TerrainRadiation::InitializeTerrainSwSplitting(size_t& i_max_unshoot, size_t& j_max_unshoot, mio::Array2D<double>& sw_Q_out)
{
	double diffmax_sw = 0.; // reference product for detecting the grid cell with most unshot shortwave radiation

	for ( size_t i = 0; i < dimx; i++ ) {
		for ( size_t j = 0; j < dimy; j++ ) {
			if (dem.grid2D(i,j)==mio::IOUtils::nodata) continue;
			if (albedo_grid(i,j)==mio::IOUtils::nodata) continue;

			// the possibly reflected shortwave radiation at the beginning
			sw_t(i,j) = sw_Q_out(i,j);

			// the grid cell with the most unshot radiation is selected taking into account albedo,
			// actual (slope) area size, reflected shortwave radiation and sum of total terrain view factor
			if ( sw_t(i,j) * viewFactorsObj.getSymetricTerrainViewFactor(i,j)  > diffmax_sw ) {
				diffmax_sw = sw_t(i,j) * viewFactorsObj.getSymetricTerrainViewFactor(i,j);
				i_max_unshoot = i;
				j_max_unshoot = j;
			}
		}
	}
}

void TerrainRadiation::InitializeTerrainRadiation(const bool& day)
{
	//int horizon_x, horizon_y; // corresponding horizon coordinates in solar azimuth direction
	lw_eps_stern = 0.; //Stop criterion for the LW terrain radiation

	for ( size_t i = 0; i < dimx; i++ ) {
		for ( size_t j = 0; j < dimy; j++ ) {
			if (dem.grid2D(i,j)==mio::IOUtils::nodata) continue;
			if (albedo_grid(i,j)==mio::IOUtils::nodata) continue;

			if (day==true) { //then, compute short wave radiation
				//set initial values for each cell for SW dir and diff

				//calculating incident direct swr for the grid point
				double& direct=tdir(i,j); //currently, horizontal, we reproject to the slope
				double& diffuse=tdiff(i,j);

				diffuse = diffuse*viewFactorsObj.getSkyViewFactor(i,j); //no reprojection for diffuse rad

				// local *reflected* direct and diffuse sky radiation (later glob_h will include terrain reflected radiation)
				glob_h(i,j) = albedo_grid(i,j) * (direct + diffuse); //HACK: despite its name, it is on the slope, not horizontal!

				// storing the start shortwave radiation distribution
				glob_start(i,j) = glob_h(i,j);

				// computing the mean reflectable start shortwave radiation distribution
				mean_glob_start += fabs( glob_start(i,j) )/ (dimx * dimy);

				// lower bound for iteration start:
				// no reflected SW terrain radiation needs to be accounted for if not a patch exists
				// with max_glob_start > 40 Wm2 and a view factor sum of > 0.1 (prevents to take mountain peaks)
				//std::cout << sw_t(i,j) << " " <<  max_glob_start << " " << viewFactorsNoraObj.getSymetricTerrainViewFactor(i,j) << std::endl;

				if ( glob_h(i,j) > max_glob_start && viewFactorsObj.getSymetricTerrainViewFactor(i,j) > 0.1 ) {
					max_glob_start = glob_h(i,j);
				}
			}

			//std::cout << "SW_T MAX=" << sw_t.getMax() << std::endl;
			//compute long wave radiation
			// saturation vapor pressure in Pa (routine in snowpack)
			//const double e_stern = mio::Atmosphere::vaporSaturationPressure( meteo2d_ta(i,j) );
			// wvp = water vapor pressure
			//const double wvp = meteo2d_rh(i,j) * e_stern;
			//Long Wave initialization
			InitializeLW(i, j, lw_eps_stern);
		}
	}

	MPIControl::instance().allreduce_sum(mean_glob_start);
	MPIControl::instance().allreduce_max(max_glob_start);

	double min_area = viewFactorsObj.min_area;
	double min_vterr = viewFactorsObj.min_vterr;

	MPIControl::instance().allreduce_min(min_area);
	MPIControl::instance().allreduce_min(min_vterr);

	// itEps_SW: SW radiosity stopping tolerance
	// epsilon = itEps_SW * |unshot radiosity|_1 sum / faktor y
	threshold_itEps_SW = itEps_SW * (mean_glob_start * dimx * dimy) * min_area	* min_vterr * (1. - max_alb) / max_alb;

	//Optimisation by GS #6 : Sort array by emitted LW radiation
	sort(lwt_byCell.begin(), lwt_byCell.end(), operator_greater);
}

//Function to calculate the min/max index around a given cell, knowing a given radius in cell units
void TerrainRadiation::CalculateIndex(const int indice, const int distance_max, int dim, int * min, int * max)
{
	int i_max = indice+distance_max;
	int i_min = indice-distance_max;

	if (i_max >= dim) {
		*max = dim-1;
	} else {
		*max = i_max;
	}
	if (i_min < 0) {
		*min = 0;
	} else {
		*min = i_min;
	}
}

//Function to calculate radiation transmit by base cell to other cell and add them in the radiation received by the other cell
//Calculate also gathering patches
void TerrainRadiation::LWTerrainRadiationCore(const double bx2, const int j_shoot, const double z_shoot, const int j, const double z, const double cellsize, const double t_snow_shoot,const double t_snow_shoot_value ,const double t_a, const double vf, double * lwi, int * s)
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
	if (vf != 0.0){
		const double by = (double)(j_shoot - j) * cellsize;	// by = dy (m) going from (i_shoot,j_shoot) to (i,j)
		const double bz = z_shoot - z;				// bz = dz (m) going from (i_shoot,j_shoot) to (i,j)
		const double dist2 = ( bx2 + (by * by) + (bz * bz) );

		/*
		const double rad = ( 0.000009886 * (pow( 0.5 * log(dist2),1.07 )) * (t_a - t_snow_shoot)
		+ 0.000003456 * t_a + 0.0001452 * t_snow_shoot - 0.0304 )
		* vf * 10000. * PI;
		*/

		//Optimisation #3.3
		//Precalculation of some parts
		const double rad = vf * ( 4.7088896243311818e-06*pow(log(dist2),1.07) * (t_a - t_snow_shoot)
		+ 0.000003456 * t_a + t_snow_shoot_value ) * M_PI * 10000;

		// the received amount is added to the total radiation at ij
		(*lwi) += rad;
		(*s)++;
	}
}

void TerrainRadiation::fillSWResultsGrids(const bool& day) {
	if (day==true) {

		//std::cout << "glob_h:" << "\n" << glob_h.toString() << "\n";
		//std::cout << "glob_start:" << "\n" << glob_start.toString() << "\n";
		//double sumrad_o = 0; // total outgoing radiation to compute effective albedos
		//double sumrad_i = 0; // total incoming radiation to compute effective albedos
		for ( size_t i = 0; i < dimx; i++ ) {
			for ( size_t j = 0; j < dimy; j++ ) {
				if (dem.grid2D(i,j)==mio::IOUtils::nodata) continue;
				// Compute the effective albedo of the model domain
				//sumrad_i += glob_start(i,j) / meteo2D(i,j).alb;
				//sumrad_o += glob_h(i,j) * sky_vf(i,j);

				// Fill in the return values for SNOWPACK or for the ebalance OUTPUT
				// the INCIDENT global shortwave radiation:
				if (albedo_grid(i,j) > 0)
                                        glob_h(i,j) = glob_h(i,j) / albedo_grid(i,j); //HACK: this is BAD!
				// TOTAL incident shortwave diffuse radiation
				//std::cout << tdiff(i,j) << " " << glob_h(i,j) << " " << glob_start(i,j) << " " << albedo_grid(i,j) << std::endl;
				if (albedo_grid(i,j) > 0)
                                        total_diff(i,j) = /*tdiff(i,j) +*/ glob_h(i,j) - (glob_start(i,j) / albedo_grid(i,j));
			}
		}
	} else {
		for ( size_t i = 0; i < dimx; i++ ) {
			for ( size_t j = 0; j < dimy; j++ ) {
				if (dem.grid2D(i,j)==mio::IOUtils::nodata) continue;
				glob_h(i,j) = 0.;
				total_diff(i,j) = 0.;
			}
		}
	}
}

void TerrainRadiation::InitializeLW (const int i, const int j, double &/*lw_eps_stern*/)
{
	// longwave sky radiation calculation
	//lw_sky(i,j) = SkyLW (wvp, meteo2D(i,j).ta, sky_vf(i,j), ea_ref); //read cloud cover input from MeteoIO?!
	lw_sky(i,j) = meteo2d_ilwr(i,j)*viewFactorsObj.getSkyViewFactor(i,j);
	lwi(i,j) = lw_sky(i,j);

	// longwave radiation of ij emittable to the surrounding terrain without attenuation and view factor
	const double tss = t_snowold(i,j);
	lw_t(i,j) = Cst::stefan_boltzmann * (tss*tss*tss*tss);

	// computing the longwave emission start distribution with l1 vector norm
	lw_start_l1 += fabs( lw_t(i,j) );

	// the grid cell with the most unshot radiation is selected taking into account an air column reduction factor,
	// the actual (slope) area size, emittable longwave radiation and sum of total terrain view factor
	const double rad_t = viewFactorsObj.getSymetricTerrainViewFactor(i,j) * lw_t(i,j);

	//Optmisation #6 by GS : put the current element in the table of LW TerrainRadiation
	//if (TERRAIN_RADIATION!=0) {
		CellsList& lwtrbc = lwt_byCell[i*dimy+j]; //temporary reference to the current element's location
		lwtrbc.radiation = rad_t;
		lwtrbc.x = i;
		lwtrbc.y = j;
	//}

	lw_eps_stern += fabs(rad_t);
}
