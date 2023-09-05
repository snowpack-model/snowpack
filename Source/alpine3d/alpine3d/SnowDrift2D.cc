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
#include <alpine3d/SnowDrift2D.h>
#include <snowpack/libsnowpack.h>
#include <alpine3d/MPIControl.h>

using namespace std;
using namespace mio;

SnowDrift2D::SnowDrift2D(const mio::Config& cfg, const double in_A3DtimeStep, const mio::DEMObject& in_dem)
              : dimx(0), dimy(0), dx(0), dt(0.), bc(CONSTANTFLUX)
{
	dx = in_dem.cellsize;		// Cell size in m, assuming equal in x and y.
	dimx = in_dem.getNx();		// x-dimension
	dimy = in_dem.getNy();		// y-dimension
	dt = in_A3DtimeStep * 86400.;	// From time steps in days to seconds.
	winderosiondeposition.set(in_dem, 0.);
	init(cfg);
}

void SnowDrift2D::init(const mio::Config& cfg)
{
	std::string tmp_bc_string="CONSTANTFLUX";	// The default option
	cfg.getValue("SNOWDRIFT2D_BC", "Alpine3D", tmp_bc_string, IOUtils::nothrow);
	if (tmp_bc_string=="CONSTANTFLUX") {
		bc=CONSTANTFLUX;
	} else if (tmp_bc_string=="ZEROFLUX") {
		bc=ZEROFLUX;
	} else if (tmp_bc_string=="PERIODIC") {
		bc=PERIODIC;
	} else {
		throw mio::InvalidArgumentException("Unknown boundary condition '" + tmp_bc_string + "' for key SNOWDRIFT2D_BC. Only CONSTANTFLUX, ZEROFLUX and PERIODIC are supported.", AT);
	}
	
	if (MPIControl::instance().master()) 
		std::cout << "[i] Initialization SnowDrift2D complete.\n";
}

void SnowDrift2D::reset_output()
{
	winderosiondeposition.set(winderosiondeposition, 0.);
}

/**
 * @brief Calculates drifting snow by applying a statistical scheme
 * @author Nander Wever
 */
void SnowDrift2D::calcSimpleSnowDrift(const mio::Grid2DObject& tmp_ErodedMass, mio::Grid2DObject& psum, const mio::Grid2DObject& tmp_vw_drift)
{
	double sum_positive_exposure = 0.;
	double sum_erodedmass = 0.;

	const double max_sx = -tmp_vw_drift.grid2D.getMax(); //positive
	for (size_t ii=0; ii<tmp_vw_drift.size(); ii++) {
		if(tmp_vw_drift(ii) < 0.) {
			sum_positive_exposure += -tmp_vw_drift(ii) / max_sx;
		}
		sum_erodedmass += tmp_ErodedMass(ii);
	}
	if (sum_positive_exposure == 0) return;

	const double ratio = sum_erodedmass / sum_positive_exposure;
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			double deposition = 0.;
			if(tmp_vw_drift(ix, iy) < 0.) {
				const double val = -tmp_vw_drift(ix, iy) * ratio / max_sx;
				deposition = val;
				psum(ix, iy) += val;
			} else {
				deposition = -tmp_ErodedMass(ix,iy);
			}
			winderosiondeposition(ix, iy) = deposition;
		}
	}
	return;
}


/**
 * @brief Calculates explicit drifting snow by solving the 2-D advection equation using a first-order upwind finite difference scheme.
 * @param grid_VW grid with wind speed
 * @param grid_DW grid with wind diection
 * @param ErodedMass grid with mass of snow particles in the air (drifting snow mass)
 * @param ustar_th grid with threshold friction velocity
 * @return 2-D grid with mass of snow particles in the air (i.e., similar to ErodedMass) after solving the 2-D advection equation
 * @author Nander Wever and Eric Keenan
 */
mio::Grid2DObject SnowDrift2D::calcExplicitSnowDrift(const mio::Grid2DObject grid_VW, const mio::Grid2DObject grid_DW, const mio::Grid2DObject& ErodedMass, const mio::Grid2DObject& ustar_th)
{
	// Retrieve and initialize grids
	mio::Grid2DObject U( grid_VW, 0. );			// U component of wind
	mio::Grid2DObject V( grid_VW, 0. );			// V component of wind
	mio::Grid2DObject Q_x( ErodedMass, 0. );		// Mass flux through a vertical gate in the x direction (kg/m/s)
	mio::Grid2DObject Q_y( ErodedMass, 0. );		// Mass flux through a vertical gate in the y direction (kg/m/s) 
	mio::Grid2DObject divQ( ErodedMass, 0. );		// (Divergence of Q (kg/m^2/s)
	mio::Grid2DObject grid_snowdrift_out = ErodedMass;	// Output mass
	grid_snowdrift_out(0.);


	// Get constants
	const double L = 10.; // Fetch length (m). Note that this is currently hard coded to 10 m. This should be replaced with a variable determined by the configuration file.

	// If there is no wind, then there is no transport of eroded snow between grid cells.
	if (grid_VW.grid2D.getMax() == 0.) {
		return ErodedMass;
	}

	// Loop over all grid cells to determine U and V
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			U(ix, iy) = IOUtils::VWDW_TO_U(std::max(0., grid_VW(ix, iy)), grid_DW(ix, iy));
			V(ix, iy) = IOUtils::VWDW_TO_V(std::max(0., grid_VW(ix, iy)), grid_DW(ix, iy));
		}
	}

	// Loop over all grid cells to determine Q_x and Q_y
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			if ( grid_VW(ix, iy) == IOUtils::nodata || grid_VW(ix, iy) <= 0.) {
				Q_x(ix, iy) = 0.;
				Q_y(ix, iy) = 0.;
			} else {
				Q_x(ix, iy) = (ErodedMass(ix, iy) * L * U(ix, iy)) / (dt * grid_VW(ix, iy));
				Q_y(ix, iy) = (ErodedMass(ix, iy) * L * V(ix, iy)) / (dt * grid_VW(ix, iy));
			}
		}
	}

	// Loop over all interior grid cells to calculate divQ using central difference.
	for (size_t iy=1; iy<dimy-1; iy++) {
		for (size_t ix=1; ix<dimx-1; ix++) {
			divQ(ix, iy) = (Q_x(ix + 1, iy) - Q_x(ix - 1, iy)) / (2. * dx) + (Q_y(ix, iy + 1) - Q_y(ix, iy - 1)) / (2. * dx);
		}
	}

	// Loop over grid cells to calculate ErodedMass perturbation
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			double dM = -divQ(ix, iy) * dt;
			winderosiondeposition(ix, iy) += dM;
			grid_snowdrift_out(ix, iy) += ErodedMass(ix, iy) + dM;
			if(grid_snowdrift_out(ix, iy) < Constants::eps) grid_snowdrift_out(ix, iy) = 0.;
		}
	}

	return grid_snowdrift_out;
}
