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
	mio::Grid2DObject Um( grid_VW, 0. );			// U component of wind on staggered grid
	mio::Grid2DObject Vm( grid_VW, 0. );			// V component of wind on staggered grid
	mio::Grid2DObject dM( ErodedMass, 0. );			// Local mass perturbation due to drifting snow redistribution.
	mio::Grid2DObject tmp_ErodedMass( ErodedMass );		// Copy of initially eroded mass.
	mio::Grid2DObject grid_snowdrift_out = tmp_ErodedMass;	// Output mass
	grid_snowdrift_out(0.);

	// If there is no wind, then there is no transport of eroded snow between grid cells.
	if (grid_VW.grid2D.getMax() == 0.) {
		return ErodedMass;
	}

	// Eroded mass must be greater than or equal to zero.
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			if (tmp_ErodedMass(ix, iy) == IOUtils::nodata || tmp_ErodedMass(ix, iy) < 0.) {
				tmp_ErodedMass(ix, iy) = 0.;
			}

			winderosiondeposition(ix, iy) = tmp_ErodedMass(ix, iy);

			// Add a determination of u and v componets for each grid cell.
			// Following Pomeroy and Gray 1990, we set the saltation velocity <U, V> to be equal to 2.8 times
			// the friction threshold velocity (ustar_th) in the direction parallel to the 10m wind.
			const double tmp_ustar_th = (ustar_th(ix, iy) == Constants::undefined) ? (0.) : (ustar_th(ix, iy));
			U(ix, iy) = IOUtils::VWDW_TO_U(std::max(0., 2.8 * tmp_ustar_th ), grid_DW(ix, iy));
			V(ix, iy) = IOUtils::VWDW_TO_V(std::max(0., 2.8 * tmp_ustar_th ), grid_DW(ix, iy));
		}
	}

	// Calculate wind speed on "staggered" grid, and determined the max wind speed components
	// --------------------------------------------------------------------------------------
	// The staggered grid is shifted left and down, w.r.t. the base grid. This means that ix=1 in the staggered gris is half-way ix=0 and ix=1 in the base grid. Similarly
	// iy=1 is half-way between iy=1 and iy=0 in the base grid.
	// In the staggered grid, ix=0 and iy=0 is only used for the periodice boundary conditions, i.e., ix=0 in the staggered grid is half-way ix=dimx-1 and ix=0.
	double Umax = 0.;
	double Vmax = 0.;
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			// U-component
			if (ix == 0) {
				// Periodic boundary conditions only
				if (U(ix, iy) == IOUtils::nodata && U(dimx-1, iy) == IOUtils::nodata) {
					Um(ix, iy) = IOUtils::nodata;
				} else if (U(dimx-1, iy) == IOUtils::nodata) {
					Um(ix, iy) = U(ix, iy);
				} else if (U(ix, iy) == IOUtils::nodata) {
					Um(ix, iy) = U(dimx-1, iy);
				} else {
					Um(ix, iy) = 0.5 * (U(dimx-1, iy) + U(ix, iy));
				}
			} else {
				if (U(ix, iy) == IOUtils::nodata && U(ix-1, iy) == IOUtils::nodata) {
					Um(ix, iy) = IOUtils::nodata;
				} else if (U(ix-1, iy) == IOUtils::nodata) {
					Um(ix, iy) = U(ix, iy);
				} else if (U(ix, iy) == IOUtils::nodata) {
					Um(ix, iy) = U(ix-1, iy);
				} else {
					Um(ix, iy) = 0.5 * (U(ix-1, iy) + U(ix, iy));
				}
			}
			// V-component
			if (iy == 0) {
				// Periodic boundary conditions only
				if (V(ix, iy) == IOUtils::nodata && V(ix, dimy-1) == IOUtils::nodata) {
					Vm(ix, iy) = IOUtils::nodata;
				} else if (V(ix, dimy-1) == IOUtils::nodata) {
					Vm(ix, iy) = V(ix, iy);
				} else if (V(ix, iy) == IOUtils::nodata) {
					Vm(ix, iy) = V(ix, dimy-1);
				} else {
					Vm(ix, iy) = 0.5 * (V(ix, dimy-1) + V(ix, iy));
				}
			} else {
				if (V(ix, iy) == IOUtils::nodata && V(ix, iy-1) == IOUtils::nodata) {
					Vm(ix, iy) = IOUtils::nodata;
				} else if (V(ix, iy-1) == IOUtils::nodata) {
					Vm(ix, iy) = V(ix, iy);
				} else if (V(ix, iy) == IOUtils::nodata) {
					Vm(ix, iy) = V(ix, iy-1);
				} else {
					Vm(ix, iy) = 0.5 * (V(ix, iy-1) + V(ix, iy));
				}
			}
			// For CFL:
			if (fabs(Vm(ix,iy)) > Vmax) Vmax = fabs(Vm(ix,iy));
			if (fabs(Um(ix,iy)) > Umax) Umax = fabs(Um(ix,iy));
		}
	}
	
	// Calculate largest possible drifting snow sub time step such that the scheme is still stable and satisfies
	// the CFL criterion.
	const double C_max = 0.999;					// Courant number used to calculate sub time step.
	double sub_dt = dt;
	if (Umax + Vmax == 0) {
		return ErodedMass;
	} else {
		sub_dt = std::min(C_max * dx / (Umax + Vmax), dt);	// Sub time step
		std::cout << "[i] Explicit snow drift sub time step = " << sub_dt << " seconds\n";
	}

	// Fill grid_snowdrift_out
	for (double time_advance = 0.; time_advance < dt; time_advance += sub_dt) {
		if (time_advance + sub_dt > dt) {
			sub_dt = dt - time_advance;
		}

		// Calculate change of suspended mass
		for (size_t iy=0; iy<dimy; iy++) {
			for (size_t ix=0; ix<dimx; ix++) {
				if (Um(ix, iy) != IOUtils::nodata && Vm(ix, iy) != IOUtils::nodata) {
					// Periodic boundary conditions are applied when set and treating ix==0 or iy==0
					// Otherwise, zero flux boundary conditions are assumed. When constant flux is set,
					// this is corrected for afterwards.
					if (ix == 0) {
						if (bc == PERIODIC) {
							if(Um(ix, iy)>0) {
								const double deltaM = tmp_ErodedMass(dimx-1, iy) * fabs(Um(ix, iy)) * (sub_dt / dx);
								dM(dimx-1, iy) -= deltaM;
								dM(ix, iy)   += deltaM;
							} else {
								const double deltaM = tmp_ErodedMass(ix, iy) * fabs(Um(ix, iy)) * (sub_dt / dx);
								dM(dimx-1, iy) += deltaM;
								dM(ix, iy)   -= deltaM;
							}
						}
					} else {
						if(Um(ix, iy)>0) {
							const double deltaM = tmp_ErodedMass(ix-1, iy) * fabs(Um(ix, iy)) * (sub_dt / dx);
							dM(ix-1, iy) -= deltaM;
							dM(ix, iy)   += deltaM;
						} else {
							const double deltaM = tmp_ErodedMass(ix, iy) * fabs(Um(ix, iy)) * (sub_dt / dx);
							dM(ix-1, iy) += deltaM;
							dM(ix, iy)   -= deltaM;
						}
					}
					if (iy == 0) {
						if (bc == PERIODIC) {
							if(Vm(ix, iy)>0) {
								const double deltaM = tmp_ErodedMass(ix, dimy-1) * fabs(Vm(ix, iy)) * (sub_dt / dx);
								dM(ix, dimy-1) -= deltaM;
								dM(ix, iy)   += deltaM;
							} else {
								const double deltaM = tmp_ErodedMass(ix, iy) * fabs(Vm(ix, iy)) * (sub_dt / dx);
								dM(ix, dimy-1) += deltaM;
								dM(ix, iy)   -= deltaM;
							}
						}
					} else {
						if(Vm(ix, iy)>0) {
							const double deltaM = tmp_ErodedMass(ix, iy-1) * fabs(Vm(ix, iy)) * (sub_dt / dx);
							dM(ix, iy-1) -= deltaM;
							dM(ix, iy)   += deltaM;
						} else {
							const double deltaM = tmp_ErodedMass(ix, iy) * fabs(Vm(ix, iy)) * (sub_dt / dx);
							dM(ix, iy-1) += deltaM;
							dM(ix, iy)   -= deltaM;
						}
					}
				}

				// Set boundary condition
				if (bc == CONSTANTFLUX) {
					// Constant-flux (i.e., dM at boundaries == 0):
					dM(0,iy) = 0.;
					dM(dimx-1,iy) = 0.;
					dM(ix,0) = 0.;
					dM(ix,dimy-1) = 0.;
				}
			}
		}

		// Apply calculated change of suspended mass
		for (size_t iy=0; iy<dimy; iy++) {
			for (size_t ix=0; ix<dimx; ix++) {
				// Negative suspended mass (tmp_ErodedMass) is impossible. If this field becomes negative, set it to zero.
				if (ErodedMass(ix, iy) != IOUtils::nodata) {
					tmp_ErodedMass(ix, iy) = ErodedMass(ix, iy) + dM(ix, iy);
					// This if statement could introduce a mass balance error.
					if (tmp_ErodedMass(ix, iy) < 0.) {
						//std::cout << "[W] Explicit Snow Drift: Potential mass balance violation #1 ";
						//std::cout << ix << " " << iy << " ";
						//std::cout << tmp_ErodedMass(ix, iy) << "\n";
						tmp_ErodedMass(ix, iy) = 0.;
					}
				}
			}
		}
	}

	double s1 = 0, s2 = 0., s3 = 0.;
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			if (dM(ix, iy) > 0.) {s1 += dM(ix, iy);}
			if (ErodedMass(ix, iy) != IOUtils::nodata) {
				// This line could introduce a mass balance error.
				// The change in mass cannot be more than was originally eroded!
				//    -95        - 10               -100
				if (dM(ix, iy) - Constants::eps < -ErodedMass(ix, iy)) {
					// std::cout << "[W] Explicit Snow Drift: Potential mass balance violation #2 ";
					// std::cout << ix << " " << iy << " ";
					// std::cout << dM(ix, iy) << -ErodedMass(ix, iy) << "\n";
					s2 += -ErodedMass(ix, iy) - dM(ix, iy);
					dM(ix, iy) = -ErodedMass(ix, iy);
				}
			}
			s3+=dM(ix, iy);
		}
	}
	if(s2!=0. || s3!=0.) {
		std::cout << "[W] Explicit Snow Drift: Potential mass balance violation #2 ";
		printf("%.10f (total pos: %.10f). Sum dM==%.10f\n", s2, s1, s3);
		//s2+=s3;
		s3=0.;
	}
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			if (ErodedMass(ix, iy) != IOUtils::nodata) {
				if (dM(ix, iy) > 0.) {dM(ix, iy) += (-s2*(dM(ix, iy)/s1));}
				winderosiondeposition(ix, iy) = dM(ix, iy);
				grid_snowdrift_out(ix, iy) = ErodedMass(ix, iy) + dM(ix, iy);
				if(grid_snowdrift_out(ix, iy) < Constants::eps) grid_snowdrift_out(ix, iy) = 0.;
			}
		}
	}

	return grid_snowdrift_out;
}