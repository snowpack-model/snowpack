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
/*------------------------------------------------------------------------------------------+
 |  This module contains the saltation model of Judith                                      |
 +------------------------------------------------------------------------------------------*/
/********************************************************************************************/
/*------------------------------------------------------------------------------------------+
 | 24.08.2000: The new very complicated model. Judith says that it is my fault, if it is    |
 | wrong. I hope not .....                                                                  |
 | Make a separate routine to model saltation only.                                         |
 | 26.04.2001.                                                                              |
 | Finally, on the Friday evening before the Swiss Bike Masters event, where Michael was    |
 | supposed to start, he also started to implement the last version of Judith's saltation   |
 | model. The GRID man Tuan Anh (Nguyen) had arrived and smiled.                            |
 | The Gaudergrat experiment GAUDEX was in good shape and almost everything was up and      |
 | running.                                                                                 |
 | 18.07.2003.                                                                              |
 +------------------------------------------------------------------------------------------*/
/********************************************************************************************/

/*------------------------------------------------------------------------------------------+
 | Includes and Defines                                                                     |
 +------------------------------------------------------------------------------------------*/
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/snowdrift/checksum.h>
#include <snowpack/libsnowpack.h>

#include <math.h>

using namespace std;

/**
 * @brief Calculate Saltation Fluxes for all Bottom Elements
 * @param setbound bool
 */
void SnowDriftA3D::compSaltation(bool setbound)
{

	unsigned int ix,iy;
	double tauS, tau_th, Ubar2;
	double weight, binding, sigma=300.;
	double dg;
	double flux_mean, cs_mean;  /* What we finally want */

	const mio::Grid2DObject store( snowpack->getGrid(SnGrids::STORE) );
	const mio::Grid2DObject swe( snowpack->getGrid(SnGrids::SWE) );

	/* Calculate the Fluxes for all Bottom Elements */
	for (iy=1; iy<ny-1; iy++) {
		for (ix=1; ix<nx-1; ix++){
			if (cH.grid2D(ix,iy)==mio::IOUtils::nodata) {
				saltation(ix,iy) = 0.0; c_salt(ix,iy) = 0.0;
				continue;
			}
			//Michi: Please CHECK these two lines!
			flux_mean=0;
			cs_mean=0;
			
			if ( cH.grid2D(ix, iy) <= 0.05 || (store(ix,iy)+swe(ix,iy))<5.) { 
				saltation(ix,iy) = 0.0; c_salt(ix,iy) = 0.0;
				continue;
			}

			/* Determine shear stress from the wind */
			/*      tauS = DENSITY_AIR*DSQR(nodes[n].Km / nodes[n].lm); */ /* First try - Nice try */
			
			if (TKE) {
				/* From TKE using similarity */
				/* Attention FUDGE */
				tauS = std::min(Constants::density_air*nodes_e.grid3D(ix,iy,1)/5.,0.2);

			} else {
				/* From wind speed using log-profile */
			    Ubar2 = mio::Optim::pow2(nodes_u.grid3D(ix,iy,2)) + mio::Optim::pow2(nodes_v.grid3D(ix,iy,2)) + mio::Optim::pow2(nodes_w.grid3D(ix,iy,2));
			    tauS =  Constants::density_air*Ubar2*mio::Optim::pow2(Saltation::karman/log((nodes_z.grid3D(ix,iy,2) - nodes_z.grid3D(ix,iy,1))/Saltation::z0_salt));
			    if (tauS == 0.) {
				printf("\n tauS zero: ix:%d, iy:%d, u:%f, v:%f"
					,ix,iy,nodes_u.grid3D(ix,iy,2),nodes_v.grid3D(ix,iy,2));
			    }
			}
			if (thresh_snow) {
				weight = 0.02*Constants::density_ice*(sp.grid2D(ix, iy) + 1.)*mio::Cst::gravity*MM_TO_M(rg.grid2D(ix, iy));
				binding = 0.0015*sigma*N3.grid2D(ix, iy)*rb.grid2D(ix, iy)*rb.grid2D(ix, iy)/rg.grid2D(ix, iy)/rg.grid2D(ix, iy);
				tau_th = std::max(tau_thresh, SnowDrift::schmidt_drift_fudge*(weight + binding));
				dg = std::min(grain_size,std::max(0.3*grain_size,2.*rg.grid2D(ix, iy)));
			} else {
				tau_th = tau_thresh; /* Pa */
				dg = grain_size;
			}/* m */

			/* Calculate Flux and Lower Concentration Boundary Condition*/
			if (!saltation_obj.compSaltation(tauS, tau_th, nodes_slope.grid3D(ix,iy,1)*(180./Constants::pi), dg, flux_mean, cs_mean)) {
				cout<<" Could not calculate Saltation"<<endl;
			    return;
			}
			saltation(ix,iy) = flux_mean; c_salt(ix,iy) = c_red*cs_mean;
		} /* for ix */
	}

	/* First set Zero Gradient Boundary Condition */
	if (setbound) {
		for (iy=0; iy<ny; iy++) {
			saltation(0,iy) =  saltation(1,iy); saltation(nx-1,iy) = saltation(nx-2,iy);
			c_salt(0,iy) = c_salt(1,iy); c_salt(nx-1,iy) = c_salt(nx-2,iy);
		}

		//Michi: the old stuffs do not set the boundary correctly!
		for (ix=0; ix<nx; ix++) {
			saltation(ix,0) =  saltation(ix,1); saltation(ix,ny-1) = saltation(ix,ny-2);
			c_salt(ix,0) = c_salt(ix,1); c_salt(ix,ny-1) = c_salt(ix,ny-2);
		}
	}
	/* End of Saltation Flux */
	//  DEBUG("Checksum saltation=%lf, c_salt=%lf",checksum(saltation),checksum(c_salt));
	return;
}

