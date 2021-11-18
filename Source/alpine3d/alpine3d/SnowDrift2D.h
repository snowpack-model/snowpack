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
#ifndef SNOWDRIFT2D_H
#define SNOWDRIFT2D_H

#include <meteoio/MeteoIO.h>

/**
 * @page snowdrift2D SnowDrift2D
 * This module calculates 2D drifting snow using the advection equation, or with a simpler version using a statistical terrain analysis.
 */
 class SnowDrift2D
  {
	public:
		SnowDrift2D(const mio::Config& cfg, const double in_A3DtimeStep, const mio::DEMObject& in_dem);
		~SnowDrift2D();

		SnowDrift2D& operator=(const SnowDrift2D&); ///<Assignement operator, required because of pointer member
		void reset_output();
		void calcSimpleSnowDrift(const mio::Grid2DObject& ErodedMass, mio::Grid2DObject& psum, const mio::Grid2DObject& tmp_vw_drift);
		mio::Grid2DObject calcExplicitSnowDrift(const mio::Grid2DObject grid_VW, const mio::Grid2DObject grid_DW, const mio::Grid2DObject& ErodedMass, const mio::Grid2DObject& ustar_th);
		mio::Grid2DObject winderosiondeposition;

	private:
		void init(const mio::Config& cfg);
		size_t dimx, dimy;	// dimensions, will be derived from a DEM provided to the class constructor
		double dx, dt;		// cell size, and time step, set in the class constructor
		enum BoundaryConditions {CONSTANTFLUX, ZEROFLUX, PERIODIC};
		BoundaryConditions bc;	// Boundary condition to use for advection equation
};

#endif
