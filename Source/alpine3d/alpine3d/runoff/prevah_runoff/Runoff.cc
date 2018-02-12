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
#include "Runoff.h"

#include "f77-uscore.h"
#include "fortran_and_c.h"

extern "C"
{
	void F77_FCN(initrunoff)(int *,int *,const char *);
	void F77_FCN(corerunoff)(int *,int *,const int *,double *);
}

using namespace mio;
using namespace std;

/**
 * @brief Constructor of Runoff instance.
 * @param sn_io is the ref to the IOManager instance from SnowPack. With this, Runoff don't
 * need to have a second IOManager
 */
Runoff::Runoff(mio::IOManager& sn_io, const double& /*in_thresh_rain*/)
{
	io = &sn_io;
	nx = ny = 0;
}

/**
 * @brief Destructor of Runoff which close the output stream which write the result files
 */
Runoff::~Runoff() {

}


/**
 * @brief This method initialize the importent parameters of Runoff.
 * @param in_dem is the Deographic object, which is needed for to take/controll the dimensions of the catchments
 * @param cfg is the configuration file, so Runoff know where to write the results maxexp
 */
bool Runoff::initialize(const mio::DEMObject& in_dem, const mio::Config& cfg)
{
	dem = in_dem;
	nx = in_dem.getNx();
	ny = in_dem.getNy();
	string config = "../setup/interpol_runoff.ini";
	cfg.getValue("RUNOFF_CFG", "Alpine3d", config, IOUtils::nothrow);
	F77_FCN(initrunoff)(&nx, &ny, config.c_str());

	runoff_total.resize(nx*ny+1, 0.); //1D array
	return true;
}

/**
 * @brief Fill runoff tables for one pixel
 * @param ix X coordinates of pixel for wich the runoff should be calculated
 * @param iy Y coordinates of pixel for wich the runoff should be calculated
 * @param soil_runoff soil runoff (ie: sub-surface)
 */
void Runoff::fillHydroTables(const unsigned int& ix, const unsigned int& iy, const double& soil_runoff)
{
	const double slope2horiz = 1. / cos( dem.slope(ix,iy)*mio::Cst::to_rad );

	runoff_total(iy*nx+ix) = soil_runoff*slope2horiz; // total runoff from snow and soil layers [mm hr-1]
}

void Runoff::setRunoff(const mio::Array2D<double>& /*runoff_surface*/, const mio::Array2D<double>& runoff_soil, const mio::Array2D<double>& /*glacier*/, const mio::Array2D<double>& /*psum*/, const mio::Array2D<double>& /*ta*/)
{
	size_t ncols, nrows;
	runoff_soil.size(ncols, nrows);
	if (ncols!=(unsigned)nx || nrows!=(unsigned)ny) {
		stringstream ss;
		ss << "DEM (" << nx << "," << ny << ") and runoff_grids (" << ncols << "," << nrows << ") have incompatible dimensions!";
		throw mio::InvalidArgumentException(ss.str(), AT);
	}

	for (unsigned int jj=0; jj<(unsigned)ny; jj++) {
		for (unsigned int ii=0; ii<(unsigned)nx; ii++) {
			fillHydroTables(ii, jj, runoff_soil(ii,jj));
		}
	}
}

/**
 * @brief Writes the Results for a specific day
 * @param i_date the Date to write the results
 */
void Runoff::output(const mio::Date& i_date)
{
	int yyyy,mm,dd,hh,mi;
	i_date.getDate(yyyy,mm,dd,hh,mi);
	const int date = yyyy*1000000+mm*10000+dd*100+hh;
	if (mi < 1e-5) {
		F77_FCN(corerunoff)(&nx, &ny, &date, &runoff_total[0]);
	}
}

