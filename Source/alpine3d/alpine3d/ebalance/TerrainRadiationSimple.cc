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
#include <alpine3d/ebalance/TerrainRadiationSimple.h>
#include <alpine3d/ebalance/ViewFactors.h>

using namespace mio;

TerrainRadiationSimple::TerrainRadiationSimple(const mio::Config& /*i_cfg*/, const mio::DEMObject& dem_in, const std::string& method)
                       : TerrainRadiationAlgorithm(method),
                         albedo_grid(dem_in.getNx(), dem_in.getNy(), IOUtils::nodata), sky_vf(dem_in.getNx(), dem_in.getNy(), 0),
                         dimx(dem_in.getNx()), dimy(dem_in.getNy()), startx(0), endx(dimx)
{
	// In case we're running with MPI enabled, calculate the slice this process is responsible for
	size_t nx = dimx;
	MPIControl::instance().getArraySliceParams(dimx, startx, nx);
	endx = startx + nx;

	initSkyViewFactor(dem_in);

}

void TerrainRadiationSimple::getRadiation(mio::Array2D<double>& /*direct*/,
                                          mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain,
                                          const mio::Array2D<double>& direct_unshaded_horizontal,
                                          const mio::Array2D<double>& /*total_ilwr*/, mio::Array2D<double>& /*sky_ilwr*/,
                                          mio::Array2D<double>& /*terrain_ilwr*/,
                                          double /*solarAzimuth*/, double /*solarElevation*/)
{
	MPIControl& mpicontrol = MPIControl::instance();
	terrain.resize(dimx, dimy, 0.);  //so reduce_sum works properly when it sums full grids
	Array2D<double> diff_corr(dimx, dimy, 0.); //so reduce_sum works properly when it sums full grids

	if (mpicontrol.master()) {
		std::cout << "[i] Calculating terrain radiation with simple method, using " << mpicontrol.size();
		if (mpicontrol.size()>1) std::cout << " processes\n";
		else std::cout << " process\n";
	}

	#pragma omp parallel for collapse(2)
	for (size_t jj=0; jj<dimy; jj++) {
		for (size_t ii=startx; ii<endx; ii++) {
			if (sky_vf(ii,jj) == IOUtils::nodata || getAlbedo(ii,jj) == IOUtils::nodata) {
				terrain(ii,jj) = IOUtils::nodata;
				diff_corr(ii,jj) = IOUtils::nodata;
			} else {
				const double terrain_reflected = getAlbedo(ii,jj) * (direct_unshaded_horizontal(ii,jj)+diffuse(ii,jj));
				const double terrain_viewFactor = 1. - sky_vf(ii,jj);
				terrain(ii,jj) = terrain_viewFactor * terrain_reflected;
				diff_corr(ii,jj) = diffuse(ii,jj) * sky_vf(ii,jj);
			}
		}
	}

	mpicontrol.reduce_sum(terrain);
	mpicontrol.reduce_sum(diff_corr);
	diffuse = diff_corr; //return the corrected diffuse radiation
}

//retrieve the albedo to use for terrain reflections
double TerrainRadiationSimple::getAlbedo(const size_t& ii, const size_t& jj)
{
	//return albedo_grid(ii,jj); //the easiest variant: the local albedo of the cell
	if (albedo_grid(ii,jj)==IOUtils::nodata) return IOUtils::nodata;

	static const size_t nr_cells_around = 2;
	//without considering unsigned value: const size_t ll_min = std::max(0, static_cast<int>(jj-nr_cells_around));
	const size_t ll_min = (jj-nr_cells_around)>dimy? 0 : jj-nr_cells_around; //since negative values will wrap around
	const size_t ll_max = std::min(dimy, jj+nr_cells_around+1); //+1 so we can compare with < instead of <=
	const size_t kk_min = (ii-nr_cells_around)>dimx? 0 : ii-nr_cells_around; //since negative values will wrap around
	const size_t kk_max = std::min(dimx, ii+nr_cells_around+1); //+1 so we can compare with < instead of <=

	unsigned short int count = 0;
	double sum = 0.;

	for (size_t ll=ll_min; ll<ll_max; ll++) {
		for (size_t kk=kk_min; kk<kk_max; kk++) {
			if (albedo_grid(kk,ll)==IOUtils::nodata) continue;
			sum += albedo_grid(kk,ll);
			count++;
		}
	}

	const double albedo = (count!=0)? sum / (double)count : IOUtils::nodata;
	return albedo;
}

void TerrainRadiationSimple::setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& /*ta*/)
{
	albedo_grid = albedo;
}

void TerrainRadiationSimple::getSkyViewFactor(mio::Array2D<double> &o_sky_vf) {
	o_sky_vf = sky_vf;
	MPIControl::instance().reduce_sum(o_sky_vf);
}

void TerrainRadiationSimple::initSkyViewFactor(const mio::DEMObject &dem)
{
	#pragma omp parallel for collapse(2)
	for (size_t jj=0; jj<dimy; jj++) {
		for (size_t ii=startx; ii<endx; ii++) {
			if (dem(ii,jj) != IOUtils::nodata)
				sky_vf(ii,jj) = mio::DEMAlgorithms::getCellSkyViewFactor(dem, ii, jj);
		}
	}
}
