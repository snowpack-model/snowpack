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
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/OMPControl.h>

using namespace mio;
using namespace std;

EnergyBalance::EnergyBalance(const unsigned int& i_nbworkers, const mio::Config& cfg, const mio::DEMObject &dem_in)
              : snowpack(NULL), terrain_radiation(NULL), radfields(i_nbworkers), dem(dem_in), vecMeteo(),
                albedo(dem_in, 0.), direct(), diffuse(), reflected(),
                timer(), dimx(dem_in.getNx()), dimy(dem_in.getNy()), nbworkers(i_nbworkers)
{
	MPIControl& instance = MPIControl::instance();

	size_t startx = 0, nx = dimx;
	instance.getArraySliceParams(dimx, startx, nx);

	#pragma omp parallel for schedule(static)
	for (size_t ii=0; ii<nbworkers; ii++) {
		size_t thread_startx, thread_nx;
		OMPControl::getArraySliceParams(nx, nbworkers, ii, thread_startx, thread_nx);
		const size_t offset = startx + thread_startx;

		#pragma omp critical(ebWorkers_status)
		std::cout << "[i] EnergyBalance worker " << ii << " on process " << instance.rank() << " will start at offset " << offset << " with nx " << thread_nx << "\n";
		radfields[ii] = new RadiationField(dem_in, offset, thread_nx);
	}

	if (instance.master())
		std::cout << "[i] EnergyBalance initialized a total of " << instance.size() << " process(es) with " << nbworkers << " worker(s) each\n";

	// Every MPI process will have its own copy of terrain_radiation object with full DEM
	const bool enable_terrain_radiation = cfg.get("Terrain_Radiation", "EBalance");
	if (enable_terrain_radiation) {
		terrain_radiation = TerrainRadiationFactory::getAlgorithm(cfg, dem, nbworkers);

		const std::string algo = terrain_radiation->algo;
		if (instance.master())
			std::cout << "[i] Using terrain radiation with model: " << algo << "\n";
	}
}

EnergyBalance::~EnergyBalance() {
	Destroy( );
}

EnergyBalance& EnergyBalance::operator=(const EnergyBalance& source) {
	if (this != &source) {
		snowpack = source.snowpack;
		terrain_radiation = source.terrain_radiation;
		radfields = source.radfields;
		dem = source.dem;
		vecMeteo = source.vecMeteo;
		albedo = source.albedo;
		direct = source.direct;
		diffuse = source.diffuse;
		reflected = source.reflected;
		timer = source.timer;
		dimx = source.dimx;
		dimy = source.dimy;
		nbworkers = source.nbworkers;
	}
	return *this;
}

std::string EnergyBalance::getGridsRequirements() const
{
	return "TOP_ALB";
}

void EnergyBalance::Destroy()
{
	while (!radfields.empty()) {
		delete radfields.back();
		radfields.pop_back();
	}

	if (terrain_radiation) {
		delete terrain_radiation;
		terrain_radiation = NULL;
	}
}

void EnergyBalance::setSnowPack(SnowpackInterface& mysnowpack)
{
	snowpack = &mysnowpack;
}

void EnergyBalance::setAlbedo(const mio::Grid2DObject& in_albedo)
{
	albedo = in_albedo;

	direct.resize(0, 0); //resetting these grids that are not valid anymore
	diffuse.resize(0, 0);
	reflected.resize(0, 0);
}

void EnergyBalance::setStations(const std::vector<mio::MeteoData>& in_vecMeteo)
{
	vecMeteo = in_vecMeteo;

	direct.resize(0, 0); //resetting these grids that are not valid anymore
	diffuse.resize(0, 0);
	reflected.resize(0, 0);
}

void EnergyBalance::setMeteo(const mio::Grid2DObject& in_ilwr,
                             const mio::Grid2DObject& in_ta, const mio::Grid2DObject& in_rh, const mio::Grid2DObject& in_p, const mio::Date timestamp)
{
	timer.restart();
	direct.resize(dimx, dimy);
	diffuse.resize(dimx, dimy);

	#pragma omp parallel for schedule(dynamic)
	for (size_t ii=0; ii<nbworkers; ii++) {
		radfields[ii]->setStations(vecMeteo, albedo); //calculate the parameters at the radiation stations
		size_t startx, nx;
		radfields[ii]->getBandOffsets(startx, nx);
		radfields[ii]->setMeteo(mio::Grid2DObject(in_ta, startx, 0, nx, dimy),
		                       mio::Grid2DObject(in_rh, startx, 0, nx, dimy),
		                       mio::Grid2DObject(in_p, startx, 0, nx, dimy),
		                       mio::Grid2DObject(albedo, startx, 0, nx, dimy));

		mio::Array2D<double> band_direct, band_diffuse;
		radfields[ii]->getRadiation(band_direct, band_diffuse);
		direct.fill(band_direct, startx, 0, nx, dimy);
		diffuse.fill(band_diffuse, startx, 0, nx, dimy);
	}
	MPIControl::instance().allreduce_sum(direct);
	MPIControl::instance().allreduce_sum(diffuse);

	if (terrain_radiation) {
		// note: parallelization has to take place inside the TerrainRadiationAlgorithm implementations
		terrain_radiation->setMeteo(albedo.grid2D, in_ta.grid2D, in_rh.grid2D, in_ilwr.grid2D);
		terrain_radiation->getRadiation(direct, diffuse, reflected);
	}

	if (MPIControl::instance().master())
		cout << "[i] Ebalance simulation done for " << timestamp.toString(Date::ISO) << "\n";

	if (snowpack) {
		double solarAzimuth, solarElevation;
		radfields[0]->getPositionSun(solarAzimuth, solarElevation); //we need it only for handing over to snowpack

		mio::Array2D<double> ilwr = in_ilwr.grid2D;
		mio::Array2D<double> global = direct+diffuse; //otherwise the compiler does not match the types

		if (!reflected.empty()) global += reflected;

		timer.stop();
		try {
			snowpack->setRadiationComponents(global, ilwr, diffuse, solarElevation, timestamp); //this triggers Snowpack calculation
		} catch(std::exception& e) {
			std::cout << "[E] Exception in snowpack->setRadiationComponents()\n";
			cout << e.what() << endl;
			std::abort(); //force core dump
		}
	}
	timer.stop();
}

double EnergyBalance::getTiming() const
{
	return timer.getElapsed();
}
