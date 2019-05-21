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
#ifndef ENERGYBALANCE_H
#define ENERGYBALANCE_H

#include <string>
#include <meteoio/MeteoIO.h>

class SnowpackInterface;

#include <alpine3d/ebalance/RadiationField.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/SnowpackInterface.h>

/**
 * @page radiation_balance Radiation balance
 * First, the measured radiation (as provided by each station that measured both ISWR, TA and RH) is interpolated over the domain:
 *      -# At each station that provides radiation the splitting coefficient (into direct and diffuse components)
 *           and an atmospheric loss factor are computed by comparing the measured radiation with the potential radiation;
 *      -# The splitting coefficient Md and the atmospheric loss factor Corr are spatially interpolated with an Inverser Distance
 *           weighting algorithm;
 *      -# At each cell, the clear sky potential radiation is computed by (Iqbal, 1983), split into direct and diffuse components using its
 *           spatially interpolated value and corrected by the spatially interpolated atmopsheric loss factor. The direct component is
 *           also projected on the local slope.
 *
 * Then the effects of the terrain on the radiation field are accounted for. These are twofold: first the topography can cast shade
 * on some parts of the domain and second some radiation can be reflected by the terrain on some other cells (see \ref principles_ebalance).
 *
 * @section topo_shading Topographical shading
 * This handles the shading of cells by the topography. The position of the sun is computed and the cells that don't have a direct view of the sun
 * only receive the diffuse fraction of the global radiation while the cells with direct view of the sun receive both the direct and diffuse
 * components. Please keep in mind that when looking for the horizon (in order to compute shading), the search will <b>stop</b> at cells that are
 * set to nodata in the dem. This means that nodata cells on the border of the domain of interest will prevent searching for shading outside of the
 * domain while nodata cells \em within the domain should be as much as possible avoided.
 *
 * @section terrain_radiation Terrain radiation
 * The second effect is computed when the key Terrain_Radiation is set to true (by default it is set to false)
 * in the [EBalance] section of the configuration file and is handled by a choice of various algorithms,
 * as chosen in the configuration file with the key Terrain_Radiation_Method that
 * can take any of the following choices:
 *     - SIMPLE : a very basic (but fast) guess, see TerrainRadiationSimple
 *     - FULL : a parallelized implementation of a radiosity algorithm, see TerrainRadiation
 *     - HELBIG : the original radiosity implementation, not parallelized, see TerrainRadiationHelbig
 *     - PETSC : a parallel implementation relying on PETSC for added performances, see TerrainRadiationPETSc
 *
 * @code
 * [EBalance]
 * Terrain_Radiation = true
 * Terrain_Radiation_Method = SIMPLE
 * @endcode
 */

class EnergyBalance
{
	public:
		EnergyBalance( const unsigned int& i_nbworkers, const mio::Config& cfg, const mio::DEMObject &dem_in);
		EnergyBalance(const EnergyBalance&);
		~EnergyBalance( );

		EnergyBalance& operator=(const EnergyBalance&); ///<Assignement operator, required because of pointer member

		void setSnowPack(SnowpackInterface &mysnowpack );
		void setAlbedo( const mio::Grid2DObject &in_albedo );

		void setMeteo(const mio::Grid2DObject& in_ilwr,
		              const mio::Grid2DObject& in_ta, const mio::Grid2DObject& in_rh, const mio::Grid2DObject& in_p, const mio::Date timestamp);

		void setStations(const std::vector<mio::MeteoData>& in_vecMeteo);
		double getTiming() const;
		void Destroy();
		std::string getGridsRequirements() const;

	private:
		SnowpackInterface *snowpack;
		TerrainRadiationAlgorithm *terrain_radiation;
		std::vector<RadiationField*> radfields;
		mio::DEMObject dem;
		std::vector<mio::MeteoData> vecMeteo;
		mio::Grid2DObject albedo;
		mio::Array2D<double> direct, diffuse, reflected;
		mio::Timer timer;

		size_t dimx, dimy;
		unsigned int nbworkers;
};

#endif
