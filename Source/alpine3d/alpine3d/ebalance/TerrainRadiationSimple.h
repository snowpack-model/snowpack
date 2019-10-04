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
#ifndef TERRAINRADIATIONSIMPLE_H
#define TERRAINRADIATIONSIMPLE_H

#include <alpine3d/MPIControl.h>
#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/ViewFactorsAlgorithm.h>

/**
 * @class  TerrainRadiationSimple
 * @brief Simple guess of terrain reflected radiation.
 * For each cell of the domain, a sky view factor is computed over 32 sectors and a 500m distance (see for example
 * Faron S. Anslow, Steven Hostetler, William R. Bidlake, and Peter U. Clark, <i>"Distributed energy balance modeling
 * of South Cascade Glacier, Washington and assessment of model uncertainty"</i>, Journal of Geophysical Research, vol. 113, 2008).
 * Then, for each cell of the domain, this view factor is transformed into a terrain view factor, multiplied by the albedo of the current
 * cell and multiplied by the sum of the direct and diffuse radiation for the current cell. This is considered to be an approximation
 * of the short wave radiation rewflected by the surroundings of the current cell.
 *
 */
class TerrainRadiationSimple : public TerrainRadiationAlgorithm {

	public:
		TerrainRadiationSimple(const mio::DEMObject &dem_in, const std::string& method);
		~TerrainRadiationSimple();

		void getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain);
		void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta, const mio::Array2D<double>& rh,const mio::Array2D<double>& ilwr);
		void getSkyViewFactor(mio::Array2D<double> &o_sky_vf) const;

	private:
		double getAlbedo(const size_t& ii, const size_t& jj);
		void initSkyViewFactors(const mio::DEMObject &dem);

		mio::Array2D<double> albedo_grid, sky_vf;

		const size_t dimx, dimy;
		size_t startx, endx;
};

#endif
