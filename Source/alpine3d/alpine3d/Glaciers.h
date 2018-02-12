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
#ifndef GLACIERS_H
#define GLACIERS_H

#include <meteoio/MeteoIO.h>

/**
 * @page glaciers Glaciers
 * This module corrects the air temperature for the effect of the katabatic flows on the glaciated pixels.
 * The flows are computed automatically for each glaciated pixel following (Quinn, 1991) (see below).
 * Then the air temperature correction is computed according to (Greuell and Bohm, 1998) with the improvements 
 * by (Ayala, Pellicciotti and Shea, 2015) (see Glaciers::correctTemperatures).
 * This model is enabled by setting to ON the GLACIER_KATABATIC_FLOW key in the SNOWPACK section.
 * 
 * Then, two calibration parameters must be set: 
 *     - KATABATIC_LAYER_HEIGHT which represents the katabatic layer height in meters. It was around 17m in the original paper and 
 * around 5m in (Ayala et al.);
 *     - KATABATIC_K_COEFFICIENT which takes into account the lateral energy fluxes on the glacier tongue in (Ayala et al.). This was
 * around 7°C but can be set to 0 (default value) to use a pure (Greuell and Bohm, 1998) model.
 *     - KATABATIC_SCALING might also be provided. This factor reduces the katabatic flow reach. A value of .25 seems to work quite well
 *
 * @section multiple_flow
 * The computation of the katabatic flows paths is based on the multiple flow direction algorithm and
 * associates to each cell an upslope area that is based on the number of cells
 * that flow into it. This follows <i>"The prediction of hillslope flow paths for distributed
 * hydrological modelling using digital terrain models"</i>, Quinn P., Chevallier P., Planchon O.,
 * hydrological processes, <b>5</b>, 1991, pp 59-79.
 *
 * @remarks This code is still experimental and awaits proper validation!!
 */
class Glaciers
{
	public:
		Glaciers(const mio::Config& cfg);
		Glaciers(const mio::Config& cfg, const mio::DEMObject& in_dem);
		~Glaciers() {}

		void setDEM(const mio::DEMObject& in_dem);
		void setGlacierMap(const mio::Grid2DObject& glacierMask);
		const mio::Grid2DObject correctTemperatures(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, const mio::Grid2DObject& ta) const;
		void correctTemperatures(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, mio::Grid2DObject& ta) const;
		void correctTemperatures(mio::Grid2DObject& ta) const;

		void getGrids(mio::Grid2DObject &alt, mio::Grid2DObject &dist) const;

	private:
		void init(const mio::Config& cfg);
		static bool enableKatabatikFlows(const mio::Grid2DObject& hs, const mio::Grid2DObject& tss, const mio::Grid2DObject& ta, const mio::Grid2DObject& isGlacier);
		static bool hillslope_distribute_cell(const mio::Grid2DObject& dem, const mio::Grid2DObject& mask, const double& A, const size_t ii, const size_t jj, mio::Grid2DObject &flow, mio::Grid2DObject &src_altitude, mio::Grid2DObject &src_distance);
		void hillslope_flow(mio::Grid2DObject glacier_mask);

		mio::DEMObject dem;
		mio::Grid2DObject isGlacier; ///< glacier pixels tagged with 1, non glacier tagged with 0 and out of domain tagged with nodata
		mio::Grid2DObject flowpath, src_altitude, src_distance;
		double KBL, K, scale;
};

#endif
