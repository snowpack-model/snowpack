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
#ifndef RUNOFF_H
#define RUNOFF_H

#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <sstream>
#include <meteoio/MeteoIO.h>
#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/SnowpackInterfaceWorker.h>

typedef unsigned long long longuint;

class SnowpackInterface; // forward declaration, cyclic header include

/**
 * @page runoff Runoff
 * This module computes the pure runoff from each cell or group of cells or a hydrological discharge (see \ref principles_runoff).
 * It is also possible to perform hydrological modeling with the PREVAH system forced with the runoff grids produced by Alpine3D
 * (see \ref prevah_runoff "PREVAH hydrological modeling").
 *
 * @section runoff_grid Runoff grid
 * In order to use the Alpine3D generated runoff as input for a hydrological model, it is often necessary to write out the generated
 * runoff at each cell and timesteps. This is done by setting the WRITE_RUNOFF_GRIDS key to true in the [Output] section. The grids will be
 * in the ARC ascii format by default but this can be configured with the RUNOFF_GRID2D key. The path where to write the grids <b>must</b> be
 * provided thanks to the RUNOFF_GRID2DPATH key. For example, to write per-pixel generated runoff grids in the ARC ascii format:
 *
 * @code
 * [Output]
 * WRITE_RUNOFF_GRIDS = true
 * RUNOFF_GRID2D = ARC
 * RUNOFF_GRID2DPATH = ../output/runoff
 * @endcode
 *
 * @section runoff_sums Sub-catchments runoff sums
 * In order to work with multiple catchments, it is possible to write out the sum of all the runoff generated within a given sub-catchment.
 * The sub-catchments must be properly defined before as shown in the \ref sub_catch_input "sub-catchments inputs" section. Then these sums
 * will be written out in a SMET file that contains both the total runoff (from the bottom of the snow/soil column)
 * as well as the precipitation, snow melt and glacier melt components (these ones at the bottom of the snow column only, ie from the surface runoff).
 * Please note that the total runoff is at the bottom of the snow/soil column while the splitting (precipitation, snowmelt and glacier melt) are at the snow/soil interface
 * (ie surface runoff), therefore the sum of the components is smaller than the total runoff (ie some of the runoff could not be classified).
 *
 * In order to write the runoff sums at each time step, it is necessary to define the following keys:
 *    - CATCHMENT, in the [Input] section that gives the path and filename of the catchment definition grid;
 *    - CATCHMENTS_PATH, in the [Ouput] section that gives the path where to write the sub-catchments sums (one SMET file per sub-catchment)
 *    - CATCHMENT_NUMBERING (optional), in the [Input] section, which specifies which numbering scheme is used by the catchment definition grid
 *      to identify the individual sub-catchments. Two values are possible for this key: ALPINE3D_OLD (if the catchments are numbered with series
 *      of powers of 2, as described in the \ref sub_catch_input "sub-catchments inputs" section), or TAUDEM (if the catchments are identified
 *      with unique numbers, which are not required to be continuous)
 *
 * On top of the runoff components, additional variables can optionally be averaged over each sub-catchment area and written in the output files
 * at each time step. These variables must be specified using key RUNOFF_FILES_EXTRA_DATA, which accepts the standard meteorological variables
 * (TA, RH, VW, ILWR, etc.) along with T_SOIL, which corresponds to the soil temperature at a certain depth. In case T_SOIL is present, the depth
 * (in meters) at which soil temperature should be returned must be specified using key SOIL_TEMPERATURE_DEPTH:
 *
 * @code
 * [Output]
 * RUNOFF_FILES_EXTRA_DATA = TA RH ISWR T_SOIL
 * SOIL_TEMPERATURE_DEPTH  = 5 ; meters
 * @endcode
 *
 */

class Runoff
{
	public:
		Runoff(const mio::Config& in_cfg, const mio::DEMObject& in_dem,
				const double& in_thresh_rain);
		Runoff(const Runoff& copy);
		
		// Methods to set other modules
		void setSnowPack(SnowpackInterface &sn_interface);
		
		virtual void output(const mio::Date& i_date, const mio::Grid2DObject& psum,
				const mio::Grid2DObject& ta);
		virtual ~Runoff();
		
		std::string getGridsRequirements() const;
		double getTiming() const;

	protected:
		mio::IOManager *io;
		SnowpackInterface *snowpack;  ///< Reference to SnowpackInterface object, used for callbacks, initialized during construction

		const double thresh_rain;
		const double tz_out;
		bool output_grids; //< output grid of runoff data? (for taking it over with an external hydro model)
		bool output_sums;  //< output total runoff per subcatchment?
		std::string catchment_out_path; //< folder in which the extracted catchment sums should be written
		double resampling_cell_size; //< [m] cell size to which the runoff grids should be resampled
		double grid_size_factor;
		const size_t n_grid_cells;
		enum NumberingType {Alpine3DOld, TauDEM} const catchment_numbering;

		mio::Timer timer;
		const mio::Grid2DObject slope_correction;
		const bool is_glacier_mask_dynamic;
		bool is_glacier_mask_set;
		mio::Grid2DObject total_runoff, glacier_mask;
		std::vector<SnGrids::Parameters> extra_meteo_variables;
		size_t n_extra_meteo_variables;
		std::map<size_t, mio::Grid2DObject> catchment_masks;

		static const double MIN_CELL_SIZE; //< [m] two points closer to each other than this value will be assumed to overlap
		static const double DISTANCE_ABSOLUTE_PRECISION; //< [m] minimum size of a grid cell

		virtual void constructCatchmentMasks(mio::Grid2DObject catchmentGrid);
		virtual void updateTotalRunoffGrid();
		virtual void updateGlacierMask();
		virtual mio::Grid2DObject computePrecipRunoff(const mio::Grid2DObject& psum,
				const mio::Grid2DObject& ta) const;
		virtual void getExtraMeteoGrids(std::vector<mio::Grid2DObject>& grids) const;
		virtual void initializeOutputFiles(const mio::Grid2DObject& dem) const;
		virtual void updateOutputFile(const size_t& catchId, const mio::Date& currTime,
				const double& totalRunoff, const double& precipRunoff,
				const double& snowRunoff, const double& glacierRunoff,
				const std::vector<double>& meteoVars) const;

		static mio::IOManager* getIOManager(mio::Config in_cfg, const bool& outputGrids);
		static mio::Grid2DObject getSlopeCorrection(const mio::DEMObject& dem);
		static NumberingType getCatchmentNumberingScheme(const mio::Config& in_cfg);
		static bool getIsGlacierDynamic(const mio::Config& cfg);
		static std::vector<SnGrids::Parameters> getExtraMeteoVariables(const mio::Config& cfg);
		static double getResamplingCellSize(const mio::DEMObject& in_dem,
				const mio::Grid2DObject& catchments);
		static bool isMultiple(const double& a, const double& b);
		static double estimateResamplingCellSize(const double& llOffset,
				const double& currSizeEstimate);
		static std::vector<size_t> factorizeCatchmentNumber(longuint value);
		static void cropMask(mio::Grid2DObject& mask);
		static double sumOverMask(const mio::Grid2DObject& grid, const mio::Grid2DObject& mask);
		static double averageOverMask(const mio::Grid2DObject& grid, const mio::Grid2DObject& mask);

	private:
		Runoff& operator=(const Runoff&) {return *this;} //< private in order to avoid being used and suppress compiler warning
};

#endif
