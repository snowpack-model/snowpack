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
#ifndef SNOWPACKINTERFACE_H
#define SNOWPACKINTERFACE_H

#include <iostream>
#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>
#include <alpine3d/MeteoObj.h>

class SnowpackInterfaceWorker;
class SnowDriftA3D;
class Runoff; // forward declaration, cyclic header include

#include <alpine3d/DataAssimilation.h>
#include <alpine3d/runoff/Runoff.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/SnowpackInterfaceWorker.h>
#include <alpine3d/Glaciers.h>
#include <alpine3d/TechSnowA3D.h>

/**
 * @page snowpack Snowpack
 * This module call the Snowpack energy balance model on each relevant cell of the domain (see \ref principles_snowpack). The cells that are excluded
 * are the following:
 *     - cells those altitude is nodata (in the dem);
 *     - cells those land cover is undefined/nodata (in the lus);
 *     - cells marked as "water body".
 *
 * During a parallel simulation, it is often more efficient (but not always possible) to let each node perform its I/O on its local disk instead
 * of requesting the data from another single node. This is enabled with the LOCAL_IO key in the [General] section (default: true).
 *
 * Several corrections are applied on glaciated pixels: the albedo of such pixels is reset to a fixed glacier albedo when no snow is present
 * and when the air temperature is higher than 5°C the atmospheric stability is forced to MONIN_OBUKHOV. It is also possible to apply an air
 * temperature correction simulating the effect of katabatic flows by setting GLACIER_KATABATIC_FLOW to true in the [Snowpack] section (this is
 * still an experimental feature).
 *
 * Several types of outputs are supported: gridded outputs and full snowpack stratigraphy.
 *
 * @section gridded_outputs Gridded outputs
 * The gridded outputs are written out for all requested the parameters every GRIDS_DAYS_BETWEEN (in days), starting at GRIDS_START
 * (in days, <i>0</i> being noon). The available parameters are provided in SnGrids::Parameters and should be written as a space delimited
 * list of parameters as GRIDS_PARAMETERS key (in the [Output] section).
 *
 * It is possible to mask the glaciated areas (in order to keep a relevant scaling
 * on snow water equivalent or similar parameters) either statically (performed once and for all at the begining of the run) or dynamically
 * (performed at every time step). The criteria to decide if a pixel is glaciated is defined in Snowpack, in SnowStation::isGlacier(). This is
 * currently defined as more than 2 meters of pure ice in the whole snowpack. The keys to enable glaciers masking and dynamic masking are
 * MASK_GLACIERS and MASK_DYNAMIC in the [Output] section.
 *
 * Finally, the soil temperature at multiple given depths can be written out simply by setting the SOIL_TEMPERATURE_DEPTHS key in the [Output] section
 * to the chosen depths (in meters). In this case, there is no need to declare a TSOILx parameter for the GRIDS_PARAMETERS key (it is currently limited
 * to at most 5 different depths but could be increased in the future). It is possible to do the same in the snow at a given depth by setting
 * the SNOW_TEMPERATURE_DEPTH key, or the average snow temperature from the surface until a given depth (SNOW_AVG_TEMPERATURE_DEPTH key) or the snow density from
 * the surface until a given depth (SNOW_AVG_DENSITY_DEPTH key). When averaging either temperature or density, if the snow height is less
 * than the requested depth, the average is performed on the whole snow pack.
 *
 * @code
 * [Output]
 * GRID2D        = ARC
 * A3D_VIEW      = TRUE                ;naming scheme compatible with the viewer of Alpine3D
 * GRID2DPATH    = ../output
 * MASK_GLACIERS = TRUE
 * MASK_DYNAMIC  = FALSE
 * SOIL_TEMPERATURE_DEPTHS = 0.1 0.5 1.5       ;output soil temperatures at 0.1m, 0.5m and 1.5m - leave this key out to skip this output
 *
 * GRIDS_WRITE        = TRUE           ;default
 * GRIDS_START        = 0.
 * GRIDS_DAYS_BETWEEN = 1.
 * GRIDS_PARAMETERS = HS SWE ISWR ILWR
 * @endcode
 *
 * @subsection output_new_grids Writing new grids out
 * If the parameters you want to write out do not already exist in  SnGrids::Parameters, then you have a few extra steps to perform:
 *     -# add the parameter in SnGrids::Parameters as well as its string representation in SnGrids::initStaticData()
 *     -# add the proper code (similarly as for the already existing parameters) in SnowpackInterfaceWorker::fillGrids();
 *     -# after recompiling, it will be possible to specify the newly created grid as explained in \ref gridded_outputs.
 *
 * Several Snowpack objects are available for the current point (and documented in Snowpack itself in the Dataclasses.h file):
 *     - meteoPixel is an object of type CurrentMeteo;
 *     - snowPixel is an object of type SnowStation (this also contains CanopyData, NodeData and ElementData);
 *     - surfaceFlux is an object of type SurfaceFluxes;
 *
 *
 * @section poi_outputs Snow stratigraphy
 * Full snowpack stratigraphy outputs are provided at specific Points Of Interest, as defined with the POI key in the [Input] section. These outputs
 * contain time series of snow profiles and fluxes as well as meteorological forcing at these points, allowing to re-run these points manually in the
 * Snowpack model. These outputs are written in the path pointed to by METEOPATH in the [Ouput] section. The time resolution is controlled similarly to
 * standard Snowpack runs with the TS_WRITE and PROF_WRITE groups of keys (see Snowpack documentation).
 *
 * @code
 * [Input]
 * POI        = A3D
 * POIFILE    =  ../input/surface-grids/poi.pts
 *
 * [Output]
 * PROF_WRITE        = TRUE
 * PROFILE_FORMAT    = PRO
 * PROF_START        = 0.0
 * PROF_DAYS_BETWEEN = 1.
 *
 * TS_WRITE          = TRUE
 * TS_START          = 0.0
 * TS_DAYS_BETWEEN   = 1.
 *
 * METEO             = SMET
 * METEOPATH         = ../output
 * @endcode
 *
 *
 */
 class SnowpackInterface
  {
	public:
		// Methods for accessing Snopack Interface Manager
		SnowpackInterface(const mio::Config &io_cfg, const size_t& nbworkers,
		                  const mio::DEMObject& dem_in,
		                  const mio::Grid2DObject& landuse_in, const std::vector<mio::Coords>& vec_pts, const mio::Date& startTime, const std::string& grids_requirements, const bool is_restart_in);
		SnowpackInterface(const SnowpackInterface&);
		~SnowpackInterface();

		SnowpackInterface& operator=(const SnowpackInterface&); ///<Assignement operator, required because of pointer member

		double getTiming() const;
		void writeOutput(const mio::Date& julian);
		void writeOutputSNO(const mio::Date& julian);

		// Methods to set other modules
		void setSnowDrift(SnowDriftA3D& drift);
		void setEnergyBalance(EnergyBalance& myeb);
		void setDataAssimilation(DataAssimilation& init_da);
		void setRunoff(Runoff& init_runoff);

		// Methods to communicate with other modules
		void assimilate(const mio::Grid2DObject& daData, const mio::Date& timestamp);
		void setSnowMassChange(const mio::Grid2DObject& new_mns, const mio::Date& timestamp);
		void setMeteo(const mio::Grid2DObject& new_psum,
		                const mio::Grid2DObject& new_psum_ph,
		                const mio::Grid2DObject& new_vw,
		                const mio::Grid2DObject& new_dw,
		                const mio::Grid2DObject& new_rh,
		                const mio::Grid2DObject& new_ta,
		                const mio::Grid2DObject& new_tsg,
		                const mio::Date& timestamp);
		void setVwDrift(const mio::Grid2DObject& new_vw_drift,
		                const mio::Date& timestamp);
		void setRadiationComponents(const mio::Array2D<double>& shortw,
		                            const mio::Array2D<double>& longwave,
		                            const mio::Array2D<double>& diff,
		                            const double& solarElevation,
		                            const mio::Date& timestamp);

		mio::Grid2DObject getGrid(const SnGrids::Parameters& param) const;

	private:
		std::string getGridsRequirements() const;
		mio::Config readAndTweakConfig(const mio::Config& io_cfg,const bool have_pts);
		bool do_grid_output(const mio::Date &date) const;
		void calcNextStep();

		void readInitalSnowCover(std::vector<SnowStation*>& snow_stations,
                             std::vector<std::pair<size_t,size_t> >& snow_stations_coord);
		void readSnowCover(const std::string& GRID_sno, const std::string& LUS_sno, const bool& is_special_point,
											 SN_SNOWSOIL_DATA &sno, ZwischenData &zwischenData, const bool& read_seaice);
		void writeSnowCover(const mio::Date& date, const std::vector<SnowStation*>& snow_station);

		void write_SMET_header(const mio::StationData& meta, const double& landuse_code) const;
		void write_SMET(const CurrentMeteo& met, const mio::StationData& meta, const SurfaceFluxes& surf) const;
		void writeOutputSpecialPoints(const mio::Date& date, const std::vector<SnowStation*>& snow_pixel, const std::vector<CurrentMeteo*>& meteo_pixel,
		                              const std::vector<SurfaceFluxes*>& surface_flux);
		void write_special_points();
		void calcLateralFlow();
		void calcSimpleSnowDrift(const mio::Grid2DObject& ErodedMass, mio::Grid2DObject& psum);

		RunInfo run_info;

		// MeteoIO, used to output grids
		mio::IOManager io;

		std::vector< std::pair<size_t,size_t> > pts; //special points

		const mio::DEMObject dem;

		// Config dependent information
		bool is_restart, useCanopy, enable_simple_snow_drift, enable_lateral_flow, a3d_view;
		bool do_io_locally; // if false all I/O will only be done on the master process
		std::string station_name; // value for the key OUTPUT::EXPERIMENT
		bool glacier_katabatic_flow, snow_preparation;
		// Output
		std::vector<std::string> Tsoil_idx; //TSOIL names in order to build the "field" header of the smet output
		double grids_start, grids_days_between; //gridded outputs
		double ts_start, ts_days_between; //time series outputs
		double prof_start, prof_days_between; //profiles outputs
		bool grids_write, ts_write, prof_write, snow_write, snow_poi_written;
		std::string meteo_outpath;
		std::string outpath;
		bool mask_glaciers; //mask glaciers in outputs?
		bool mask_dynamic; //mask glaciers in outputs changes over time?
		mio::Grid2DObject maskGlacier; // save the mask
		double tz_out;

		SnowpackConfig sn_cfg;
		// SnowpackIO, used to output non grids data
		SnowpackIO snowpackIO;

		size_t dimx, dimy;
		size_t mpi_offset, mpi_nx;
		mio::Grid2DObject landuse;
		// meteo forcing variables
		mio::Grid2DObject mns, shortwave, longwave, diffuse;
		mio::Grid2DObject psum, psum_ph, psum_tech, grooming, vw, vw_drift, dw, rh, ta, tsg;
		mio::Grid2DObject winderosiondeposition;
		double solarElevation;

		std::vector<std::string> output_grids; //which grids should be written out
		std::vector<SnowpackInterfaceWorker*> workers;
		std::vector<size_t> worker_startx; // stores offset for each workers slice
		std::vector<size_t> worker_deltax; // stores size for each workers slize
		std::vector<std::vector<std::pair<size_t,size_t> > > worker_stations_coord; // ttores te grid coordiante of each worker
		// time relevant
		mio::Timer timer; // used to mesure calc time of one step
		mio::Date nextStepTimestamp;
		double timeStep; // size of timestep

		bool dataMeteo2D, dataDa, dataSnowDrift, dataRadiation; // say, if data are present for the actuall step

		// ==== other Modules ====
		SnowDriftA3D *drift;
		EnergyBalance *eb;
		DataAssimilation *da;
		Runoff *runoff;

		Glaciers *glaciers;
		TechSnowA3D *techSnow;
};

#endif
