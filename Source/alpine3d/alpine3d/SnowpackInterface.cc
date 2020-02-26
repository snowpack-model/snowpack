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
#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/AlpineMain.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/OMPControl.h>

#include <errno.h>
#include <algorithm>

using namespace std;
using namespace mio;

//sort the by increasing y and increasing x as a second key
inline bool pair_comparator(const std::pair<double, double>& l, const std::pair<double, double>& r)
{
	if (l.second == r.second)
		return l.first < r.first;

	return l.second < r.second;
}

//convert the POI to grid index representation and warn of duplicates
std::vector< std::pair<size_t,size_t> > prepare_pts(const std::vector<Coords>& vec_pts)
{
	std::vector< std::pair<size_t,size_t> > pts;
	std::vector<size_t> vec_idx;
	for (size_t ii=0; ii<vec_pts.size(); ii++) {
		const size_t ix = vec_pts[ii].getGridI();
		const size_t iy = vec_pts[ii].getGridJ();
		std::pair<size_t,size_t> tmp(ix,iy);

		const std::vector< std::pair<size_t,size_t> >::const_iterator it = std::find(pts.begin(), pts.end(), tmp);
		if (it != pts.end()) {
			if (MPIControl::instance().master()) {
				const size_t orig_idx = vec_idx[ it - pts.begin() ];
				std::cout << "[W] POI #" << ii << " " << vec_pts[ii].toString(Coords::CARTESIAN);
				std::cout << " is a duplicate of POI #" <<  orig_idx << " " << vec_pts[orig_idx].toString(Coords::CARTESIAN) << std::endl; //flush for openmp
			}
		} else {
			pts.push_back( tmp );
			vec_idx.push_back( ii ); //in order to be able to find out where a duplicate is originating from
		}
	}
	sort(pts.begin(), pts.end(), pair_comparator);
	return pts;
}

/**
 * @brief Constructs and initialise Snowpack Interface Master. He creates the Worker and
 * Distributes the Data from the other modules to the Worker. Is the acces interface A3D side
 * @param io_cfg is used to init Runoff and to create IOManager, which is used to write the standart output
 * @param nbworkers gives the new values for the wind speed
 * @param dem_in gives the demographic Data. Also tetermines size and position of the geographical modeling scope
 * @param landuse_in gives the landuse Data. Also tetermines size and position of the landuse for modeling scope
 * @param vec_pts gives the spezial points. For this points more output is done then for the others. Calcualtion is the same.
 * @param startTime is the time and date the first simulation step is done
 * @param grids_requirements list of grids that must be prepared for other modules (similar to the Output::GRIDS_PARAMETERS configuration key)
 * @param is_restart_in is used to know which files the worker needs to read to init the pixel
 */
SnowpackInterface::SnowpackInterface(const mio::Config& io_cfg, const size_t& nbworkers,
                                     const mio::DEMObject& dem_in,
                                     const mio::Grid2DObject& landuse_in,
                                     const std::vector<mio::Coords>& vec_pts,
                                     const mio::Date &startTime,
                                     const std::string& grids_requirements,
                                     const bool is_restart_in)
                : run_info(), io(io_cfg), pts(prepare_pts(vec_pts)),dem(dem_in),
                  is_restart(is_restart_in), useCanopy(false), enable_simple_snow_drift(false), enable_lateral_flow(false), a3d_view(false),
                  do_io_locally(true), station_name(),glacier_katabatic_flow(false), snow_preparation(false),
                  Tsoil_idx(), grids_start(0), grids_days_between(0), ts_start(0.), ts_days_between(0.), prof_start(0.), prof_days_between(0.),
                  grids_write(true), ts_write(false), prof_write(false), snow_write(false), snow_poi_written(false),
                  meteo_outpath(), outpath(), mask_glaciers(false), mask_dynamic(false), maskGlacier(), tz_out(0.),
                  sn_cfg(readAndTweakConfig(io_cfg,!pts.empty())), snowpackIO(sn_cfg), dimx(dem_in.getNx()), dimy(dem_in.getNy()), mpi_offset(0), mpi_nx(dimx),
                  landuse(landuse_in), mns(dem_in, IOUtils::nodata), shortwave(dem_in, IOUtils::nodata), longwave(dem_in, IOUtils::nodata), diffuse(dem_in, IOUtils::nodata),
                  psum(dem_in, IOUtils::nodata), psum_ph(dem_in, IOUtils::nodata), psum_tech(dem_in, IOUtils::nodata), grooming(dem_in, IOUtils::nodata),
                  vw(dem_in, IOUtils::nodata), vw_drift(dem_in, IOUtils::nodata), dw(dem_in, IOUtils::nodata), rh(dem_in, IOUtils::nodata),
                  ta(dem_in, IOUtils::nodata), tsg(dem_in, IOUtils::nodata), winderosiondeposition(dem_in, 0),
                  solarElevation(0.), output_grids(), workers(nbworkers), worker_startx(nbworkers), worker_deltax(nbworkers), worker_stations_coord(nbworkers),
                  timer(), nextStepTimestamp(startTime), timeStep(dt_main/86400.), dataMeteo2D(false), dataDa(false), dataSnowDrift(false), dataRadiation(false),
                  drift(NULL), eb(NULL), da(NULL), runoff(NULL), glaciers(NULL), techSnow(NULL)
{
	MPIControl& mpicontrol = MPIControl::instance();

	std::vector<SnowStation*> snow_stations;
	std::vector<std::pair<size_t,size_t> > snow_stations_coord;
	readInitalSnowCover(snow_stations,snow_stations_coord);

	if (mpicontrol.master()) {
		std::cout << "[i] SnowpackInterface initializing a total of " << mpicontrol.size();
		if (mpicontrol.size()>1) std::cout << " processes with " << nbworkers;
		else std::cout << " process with " << nbworkers;
		if (nbworkers>1) std::cout << " workers";
		else  std::cout << " worker";
		std::cout << " each using Snowpack " << snowpack::getLibVersion() << "\n";
	}

	//check if simple snow drift is enabled (needs to be determined before grid requirements check!)
	enable_simple_snow_drift = false;
	sn_cfg.getValue("SIMPLE_SNOW_DRIFT", "Alpine3D", enable_simple_snow_drift, IOUtils::nothrow);

	//create and prepare  the vector of output grids
	if (grids_write) {
		sn_cfg.getValue("GRIDS_PARAMETERS", "output", output_grids);
		std::vector<double> soil_temp_depths;
		sn_cfg.getValue("SOIL_TEMPERATURE_DEPTHS", "Output", soil_temp_depths, IOUtils::nothrow);
		const unsigned short max_Tsoil( SnGrids::lastparam - SnGrids::TSOIL1 + 1 );
		if (soil_temp_depths.size()>max_Tsoil)
			throw InvalidArgumentException("Too many soil temperatures requested", AT);
		for (size_t ii=0; ii<soil_temp_depths.size(); ii++) {
			std::stringstream ss;
			ss << (ii+1);
			const std::string ii_str(ss.str());
			Tsoil_idx.push_back( ii_str );
			output_grids.push_back( "TSOIL"+ii_str );
		}
		SnowpackInterfaceWorker::uniqueOutputGrids(output_grids);
	}
	//add the grids that are necessary for the other modules
	const std::string all_grids = sn_cfg.get("GRIDS_PARAMETERS", "output", "");
	sn_cfg.addKey("GRIDS_PARAMETERS", "output", all_grids + " " + grids_requirements + " " + getGridsRequirements()); //also consider own requirements

	//check if lateral flow is enabled
	sn_cfg.getValue("LATERAL_FLOW", "Alpine3D", enable_lateral_flow, IOUtils::nothrow);

	//check if A3D viez should be used for grids
	sn_cfg.getValue("A3d_VIEW", "Output", a3d_view, IOUtils::nothrow);

	//If MPI is active, every node gets a slice of the DEM to work on
	mpicontrol.getArraySliceParamsOptim(dimx, mpi_offset, mpi_nx,dem,landuse);
	std::cout << "[i] MPI instance "<< mpicontrol.rank() <<" for solving snowpack : grid range = ["
	<< mpi_offset << " to " << mpi_offset+mpi_nx-1 << "] " << mpi_nx << " columns\n";
	//Cut DEM and landuse in MPI domain, MPI domain are computed
	//by trying to havve the same number of cell to compute per domain
	const DEMObject mpi_sub_dem(dem_in, mpi_offset, 0, mpi_nx, dimy, false);
	const Grid2DObject mpi_sub_landuse(landuse_in, mpi_offset, 0, mpi_nx, dimy);
	std::vector<std::vector<size_t> > omp_snow_stations_ind;
	//The OMP slicing for snowpack computation is not based on rectangular domain.
	//It return a table of pixel to compute and coordiantes to have the same
	//number of pixel per slice.
	OMPControl::getArraySliceParamsOptim(nbworkers,snow_stations,mpi_sub_dem, mpi_sub_landuse,omp_snow_stations_ind);
	// construct slices and workers
	#pragma omp parallel for schedule(static)
	for (size_t ii=0; ii<nbworkers; ii++) {
		// Could be optimised, but not really big gain.
		// Each worker will check again the point to be sure they belong
		// to it, so no need to double check here
		std::vector< std::pair<size_t,size_t> > sub_pts;
		const size_t n_pts = pts.size();
		for (size_t kk=0; kk<n_pts; kk++) { // could be optimised... but not really big gain
			if (pts[kk].first >= mpi_offset && pts[kk].first < mpi_offset + mpi_nx) {
				sub_pts.push_back( pts[kk] );
				sub_pts.back().first -= mpi_offset;
			}
		}
		// In this  implementation, all OMP workers see the whole MPI grid
		// This could be change but with this when passing the meteo grids
		// Only one slicing and copy is necessary per MPI, insetead of one per OMP
		const DEMObject sub_dem(mpi_sub_dem);
		const Grid2DObject sub_landuse(mpi_sub_landuse);

		// Generate workers
		std::vector<SnowStation*> thread_stations;
		std::vector<std::pair<size_t,size_t> > thread_stations_coord;
		for (std::vector<size_t>::iterator it = omp_snow_stations_ind.at(ii).begin(); it != omp_snow_stations_ind.at(ii).end(); ++it){
			thread_stations.push_back (snow_stations.at(*it));
			thread_stations_coord.push_back(snow_stations_coord.at(*it));
		}

		// The OMP slicing into rectangle is only used for post computation
		// Over the grid (i.e. lateral flow and snow preparation)
		size_t omp_offset, omp_nx;
		OMPControl::getArraySliceParams(mpi_nx, nbworkers, ii, omp_offset, omp_nx);
		const size_t offset = mpi_offset + omp_offset;

		workers[ii] = new SnowpackInterfaceWorker(sn_cfg, sub_dem, sub_landuse, sub_pts, thread_stations, thread_stations_coord, offset);

		worker_startx[ii] = offset;
		worker_deltax[ii] = omp_nx;
		#pragma omp critical(snowpackWorkers_status)
		{
			const std::pair<size_t,size_t> coord_start(snow_stations_coord.at(omp_snow_stations_ind.at(ii).front()));
			const std::pair<size_t,size_t> coord_end(snow_stations_coord.at(omp_snow_stations_ind.at(ii).back()));

			std::cout << "[i] SnowpackInterface worker for solving snowpack " << ii
								<< " on process " << mpicontrol.rank() << ": coord. x-y range = [" << coord_start.first + mpi_offset << "-"
								<< coord_start.second << " to " << coord_end.first + mpi_offset << "-"
								<< coord_end.second<<"] " << omp_snow_stations_ind.at(ii).size() << " cells\n";
								std::cout << "[i] SnowpackInterface worker for grid computation " << ii
								<< " on process " << mpicontrol.rank() << ": grid range = [" << offset
								<< " to " << omp_nx+offset-1 << "]  " << omp_nx << " columns\n";
		}
	}

	// init glacier map (after creating and init workers) for output
	if (mask_glaciers || glacier_katabatic_flow) {
		maskGlacier = getGrid(SnGrids::GLACIER);
		if (glacier_katabatic_flow) {
			glaciers = new Glaciers(io_cfg, dem);
			glaciers->setGlacierMap(maskGlacier);
		}
	}

	//init snow preparation
	if (snow_preparation) {
		techSnow = new TechSnowA3D(io_cfg, dem);
	}
}

SnowpackInterface& SnowpackInterface::operator=(const SnowpackInterface& source) {
	if (this != &source) {
		run_info = source.run_info;
		snowpackIO = source.snowpackIO;
		dimx = source.dimx;
		dimy = source.dimy;
		landuse = source.landuse;
		mns = source.mns;
		shortwave = source.shortwave;
		longwave = source.longwave;
		diffuse = source.diffuse;
		psum = source.psum;
		psum_ph = source.psum_ph;
		psum_tech = source.psum_tech;
		grooming = source.grooming;
		vw = source.vw;
		vw_drift = source.vw_drift;
		dw = source.dw;
		rh = source.rh;
		ta = source.ta;
		tsg = source.tsg;
		solarElevation = source.solarElevation;
		output_grids = source.output_grids;
		workers = source.workers;
		worker_startx = source.worker_startx;
		worker_deltax = source.worker_deltax;
		timer = source.timer;
		nextStepTimestamp = source.nextStepTimestamp;
		timeStep = source.timeStep;

		drift = source.drift;
		eb = source.eb;
		da = source.da;
		runoff = source.runoff;
		dataMeteo2D = source.dataMeteo2D;
		dataDa = source.dataDa;
		dataSnowDrift = source.dataSnowDrift;
		dataRadiation = source.dataRadiation;

		//io = source.io;

		outpath = source.outpath;
		mask_glaciers = source.mask_glaciers;
		mask_dynamic = source.mask_dynamic;
		maskGlacier = source.maskGlacier;

		glacier_katabatic_flow = source.glacier_katabatic_flow;
		snow_preparation = source.snow_preparation;
		glaciers = source.glaciers;
		techSnow = source.techSnow;
		enable_lateral_flow = source.enable_lateral_flow;
		a3d_view = source.a3d_view;

		sn_cfg = source.sn_cfg;
		//dem = source.dem;
		is_restart = source.is_restart;
		useCanopy = source.useCanopy;
		do_io_locally = source.do_io_locally;
		station_name = source.station_name;

		Tsoil_idx = source.Tsoil_idx;
		grids_start = source.grids_start;
		grids_days_between = source.grids_days_between;
		ts_start = source.ts_start;
		ts_days_between = source.ts_days_between;
		prof_start = source.prof_start;
		prof_days_between = source.prof_days_between;
		grids_write = source.grids_write;
		ts_write = source.ts_write;
		prof_write = source.prof_write;
		snow_write = source.snow_write;
		snow_poi_written = source.snow_poi_written;
		meteo_outpath = source.meteo_outpath;
		tz_out = source.tz_out;
		pts = source.pts;
	}
	return *this;
}

std::string SnowpackInterface::getGridsRequirements() const
{
	std::string ret = "";
	if (glacier_katabatic_flow) {
		ret += " GLACIER TSS HS";
	}
	if (enable_simple_snow_drift) {
		 ret += " ERODEDMASS";
	}
	return ret;
}

mio::Config SnowpackInterface::readAndTweakConfig(const mio::Config& io_cfg, const bool have_pts)
{
	SnowpackConfig tmp_cfg(io_cfg);
	//force some keys
	double calculation_step_length;
	tmp_cfg.getValue("CALCULATION_STEP_LENGTH", "Snowpack", calculation_step_length);
	std::stringstream ss;
	ss << calculation_step_length;
	tmp_cfg.addKey("METEO_STEP_LENGTH", "Snowpack", ss.str());
	tmp_cfg.addKey("ALPINE3D", "SnowpackAdvanced", "true");
	tmp_cfg.addKey("ALPINE3D_PTS", "SnowpackAdvanced",have_pts?"true":"false");
	tmp_cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "true");
	tmp_cfg.getValue("LOCAL_IO", "General", do_io_locally, IOUtils::nothrow);
	tmp_cfg.getValue("GRID2DPATH", "Output", outpath);
	tmp_cfg.getValue("MASK_GLACIERS", "Output", mask_glaciers, IOUtils::nothrow);
	tmp_cfg.getValue("MASK_DYNAMIC", "Output", mask_dynamic, IOUtils::nothrow);
	tmp_cfg.getValue("GLACIER_KATABATIC_FLOW", "Snowpack", glacier_katabatic_flow, IOUtils::nothrow);
	tmp_cfg.getValue("SNOW_PREPARATION", "Input", snow_preparation, IOUtils::nothrow);
	tmp_cfg.getValue("GRIDS_WRITE", "Output", grids_write);
	tmp_cfg.getValue("GRIDS_START", "Output", grids_start);
	tmp_cfg.getValue("GRIDS_DAYS_BETWEEN", "Output", grids_days_between);
	tmp_cfg.getValue("TS_WRITE", "Output", ts_write);
	tmp_cfg.getValue("TS_START", "Output", ts_start);
	tmp_cfg.getValue("TS_DAYS_BETWEEN", "Output", ts_days_between);
	tmp_cfg.getValue("PROF_WRITE", "Output", prof_write);
	tmp_cfg.getValue("PROF_START", "Output", prof_start);
	tmp_cfg.getValue("PROF_DAYS_BETWEEN", "Output", prof_days_between);

	tmp_cfg.getValue("METEOPATH", "Output", meteo_outpath);
	tmp_cfg.getValue("TIME_ZONE", "Output", tz_out, IOUtils::nothrow);
	tmp_cfg.getValue("EXPERIMENT", "Output", station_name, IOUtils::dothrow);

	tmp_cfg.getValue("SNOW_WRITE", "Output", snow_write);
	tmp_cfg.getValue("CANOPY", "Snowpack", useCanopy);

	return tmp_cfg;
}

/**
 * @brief Destructor of SnowpackInterface Master. Handels special cases with POP-C++ and
 * also free correctly runoff and workers.
 */
SnowpackInterface::~SnowpackInterface()
{
	if (glacier_katabatic_flow) delete glaciers;
	//if (runoff) delete runoff;
	while (!workers.empty()) delete workers.back(), workers.pop_back();
}


/**
 * @brief get time who was used to exchange Data with Workers and run on each Pixel
 * the Snowpack Model throught workers.
 */
double SnowpackInterface::getTiming() const
{
	return timer.getElapsed();
}


/**
 * @brief Internal in Snowpack Interface Master used method to write standard
 * results in output files.
 * Attantion: to have old format files output, set in ini file following key:
 * **A3D_VIEW	= true**
 *
 * @param date is the date for which the output is done
 */
void SnowpackInterface::writeOutput(const mio::Date& date)
{
	MPIControl& mpicontrol = MPIControl::instance();
	const bool isMaster = mpicontrol.master();

	if (do_grid_output(date)) {
		//no OpenMP pragma here, otherwise multiple threads might call an MPI allreduce_sum()
		for (size_t ii=0; ii<output_grids.size(); ii++) {
			const size_t SnGrids_idx = SnGrids::getParameterIndex( output_grids[ii] );
			mio::Grid2DObject grid( getGrid( static_cast<SnGrids::Parameters>(SnGrids_idx)) );

			if (isMaster) {
				if (mask_glaciers) grid *= maskGlacier;
				const size_t meteoGrids_idx = MeteoGrids::getParameterIndex( output_grids[ii] );
				if (meteoGrids_idx!=IOUtils::npos) { //for this, the grid plugins must be thread-safe!
					io.write2DGrid(grid, static_cast<MeteoGrids::Parameters>(meteoGrids_idx), date);
				} else {
					std::string grid_type;
					sn_cfg.getValue("GRID2D", "output",grid_type);
					if (grid_type == "ARC") {
						std::string file_name;
						if (a3d_view) {
							std::string dateStr( date.toString(Date::NUM) );
							dateStr.erase( dateStr.size()-2, string::npos); //remove the seconds
							file_name =  dateStr + "_" + output_grids[ii] + ".asc" ;
						} else {
							std::string date_str(date.toString(mio::Date::ISO));
							std::replace( date_str.begin(), date_str.end(), ':', '.');
							file_name =  date_str + "_" + output_grids[ii] + ".asc" ;
						}
						io.write2DGrid(grid, file_name);
					} else if (grid_type=="NETCDF") {
						std::string file;
						sn_cfg.getValue("GRID2DFILE", "output", file);
						io.write2DGrid(grid, output_grids[ii]+"@"+date.toString(mio::Date::ISO));
					} else {
						throw InvalidFormatException("[E] Only ARC and NetCDF allow writing out non-standard grids such as "+output_grids[ii], AT);
					}

					// Reset WINDEROSIONDEPOSITION, which is cumulative since the previous output grid
					if (output_grids[ii] == "WINDEROSIONDEPOSITION") winderosiondeposition.set(winderosiondeposition, 0.);
				}
			}
		}
	}
	// Output Runoff: at each time step
	if (runoff) runoff->output(date, psum, ta);
}

/**
 * @brief Method tells if on given date, gridded output should be done (read this out of snowpack ini)
 * @param date is the date object which is controlled, if output needs to be done
 */
bool SnowpackInterface::do_grid_output(const mio::Date &date) const
{
	return (grids_write && booleanTime(date.getJulian(), grids_days_between, grids_start, dt_main/60.));
}

/**
 * @brief commands worker to write .sno files. Is triggered by Alpine Control
 *
 * Hack --> Find better software architecture to do this then SnowpackInterface also
 * does output for other modules of A3D here..
 * @param date is the date witch the output is done
 */
void SnowpackInterface::writeOutputSNO(const mio::Date& date)
{
	MPIControl& mpicontrol = MPIControl::instance();

	vector<SnowStation*> snow_station;

	for (size_t ii=0; ii<workers.size(); ii++)
		workers[ii]->getOutputSNO(snow_station);

	if (mpicontrol.master()) {
		std::cout << "[i] Writing SNO output for process " << mpicontrol.master_rank() << "\n";
		writeSnowCover(date, snow_station); //local data

		//Now gather all elements on the master node
		for (size_t ii=0; ii<mpicontrol.size(); ii++) {
			if (ii == mpicontrol.master_rank() || do_io_locally) continue;
			std::cout << "[i] Writing SNO output for process " << ii << "\n";
			vector<SnowStation*> snow_station_tmp;

			mpicontrol.receive(snow_station_tmp, ii);
			writeSnowCover(date, snow_station_tmp);
			for (size_t jj=0; jj<snow_station_tmp.size(); jj++) delete snow_station_tmp[jj];
		}
	} else {
		if (do_io_locally) {
			std::cout << "[i] Writing SNO output for process " << mpicontrol.rank() << "\n";
			writeSnowCover(date, snow_station); //local data
		} else {
			mpicontrol.send(snow_station, mpicontrol.master_rank());
		}
	}
}

void SnowpackInterface::writeSnowCover(const mio::Date& date, const std::vector<SnowStation*>& snow_station)
{
	for (size_t jj=0; jj<snow_station.size(); jj++)
		snowpackIO.writeSnowCover(date, *(snow_station[jj]), ZwischenData());
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Methods to set references to other methodes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/**
 * @brief Set reference to SnowDrift module, to comunicate with it.
 * @param mydrift is the reference to the SnowDrift module
 */
void SnowpackInterface::setSnowDrift(SnowDriftA3D& mydrift)
{
	drift = &mydrift;

	if (drift) {
		for (size_t i = 0; i < workers.size(); i++) workers[i]->setUseDrift(true);

		// Provide initial snow parameters to SnowDrift
		const Grid2DObject cH( getGrid(SnGrids::HS) );
		const Grid2DObject sp( getGrid(SnGrids::SP) );
		const Grid2DObject rg( getGrid(SnGrids::RG) );
		const Grid2DObject N3( getGrid(SnGrids::N3) );
		const Grid2DObject rb( getGrid(SnGrids::RB) );
		drift->setSnowSurfaceData(cH, sp, rg, N3, rb);
	}
}

/**
 * @brief Set reference to EnergyBalance module, to comunicate with it.
 * @param myeb is the reference to the EnergyBalance module
 */
void SnowpackInterface::setEnergyBalance(EnergyBalance& myeb)
{

	eb = &myeb;
	if (eb) {
		for (size_t i = 0; i < workers.size(); i++){
			workers[i]->setUseEBalance(true);
		}
		// Provide initial albedo to EnergyBalance
		const Grid2DObject alb( getGrid(SnGrids::TOP_ALB) );
		eb->setAlbedo(alb);
	}
}

/**
 * @brief Set reference to DataAssimilation module, to comunicate with it.
 * @param init_da is the reference to the DataAssimilation module
 */
void SnowpackInterface::setDataAssimilation(DataAssimilation& init_da)
{
	da = &init_da;
}

void SnowpackInterface::setRunoff(Runoff& init_runoff)
{
	runoff = &init_runoff;
}

/**
 * @brief Interface that DataAssimilation can push the data to the SnowpackInterface
 * This is currently never used...
 * @param daData are the data from the DataAssimilation module
 * @param timestamp for controlle if data from DataAssimilation are form the correct simulation step
 */
void SnowpackInterface::assimilate(const Grid2DObject& /*daData*/, const mio::Date& timestamp)
{

	if (nextStepTimestamp != timestamp) {
		throw InvalidArgumentException("Assimilation and snowpack time steps don't match", AT);
	}

	cout << "Updating state variables...\n";
	/*for (size_t iy = 0; iy < dimy; iy++) {
		for (size_t ix = 0; ix < dimx; ix++) {
			const size_t i = dimy*iy + ix;
			if ( daData.grid2D(ix,iy) == 1.0 ) {
				//if DA-data DOES NOT have snow
				if (sn_Xdata[i].cH-sn_Xdata[i].Ground > 0.) {
					//and snowpack DOES have snow
					//clear the pixel
					store(ix,iy) = 0.;
					sn_Xdata[i].cH = sn_Xdata[i].Ground;
					sn_Xdata[i].mH = sn_Xdata[i].Ground;
					sn_Xdata[i].nElems = sn_Xdata[i].SoilNode;
					sn_Xdata[i].nNodes = sn_Xdata[i].nElems + 1;
				}
			} else {
				if ( daData.grid2D(ix,iy) == 4 ) {
					//if DA-data DOES have snow
					if (sn_Xdata[i].cH-sn_Xdata[i].Ground < 1e-15) {
						//if snowpack DOES NOT have snow -> add some snow
						store(ix,iy) += 40*M_TO_H(calculation_step_length);
					}
				}
			}
		}
	}*/

	dataDa = true;
	calcNextStep();
}

/**
 * @brief Interface that SnowDrift can push the data to the SnowpackInterface
 * @param new_mns are the data about the new Snow masses from the DataAssimilation module
 * @param timestamp for controlle if data from DataAssimilation are form the correct simulation step
 */
void SnowpackInterface::setSnowMassChange(const mio::Grid2DObject& new_mns, const mio::Date& timestamp)
{
	if (nextStepTimestamp != timestamp) {
		if (MPIControl::instance().master()) {
			std::cerr << "Providing drift snow mass at " << timestamp.toString(Date::ISO);
			std::cerr << " for Snowpack timestamp " << nextStepTimestamp.toString(Date::ISO) << "\n";
		}
		throw InvalidArgumentException("Snowdrift and snowpack steps don't match", AT);
	}

	if (!new_mns.isSameGeolocalization(dem)) {
		std::ostringstream ss;
		ss << "Trying to set snow mass changes from a (" << new_mns.getNx() << "," << new_mns.getNy() << ") grid ";
		ss << "when the dem is (" << dem.getNx() << "," << dem.getNy() << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}

	mns = new_mns;
	dataSnowDrift = true;
	calcNextStep();
}

/**
 * @brief get Meteo changes from AlpineControl or SnowDrift module
 * @param new_psum gives the new values for new Water High
 * @param new_psum_ph gives the new values for the precipitation phase
 * @param new_vw gives the new values for the wind speed
 * @param new_dw gives the new values for the wind direction
 * @param new_rh gives the new values for the realtiv Humidity
 * @param new_ta gives the new values for the Aire temperature
 * @param timestamp is the time of the calculation step from which this new values are comming
 */
void SnowpackInterface::setMeteo(const Grid2DObject& new_psum, const Grid2DObject& new_psum_ph, const Grid2DObject& new_vw, const Grid2DObject& new_dw, const Grid2DObject& new_rh, const Grid2DObject& new_ta, const Grid2DObject& new_tsg, const mio::Date& timestamp)
{
	if (nextStepTimestamp != timestamp) {
		if (MPIControl::instance().master()) {
			std::cerr << "Providing meteo fields at " << timestamp.toString(Date::ISO);
			std::cerr << " for Snowpack timestamp " << nextStepTimestamp.toString(Date::ISO) << "\n";
		}
		throw InvalidArgumentException("Meteo and snowpack time steps don't match", AT);
	}

	psum = new_psum;
	psum_ph = new_psum_ph;
	vw = new_vw;
	dw = new_dw;
	rh = new_rh;
	if (mask_dynamic) maskGlacier = getGrid(SnGrids::GLACIER); //so the updated glacier map is available for all

	if (!glacier_katabatic_flow) {
		ta = new_ta;
	} else {
		if (mask_dynamic) glaciers->setGlacierMap(maskGlacier);
		const Grid2DObject TSS( getGrid(SnGrids::TSS) );
		const Grid2DObject cH( getGrid(SnGrids::HS) );
		ta = glaciers->correctTemperatures(cH, TSS, new_ta);
	}
	tsg = new_tsg;

	if (snow_preparation) {
		const Grid2DObject cH( getGrid(SnGrids::HS) );
		techSnow->setMeteo(new_ta, new_rh, cH, timestamp);
		psum_tech = techSnow->getGrid(SnGrids::PSUM_TECH);
		grooming = techSnow->getGrid(SnGrids::GROOMING);
	}

	dataMeteo2D = true;
	calcNextStep();
}

void SnowpackInterface::setVwDrift(const Grid2DObject& new_vw_drift, const mio::Date& timestamp)
{
	if (nextStepTimestamp != timestamp) {
		if (MPIControl::instance().master()) {
			std::cerr << "Providing meteo fields at " << timestamp.toString(Date::ISO);
			std::cerr << " for Snowpack timestamp " << nextStepTimestamp.toString(Date::ISO) << "\n";
		}
		throw InvalidArgumentException("Meteo and snowpack time steps don't match", AT);
	}

	vw_drift = new_vw_drift;
}

/**
 * @brief get values from Energy Balance
 * @param shortwave_in is the 2D double Array map of ISWR [W m-2]
 * @param longwave_in is the 2D double Array map of ILWR [W m-2]
 * @param diff_in is the 2D double Array map of Diffuse radiation from the sky [W m-2]
 * @param solarElevation_in double of Solar elevation to be used for Canopy (in dec)
 * @param timestamp is the time of the calculation step from which this new values are comming
 */
void SnowpackInterface::setRadiationComponents(const mio::Array2D<double>& shortwave_in, const mio::Array2D<double>& longwave_in, const mio::Array2D<double>& diff_in, const double& solarElevation_in, const mio::Date& timestamp)
{
	if (nextStepTimestamp != timestamp) {
		if (MPIControl::instance().master()) {
			std::cerr << "Providing radiation fields at " << timestamp.toString(Date::ISO);
			std::cerr << " for Snowpack timestamp " << nextStepTimestamp.toString(Date::ISO) << "\n";
		}
		throw InvalidArgumentException("Radiation and snowpack time steps don't match", AT);
	}

	shortwave.grid2D = shortwave_in;
	longwave.grid2D = longwave_in;
	diffuse.grid2D = diff_in;
	solarElevation = solarElevation_in;

	dataRadiation = true;
	calcNextStep();
}

/**
 * @brief Request specific grid by parameter type
 * @param param parameter
 * @return 2D output grid (empty if the requested parameter was not available)
 */
mio::Grid2DObject SnowpackInterface::getGrid(const SnGrids::Parameters& param) const
{
	//special case for the meteo forcing grids
	switch (param) {
		case SnGrids::TA:
			return ta;
		case SnGrids::RH:
			return rh;
		case SnGrids::VW:
			return vw;
		case SnGrids::VW_DRIFT:
			return vw_drift;
		case SnGrids::DW:
			return dw;
		case SnGrids::PSUM:
			return psum;
		case SnGrids::PSUM_PH:
			return psum_ph;
		case SnGrids::PSUM_TECH:
			return psum_tech;
		case SnGrids::ISWR:
			return shortwave;
		case SnGrids::ILWR:
			return longwave;
		case SnGrids::WINDEROSIONDEPOSITION:
			return winderosiondeposition;
		default: ; //so compilers do not complain about missing conditions
	}

	mio::Grid2DObject o_grid2D(dem, 0.); //so the allreduce_sum works
	size_t errCount = 0;
	mio::Grid2DObject tmp_grid2D_(dem,mpi_offset, 0, mpi_nx, dimy);
	mio::Grid2DObject tmp_grid2D(tmp_grid2D_,mio::IOUtils::nodata);
	#pragma omp parallel for schedule(dynamic) reduction(+: errCount)
	for (size_t ii = 0; ii < workers.size(); ii++) {
		const mio::Grid2DObject tmp( workers[ii]->getGrid(param) );
		if (!tmp.empty()) {
			for (size_t i=0; i<tmp_grid2D.getNx(); ++i){
				for (size_t j=0; j<tmp_grid2D.getNy(); ++j){
					if (tmp(i,j)!=mio::IOUtils::nodata)
					{
							tmp_grid2D(i,j)=tmp(i,j);
					}
				}
			}
		} else {
			errCount++;
		}
	}
	o_grid2D.grid2D.fill(tmp_grid2D.grid2D, mpi_offset, 0, mpi_nx, dimy);

	MPIControl::instance().allreduce_sum(o_grid2D);
	//with some MPI implementations, when transfering large amounts of data, the buffers might get full and lead to a crash

	if (errCount>0) {
		std::cerr << "[W] Requested " << SnGrids::getParameterName( param ) << " but this was not available in the workers\n";
		o_grid2D.clear(); //the requested parameter was not available
	}
	return o_grid2D;
}

/**
 * @brief Get data from other modules and run one simulation step.
 * Once the simulation step has been performed, the data are pushed to ther other modules.
 */
void SnowpackInterface::calcNextStep()
{
	//Control if all data are present
	if (!dataMeteo2D) {
		return;
	} else {
		if (drift != NULL && !dataSnowDrift) return;
		if (da != NULL && !dataDa) return;
		if (eb != NULL && !dataRadiation) return;
	}
	// control if the necessary data are available
	if (!dataRadiation) {
		throw NoDataException("Radiation data not available", AT);
	}

	dataDa = dataMeteo2D = dataSnowDrift = dataRadiation = false; //the external modules will turn them back to true when pushing their data

	// timing
	timer.restart();

	if (enable_simple_snow_drift) {
		// calc simple snow drift, by using eroded snow from previous time step
		const Grid2DObject erodedmass( getGrid(SnGrids::ERODEDMASS) );
		calcSimpleSnowDrift(erodedmass, psum);
	}

	size_t errCount = 0;
	const mio::Grid2DObject tmp_psum(psum, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_psum_ph(psum_ph,  mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_psum_tech(psum_tech,  mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_rh(rh, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_ta(ta, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_tsg(tsg, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_vw(vw, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_vw_drift(vw_drift, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_dw(dw, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_mns(mns, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_shortwave(shortwave, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_diffuse(diffuse, mpi_offset, 0, mpi_nx, dimy);
	const mio::Grid2DObject tmp_longwave(longwave, mpi_offset, 0, mpi_nx, dimy);
	#pragma omp parallel for schedule(dynamic, 1) reduction(+: errCount)
	for (size_t ii = 0; ii < workers.size(); ii++) { // make slices
		// run model, process exceptions in a way that is compatible with openmp
		try {
			workers[ii]->runModel(nextStepTimestamp, tmp_psum, tmp_psum_ph, tmp_psum_tech, tmp_rh, tmp_ta, tmp_tsg, tmp_vw, tmp_vw_drift, tmp_dw, tmp_mns, tmp_shortwave, tmp_diffuse, tmp_longwave, solarElevation);
			if (snow_preparation) {
				const mio::Grid2DObject tmp_grooming(grooming, mpi_offset, 0, mpi_nx, dimy);
				workers[ii]->grooming( nextStepTimestamp, tmp_grooming );
			}
		} catch(const std::exception& e) {
			++errCount;
			cout << e.what() << std::endl;
		}
	}

	//Lateral flow
	if (enable_lateral_flow) {
		calcLateralFlow();
	}

	//Retrieve special points data and write files
	if (!pts.empty()) write_special_points();

	if (errCount>0) {
		//something wrong took place, quitting. At least we tried writing the special points out
		std::abort(); //force core dump
	}

	// Gather data if needed and make exchange for SnowDrift
	if (drift) {
		const Grid2DObject cH( getGrid(SnGrids::HS) );
		const Grid2DObject sp( getGrid(SnGrids::SP) );
		const Grid2DObject rg( getGrid(SnGrids::RG) );
		const Grid2DObject N3( getGrid(SnGrids::N3) );
		const Grid2DObject rb( getGrid(SnGrids::RB) );
		drift->setSnowSurfaceData(cH, sp, rg, N3, rb);
	}

	// Gather data if needed and make exchange for EnergyBalance
	if (eb) {
		const Grid2DObject alb( getGrid(SnGrids::TOP_ALB) );
		eb->setAlbedo(alb);
	}

	//make output
	writeOutput(nextStepTimestamp);

	timer.stop();
	if (MPIControl::instance().master())
		cout << "[i] Snowpack simulations done for " << nextStepTimestamp.toString(Date::ISO) << "\n";
	nextStepTimestamp = nextStepTimestamp + timeStep;
	MPIControl::instance().barrier();

}

void SnowpackInterface::write_special_points()
{
	MPIControl& mpicontrol = MPIControl::instance();

	std::vector<SnowStation*> snow_pixel;
	std::vector<CurrentMeteo*> meteo_pixel;
	std::vector<SurfaceFluxes*> surface_flux;

	// note: do not parallelize this with OpenMP
	for (size_t ii=0; ii<workers.size(); ii++)
		workers[ii]->getOutputSpecialPoints(snow_pixel, meteo_pixel, surface_flux);

	if (do_io_locally) {
		writeOutputSpecialPoints(nextStepTimestamp, snow_pixel, meteo_pixel, surface_flux);
		if (!snow_write && !snow_poi_written) {
			writeSnowCover(nextStepTimestamp, snow_pixel); //also write the .sno of the special points
			snow_poi_written = true;
		}
	} else { // data has to be sent to the master process
		if (mpicontrol.master()) {
			// Write out local data first and then gather data from all processes
			writeOutputSpecialPoints(nextStepTimestamp, snow_pixel, meteo_pixel, surface_flux);
			if (!snow_write && !snow_poi_written) {
				writeSnowCover(nextStepTimestamp, snow_pixel); //also write the .sno of the special points
			}

			for (size_t ii=0; ii<mpicontrol.size(); ii++) {
				if (ii == mpicontrol.master_rank()) continue;
				snow_pixel.clear(); meteo_pixel.clear(); surface_flux.clear();

				mpicontrol.receive(snow_pixel, ii);
				mpicontrol.receive(meteo_pixel, ii);
				mpicontrol.receive(surface_flux, ii);

				writeOutputSpecialPoints(nextStepTimestamp, snow_pixel, meteo_pixel, surface_flux);
				if (!snow_write && !snow_poi_written) {
					writeSnowCover(nextStepTimestamp, snow_pixel); //also write the .sno of the special points
				}
			}
			snow_poi_written = true;
		} else {
			mpicontrol.send(snow_pixel, mpicontrol.master_rank());
			mpicontrol.send(meteo_pixel, mpicontrol.master_rank());
			mpicontrol.send(surface_flux, mpicontrol.master_rank());
		}
	}

	#pragma omp parallel for schedule(static)
	for (size_t ii=0; ii<workers.size(); ii++) workers[ii]->clearSpecialPointsData();
}

/**
 * @brief Write the output which is asked to have more for the special points
 * @param date Output date
 * @param snow_pixel The SnowStation data for all the special points
 * @param meteoPixel The CurrentMeteo data for all the special points
 * @param surface_flux The SurfaceFlux data for all the special points
 */
void SnowpackInterface::writeOutputSpecialPoints(const mio::Date& date, const std::vector<SnowStation*>& snow_pixel, const std::vector<CurrentMeteo*>& meteo_pixel,
                                                 const std::vector<SurfaceFluxes*>& surface_flux)
{
	const bool TS = (ts_write && booleanTime(date.getJulian(), ts_days_between, ts_start, dt_main/60.));
	const bool PR = (prof_write && booleanTime(date.getJulian(), prof_days_between, prof_days_between, dt_main/60));

	const ProcessDat Hdata; // empty ProcessDat, get it from where ??
	for (size_t ii=0; ii<snow_pixel.size(); ii++) {
		write_SMET(*meteo_pixel[ii], snow_pixel[ii]->meta, *surface_flux[ii]);
		if (TS) snowpackIO.writeTimeSeries(*snow_pixel[ii], *surface_flux[ii], *meteo_pixel[ii], Hdata, 0.);
		if (PR) snowpackIO.writeProfile(date, *snow_pixel[ii]);
	}
}

/**
 * @brief Write header of SMET file for specific point
 * @param meta StationData for the SMET header initialization
 */
void SnowpackInterface::write_SMET_header(const mio::StationData& meta, const double& landuse_code) const
{
	const std::string filename( meteo_outpath + "/" + meta.stationName + "_meteo.smet" );
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename,AT);
	errno = 0;

	std::ofstream smet_out; //Output file streams
	smet_out.open(filename.c_str());
	if (smet_out.fail()) throw AccessException(filename.c_str(), AT);

	smet_out << "SMET 1.1 ASCII\n";
	smet_out << "[HEADER]\n";
	smet_out << "station_name = " << meta.stationName << "\n";
	smet_out << "station_id   = " << meta.stationID << "\n";
	smet_out << std::right;
	smet_out << std::fixed;

	smet_out << "altitude     = " << std::setw(11)  << std::setprecision(1) << meta.position.getAltitude() << "\n";
	smet_out << "latitude     = " << std::setw(11) << std::setprecision(8) << meta.position.getLat() << "\n";
	smet_out << "longitude    = " << std::setw(11) << std::setprecision(8) << meta.position.getLon() << "\n";
	smet_out << "easting      = " << std::setw(11) << std::setprecision(1) << meta.position.getEasting() << "\n";
	smet_out << "northing     = " << std::setw(11) << std::setprecision(1) << meta.position.getNorthing() << "\n";
	smet_out << "epsg         = " << std::setw(11) << std::setprecision(0) << meta.position.getEPSG() << "\n";
	smet_out << "slope        = " << std::setw(11)  << std::setprecision(1) << meta.getSlopeAngle() << "\n";
	smet_out << "azimuth      = " << std::setw(11)  << std::setprecision(1) << meta.getAzimuth() << "\n";
	smet_out << "landuse      = " << std::setw(11)  << std::setprecision(0) << SnowpackInterfaceWorker::round_landuse( landuse_code ) << "\n";
	smet_out << "nodata       = " << std::setw(11)  << std::setprecision(0) << mio::IOUtils::nodata << "\n";
	smet_out << "tz           = " << std::setw(11)  << std::setprecision(0) << tz_out << "\n";
	smet_out << "source       = " <<  "Alpine3D version " << A3D_VERSION << " run by " << run_info.user << "\n";
	smet_out << "creation     = " << run_info.computation_date.toString(Date::ISO) << "\n";
	if (useCanopy) smet_out << "comment      = " << "ISWR/RSWR are above the canopy, ISWR_can/RSWR_can and PSUM/PSUM_PH are below the canopy\n";

	smet_out << "fields       = timestamp TA TSS TSG VW DW VW_MAX ISWR OSWR ILWR PSUM PSUM_PH HS RH";
	for (size_t ii=0; ii<Tsoil_idx.size(); ii++)
		smet_out << " TSOIL" << Tsoil_idx[ii];

	if (useCanopy) smet_out << " ISWR_can RSWR_can";
	smet_out << "\n[DATA]\n";

	smet_out.close();
}

/**
 * @brief Write the meteorogical data for the current step into the SMET file for the respective point
 * @param met The CurrentMeteo data for one special point
 * @param ix is the x-coordiante of the special point
 * @param iy is the y-coordiante of the special point
 */
void SnowpackInterface::write_SMET(const CurrentMeteo& met, const mio::StationData& meta, const SurfaceFluxes& surf) const
{
	const std::string filename( meteo_outpath + "/" + meta.stationName + "_meteo.smet" );
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename,AT);

	std::ofstream smet_out; //Output file streams
	smet_out.open(filename.c_str(), std::ios::out | std::ios::app );
	if (smet_out.fail()) throw AccessException(filename.c_str(), AT);

	// write line
	smet_out.fill(' ');
	smet_out << std::right;
	smet_out << std::fixed;
	smet_out << met.date.toString(mio::Date::ISO) << " ";
	smet_out << std::setw(8) << std::setprecision(2) << met.ta << " ";
	smet_out << std::setw(8) << std::setprecision(2) << met.tss << " ";
	smet_out << std::setw(8) << std::setprecision(2) << met.ts0 << " ";
	smet_out << std::setw(6) << std::setprecision(1) << met.vw << " ";
	smet_out << std::setw(5) << std::setprecision(0) << met.dw << " ";
	smet_out << std::setw(6) << std::setprecision(1) << met.vw_max << " ";
	smet_out << std::setw(6) << std::setprecision(0) << met.iswr << " ";
	smet_out << std::setw(6) << std::setprecision(0) << met.rswr << " ";
	smet_out << std::setw(6) << std::setprecision(3) << mio::Atmosphere::blkBody_Radiation(met.ea, met.ta) << " ";
	smet_out << std::setw(6) << std::setprecision(3) << met.psum << " ";
	smet_out << std::setw(6) << std::setprecision(3) << met.psum_ph << " ";
	smet_out << std::setw(8) << std::setprecision(3) << met.hs / cos(meta.getSlopeAngle()*Cst::to_rad) << " ";
	smet_out << std::setw(7) << std::setprecision(3) << met.rh << " ";
	if (!met.ts.empty())
		smet_out << std::setw(8) << std::setprecision(2) << met.ts[0] << " ";
	if (useCanopy) {
		smet_out << std::setw(6) << std::setprecision(0) << surf.sw_in << " ";
		smet_out << std::setw(6) << std::setprecision(0) << surf.sw_out << " ";
	}
	smet_out << "\n";

	smet_out.close();
}

/**
 * @page reading_snow_files Reading initial snow cover
 * The initial snow cover consist of an instantaneous snow/soil profile from which the time evolution will be computed.
 * When this is for a normal "cold" start, the file names are built based on the landuse code. For restarts, the
 * file names are built based on the cell (ii,jj) indices, for example:
 * 	+ {station_name}_{landuse_code}.{ext} for a "cold" start;
 * 	+ {ii}_{jj}_{station_name}.{ext} for a restart;
 *
 * The station name is given in the [Output] section as "EXPERIMENT" key. The other keys controlling the process (including
 * the file extension) are:
 * 	+ in the [Snowpack] section:
 * 		+ CANOPY: should the pixels enable the canopy module?
 * 		+ SNP_SOIL: should the pixels use soil layers?
 * 		+ in the [Input] section:
 * 		+ SNOW: file format of the "sno" files, either SMET or SNOOLD (default: SMET);
 * 		+ COORDSYS, COORDPARAM: in order to convert (ii,jj) coordinates to geographic coordinates so each pixel's metadata
 * can be reused (for example in order to rerun a \ref poi_outputs "Point Of Interest" offline in the SNOWPACK standalone model).
 */
 void SnowpackInterface::readInitalSnowCover(std::vector<SnowStation*>& snow_stations,
                                             std::vector<std::pair<size_t,size_t> >& snow_stations_coord){
  //HACK: with nextStepTimestamp, check that the snow cover is older than the start timestep!

	if (MPIControl::instance().master() || do_io_locally) {
		const bool useSoil = sn_cfg.get("SNP_SOIL", "Snowpack");
		const std::string coordsys = sn_cfg.get("COORDSYS", "Input");
		const std::string coordparam = sn_cfg.get("COORDPARAM", "Input", "");
		Coords llcorner_out( dem.llcorner );
		llcorner_out.setProj(coordsys, coordparam);
		const double refX = llcorner_out.getEasting();
		const double refY = llcorner_out.getNorthing();
		const double cellsize = dem.cellsize;

		const size_t nrWorkers = MPIControl::instance().size();
		for (size_t ii=0; ii<nrWorkers; ii++) {
			if (do_io_locally && (ii != MPIControl::instance().rank())) continue; // only read/write points managed by this process

			SN_SNOWSOIL_DATA snow_soil;
			size_t startx, deltax;
			MPIControl::instance().getArraySliceParamsOptim(dimx, ii, startx, deltax,dem,landuse);
			vector<SnowStation*> snow_stations_tmp;
			vector<pair<size_t,size_t> > snow_stations_coord_tmp;
			snow_stations_tmp.reserve( dimy*deltax );

			// read snow cover for all points that are dealt with on this process
			for (size_t iy = 0; iy < dimy; iy++) {
				for (size_t ix = startx; ix < (startx+deltax); ix++) {
					snow_stations_coord_tmp.push_back(std::pair<size_t,size_t>(ix - startx,iy));
					if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) { //skip nodata cells as well as water bodies, etc
						snow_stations_tmp.push_back( NULL );
						continue;
					}
					snow_stations_tmp.push_back( new SnowStation(useCanopy, useSoil) );

					SnowStation& snowPixel = *(snow_stations_tmp.back());
					const bool is_special_point = SnowpackInterfaceWorker::is_special(pts, ix, iy);

					// get potential filenames for initial snow pixel values
					std::stringstream LUS_sno, GRID_sno;
					LUS_sno << station_name << "_" << SnowpackInterfaceWorker::round_landuse(landuse.grid2D(ix,iy));
					GRID_sno << ix << "_" << iy << "_" << station_name;

					// read standard values of pixel
					try {
						ZwischenData zwischenData; //not used by Alpine3D but necessary for Snowpack
						readSnowCover(GRID_sno.str(), LUS_sno.str(), is_special_point, snow_soil, zwischenData, (snowPixel.Seaice!=NULL));
					} catch (exception& e) {
						cout << e.what()<<"\n";
						throw IOException("Can not read snow files", AT);
					}

					// Copy standard values to specific pixel (station) data and init it
					try {
						snowPixel.initialize(snow_soil, 0); //force sector 0
						snowPixel.mH = Constants::undefined;
					} catch (exception&) {
						cout << "[E] Could not intialize cell (" << ix << "," << iy << ")!\n";
						throw IOException("Can not initialize snow pixel", AT);
					}

					// Set proper pixel metadata
					snowPixel.meta.position.setProj(coordsys, coordparam);
					snowPixel.meta.position.setXY(refX+double(ix)*cellsize, refY+double(iy)*cellsize, dem.grid2D(ix,iy));
					snowPixel.meta.position.setGridIndex((int)ix, (int)iy, 0, true);
					// Note that pixels without valid azimuth will be asigned azimuth = 0 and slope angle = 0.
					snowPixel.meta.setSlope( (dem.azi(ix,iy) == mio::IOUtils::nodata) ? (0.) : (dem.slope(ix,iy)),
							         (dem.azi(ix,iy) == mio::IOUtils::nodata) ? (0.) : (dem.azi(ix,iy)) );
					snowPixel.cos_sl = cos( snowPixel.meta.getSlopeAngle()*mio::Cst::to_rad );

					// Initialize the station name for the pixel
					stringstream station_idx;
					station_idx << ix << "_" << iy;
					snowPixel.meta.stationName = station_idx.str() + "_" + station_name;
					snowPixel.meta.stationID = station_idx.str();
					if (is_special_point) { //create SMET files for special points
						write_SMET_header(snowPixel.meta, landuse(ix, iy));
					}
				}
			}
			if (ii == MPIControl::instance().master_rank() || do_io_locally) {
				snow_stations = snow_stations_tmp; //simply copy the pointers
				snow_stations_coord = snow_stations_coord_tmp;
			} else {
				MPIControl::instance().send(snow_stations_tmp, ii);
				MPIControl::instance().send(snow_stations_coord_tmp, ii);
				while (!snow_stations_tmp.empty()) delete snow_stations_tmp.back(), snow_stations_tmp.pop_back();
			}
		}
		std::cout << "[i] Read initial snow cover for process " << MPIControl::instance().rank() << "\n";
	} else {
		MPIControl::instance().receive(snow_stations, MPIControl::instance().master_rank());
		MPIControl::instance().receive(snow_stations_coord, MPIControl::instance().master_rank());
	}
}

void SnowpackInterface::readSnowCover(const std::string& GRID_sno, const std::string& LUS_sno, const bool& is_special_point,
																			SN_SNOWSOIL_DATA &sno, ZwischenData &zwischenData, const bool& read_seaice)
{
	// read standard values of pixel
	if (is_special_point && !is_restart) {
		//special points can come either from LUS snow files or GRID snow files
		if (snowpackIO.snowCoverExists(GRID_sno, station_name)) {
			snowpackIO.readSnowCover(GRID_sno, station_name, sno, zwischenData, read_seaice);
		} else {
			snowpackIO.readSnowCover(LUS_sno, station_name, sno, zwischenData, read_seaice);
		}
	} else {
		if (is_restart) {
			snowpackIO.readSnowCover(GRID_sno, station_name, sno, zwischenData, read_seaice);
		} else {
			snowpackIO.readSnowCover(LUS_sno, station_name, sno, zwischenData, read_seaice);
		}
	}

	//check that the layers are older than the start date
	if (sno.nLayers>0 && sno.Ldata.front().depositionDate>nextStepTimestamp) {
		ostringstream ss;
		ss <<  "A layer can not be younger than the start date!";
		if (snowpackIO.snowCoverExists(GRID_sno, station_name))
			ss << " Please check profile '" << GRID_sno << "'";
		else
			ss << " Please check profile '" << LUS_sno << "'";
		throw IOException(ss.str(), AT);
	}
}

/**
 * @brief Calculates lateral flow
 * @author Nander Wever
 */
void SnowpackInterface::calcLateralFlow()
{
	std::vector<SnowStation*> snow_pixel;
	// Retrieve snow stations
	for (size_t ii = 0; ii < workers.size(); ii++) {
		workers[ii]->getLateralFlow(snow_pixel);
	}
	// Translate lateral flow in source/sink term
	size_t ix=0;											// The source cell x coordinate
	int ixd=-1, iyd=-1;										// The destination cell x and y coordinates
	size_t errCount = 0;
	#pragma omp parallel for schedule(dynamic, 1) reduction(+: errCount)
	for (size_t ii = 0; ii < workers.size(); ii++) {						// Cycle over all workers
		try {
			for (size_t jj = 0; jj < worker_deltax[ii]; jj++) {				// Cycle over x range per worker
				for (size_t iy = 0; iy < dimy; iy++) {					// Cycle over y
					ix = worker_startx[ii] + jj;
					const size_t index_SnowStation_src = ix * dimy + iy;		// Index of source cell of water
					double tmp_dist = -1;						// Cell distance
					if (snow_pixel[index_SnowStation_src] != NULL) {		// Make sure it is not a NULL pointer (in case of skipped cells)
						// Now determine destination cell for the water, based on azimuth
						if ((snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 337.5 && snow_pixel[index_SnowStation_src]->meta.getAzimuth() <= 360.) || (snow_pixel[index_SnowStation_src]->meta.getAzimuth() >=0. && snow_pixel[index_SnowStation_src]->meta.getAzimuth() <= 22.5)) {
							ixd = ix;
							iyd = iy-1;
							tmp_dist = dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 22.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 67.5) {
							ixd = ix+1;
							iyd = iy-1;
							tmp_dist = mio::Cst::Sqrt2 * dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 67.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 112.5) {
							ixd = ix+1;
							iyd = iy;
							tmp_dist = dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 112.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 157.5) {
							ixd = ix+1;
							iyd = iy+1;
							tmp_dist = mio::Cst::Sqrt2 * dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 157.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 202.5) {
							ixd = ix;
							iyd = iy+1;
							tmp_dist = dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 202.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 247.5) {
							ixd = ix-1;
							iyd = iy+1;
							tmp_dist = mio::Cst::Sqrt2 * dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 247.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 292.5) {
							ixd = ix-1;
							iyd = iy;
							tmp_dist = dem.cellsize;
						} else if (snow_pixel[index_SnowStation_src]->meta.getAzimuth() > 292.5 || snow_pixel[index_SnowStation_src]->meta.getAzimuth() < 337.5) {
							ixd = ix-1;
							iyd = iy-1;
							tmp_dist = mio::Cst::Sqrt2 * dem.cellsize;
						} else {
							// Undefined aspect, don't route the lateral flow by setting destination cell outside domain
							ixd=-1;
							iyd=-1;
						}
						if (ixd >= 0 && iyd >= 0 && ixd < int(dimx) && iyd < int(dimy)) {								// Check if destination cell is inside domain
							const size_t index_SnowStation_dst = dimy * ixd + iyd;							// Destination cell of water
							if (snow_pixel[index_SnowStation_dst] != NULL) {							// Make sure destination cell is not a NULL pointer (in case of skipped cells)
								for (size_t n=0; n < snow_pixel[index_SnowStation_src]->getNumberOfElements(); n++) {			// Loop over all layers in source cell
									for (size_t nn=0; nn < snow_pixel[index_SnowStation_dst]->getNumberOfElements(); nn++) {	// Loop over all layers in destination cell
										// Now check if deposition date is equal or newer (i.e., never put lateral water in an older layer, except when we cycled over all layers and we are at the top layer)
										if (snow_pixel[index_SnowStation_dst]->Edata[nn].depositionDate >= snow_pixel[index_SnowStation_src]->Edata[n].depositionDate
											|| nn == snow_pixel[index_SnowStation_dst]->getNumberOfElements()-1) {
											// The flux into the pixel is a source term for the destination cell
											snow_pixel[index_SnowStation_dst]->Edata[nn].lwc_source += snow_pixel[index_SnowStation_src]->Edata[n].SlopeParFlux / tmp_dist * (snow_pixel[index_SnowStation_dst]->Edata[nn].L / snow_pixel[index_SnowStation_src]->Edata[n].L);
											// The flux out of the pixel is a sink term for the source cell
											snow_pixel[index_SnowStation_src]->Edata[n].lwc_source -= snow_pixel[index_SnowStation_src]->Edata[n].SlopeParFlux / tmp_dist * (snow_pixel[index_SnowStation_dst]->Edata[nn].L / snow_pixel[index_SnowStation_src]->Edata[n].L);
											// Set the SlopeParFlux to zero, now that we have redistributed it.
											snow_pixel[index_SnowStation_src]->Edata[n].SlopeParFlux = 0.;
											break;
										}
									}
								}
							}
						}
					}
				}
				ix++;
			}
		} catch(const std::exception& e) {
			++errCount;
			cout << e.what() << std::endl;
		}
	}

	if (errCount>0) {
		//something wrong took place, quitting. At least we tried writing the special points out
		std::abort(); //force core dump
	}

	// Send back SnowStations to the workers
	#pragma omp parallel for schedule(dynamic, 1) reduction(+: errCount)
	for (size_t ii = 0; ii < workers.size(); ii++) {
		std::vector<SnowStation*> snow_pixel_out;					// Construct vector to send snow stations to the associated worker
		for (size_t jj = 0; jj < worker_deltax[ii]; jj++) {				// Cycle over x range per worker
			for (size_t iy = 0; iy < dimy; iy++) {					// Cycle over y
				ix = worker_startx[ii] + jj;
				const size_t idx = ix * dimy + iy;				// Index of snow pixel
				snow_pixel_out.push_back( snow_pixel[idx] );
			}
		}
		workers[ii]->setLateralFlow(snow_pixel_out);
	}
	return;
}


/**
 * @brief Get data from other modules and run one simulation step.
 * Once the simulation step has been performed, the data are pushed to ther other modules.
 */
void SnowpackInterface::calcSimpleSnowDrift(const mio::Grid2DObject& tmp_ErodedMass, mio::Grid2DObject& tmp_psum)
{
	double sum_positive_exposure = 0.;
	double sum_erodedmass = 0.;

	const double max_sx = -vw_drift.grid2D.getMax(); //positive
	for (size_t ii=0; ii<vw_drift.size(); ii++) {
		if(vw_drift(ii) < 0.) {
			sum_positive_exposure += -vw_drift(ii) / max_sx;
		}
		sum_erodedmass += tmp_ErodedMass(ii);
	}
	if (sum_positive_exposure == 0) return;

	const double ratio = sum_erodedmass / sum_positive_exposure;
	for (size_t iy=0; iy<dimy; iy++) {
		for (size_t ix=0; ix<dimx; ix++) {
			double deposition = 0.;
			if(vw_drift(ix, iy) < 0.) {
				const double val = -vw_drift(ix, iy) * ratio / max_sx;
				deposition = val;
				tmp_psum(ix, iy) += val;
			} else {
				deposition = -tmp_ErodedMass(ix,iy);
			}
			winderosiondeposition(ix, iy) = deposition;
		}
	}
	return;
}
