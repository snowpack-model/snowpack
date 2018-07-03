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
#include <alpine3d/runoff/Runoff.h>

const double Runoff::DISTANCE_ABSOLUTE_PRECISION = 1e-3; //< [m]
const double Runoff::MIN_CELL_SIZE = 5.0; //< [m]

/**
 * @brief Constructor of Runoff instance.
 * @param in_cfg reference to the config object used by SnowpackInterfaceWorker in order to build our own IOManager in Runoff
 * @param in_dem reference to the DEM used by Alpine3D
 * @param in_thresh_rain rain/snow temperature threshold. This will be used to correctly split the output between melt/precipitation
 */
Runoff::Runoff(const mio::Config& in_cfg, const mio::DEMObject& in_dem, const double& in_thresh_rain) :
		io(NULL), snowpack(NULL), thresh_rain(in_thresh_rain),
		tz_out(in_cfg.get("TIME_ZONE", "Output")),  output_grids(false),
		output_sums(false), catchment_out_path(), resampling_cell_size(),
		grid_size_factor(), n_grid_cells(in_dem.getNx()*in_dem.getNy()),
		catchment_numbering(getCatchmentNumberingScheme(in_cfg)),
		timer(), slope_correction(getSlopeCorrection(in_dem)),
		is_glacier_mask_dynamic(getIsGlacierDynamic(in_cfg)),
		is_glacier_mask_set(false), total_runoff(), glacier_mask(),
		extra_meteo_variables(getExtraMeteoVariables(in_cfg)),
		n_extra_meteo_variables(extra_meteo_variables.size()), catchment_masks()
{
	in_cfg.getValue("WRITE_RUNOFF_GRIDS", "OUTPUT", output_grids, mio::IOUtils::nothrow);

	std::string catchmentInFile;
	in_cfg.getValue("CATCHMENT", "INPUT", catchmentInFile, mio::IOUtils::nothrow);
	output_sums = !catchmentInFile.empty();

	if (output_grids || output_sums) {
		io = getIOManager(in_cfg, output_grids);
		total_runoff.set(in_dem, mio::IOUtils::nodata);
	}

	if (output_sums) {
		in_cfg.getValue("CATCHMENTS_PATH", "OUTPUT", catchment_out_path);
		mio::Grid2DObject catchmentGrid;
		io->read2DGrid(catchmentGrid, catchmentInFile);
		resampling_cell_size = getResamplingCellSize(in_dem, catchmentGrid);
		grid_size_factor = in_dem.cellsize/resampling_cell_size;
		constructCatchmentMasks(catchmentGrid);
		initializeOutputFiles(in_dem);
		glacier_mask.set(in_dem, mio::IOUtils::nodata);
	}
	if (MPIControl::instance().master()) {
		std::cout << "[i] Runoff initialised";
		if (output_sums) std::cout << " - Number of catchments: " << catchment_masks.size() << "\n";
		else std::cout << "\n";
	}
}

/*
 * @brief Copy constructor
 * @param copy Runoff object to copy
 */
Runoff::Runoff(const Runoff& copy) :
		io(new mio::IOManager(*(copy.io))), snowpack(copy.snowpack),
		thresh_rain(copy.thresh_rain), tz_out(copy.tz_out),
		output_grids(copy.output_grids), output_sums(copy.output_sums),
		catchment_out_path(copy.catchment_out_path),
		resampling_cell_size(copy.resampling_cell_size),
		grid_size_factor(copy.grid_size_factor), n_grid_cells(copy.n_grid_cells),
		catchment_numbering(copy.catchment_numbering),
		timer(copy.timer), slope_correction(copy.slope_correction),
		is_glacier_mask_dynamic(copy.is_glacier_mask_dynamic),
		is_glacier_mask_set(copy.is_glacier_mask_set),
		total_runoff(copy.total_runoff), glacier_mask(copy.glacier_mask),
		extra_meteo_variables(copy.extra_meteo_variables),
		n_extra_meteo_variables(copy.n_extra_meteo_variables),
		catchment_masks(copy.catchment_masks)
{}

std::string Runoff::getGridsRequirements() const
{
	return "PSUM TA MS_SOIL_RUNOFF GLACIER MS_SNOWPACK_RUNOFF";
}


/**
 * @brief Sets the internal reference to SnowpackInterface object
 * @param sn_interface Reference to the SnowpackInterface object
 */
void Runoff::setSnowPack(SnowpackInterface &sn_interface) {
	snowpack = &sn_interface;
}


/**
 * @brief Writes the results for a specific day
 * @param i_date the date at which the results should be written
 * @param psum grid of precipitation in mm/h
 * @param ta grid of air temperature
 */
void Runoff::output(const mio::Date& i_date, const mio::Grid2DObject& psum, const mio::Grid2DObject& ta)
{
	if (!output_grids && !output_sums) return;

	timer.restart();
	updateTotalRunoffGrid();

	if(MPIControl::instance().master() && output_grids)
		io->write2DGrid(total_runoff, mio::MeteoGrids::ROT, i_date);

	if (output_sums) {
		if(!is_glacier_mask_set || is_glacier_mask_dynamic) {
			updateGlacierMask();
			is_glacier_mask_set = true;
		}

		//Compute precip, glacier and melt runoffs
		mio::Grid2DObject precipRunoff(computePrecipRunoff(psum, ta));
		mio::Grid2DObject meltRunoff(total_runoff - precipRunoff);
		mio::Grid2DObject glacierRunoff(meltRunoff*glacier_mask);
		mio::Grid2DObject snowRunoff(meltRunoff - glacierRunoff);
		mio::Grid2DObject totalRunoff(total_runoff);

		//Get the grids of the additional meteo variables
		std::vector<mio::Grid2DObject> extraGrids;
		getExtraMeteoGrids(extraGrids);


		if(MPIControl::instance().master()) {
		//Resample runoff grids to match the cell size of the catchment masks
			if (fabs(grid_size_factor - 1.0) > 1e-5) {
				totalRunoff   = mio::LibResampling2D::Nearest(totalRunoff, grid_size_factor);
				precipRunoff  = mio::LibResampling2D::Nearest(precipRunoff, grid_size_factor);
				snowRunoff    = mio::LibResampling2D::Nearest(snowRunoff, grid_size_factor);
				glacierRunoff = mio::LibResampling2D::Nearest(glacierRunoff, grid_size_factor);
				for (size_t iVar(0); iVar < n_extra_meteo_variables; ++iVar)
					extraGrids[iVar] = mio::LibResampling2D::Nearest(extraGrids[iVar], grid_size_factor);
			}

			//Sum the runoffs over each mask and write them in the output files
			double currTotalRunoff, currPrecipRunoff, currSnowRunoff, currGlacierRunoff;
			std::vector<double> currMeteoVars(n_extra_meteo_variables);
			std::map<size_t, mio::Grid2DObject>::const_iterator itMask;
			for (itMask = catchment_masks.begin(); itMask != catchment_masks.end(); ++itMask) {
				currTotalRunoff   = sumOverMask(totalRunoff,   itMask->second);
				currPrecipRunoff  = sumOverMask(precipRunoff,  itMask->second);
				currSnowRunoff    = sumOverMask(snowRunoff,    itMask->second);
				currGlacierRunoff = sumOverMask(glacierRunoff, itMask->second);
				for (size_t iVar(0); iVar < n_extra_meteo_variables; ++iVar)
					currMeteoVars[iVar] = averageOverMask(extraGrids[iVar], itMask->second);

				updateOutputFile(itMask->first, i_date, currTotalRunoff,
						currPrecipRunoff, currSnowRunoff, currGlacierRunoff,
						currMeteoVars);
			}
		}
	}
	timer.stop();
}


/**
 * @brief Destructor of class Runoff
 */
Runoff::~Runoff()
{
	delete io;
}


/**
 * @brief Initializes private attribute catchment_masks
 * @param catchmentGrid grid defining the catchments. The catchment numbering
 * scheme must be specified in the ini file using the key CATCHMENT_NUMBERING
 * in section INPUT. This scheme can be either ALPINE3D_OLD (for catchments
 * numbered with powers of 2), or TAUDEM (for standard numbering).
 */
void Runoff::constructCatchmentMasks(mio::Grid2DObject catchmentGrid)
{
	const double factor = catchmentGrid.cellsize/resampling_cell_size;
	if (fabs(factor - 1.0) > 1e-5) {
		catchmentGrid = mio::LibResampling2D::Nearest(catchmentGrid, factor);
	}
	mio::Grid2DObject emptyGrid(catchmentGrid, mio::IOUtils::nodata);

	std::vector<size_t> currIndices;
	for (size_t iCell = 0; iCell < catchmentGrid.size(); ++iCell) {
		if (catchmentGrid(iCell) == mio::IOUtils::nodata) continue;
		const longuint currValue = static_cast<longuint>( round(catchmentGrid(iCell)) );

		if (catchment_numbering == Alpine3DOld) {
			currIndices = factorizeCatchmentNumber(currValue);
		} else {
			if (currValue > std::numeric_limits<size_t>::max()) {
				std::ostringstream os;
				os << "The ID numbers of some of the subwatersheds defined in "
				   << "the catchment file exceed the maximum index value ("
				   << std::fixed << std::numeric_limits<size_t>::max() << ")."
				   << " Did you forget to add \"CATCHMENT_NUMBERING = ALPINE3D_OLD\""
				   << " in section [INPUT] of your configuration file?";
				throw mio::IndexOutOfBoundsException(os.str(), AT);
			}
			currIndices.assign(1, static_cast<size_t>( currValue ));
		}

		for (std::vector<size_t>::const_iterator it = currIndices.begin(); it != currIndices.end(); ++it) {
			if (catchment_masks.find(*it) == catchment_masks.end())
				catchment_masks[*it] = emptyGrid;
			catchment_masks[*it](iCell) = 1.0;
		}
	}

	for (std::map<size_t, mio::Grid2DObject>::iterator it = catchment_masks.begin(); it != catchment_masks.end(); ++it) {
		cropMask(it->second);
	}
}


/**
 * @brief Updates protected attribute total_runoff.
 */
void Runoff::updateTotalRunoffGrid()
{
	total_runoff = snowpack->getGrid(SnGrids::MS_SOIL_RUNOFF);

	if ( (total_runoff.getNx() != slope_correction.getNx()) ||
	     (total_runoff.getNy() != slope_correction.getNy()) )
		throw mio::InvalidArgumentException("DEM and soil runoff grid have "
				"incompatible dimensions!", AT);

	total_runoff *= slope_correction;
}


/**
 * @brief Updates protected attribute glacier_mask.
 */
void Runoff::updateGlacierMask() {
	glacier_mask = snowpack->getGrid(SnGrids::GLACIER);

	if ( (glacier_mask.getNx() != slope_correction.getNx()) ||
	     (glacier_mask.getNy() != slope_correction.getNy()) )
		throw mio::InvalidArgumentException("DEM and glacier mask grid have "
				"incompatible dimensions!", AT);

	for (size_t iCell = 0; iCell < n_grid_cells; ++iCell)
		glacier_mask(iCell) = static_cast<double>(glacier_mask(iCell) == mio::IOUtils::nodata);
}


/**
 * @brief Computes the grid storing the runoff which originates from liquid
 * precipitation
 * @param psum grid of precipitation in mm/h
 * @param ta grid of air temperature
 * @return Precipitation runoff grid
 */
mio::Grid2DObject Runoff::computePrecipRunoff(const mio::Grid2DObject& psum, const mio::Grid2DObject& ta) const
{
	const mio::Grid2DObject surfRunoff( snowpack->getGrid(SnGrids::MS_SNOWPACK_RUNOFF) );
	if ( (psum.getNx()            != slope_correction.getNx()) ||
	     (ta.getNx()              != slope_correction.getNx()) ||
	     (surfRunoff.getNx()      != slope_correction.getNx()) ||
	     (psum.getNy()            != slope_correction.getNy()) ||
	     (ta.getNy()              != slope_correction.getNy()) ||
	     (surfRunoff.getNy()      != slope_correction.getNy()) )
		throw mio::InvalidArgumentException("DEM and input grids have incompatible "
				"dimensions!", AT);

	mio::Grid2DObject precipRunoff(total_runoff, mio::IOUtils::nodata);
	for (size_t iCell = 0; iCell < n_grid_cells; ++iCell) {
		const double& precip = psum(iCell);
		const double& runoff = total_runoff(iCell);
		if (runoff == mio::IOUtils::nodata || precip == mio::IOUtils::nodata ||
				ta(iCell) == mio::IOUtils::nodata || surfRunoff(iCell) == mio::IOUtils::nodata)
			precipRunoff(iCell) = mio::IOUtils::nodata;
		else if (ta(iCell) > thresh_rain)
			precipRunoff(iCell) = std::min(precip, runoff);
		else
			precipRunoff(iCell) = 0.0;
	}

	return precipRunoff;
}


/**
 * @brief Returns the grids corresponding to the additional meteo variables
 * which have to be averaged over the subcatchments and written in the output
 * files.
 * @param[out] grids vector containing the meteo grids
 */
void Runoff::getExtraMeteoGrids(std::vector<mio::Grid2DObject>& grids) const {
	grids.clear();
	grids.reserve(n_extra_meteo_variables);

	for(size_t iVar(0); iVar < n_extra_meteo_variables; ++iVar) {
		const SnGrids::Parameters& currParam(extra_meteo_variables.at(iVar));
		const mio::Grid2DObject tmp( snowpack->getGrid(currParam) );

		if(tmp.grid2D.getCount() == 0)
			throw mio::InvalidArgumentException("Cannot average parameter " +
					SnGrids::getParameterName(currParam) + " over the "
					"subcatchments: cannot retrieve the parameter value over "
					"the entire watershed", AT);
		else if ( (tmp.getNx() != slope_correction.getNx()) ||
		          (tmp.getNy() != slope_correction.getNy()) )
			throw mio::InvalidArgumentException("DEM and " +
					SnGrids::getParameterName(currParam) + " grid have "
					"incompatible dimensions!", AT);

		grids.push_back(tmp);
	}
}


/**
 * @brief Initializes the SMET files in which the per-catchment-aggregated
 * runoff values will be written
 * @param dem Reference to the DEM object
 */
void Runoff::initializeOutputFiles(const mio::Grid2DObject& dem) const
{
	std::stringstream ss;
	const double cellArea = pow(catchment_masks.begin()->second.cellsize, 2); // in m^2
	double catchArea;

	for (std::map<size_t, mio::Grid2DObject>::const_iterator itMask = catchment_masks.begin();
			itMask != catchment_masks.end(); ++itMask, ss.str(""), ss.clear()) {
		catchArea = double(itMask->second.grid2D.getCount())*cellArea*1e-6; //< in km^2

		ss << "catch" << std::setfill('0') << std::setw(2) << itMask->first;
		const std::string id = ss.str();

		const std::string filename = catchment_out_path + "/" + id + ".smet";
		std::ofstream smet_out; //Output file streams
		smet_out.open(filename.c_str());
		if (smet_out.fail()) throw mio::AccessException(filename.c_str(), AT);

		smet_out << "SMET 1.1 ASCII\n";
		smet_out << "[HEADER]\n";
		smet_out << "station_id   = " << id << "\n";
		smet_out << "station_name = " << id << "\n";
		smet_out << std::right;
		smet_out << std::fixed;

		smet_out << "latitude     = " << std::setw(11) << std::setprecision(8) << dem.llcorner.getLat() << "\n";
		smet_out << "longitude    = " << std::setw(11) << std::setprecision(8) << dem.llcorner.getLon() << "\n";
		smet_out << "altitude     = " << std::setw(7)  << std::setprecision(1) << dem.llcorner.getAltitude() << "\n";
		smet_out << "catchment_surface = " << std::setw(7)  << std::setprecision(3) << catchArea << "\n";
		smet_out << "comment = surface is in km^2, runoff (Ro) in m^3/h\n";
		smet_out << "nodata       = " << std::setw(7)  << std::setprecision(0) << mio::IOUtils::nodata << "\n";
		smet_out << "tz           = " << std::setw(7)  << std::setprecision(0) << tz_out << "\n";

		smet_out << "fields       = timestamp precip_Ro snowmelt_Ro glacier_melt_Ro total_Ro";
		for (size_t iVar(0); iVar < n_extra_meteo_variables; ++iVar)
			smet_out << " " << SnGrids::getParameterName(extra_meteo_variables.at(iVar));
		smet_out << "\n[DATA]\n";

		smet_out.close();
	}
}


/**
 * @brief Writes the runoff values aggregated over a given catchment in the
 * corresponding SMET file
 * @param catchId catchment id number
 * @param currTime time corresponding to the runoff values
 * @param totalRunoff total runoff (corrected for slope) over the catchment,
 * in mm/h
 * @param precipRunoff part of the total runoff attributable to liquid
 * precipitation, in mm/h
 * @param snowRunoff part of the total runoff attributable to snow melt, in
 * mm/h
 * @param glacierRunoff part of the total runoff attributable to glacier melt,
 * in mm/h
 * @param meteoVars additional meteorological variables averaged over the
 * catchment area. The units of each variable correspond to those used in
 * Alpine3D or Snowpack.
 */
void Runoff::updateOutputFile(const size_t& catchId, const mio::Date& currTime,
		const double& totalRunoff, const double& precipRunoff,
		const double& snowRunoff, const double& glacierRunoff,
		const std::vector<double>& meteoVars) const
{
	const double cellArea = resampling_cell_size*resampling_cell_size;

	std::stringstream ss;
	ss << "catch" << std::setfill('0') << std::setw(2) << catchId;
	const std::string filename = catchment_out_path + "/" + ss.str() + ".smet";
	std::ofstream smet_out; //Output file streams
	smet_out.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	if (smet_out.fail()) throw mio::AccessException(filename.c_str(), AT);

	smet_out.fill(' ');
	smet_out << std::right;
	smet_out << std::fixed;
	smet_out << currTime.toString(mio::Date::ISO) << " ";
	smet_out << std::setw(10) << std::setprecision(3) << cellArea*precipRunoff*1e-3  << " ";
	smet_out << std::setw(10) << std::setprecision(3) << cellArea*snowRunoff*1e-3    << " ";
	smet_out << std::setw(10) << std::setprecision(3) << cellArea*glacierRunoff*1e-3 << " ";
	smet_out << std::setw(10) << std::setprecision(3) << cellArea*totalRunoff*1e-3;
	for (size_t iVar(0); iVar < n_extra_meteo_variables; ++iVar)
		smet_out << " " << std::setw(10) << std::setprecision(2) << meteoVars.at(iVar);

	smet_out << "\n";
	smet_out.close();
}


/**
 * @brief Returns a IOManager object which will write the runoff grids in the
 * correct output folder
 * @param cfg Config object used in Alpine3D
 * @param outputGrids should the runoff grids be written out?
 */
mio::IOManager* Runoff::getIOManager(mio::Config cfg, const bool& outputGrids)
{
	if (outputGrids) {
		std::string runoff_grid2d("ARC");
		cfg.getValue("RUNOFF_GRID2D", "Output", runoff_grid2d, mio::IOUtils::nothrow);
		mio::IOUtils::toUpper(runoff_grid2d);

		cfg.addKey("GRID2D", "Output", runoff_grid2d);

		if(runoff_grid2d == "NETCDF") {
			const std::string runoff_grid2d_file = cfg.get("RUNOFF_GRID2DFILE", "Output");
			cfg.addKey("GRID2DFILE", "Output", runoff_grid2d_file);
		} else {
			const std::string runoff_grid2d_path = cfg.get("RUNOFF_GRID2DPATH", "Output");
			cfg.addKey("GRID2DPATH", "Output", runoff_grid2d_path);
		}
	}

	return new mio::IOManager(cfg);
}


/**
 * @brief Computes the grid used to correct the runoff values for slope
 * (Alpine3D assumes the cells to be flat, the predicted runoffs therefore
 * have to be corrected for terrain inclination)
 * @param dem Reference to the DEM object used by Alpine3D
 */
mio::Grid2DObject Runoff::getSlopeCorrection(const mio::DEMObject& dem)
{
	const double to_rad = M_PI/180.;

	mio::Grid2DObject slope_corr(dem, mio::IOUtils::nodata);
	for (size_t iCell = 0; iCell < dem.getNx()*dem.getNy(); ++iCell) {
		if (dem(iCell) != mio::IOUtils::nodata)
			slope_corr(iCell) = 1.0/cos(dem.slope(iCell)*to_rad);
	}

	return slope_corr;
}


/**
 * @brief Returns the numbering scheme used in the grid defining the catchments
 * @param in_cfg Reference to the Config object used in Alpine3D
 * @return Catchment numbering scheme
 */
Runoff::NumberingType Runoff::getCatchmentNumberingScheme(const mio::Config& in_cfg)
{
	NumberingType type;
	std::string numbering("Alpine3D_old");
	in_cfg.getValue("CATCHMENT_NUMBERING", "INPUT", numbering, mio::IOUtils::nothrow);
	mio::IOUtils::toUpper(numbering);
	if (numbering == "TAUDEM")
		type = TauDEM;
	else if (numbering == "ALPINE3D_OLD")
		type = Alpine3DOld;
	else
		throw mio::IOException("Key CATCHMENT_NUMBERING can only be either "
				"'TauDEM' or 'Alpine3D_old'", AT);
	return type;
}


/**
 * @brief Returns whether the glaciers are dynamic, i.e. whether the glacier
 * mask should be expected to change from one time step to the other.
 * @param cfg Reference to the Config object holding the simulation parameters
 * @return Is glacier mask dynamic?
 */
bool Runoff::getIsGlacierDynamic(const mio::Config& cfg) {
	bool isDynamic(false);
	cfg.getValue("MASK_DYNAMIC", "Output", isDynamic, mio::IOUtils::nothrow);
	return isDynamic;
}


/**
 * @brief This method returns the meteo variables which have to be averaged over
 * the subwatershed areas and written in the output files on top of runoff.
 * @param cfg Reference to the Config object holding the simulation parameters
 * @return Variables which have to be averaged over the subwatershed areas
 */
std::vector<SnGrids::Parameters> Runoff::getExtraMeteoVariables(const mio::Config& cfg) {
	std::vector<std::string> extraDataNames;
	cfg.getValue("RUNOFF_FILES_EXTRA_DATA", "Output", extraDataNames, mio::IOUtils::nothrow);
	const size_t nData(extraDataNames.size());

	std::vector<SnGrids::Parameters> extraData;
	extraData.reserve(nData);
	for(size_t iData(0); iData < nData; ++iData) {
		const size_t iParam(SnGrids::getParameterIndex(extraDataNames[iData]));
		extraData.push_back((SnGrids::Parameters)iParam);
	}

	return extraData;
}


/**
 * @brief Returns the cell size which must be used to resample both the
 * catchment masks and the runoff grids
 * @param in_dem reference to the Config object used in Alpine3D
 * @param catchmentGrid grid defining the catchments
 * @return Resampling cell size in m
 */
double Runoff::getResamplingCellSize(const mio::DEMObject& in_dem,
		const mio::Grid2DObject& catchmentGrid)
{
	double resampCellSize;
	const double minSize = std::min(catchmentGrid.cellsize, in_dem.cellsize);
	const double maxSize = std::max(catchmentGrid.cellsize, in_dem.cellsize);
	if (isMultiple(maxSize, minSize)) {
		const double llxOffset = fabs(catchmentGrid.llcorner.getEasting() -
				in_dem.llcorner.getEasting());
		const double llyOffset = fabs(catchmentGrid.llcorner.getNorthing() -
				in_dem.llcorner.getNorthing());
		const double xCellSize(estimateResamplingCellSize(llxOffset, minSize));
		const double yCellSize(estimateResamplingCellSize(llyOffset, minSize));
		if (isMultiple(xCellSize, yCellSize) || isMultiple(yCellSize, xCellSize))
			resampCellSize = std::min(xCellSize, yCellSize);
		else
			resampCellSize = MIN_CELL_SIZE;
	} else {
		resampCellSize = MIN_CELL_SIZE;
	}
	if (resampCellSize <= MIN_CELL_SIZE + DISTANCE_ABSOLUTE_PRECISION)
		std::cout << "[w] The resolution of the grid defining the sub-catchments "
		          << "matches poorly with the resolution of the DEM. The runoff "
		          << "module will re-interpolate all runoff grids using a new cell "
		          << "size of " << MIN_CELL_SIZE << " meters, which might "
		          << "significantly impact simulation time in case of a large DEM. "
		          << "It is recommended that you define your sub-catchments using the "
		          << "same DEM as the one used by Alpine3D to solve this issue.\n";

	return resampCellSize;
}


/**
 * @brief Returns whether the first input is a multiple of the second one
 */
bool Runoff::isMultiple(const double& a, const double& b) {
	const double res = fmod(a, b);
	return fabs(res - b/2.0) > b/2.0 - DISTANCE_ABSOLUTE_PRECISION;
}


double Runoff::estimateResamplingCellSize(const double& llOffset,
		const double& currSizeEstimate) {
	double cellSize;
	if (isMultiple(llOffset, currSizeEstimate)) {
		cellSize = std::max(currSizeEstimate, MIN_CELL_SIZE);
	} else {
		double xShift = fmod(llOffset, currSizeEstimate);
		xShift   = std::min(xShift, currSizeEstimate - xShift);
		cellSize = currSizeEstimate/round(currSizeEstimate/xShift);
		cellSize = std::max(cellSize, MIN_CELL_SIZE);
	}
	return cellSize;
}


/**
 * @brief This function splits up a given unsigned int value into a sum of powers of 2
 * @param value The value that shall be split up into a sum of powers of two
 * @return A vector that will hold all exponents of the sum of powers
 */
std::vector<size_t> Runoff::factorizeCatchmentNumber(longuint value)
{
	std::vector<size_t> result;
	while (value != 0) {
		const double lnresult = log(double(value))/log(2);
		const size_t exp = (unsigned int)floor(lnresult);
		const longuint pot = longuint(floor(pow(double(2), double(exp)) + 0.00001)); //numerically stable
		result.push_back(exp);
		value -= pot;
	}
	return result;
}


/**
 * @brief Crops the mask so as to remove as many nodata cells as possible
 * (note: this method possibly changes the georeferencing of the mask by
 * changing its lower-left corner as well as its number of rows or columns)
 */
void Runoff::cropMask(mio::Grid2DObject& mask)
{
	size_t min_ix(mask.getNx()-1), max_ix(0), min_iy(mask.getNy()-1), max_iy(0);
	for (size_t iy = 0; iy < mask.getNy(); ++iy) {
		for (size_t ix = 0; ix < mask.getNx(); ++ix) {
			if (mask(ix,iy) != mio::IOUtils::nodata) {
				if (ix < min_ix) min_ix = ix;
				if (ix > max_ix)	max_ix = ix;
				if (iy < min_iy) min_iy = iy;
				if (iy > max_iy)	max_iy = iy;
			}
		}
	}
	mio::Grid2DObject tmp(mask, min_ix, min_iy, max_ix-min_ix+1, max_iy-min_iy+1);
	mask = tmp;
}


/**
 * @brief Sums the values of the grid cells which are located over the mask
 * @param grid grid whose cell values have to be summed
 * @param mask mask over which the grid cell values have to be summed
 * @return Sum of the masked grid cell values
 */
double Runoff::sumOverMask(const mio::Grid2DObject& grid, const mio::Grid2DObject& mask)
{
	mio::Coords llCorner(mask.llcorner);
	grid.gridify(llCorner);
	mio::Grid2DObject subgrid;
	try {
		mio::Grid2DObject tmp(grid, llCorner.getGridI(), llCorner.getGridJ(),
				mask.getNx(), mask.getNy());
		subgrid = tmp;
	} catch (...) {
		throw mio::InvalidFormatException("Catchment mask extends beyond the DEM boundaries", AT);
	}

	double sum(0.0);
	for (size_t iCell = 0; iCell < mask.getNx()*mask.getNy(); ++iCell) {
		if (mask(iCell) != mio::IOUtils::nodata && subgrid(iCell) != mio::IOUtils::nodata)
			sum += subgrid(iCell);
	}
	return sum;
}


/**
 * @brief Averages the values of the grid cells which are located over the mask
 * @param grid grid whose cell values have to be averaged
 * @param mask mask over which the grid cell values have to be averaged
 * @return Average of the masked grid cell values
 */
double Runoff::averageOverMask(const mio::Grid2DObject& grid, const mio::Grid2DObject& mask)
{
	return sumOverMask(grid, mask)/double(mask.grid2D.getCount());
}

double Runoff::getTiming() const
{
	return timer.getElapsed();
}
