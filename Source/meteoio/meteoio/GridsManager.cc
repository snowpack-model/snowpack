// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/GridsManager.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/MathOptim.h>

using namespace std;

namespace mio {

GridsManager::GridsManager(IOHandler& in_iohandler, const Config& in_cfg)
             : iohandler(in_iohandler), cfg(in_cfg), buffer(0), gridprocessor(cfg), grids2d_list(), grids2d_start(), grids2d_end(),
               grid2d_list_buffer_size(370.), processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated), dem_altimeter(false)
{
	size_t max_grids = 10;
	cfg.getValue("BUFF_GRIDS", "General", max_grids, IOUtils::nothrow);  //HACK document it!
	buffer.setMaxGrids(max_grids);
	cfg.getValue("BUFFER_SIZE", "General", grid2d_list_buffer_size, IOUtils::nothrow);
	cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow); //HACK document it! if no dem is found but local and sea level pressure grids are found, use them to rebuild a DEM; [Input] section
}

/**
* @brief Set the desired ProcessingLevel
* @details The processing level affects the way meteo data is read and processed. Three values are possible:
*    - IOUtils::raw data shall be read directly from the buffer;
*    - IOUtils::filtered data shall be filtered before returned to the user;
*    - IOUtils::resampled data shall be resampled before returned to the user;
*
* This only affects the function getMeteoData(const Date&, METEO_DATASET&). The three values
* can be combined: e.g. IOUtils::filtered | IOUtils:resampled
* @param i_level The ProcessingLevel values that shall be used to process data
*/
void GridsManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOUtils::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOUtils::raw) == IOUtils::raw)
	    && ((i_level & IOUtils::filtered) == IOUtils::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

/**
* @brief Read the requested grid, according to the configured processing level
* @details If the grid has been buffered, it will be returned from the buffer. If it is not available but can be generated, it will
* be generated transparently. If everything fails, it will throw an exception.
* @param[out] grid2D a grid filled with the requested parameter
* @param option a parameter used to figure out which filename or grid to read
*/
void GridsManager::read2DGrid(Grid2DObject& grid2D, const std::string& option)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, option);
	} else {
		if (buffer.get(grid2D, option)) return;

		iohandler.read2DGrid(grid2D, option);
		buffer.push(grid2D, option);
	}
}

/**
* @brief Read the requested grid, according to the configured processing level
* @details If the grid has been buffered, it will be returned from the buffer. If it is not available but can be generated, it will
* be generated transparently. If everything fails, it will thrown an exception.
* @param[out] grid2D a grid filled with the requested parameter
* @param[in] parameter the parameter to get
* @param[in] date the timestamp that we would like to have
* @param[in] enable_grid_resampling Completely enable or disable temporal grid resampling.
*/
void GridsManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date, const bool& enable_grid_resampling)
{
	grid2D = getGrid(parameter, date, false, enable_grid_resampling);
}

void GridsManager::readDEM(DEMObject& grid2D)
{
	//TODO: dem_altimeter; reading DEM data with no associated date (ie Date()) OR with associated date (for example, TLS data)
	if (processing_level == IOUtils::raw){
		iohandler.readDEM(grid2D);
	} else {
		if (!buffer.get(grid2D, "/:DEM")) {
			iohandler.readDEM(grid2D);
			buffer.push(grid2D, "/:DEM");
		}
	}

	//reproject grid if it is lat/lon
	if (grid2D.isLatlon())
		grid2D.reproject(); //this is currently very, very primitive
}

void GridsManager::readLanduse(Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readLanduse(grid2D);
	} else {
		if (buffer.get(grid2D, "/:LANDUSE"))
			return;

		iohandler.readLanduse(grid2D);
		buffer.push(grid2D, "/:LANDUSE");
	}
}

void GridsManager::readGlacier(Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readGlacier(grid2D);
	} else {
		if (buffer.get(grid2D, "/:GLACIER"))
			return;

		iohandler.readGlacier(grid2D);
		buffer.push(grid2D, "/:GLACIER");
	}
}

void GridsManager::readAssimilationData(const Date& date, Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readAssimilationData(date, grid2D);
	} else {
		const string grid_hash = "/:ASSIMILATIONDATA"+date.toString(Date::ISO);
		if (buffer.get(grid2D, grid_hash))
			return;

		iohandler.readAssimilationData(date, grid2D);
		buffer.push(grid2D, grid_hash);
	}
}

/** @brief Create a list of virtual stations from each grid point
 * @details Please note that the two options are mutually exclusive.
 * @param[in] dem the Digital Elevation Model to use (each grid point will be made into a virtual station)
 * @return a vector of virtual stations
 */
std::vector<StationData> GridsManager::initVirtualStationsAtAllGridPoints(const DEMObject& dem) const
{
	std::vector<StationData> v_stations;
	//get virtual stations coordinates
	std::string coordin, coordinparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	Coords llcorner( dem.llcorner );
	llcorner.setProj(coordin, coordinparam); //make sure the DEM and the VStations are in the same projection
	const double dem_easting = llcorner.getEasting();
	const double dem_northing = llcorner.getNorthing();

	Coords curr_point(coordin, coordinparam);
	size_t stat_id=0;
	for (size_t jj=0; jj<dem.getNy(); jj++) {
		for (size_t ii=0; ii<dem.getNx(); ii++) {
			const double easting = dem_easting + dem.cellsize*static_cast<double>(ii);
			const double northing = dem_northing + dem.cellsize*static_cast<double>(jj);
			curr_point.setXY(easting, northing, dem(ii,jj));
			curr_point.setGridIndex(static_cast<int>(ii), static_cast<int>(jj), IOUtils::inodata, true);

			//extract vstation number, build the station name and station ID
			stat_id++;
			const std::string id_str( IOUtils::toString(stat_id) );
			const std::string id_grid( IOUtils::toString(ii)+"-"+IOUtils::toString(jj) );
			StationData sd(curr_point, "GRID_"+id_grid, "Virtual_Station_"+id_str);
			sd.setSlope(dem.slope(ii,jj), dem.azi(ii,jj));
			v_stations.push_back( sd );
		}
	}

	return v_stations;
}

/** @brief Create a list of virtual stations from the user-provided input
 * @details Please note that the two options are mutually exclusive.
 * @param[in] dem the Digital Elevation Model to use
 * @param[in] adjust_coordinates should the coordinates be recomputed to match DEM cells?
 * @param[in] fourNeighbors pick the surrounding four nodes instead of only the exact one?
 * @return a vector of virtual stations
 */
std::vector<StationData> GridsManager::initVirtualStations(const DEMObject& dem, const bool& adjust_coordinates, const bool& fourNeighbors) const
{
	//get virtual stations coordinates
	std::string coordin, coordinparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	Coords llcorner( dem.llcorner );
	llcorner.setProj(coordin, coordinparam); //make sure the DEM and the VStations are in the same projection
	const double dem_easting = llcorner.getEasting();
	const double dem_northing = llcorner.getNorthing();

	//read the provided coordinates, remove duplicates and generate metadata
	std::vector<StationData> v_stations;
	const std::vector< std::pair<std::string, std::string> > vecStation( cfg.getValues("Vstation", "InputEditing") );
	
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		if (vecStation[ii].first.find('_') != std::string::npos) continue; //so we skip the other vstations_xxx parameters

		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		Coords curr_point(coordin, coordinparam, vecStation[ii].second);
		if (curr_point.isNodata()) continue;

		if (!dem.gridify(curr_point))
			throw NoDataException("Virtual station "+vecStation[ii].second+" is not contained in provided DEM "+dem.toString(DEMObject::SHORT), AT);

		//extract vstation number, used to build the station name and station ID
		const std::string id_num( vecStation[ii].first.substr(std::string("Vstation").length()) );
		const bool has_id = cfg.keyExists("Vid"+id_num, "InputEditing");
		const std::string vir_id = (has_id)? cfg.get("Vid"+id_num, "InputEditing"): "VIR"+id_num;
		const bool has_name = cfg.keyExists("Vname"+id_num, "InputEditing");
		const std::string vir_name = (has_name)? cfg.get("Vname"+id_num, "InputEditing"): std::string("");

		size_t i = curr_point.getGridI(), j = curr_point.getGridJ();
		if (fourNeighbors) { //pick the surrounding four nodes
			if (i==dem.getNx()) i--;
			if (j==dem.getNy()) j--;

			for (size_t mm=i; mm<=i+1; mm++) {
				for (size_t nn=j; nn<=j+1; nn++) {
					const double easting = dem_easting + dem.cellsize*static_cast<double>(mm);
					const double northing = dem_northing + dem.cellsize*static_cast<double>(nn);
					curr_point.setXY(easting, northing, dem(mm,nn));
					curr_point.setGridIndex(static_cast<int>(mm), static_cast<int>(nn), IOUtils::inodata, true);
					const std::string sub_id( IOUtils::toString((mm-i)*2+(nn-j)+1) );
					const std::string grid_pos( IOUtils::toString(mm) + "-" + IOUtils::toString(nn) );
					StationData sd(curr_point, vir_id+"_"+sub_id, "Virtual_Station_"+grid_pos);
					sd.setSlope(dem.slope(mm,nn), dem.azi(mm,nn));
					v_stations.push_back( sd );
				}
			}
		} else { //pick the exact node
			if (adjust_coordinates) { //adjust coordinates to match the chosen cell
				if (dem.isLatlon()) {
					const double dem_lat = llcorner.getLat();
					const double dem_lon = llcorner.getLon();
					const double new_lat = dem_lat + static_cast<double>(j)*(dem.ur_lat - dem_lat)/static_cast<double>(dem.getNy()-1);
					const double new_lon = dem_lon + static_cast<double>(i)*(dem.ur_lon - dem_lon)/static_cast<double>(dem.getNx()-1);
					curr_point.setLatLon(new_lat, new_lon, dem(i,j));
				} else {
					const double easting = dem_easting + static_cast<double>(i)*dem.cellsize;
					const double northing = dem_northing + static_cast<double>(j)*dem.cellsize;
					curr_point.setXY(easting, northing, dem(i,j));
				}
				curr_point.setGridIndex(static_cast<int>(i), static_cast<int>(j), IOUtils::inodata, true);
			}

			StationData sd(curr_point, vir_id, (has_name)? (vir_name): ("Virtual_Station_"+id_num) );
			sd.setSlope(dem.slope(i,j), dem.azi(i,j));
			v_stations.push_back( sd );
		}
	}

	if (v_stations.empty()) throw NoDataException("No virtual stations provided", AT);

	//removing potential duplicates
	if (StationData::unique( v_stations, true ) && !fourNeighbors)
		std::cout << "[W] Some of the input stations have been identified as duplicates and removed\n";

	return v_stations;
}

/**
* @brief Extract time series from grids at the specified points (virtual stations).
* @param[in] dem the Digital Elevation Model to use to check the geolocalization of the grids
* @param[in] v_params the MeteoGrids parameter index that have to be extracted
* @param[in] v_stations a vector of StationData where to provide the meteorological data
* @param[in] date when to extract the virtual stations
* @param[in] PtsExtract use the method to only request points instead of full grids from the plugin?
* @return a vector of meteodata for the configured virtual stations at the provided date, for the provided parameters
*/
METEO_SET GridsManager::getVirtualStationsFromGrid(const DEMObject& dem, const std::vector<size_t>& v_params, const std::vector<StationData>& v_stations, const Date& date, const bool& PtsExtract)
{
	//HACK handle extra parameters when possible
	const size_t nrStations = v_stations.size();
	METEO_SET vecMeteo( nrStations );

	//create stations without measurements
	std::vector< std::pair <size_t, size_t> > Pts;
	for (size_t ii=0; ii<nrStations; ii++) {
		MeteoData md(date, v_stations[ii]);
		vecMeteo[ii] = md;
		Pts.push_back(std::pair <size_t, size_t> (v_stations[ii].position.getGridI(), v_stations[ii].position.getGridJ())); //this should work since invalid stations have been removed in init
	}

	for (size_t param=0; param<v_params.size(); param++) { //loop over required parameters
		const MeteoGrids::Parameters grid_param = static_cast<MeteoGrids::Parameters>( v_params[param] );
		if (PtsExtract) {
			const vector <double> retVec = getPtsfromGrid(grid_param, date, Pts);
			for (size_t ii=0; ii<nrStations; ii++) { //loop over all virtual stations
				//check if this is a standard MeteoData parameter
				const size_t meteo_param = vecMeteo[ii].getParameterIndex( MeteoGrids::getParameterName(grid_param) ); //is this name also a meteoparameter?
				if (meteo_param!=IOUtils::npos)
					vecMeteo[ii]( static_cast<MeteoData::Parameters>(meteo_param) ) = retVec[ii];
				}
		} else {
			const Grid2DObject grid( getGrid(grid_param, date, false) ); //keep lat/lon grids if they are so

			if (!grid.isSameGeolocalization(dem))
				throw InvalidArgumentException("In GRID_EXTRACT, the DEM and the source grid don't match for '"+MeteoGrids::getParameterName(grid_param)+"' on "+date.toString(Date::ISO), AT);

			for (size_t ii=0; ii<nrStations; ii++) { //loop over all virtual stations
				const size_t grid_i = v_stations[ii].position.getGridI(); //this should work since invalid stations have been removed in init
				const size_t grid_j = v_stations[ii].position.getGridJ();

				//check if this is a standard MeteoData parameter
				const size_t meteo_param = vecMeteo[ii].getParameterIndex( MeteoGrids::getParameterName(grid_param) ); //is this name also a meteoparameter?
				if (meteo_param!=IOUtils::npos)
					vecMeteo[ii]( static_cast<MeteoData::Parameters>(meteo_param) ) = grid(grid_i, grid_j);
			}
		}
	}

	return vecMeteo;
}

/**
* @brief Extract time series from grids at the specified points (virtual stations).
* @param[in] dem the Digital Elevation Model to use to check the geolocalization of the grids
* @param[in] v_params the MeteoGrids parameter index that have to be extracted
* @param[in] v_stations a vector of StationData where to provide the meteorological data
* @param[in] dateStart when to start extracting the virtual stations
* @param[in] dateEnd when to stop extracting the virtual stations
* @param[in] PtsExtract use the method to only request points instead of full grids from the plugin?
* @return a vector of meteodata for the configured virtual stations at the provided date, for the provided parameters
*/
std::vector<METEO_SET> GridsManager::getVirtualStationsFromGrid(const DEMObject& dem, const std::vector<size_t>& v_params, const std::vector<StationData>& v_stations, const Date& dateStart, const Date& dateEnd, const bool& PtsExtract)
{
	const size_t nrStations = v_stations.size();
	std::vector<METEO_SET> vecvecMeteo( nrStations );

	const bool status = setGrids2d_list(dateStart, dateEnd);
	if (!status)
		throw InvalidArgumentException("The chosen plugin seems not to support the list2DGrids call that is required for gridded data extraction", AT);

	//look for the last date in grids2d_list just before dateStart
	std::map<Date, std::set<size_t> >::const_iterator it;
	for (it=grids2d_list.begin(); it!=grids2d_list.end(); ++it) {
		if (it->first>=dateStart) break;
	}
	if (it==grids2d_list.end()) return std::vector<METEO_SET>();
	if (it!=grids2d_list.begin() && it->first!=dateStart) --it; //we want to ensure the range contains the start date (for interpolations)

	//now, we read the data for each available timestep
	for (; it!=grids2d_list.end(); ++it) {
		const METEO_SET vecMeteo( getVirtualStationsFromGrid(dem, v_params, v_stations, it->first, PtsExtract) ); //the number of stations can not change
		for (size_t ii=0; ii<nrStations; ii++)
			vecvecMeteo[ii].push_back( vecMeteo[ii] );

		if (it->first>=dateEnd) break;
	}

	return vecvecMeteo;
}

////////////////////////////////////////////////////// Private members //////////////////////////////////////////

/**
* @brief Check if a given parameter is available in a grid (in buffer or in raw data)
* @param available_params list of parameters available by a direct (=raw) read from the plugin
* @param parameter the parameter to check the availability of
* @param date the data of availability
* @return true if the given parameter is available (ie a "read" will return the data for this parameter)
*/
bool GridsManager::isAvailable(const std::set<size_t>& available_params, const MeteoGrids::Parameters& parameter, const Date& date) const
{
	const bool in_buffer = buffer.has(parameter, date);
	if (in_buffer) return true;

	const bool in_raw = (available_params.find( parameter ) != available_params.end());
	return in_raw;
}

/**
* @brief Check if the grids2d_list "buffer" covers the proper range and rebuffer if not.
* @details This does not mean that the requested date is in the buffer, but that the range covered by the buffer contains this date.
* @param[in] date the timestamp that we would like to have
* @return true if the buffer is ready (it was already ready or it has been rebuffered) or false if the grid plugin could not provide the information
*/
bool GridsManager::setGrids2d_list(const Date& date)
{
	if (grids2d_list.empty() || date<grids2d_start || date>grids2d_end) {
		grids2d_start = date - 1.;
		grids2d_end = date + grid2d_list_buffer_size;

		const bool status = iohandler.list2DGrids(grids2d_start, grids2d_end, grids2d_list);
		if (status) {
			//the plugin might have returned a range larger than requested, so adjust the min/max dates if necessary
			if (!grids2d_list.empty()) {
				if (grids2d_start > grids2d_list.begin()->first) grids2d_start = grids2d_list.begin()->first;
				if (grids2d_end < grids2d_list.rbegin()->first) grids2d_end = grids2d_list.rbegin()->first;
			}
			return true;
		}
		return false;
	}
	return true;
}

bool GridsManager::setGrids2d_list(const Date& dateStart, const Date& dateEnd)
{
	if (grids2d_list.empty() || dateStart<grids2d_start || dateEnd>grids2d_end) {
		grids2d_start = dateStart;
		grids2d_end = dateEnd;
		const bool status = iohandler.list2DGrids(grids2d_start, grids2d_end, grids2d_list);
		if (status) {
			//the plugin might have returned a range larger than requested, so adjust the min/max dates if necessary
			if (!grids2d_list.empty()) {
				if (grids2d_start > grids2d_list.begin()->first) grids2d_start = grids2d_list.begin()->first;
				if (grids2d_end < grids2d_list.rbegin()->first) grids2d_end = grids2d_list.rbegin()->first;
			}
			return true;
		}
		return false;
	}
	return true;
}

/**
* @brief Read a grid, either from the buffer or from the plugin
* @param parameter the parameter to get
* @param date the data associated with the parameter
* @return 2D grid containing the data for this parameter at this date
*/
Grid2DObject GridsManager::getRawGrid(const MeteoGrids::Parameters& parameter, const Date& date)
{
	Grid2DObject grid2D;
	if (!buffer.get(grid2D, parameter, date)) {
		iohandler.read2DGrid(grid2D, parameter, date);
		buffer.push(grid2D, parameter, date);
	}

	return grid2D;
}

/**
* @brief Get the requested grid, according to the configured processing level
* @details If the grid has been buffered, it will be returned from the buffer. If it is not available but can be generated, it will
* be generated transparently. If everything fails, it will thrown an exception.
* @param[in] parameter the parameter to get
* @param[in] date the timestamp that we would like to have
* @param[in] enforce_cartesian set to true to garantee a cartesian grid (it will be reprojected if necessary)
* @return a grid filled with the requested parameter
*/
Grid2DObject GridsManager::getGrid(const MeteoGrids::Parameters& parameter, const Date& date, const bool& enforce_cartesian, const bool& enable_grid_resampling)
{
	Grid2DObject grid2D;

	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, parameter, date);
	} else {
		if (!buffer.get(grid2D, parameter, date)) {
			const bool status = setGrids2d_list( date ); //rebuffer the grid list if necessary
			if (!status) { //this means that the list2DGrids call is not implemeted in the plugin, we try to save the day...
				iohandler.read2DGrid(grid2D, parameter, date);
				buffer.push(grid2D, parameter, date);
			} else { //the list2DGrids call is implemented in the plugin
				const std::map<Date, std::set<size_t> >::const_iterator it = grids2d_list.find(date);
				if (it!=grids2d_list.end()) {
					if ( it->second.find(parameter) != it->second.end() ) {
						iohandler.read2DGrid(grid2D, parameter, date);
						buffer.push(grid2D, parameter, date);
					} else { //the right parameter could not be found, can we generate it?
						if (!generateGrid(grid2D, it->second, parameter, date))
							throw NoDataException("Could not find or generate a grid of "+MeteoGrids::getParameterName( parameter )+" at time "+date.toString(Date::ISO), AT);
					}
				} else {
					if (enable_grid_resampling) { //user requests to temporally interpolate inbetween grids
						const Date sdate( date - gridprocessor.getWindowSize() );
						const Date edate( date + gridprocessor.getWindowSize() );
						const bool status_all = setGrids2d_list(sdate, edate); //rebuffer the grid list if necessary
						if (!status_all)
							throw InvalidArgumentException("The data input plugin does not support the necessary call to query all available grids.", AT);
						const std::map<Date, Grid2DObject> all_grids( getAllGridsForParameter(parameter) );
						gridprocessor.resample(date, parameter, all_grids, grid2D);
						buffer.push(grid2D, parameter, date);
					} else {
						const std::string msg1("Could not find grid for "+MeteoGrids::getParameterName( parameter )+" at time " + date.toString(Date::ISO) + ". " );
						std::string msg2("There are no grids.");
						if (!grids2d_list.empty())
							msg2 = "There are grids from " + grids2d_list.begin()->first.toString(Date::ISO) + " until " + grids2d_list.rbegin()->first.toString(Date::ISO);
						throw NoDataException(msg1 + msg2, AT);
					} //endif enable_grid_resampling
				}
			} //status
		} //buffer
	} //processing level

	 //reproject grid if it is lat/lon
	if (enforce_cartesian && grid2D.isLatlon())
		grid2D.reproject(); //HACK this is currently very, very primitive
	
	return grid2D;
}

/**
 * @brief Build a list of grids that are available as input for a specific parameter.
 * @details This function queries the input plugin for grids that are available for a specific
 * parameter. If gridded data is available from which the parameter can losslessly be generated
 * this is treated as equal and also returned.
 * @param[in] parameter The meteo parameter to go look for.
 * @return A list of available grids and their dates.
 */
std::map<Date, Grid2DObject> GridsManager::getAllGridsForParameter(const MeteoGrids::Parameters& parameter)
{
	std::map<Date, Grid2DObject> all_grids;
	for (auto it = grids2d_list.begin(); it != grids2d_list.end(); ++it) {
		if (isAvailable(it->second, parameter, it->first)) {
			const auto pair( std::make_pair(it->first, getRawGrid(parameter, it->first)) );
			all_grids.insert(pair);
		} else {
			Grid2DObject generated;
			if (generateGrid(generated, it->second, parameter, it->first)) {
				const auto pair( std::make_pair(it->first, generated) );
				all_grids.insert(pair);
			}
		}
	} //endfor
	return all_grids;
}

/**
* @brief Generate a grid for a given parameter, based on the available parameters
* @details Even if a given parameter is not available, it might be possible to generate it
* on the fly based on the available data (for example, U and V wind components can be used to generate
* the VW and DW vector wind components).
*
* It is assumed that the meteo parameters are coming out of models, so the available_params are
* all available at all the timesteps, so we don't need to search a combination of parameters and timesteps
*
* @param[out] grid2D a grid filled with the requested parameter or empty if it could not be generated
* @param available_params list of parameters available by a direct (=raw) read from the plugin
* @param parameter the parameter to get
* @param date the data associated with the parameter
* @return true if the requested grid could be generated
*
* @note HACK missing: checking that all grids used together have the same geolocalization
*/
bool GridsManager::generateGrid(Grid2DObject& grid2D, const std::set<size_t>& available_params, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (parameter==MeteoGrids::DEM) {
		if (!dem_altimeter) return false;
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		const bool hasP = isAvailable(available_params, MeteoGrids::P, date);
		const bool hasP_sea = isAvailable(available_params, MeteoGrids::P_SEA, date);
		if (!hasTA || !hasP || !hasP_sea) return false;

		const Grid2DObject ta( getRawGrid(MeteoGrids::TA, date) );
		const Grid2DObject p( getRawGrid(MeteoGrids::P, date) );
		const Grid2DObject p_sea( getRawGrid(MeteoGrids::P_SEA, date) );

		static const double k = Cst::gravity / (Cst::mean_adiabatique_lapse_rate * Cst::gaz_constant_dry_air);
		static const double k_inv = 1./k;
		grid2D.set(p, IOUtils::nodata);
		for (size_t ii=0; ii<grid2D.size(); ii++) {
			if (grid2D(ii)==IOUtils::nodata || ta(ii)==IOUtils::nodata || p_sea(ii)==IOUtils::nodata) continue;
			const double K = pow(grid2D(ii)/p_sea(ii), k_inv); //because grid2D has been initialized with P
			grid2D(ii) = ta(ii)*Cst::earth_R0*(1.-K) / (Cst::mean_adiabatique_lapse_rate * Cst::earth_R0 - ta(ii)*(1.-K));
		}
		buffer.push(grid2D, parameter, date);
		return true;
	}

	if (parameter==MeteoGrids::VW || parameter==MeteoGrids::DW) {
		const bool hasU = isAvailable(available_params, MeteoGrids::U, date);
		const bool hasV = isAvailable(available_params, MeteoGrids::V, date);
		if (!hasU || !hasV) return false;

		const Grid2DObject V( getRawGrid(MeteoGrids::V, date) );
		Grid2DObject U( getRawGrid(MeteoGrids::U, date) );

		grid2D.set(U, IOUtils::nodata);
		if (parameter==MeteoGrids::VW) {
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = sqrt( Optim::pow2(U(ii)) + Optim::pow2(V(ii)) );
		} else {
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = IOUtils::UV_TO_DW(U(ii), V(ii)); // turn into degrees [0;360)
		}

		buffer.push(grid2D, parameter, date);
		return true;
	}

	if (parameter==MeteoGrids::RH) {
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		if (!hasTA) return false;
		const Grid2DObject TA( getRawGrid(MeteoGrids::TA, date) );

		const bool hasDEM = isAvailable(available_params, MeteoGrids::DEM, date);
		const bool hasQI = isAvailable(available_params, MeteoGrids::QI, date);
		const bool hasTD = isAvailable(available_params, MeteoGrids::TD, date);

		if (hasTD) {
			grid2D = getRawGrid(MeteoGrids::TD, date);

			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = Atmosphere::DewPointtoRh(grid2D(ii), TA(ii), false);
			return true;
		} else if (hasQI && hasDEM) {
			const Grid2DObject dem( getRawGrid(MeteoGrids::DEM, date) ); //HACK use readDEM instead?
			grid2D = getRawGrid(MeteoGrids::QI, date);
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = Atmosphere::specToRelHumidity(dem(ii), TA(ii), grid2D(ii));
			return true;
		}

		return false;
	}

	if (parameter==MeteoGrids::ISWR) {
		const bool hasISWR_DIFF = isAvailable(available_params, MeteoGrids::ISWR_DIFF, date);
		const bool hasISWR_DIR = isAvailable(available_params, MeteoGrids::ISWR_DIR, date);

		if (hasISWR_DIFF && hasISWR_DIR) {
			const Grid2DObject iswr_diff( getRawGrid(MeteoGrids::ISWR_DIFF, date) );
			grid2D = getRawGrid(MeteoGrids::ISWR_DIR, date);
			grid2D += iswr_diff;
			buffer.push(grid2D, MeteoGrids::ISWR, date);
			return true;
		}

		const bool hasALB = isAvailable(available_params, MeteoGrids::ALB, date);
		const bool hasRSWR = isAvailable(available_params, MeteoGrids::ISWR, date);
		if (hasRSWR && hasALB) {
			const Grid2DObject alb( getRawGrid(MeteoGrids::ALB, date) );
			grid2D = getRawGrid(MeteoGrids::RSWR, date);
			grid2D /= alb;
			buffer.push(grid2D, MeteoGrids::ISWR, date);
			return true;
		}

		return false;
	}

	if (parameter==MeteoGrids::RSWR) {
		const bool hasALB = isAvailable(available_params, MeteoGrids::ALB, date);
		if (!hasALB) return false;
		const bool hasISWR = isAvailable(available_params, MeteoGrids::ISWR, date);

		if (!hasISWR) {
			const bool hasISWR_DIFF = isAvailable(available_params, MeteoGrids::ISWR_DIFF, date);
			const bool hasISWR_DIR = isAvailable(available_params, MeteoGrids::ISWR_DIR, date);
			if (!hasISWR_DIFF || !hasISWR_DIR) return false;

			const Grid2DObject iswr_diff( getRawGrid(MeteoGrids::ISWR_DIFF, date) );
			grid2D = getRawGrid(MeteoGrids::ISWR_DIR, date);
			grid2D += iswr_diff;
			buffer.push(grid2D, MeteoGrids::ISWR, date);
		} else {
			grid2D = getRawGrid(MeteoGrids::ISWR, date);
		}

		const Grid2DObject alb( getRawGrid(MeteoGrids::ALB, date) );

		grid2D *= alb;
		buffer.push(grid2D, MeteoGrids::RSWR, date);
		return true;
	}

	if (parameter==MeteoGrids::ILWR) {
		const bool hasOLWR = isAvailable(available_params, MeteoGrids::OLWR, date);
		const bool hasLWR_NET = isAvailable(available_params, MeteoGrids::LWR_NET, date);

		if (hasOLWR && hasLWR_NET) {
			const Grid2DObject lwr_net( getRawGrid(MeteoGrids::LWR_NET, date) );
			grid2D = getRawGrid(MeteoGrids::OLWR, date);
			grid2D += lwr_net;
			buffer.push(grid2D, MeteoGrids::ILWR, date);
			return true;
		}
	}

	if (parameter==MeteoGrids::HS) {
		const bool hasRSNO = isAvailable(available_params, MeteoGrids::RSNO, date);
		const bool hasSWE = isAvailable(available_params, MeteoGrids::SWE, date);

		if (hasRSNO && hasSWE) {
			const Grid2DObject rsno( getRawGrid(MeteoGrids::RSNO, date) );
			grid2D = getRawGrid(MeteoGrids::SWE, date);
			grid2D *= 1000.0; //convert mm=kg/m^3 into kg
			grid2D /= rsno;
			buffer.push(grid2D, MeteoGrids::HS, date);
			return true;
		}
		return false;
	}

	if (parameter==MeteoGrids::PSUM) {
		const bool hasPSUM_S = isAvailable(available_params, MeteoGrids::PSUM_S, date);
		const bool hasPSUM_L = isAvailable(available_params, MeteoGrids::PSUM_L, date);
		const bool hasPSUM_LC = isAvailable(available_params, MeteoGrids::PSUM_LC, date);

		if (hasPSUM_S && hasPSUM_L && hasPSUM_LC) {
			const Grid2DObject psum_l( getRawGrid(MeteoGrids::PSUM_L, date) );
			const Grid2DObject psum_lc( getRawGrid(MeteoGrids::PSUM_LC, date) );
			grid2D = getRawGrid(MeteoGrids::PSUM_S, date);
			grid2D += psum_l;
			grid2D += psum_lc;
			buffer.push(grid2D, MeteoGrids::PSUM, date);
			return true;
		}

		if (hasPSUM_S && hasPSUM_L) {
			const Grid2DObject psum_l( getRawGrid(MeteoGrids::PSUM_L, date) );
			grid2D = getRawGrid(MeteoGrids::PSUM_S, date);
			grid2D += psum_l;
			buffer.push(grid2D, MeteoGrids::PSUM, date);
			return true;
		}
	}

	if (parameter==MeteoGrids::PSUM_PH) {
		const bool hasPSUM_S = isAvailable(available_params, MeteoGrids::PSUM_S, date);
		const bool hasPSUM_L = isAvailable(available_params, MeteoGrids::PSUM_L, date);
		const bool hasPSUM_LC = isAvailable(available_params, MeteoGrids::PSUM_LC, date);

		if (hasPSUM_S && hasPSUM_L && hasPSUM_LC) {
			const Grid2DObject psum_l( getRawGrid(MeteoGrids::PSUM_L, date) );
			const Grid2DObject psum_lc( getRawGrid(MeteoGrids::PSUM_LC, date) );
			grid2D = getRawGrid(MeteoGrids::PSUM_S, date);
			grid2D += psum_l + psum_lc;
			buffer.push(grid2D, MeteoGrids::PSUM, date);

			for (size_t ii=0; ii<grid2D.size(); ii++) {
				const double psum = grid2D(ii);
				if (psum!=IOUtils::nodata && psum>0)
					grid2D(ii) = (psum_l(ii) + psum_lc(ii)) / psum;
			}
			return true;
		}

		if (hasPSUM_S && hasPSUM_L) {
			const Grid2DObject psum_l( getRawGrid(MeteoGrids::PSUM_L, date) );
			grid2D = getRawGrid(MeteoGrids::PSUM_S, date);
			grid2D += psum_l;
			buffer.push(grid2D, MeteoGrids::PSUM, date);

			for (size_t ii=0; ii<grid2D.size(); ii++) {
				const double psum = grid2D(ii);
				if (psum!=IOUtils::nodata && psum>0)
					grid2D(ii) = psum_l(ii) / psum;
			}
			return true;
		}

		return false;
	}

	return false;
}

/**
* @brief Get the requested grid points from the grid, according to the configured processing level
* @details If the grid has been buffered, it will be returned from the buffer. If it is not available but can be generated, it will
* be generated transparently. If everything fails, it will thrown an exception.
* @param[in] parameter the parameter to get
* @param[in] date the timestamp that we would like to have
* @param[in] enforce_cartesian set to true to garantee a cartesian grid (it will be reprojected if necessary)
* @return a grid filled with the requested parameter
*/
std::vector < double > GridsManager::getPtsfromGrid(const MeteoGrids::Parameters& parameter, const Date& date, const std::vector< std::pair<size_t, size_t> >& Pts)
{
	std::vector<double> retVec;

	if (processing_level == IOUtils::raw){
		iohandler.readPointsIn2DGrid(retVec, parameter, date, Pts);
	} else {
		Grid2DObject grid2D;
		if (!buffer.get(grid2D, parameter, date)) {
			const bool status = setGrids2d_list( date ); //rebuffer the grid list if necessary
			if (!status) { //this means that the list2DGrids call is not implemeted in the plugin, we try to save the day...
				iohandler.readPointsIn2DGrid(retVec, parameter, date, Pts);
			} else { //the list2DGrids call is implemented in the plugin
				const std::map<Date, std::set<size_t> >::const_iterator it = grids2d_list.find(date);
				if (it!=grids2d_list.end()) {
					if ( it->second.find(parameter) != it->second.end() ) {
						iohandler.readPointsIn2DGrid(retVec, parameter, date, Pts);
					} else { //the right parameter could not be found, can we generate it?
						if (!getPtsfromgenerateGrid(retVec, it->second, parameter, date, Pts))
							throw NoDataException("Could not find or generate a grid of "+MeteoGrids::getParameterName( parameter )+" at time "+date.toString(Date::ISO), AT);
					}
				} else {
					const std::string msg1("Could not find grid for "+MeteoGrids::getParameterName( parameter )+" at time " + date.toString(Date::ISO) + ". " );
					const std::string msg2("There are grids from " + grids2d_list.begin()->first.toString(Date::ISO) + " until " + grids2d_list.rbegin()->first.toString(Date::ISO));
					throw NoDataException(msg1 + msg2, AT);
				}
			}
		} else {
			// If grid exists in buffer
			retVec = grid2D.extractPoints(Pts);
		}
	}

	return retVec;
}

/**
* @brief Generate point data from a grid that could be generated for a given parameter, based on the available parameters
* @details Even if a given parameter is not available, it might be possible to generate it
* on the fly based on the available data (for example, U and V wind components can be used to generate
* the VW and DW vector wind components).
*
* It is assumed that the meteo parameters are coming out of models, so the available_params are
* all available at all the timesteps, so we don't need to search a combination of parameters and timesteps
*
* @param[out] Vec a vector of doubles filled with the requested parameter or empty if it could not be generated
* @param available_params list of parameters available by a direct (=raw) read from the plugin
* @param parameter the parameter to get
* @param date the data associated with the parameter
* @return true if the requested grid could be generated
*
* @note HACK missing: checking that all grids used together have the same geolocalization
*/
bool GridsManager::getPtsfromgenerateGrid(std::vector<double>& Vec, const std::set<size_t>& available_params, const MeteoGrids::Parameters& parameter, const Date& date, const std::vector< std::pair<size_t, size_t> >& Pts)
{
	if (parameter==MeteoGrids::DEM) {
		if (!dem_altimeter) return false;
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		const bool hasP = isAvailable(available_params, MeteoGrids::P, date);
		const bool hasP_sea = isAvailable(available_params, MeteoGrids::P_SEA, date);
		if (!hasTA || !hasP || !hasP_sea) return false;

		std::vector<double> ta;
		iohandler.readPointsIn2DGrid(ta, MeteoGrids::TA, date, Pts);
		std::vector<double> p;
		iohandler.readPointsIn2DGrid(p, MeteoGrids::P, date, Pts);
		std::vector<double> p_sea;
		iohandler.readPointsIn2DGrid(p_sea, MeteoGrids::P_SEA, date, Pts);

		static const double k = Cst::gravity / (Cst::mean_adiabatique_lapse_rate * Cst::gaz_constant_dry_air);
		static const double k_inv = 1./k;
		for (size_t ii=0; ii<ta.size(); ii++) {
			if (p[ii]==IOUtils::nodata || ta[ii]==IOUtils::nodata || p_sea[ii]==IOUtils::nodata) {
				Vec.push_back(IOUtils::nodata);
			} else {
				const double K = pow(p[ii]/p_sea[ii], k_inv);
				Vec.push_back(ta[ii]*Cst::earth_R0*(1.-K) / (Cst::mean_adiabatique_lapse_rate * Cst::earth_R0 - ta[ii]*(1.-K)));
			}
		}
		return true;
	}

	if (parameter==MeteoGrids::VW || parameter==MeteoGrids::DW) {
		const bool hasU = isAvailable(available_params, MeteoGrids::U, date);
		const bool hasV = isAvailable(available_params, MeteoGrids::V, date);
		if (!hasU || !hasV) return false;

		std::vector<double> U;
		std::vector<double> V;
		iohandler.readPointsIn2DGrid(U, MeteoGrids::U, date, Pts);
		iohandler.readPointsIn2DGrid(V, MeteoGrids::V, date, Pts);

		if (parameter==MeteoGrids::VW) {
			for (size_t ii=0; ii<U.size(); ii++)
				Vec.push_back((U[ii]!=IOUtils::nodata && V[ii]!=IOUtils::nodata) ? (sqrt( Optim::pow2(U[ii]) + Optim::pow2(V[ii]) )) : (IOUtils::nodata));
		} else {
			for (size_t ii=0; ii<U.size(); ii++)
				Vec.push_back((U[ii]!=IOUtils::nodata && V[ii]!=IOUtils::nodata) ? (IOUtils::UV_TO_DW(U[ii], V[ii])) : (IOUtils::nodata)); // turn into degrees [0;360)
		}

		return true;
	}

	if (parameter==MeteoGrids::RH) {
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		if (!hasTA) return false;
		std::vector<double> TA;
		iohandler.readPointsIn2DGrid(TA, MeteoGrids::TA, date, Pts);

		const bool hasDEM = isAvailable(available_params, MeteoGrids::DEM, date);
		const bool hasQI = isAvailable(available_params, MeteoGrids::QI, date);
		const bool hasTD = isAvailable(available_params, MeteoGrids::TD, date);

		if (hasTD) {
			std::vector<double> TD;
			iohandler.readPointsIn2DGrid(TD, MeteoGrids::TD, date, Pts);

			for (size_t ii=0; ii<TD.size(); ii++)
				Vec.push_back(Atmosphere::DewPointtoRh(TD[ii], TA[ii], false));
			return true;
		} else if (hasQI && hasDEM) {
			std::vector<double> dem_pts;
			const Grid2DObject dem( getRawGrid(MeteoGrids::DEM, date) ); //HACK use readDEM instead?
			dem_pts = dem.extractPoints(Pts);
			std::vector<double> QI;
			iohandler.readPointsIn2DGrid(QI, MeteoGrids::QI, date, Pts);
			for (size_t ii=0; ii<QI.size(); ii++)
				Vec.push_back(Atmosphere::specToRelHumidity(dem_pts[ii], TA[ii], QI[ii]));
			return true;
		}

		return false;
	}

	if (parameter==MeteoGrids::ISWR) {
		const bool hasISWR_DIFF = isAvailable(available_params, MeteoGrids::ISWR_DIFF, date);
		const bool hasISWR_DIR = isAvailable(available_params, MeteoGrids::ISWR_DIR, date);

		if (hasISWR_DIFF && hasISWR_DIR) {
			std::vector<double> iswr_diff;
			iohandler.readPointsIn2DGrid(iswr_diff, MeteoGrids::ISWR_DIFF, date, Pts);
			std::vector<double> iswr_dir;
			iohandler.readPointsIn2DGrid(iswr_dir, MeteoGrids::ISWR_DIR, date, Pts);
			for (size_t ii=0; ii<iswr_diff.size(); ii++)
				Vec.push_back((iswr_diff[ii]!=IOUtils::nodata && iswr_dir[ii]!=IOUtils::nodata) ? (iswr_diff[ii] + iswr_dir[ii]) : (IOUtils::nodata));
			return true;
		}

		const bool hasALB = isAvailable(available_params, MeteoGrids::ALB, date);
		const bool hasRSWR = isAvailable(available_params, MeteoGrids::ISWR, date);
		if (hasRSWR && hasALB) {
			std::vector<double> rswr;
			iohandler.readPointsIn2DGrid(rswr, MeteoGrids::RSWR, date, Pts);
			std::vector<double> alb;
			iohandler.readPointsIn2DGrid(alb, MeteoGrids::ALB, date, Pts);
			for (size_t ii=0; ii<rswr.size(); ii++)
				Vec.push_back((rswr[ii]!=IOUtils::nodata && alb[ii]!=IOUtils::nodata) ? (rswr[ii] / alb[ii]) : (IOUtils::nodata));
			return true;
		}

		return false;
	}

	if (parameter==MeteoGrids::RSWR) {
		const bool hasALB = isAvailable(available_params, MeteoGrids::ALB, date);
		if (!hasALB) return false;
		const bool hasISWR = isAvailable(available_params, MeteoGrids::ISWR, date);

		std::vector<double> iswr;
		if (!hasISWR) {
			const bool hasISWR_DIFF = isAvailable(available_params, MeteoGrids::ISWR_DIFF, date);
			const bool hasISWR_DIR = isAvailable(available_params, MeteoGrids::ISWR_DIR, date);
			if (!hasISWR_DIFF || !hasISWR_DIR) return false;

			std::vector<double> iswr_diff;
			iohandler.readPointsIn2DGrid(iswr_diff, MeteoGrids::ISWR_DIFF, date, Pts);
			std::vector<double> iswr_dir;
			iohandler.readPointsIn2DGrid(iswr_dir, MeteoGrids::ISWR_DIR, date, Pts);
			for (size_t ii=0; ii<iswr_diff.size(); ii++)
				iswr.push_back((iswr_diff[ii]!=IOUtils::nodata && iswr_dir[ii]!=IOUtils::nodata) ? (iswr_diff[ii] + iswr_dir[ii]) : (IOUtils::nodata));
		} else {
			iohandler.readPointsIn2DGrid(iswr, MeteoGrids::ISWR, date, Pts);
		}

		std::vector<double> alb;
		iohandler.readPointsIn2DGrid(alb, MeteoGrids::ALB, date, Pts);
		for (size_t ii=0; ii<iswr.size(); ii++)
			Vec.push_back((iswr[ii]!=IOUtils::nodata && alb[ii]!=IOUtils::nodata) ? (iswr[ii] * alb[ii]) : (IOUtils::nodata));

		return true;
	}

	if (parameter==MeteoGrids::ILWR) {
		const bool hasOLWR = isAvailable(available_params, MeteoGrids::OLWR, date);
		const bool hasLWR_NET = isAvailable(available_params, MeteoGrids::LWR_NET, date);

		if (hasOLWR && hasLWR_NET) {
			std::vector<double> lwr_net;
			iohandler.readPointsIn2DGrid(lwr_net, MeteoGrids::LWR_NET, date, Pts);
			std::vector<double> olwr;
			iohandler.readPointsIn2DGrid(olwr, MeteoGrids::OLWR, date, Pts);
			for (size_t ii=0; ii<lwr_net.size(); ii++)
				Vec.push_back((lwr_net[ii]!=IOUtils::nodata && olwr[ii]!=IOUtils::nodata) ? (lwr_net[ii] + olwr[ii]) : (IOUtils::nodata));
			return true;
		}
	}

	if (parameter==MeteoGrids::HS) {
		const bool hasRSNO = isAvailable(available_params, MeteoGrids::RSNO, date);
		const bool hasSWE = isAvailable(available_params, MeteoGrids::SWE, date);

		if (hasRSNO && hasSWE) {
			std::vector<double> rsno;
			iohandler.readPointsIn2DGrid(rsno, MeteoGrids::RSNO, date, Pts);
			std::vector<double> swe;
			iohandler.readPointsIn2DGrid(swe, MeteoGrids::SWE, date, Pts);
			for (size_t ii=0; ii<rsno.size(); ii++)
				Vec.push_back((rsno[ii]!=IOUtils::nodata && swe[ii]!=IOUtils::nodata) ? ((1000.*swe[ii])/rsno[ii]) : (IOUtils::nodata));  //convert mm=kg/m^3 into kg
			return true;
		}
		return false;
	}

	if (parameter==MeteoGrids::PSUM) {
		const bool hasPSUM_S = isAvailable(available_params, MeteoGrids::PSUM_S, date);
		const bool hasPSUM_L = isAvailable(available_params, MeteoGrids::PSUM_L, date);
		const bool hasPSUM_LC = isAvailable(available_params, MeteoGrids::PSUM_LC, date);

		if (hasPSUM_S && hasPSUM_L && hasPSUM_LC) {
			std::vector<double> psum_l;
			iohandler.readPointsIn2DGrid(psum_l, MeteoGrids::PSUM_L, date, Pts);
			std::vector<double> psum_lc;
			iohandler.readPointsIn2DGrid(psum_lc, MeteoGrids::PSUM_LC, date, Pts);
			std::vector<double> psum_s;
			iohandler.readPointsIn2DGrid(psum_s, MeteoGrids::PSUM_S, date, Pts);
			for (size_t ii=0; ii<psum_l.size(); ii++)
				Vec.push_back((psum_l[ii]!=IOUtils::nodata && psum_lc[ii]!=IOUtils::nodata && psum_s[ii]!=IOUtils::nodata) ? (psum_l[ii] + psum_lc[ii] + psum_s[ii]) : (IOUtils::nodata));
			return true;
		}

		if (hasPSUM_S && hasPSUM_L) {
			std::vector<double> psum_l;
			iohandler.readPointsIn2DGrid(psum_l, MeteoGrids::PSUM_L, date, Pts);
			std::vector<double> psum_s;
			iohandler.readPointsIn2DGrid(psum_s, MeteoGrids::PSUM_S, date, Pts);
			for (size_t ii=0; ii<psum_l.size(); ii++)
				Vec.push_back((psum_l[ii]!=IOUtils::nodata && psum_s[ii]!=IOUtils::nodata) ? (psum_l[ii] + psum_s[ii]) : (IOUtils::nodata));
			return true;
		}
	}

	if (parameter==MeteoGrids::PSUM_PH) {
		const bool hasPSUM_S = isAvailable(available_params, MeteoGrids::PSUM_S, date);
		const bool hasPSUM_L = isAvailable(available_params, MeteoGrids::PSUM_L, date);

		if (hasPSUM_S && hasPSUM_L) {
			std::vector<double> psum_l;
			iohandler.readPointsIn2DGrid(psum_l, MeteoGrids::PSUM_L, date, Pts);
			std::vector<double> psum_s;
			iohandler.readPointsIn2DGrid(psum_s, MeteoGrids::PSUM_S, date, Pts);
			for (size_t ii=0; ii<psum_l.size(); ii++) {
				const double psum = ((psum_l[ii]!=IOUtils::nodata && psum_s[ii]!=IOUtils::nodata) ? (psum_l[ii] + psum_s[ii]) : (IOUtils::nodata));
				Vec.push_back((psum != IOUtils::nodata && psum > 0.) ? (psum_l[ii]/psum) : (IOUtils::nodata));
			}

			return true;
		}

		return false;
	}

	return false;
}


const std::string GridsManager::toString() const {
	ostringstream os;
	os << "<GridsManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler& iohandler = " << hex << &iohandler << dec << "\n";
	os << "Processing level = " << processing_level << "\n";
	os << buffer.toString();
	os << "</GridsManager>\n";
	return os.str();
}

} //namespace
