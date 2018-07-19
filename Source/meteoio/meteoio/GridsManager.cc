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

#include <meteoio/GridsManager.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/MathOptim.h>

using namespace std;

namespace mio {

GridsManager::GridsManager(IOHandler& in_iohandler, const Config& in_cfg)
             : iohandler(in_iohandler), cfg(in_cfg), buffer(0), grids2d_list(), grids2d_start(), grids2d_end(),
               grid2d_list_buffer_size(370.), processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated), dem_altimeter(false)
{
	size_t max_grids = 10;
	cfg.getValue("BUFF_GRIDS", "General", max_grids, IOUtils::nothrow);
	buffer.setMaxGrids(max_grids);
	cfg.getValue("BUFFER_SIZE", "General", grid2d_list_buffer_size, IOUtils::nothrow);
	cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow); //HACK document it! if no dem is found but local and sea level pressure grids are found, use them to rebuild a DEM; [Input] section
}

/**
* @brief Set the desired ProcessingLevel
*        The processing level affects the way meteo data is read and processed
*        Three values are possible:
*        - IOUtils::raw data shall be read directly from the buffer
*        - IOUtils::filtered data shall be filtered before returned to the user
*        - IOUtils::resampled data shall be resampled before returned to the user
*          this only affects the function getMeteoData(const Date&, METEO_DATASET&);
*
*        The three values can be combined: e.g. IOUtils::filtered | IOUtils:resampled
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

void GridsManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, filename);
	} else {
		if (buffer.get(grid2D, filename)) return;

		iohandler.read2DGrid(grid2D, filename);
		buffer.push(grid2D, filename);
	}
}

//do we have the requested parameter either in buffer or raw?
bool GridsManager::isAvailable(const std::set<size_t>& available_params, const MeteoGrids::Parameters& parameter, const Date& date) const
{
	const bool in_buffer = buffer.has(parameter, date);
	if (in_buffer) return true;
	
	const bool in_raw = (available_params.find( parameter ) != available_params.end());
	return in_raw;
}

//get the requested grid either from buffer or by a direct reading
//the goal is that subsequent calls to the same grid find it in the buffer transparently
void GridsManager::getGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (buffer.get(grid2D, parameter, date)) return;
	iohandler.read2DGrid(grid2D, parameter, date);
}

//Grids of meteo params are coming out of models, so we assume that all possible parameters
//are available at each time steps. So we don't need to search a combination of parameters and timesteps
bool GridsManager::read2DGrid(Grid2DObject& grid2D, const std::set<size_t>& available_params, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (parameter==MeteoGrids::DEM) {
		if (!dem_altimeter) return false;
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		const bool hasP = isAvailable(available_params, MeteoGrids::P, date);
		const bool hasP_sea = isAvailable(available_params, MeteoGrids::P_SEA, date);
		if (!hasTA || !hasP || !hasP_sea) return false;

		Grid2DObject ta, p, p_sea;
		getGrid(ta, MeteoGrids::TA, date);
		getGrid(p, MeteoGrids::P, date);
		getGrid(p_sea, MeteoGrids::P_SEA, date);

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
		
		Grid2DObject U,V;
		getGrid(U, MeteoGrids::U, date);
		buffer.push(U, MeteoGrids::U, date);
		getGrid(V, MeteoGrids::V, date);
		buffer.push(V, MeteoGrids::V, date);
		
		grid2D.set(U, IOUtils::nodata);
		if (parameter==MeteoGrids::VW) {
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = sqrt( Optim::pow2(U(ii)) + Optim::pow2(V(ii)) );
		} else {
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) =  fmod( atan2( U(ii), V(ii) ) * Cst::to_deg + 360., 360.); // turn into degrees [0;360)
		}
		buffer.push(grid2D, parameter, date);
		return true;
	}
	
	if (parameter==MeteoGrids::RH) {
		const bool hasTA = isAvailable(available_params, MeteoGrids::TA, date);
		if (!hasTA) return false;
		Grid2DObject TA;
		getGrid(TA, MeteoGrids::TA, date);
		buffer.push(TA, MeteoGrids::TA, date);
		
		const bool hasDEM = isAvailable(available_params, MeteoGrids::DEM, date);
		const bool hasQI = isAvailable(available_params, MeteoGrids::QI, date);
		const bool hasTD = isAvailable(available_params, MeteoGrids::TD, date);
		
		if (hasTA && hasTD) {
			getGrid(grid2D, MeteoGrids::TD, date);
			buffer.push(grid2D, MeteoGrids::TD, date);
			
			for (size_t ii=0; ii<grid2D.size(); ii++)
				grid2D(ii) = Atmosphere::DewPointtoRh(grid2D(ii), TA(ii), false);
			return true;
		} else if (hasQI && hasDEM && hasTA) {
			Grid2DObject dem;
			getGrid(dem, MeteoGrids::DEM, date); //HACK use readDEM instead?
			getGrid(grid2D, MeteoGrids::QI, date);
			buffer.push(grid2D, MeteoGrids::QI, date);
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
			Grid2DObject iswr_diff;
			getGrid(iswr_diff, MeteoGrids::ISWR_DIFF, date);
			getGrid(grid2D, MeteoGrids::ISWR_DIR, date);
			grid2D += iswr_diff;
			buffer.push(grid2D, MeteoGrids::ISWR, date);
			return true;
		}
		
		const bool hasALB = isAvailable(available_params, MeteoGrids::ALB, date);
		const bool hasRSWR = isAvailable(available_params, MeteoGrids::ISWR, date);
		if (hasRSWR && hasALB) {
			Grid2DObject alb;
			getGrid(alb, MeteoGrids::ALB, date);
			getGrid(grid2D, MeteoGrids::RSWR, date);
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
			
			Grid2DObject iswr_diff;
			getGrid(iswr_diff, MeteoGrids::ISWR_DIFF, date);
			getGrid(grid2D, MeteoGrids::ISWR_DIR, date);
			grid2D += iswr_diff;
			buffer.push(grid2D, MeteoGrids::ISWR, date);
		} else {
			getGrid(grid2D, MeteoGrids::ISWR, date);
		}
		
		Grid2DObject alb;
		getGrid(alb, MeteoGrids::ALB, date);
		
		grid2D *= alb;
		buffer.push(grid2D, MeteoGrids::RSWR, date);
		return true;
	}
	
	if (parameter==MeteoGrids::HS) {
		const bool hasRSNO = isAvailable(available_params, MeteoGrids::RSNO, date);
		const bool hasSWE = isAvailable(available_params, MeteoGrids::SWE, date);
		
		if (hasRSNO && hasSWE) {
			Grid2DObject rsno;
			getGrid(rsno, MeteoGrids::RSNO, date);
			getGrid(grid2D, MeteoGrids::SWE, date);
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
		
		if (hasPSUM_S && hasPSUM_L) {
			Grid2DObject psum_l;
			getGrid(grid2D, MeteoGrids::PSUM_S, date);
			buffer.push(grid2D, MeteoGrids::PSUM_S, date);
			getGrid(psum_l, MeteoGrids::PSUM_L, date);
			buffer.push(psum_l, MeteoGrids::PSUM_L, date);
			grid2D += psum_l;
			buffer.push(grid2D, MeteoGrids::PSUM, date);
			return true;
		}
	}
	
	if (parameter==MeteoGrids::PSUM_PH) {
		const bool hasPSUM_S = isAvailable(available_params, MeteoGrids::PSUM_S, date);
		const bool hasPSUM_L = isAvailable(available_params, MeteoGrids::PSUM_L, date);
		
		if (hasPSUM_S && hasPSUM_L) {
			Grid2DObject psum_l;
			getGrid(grid2D, MeteoGrids::PSUM_S, date);
			buffer.push(grid2D, MeteoGrids::PSUM_S, date);
			getGrid(psum_l, MeteoGrids::PSUM_L, date);
			buffer.push(psum_l, MeteoGrids::PSUM_L, date);
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

void GridsManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, parameter, date);
	} else {
		if (buffer.get(grid2D, parameter, date)) return;
		
		const bool status = setGrids2d_list( date ); //rebuffer the grid list if necessary
		if (!status) { //this means that the list2DGrids call is not implemeted in the plugin, we try to save the day...
			iohandler.read2DGrid(grid2D, parameter, date);
			buffer.push(grid2D, parameter, date);
			return;
		}
		
		const std::map<Date, std::set<size_t> >::const_iterator it = grids2d_list.find(date);
		if (it!=grids2d_list.end()) {
			if ( it->second.find(parameter) != it->second.end() ) {
				iohandler.read2DGrid(grid2D, parameter, date);
				buffer.push(grid2D, parameter, date);
			} else { //the right parameter could not be found, can we generate it?
				if (!read2DGrid(grid2D, it->second, parameter, date))
					throw NoDataException("Could not find or generate a grid of "+MeteoGrids::getParameterName( parameter )+" at time "+date.toString(Date::ISO), AT);
			}
		} else {
			const std::string msg1("Could not find any grids at time " + date.toString(Date::ISO) + ". " );
			const std::string msg2("There are grids from " + grids2d_list.begin()->first.toString(Date::ISO) + " until " + grids2d_list.rbegin()->first.toString(Date::ISO));
			throw NoDataException(msg1 + msg2, AT);
		}
	}
}

//HACK buffer 3D grids!
void GridsManager::read3DGrid(Grid3DObject& grid_out, const std::string& i_filename)
{
	iohandler.read3DGrid(grid_out, i_filename);
}

void GridsManager::read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.read3DGrid(grid_out, parameter, date);
}

void GridsManager::readDEM(DEMObject& grid2D)
{
	//TODO: dem_altimeter; reading DEM data with no associated date (ie Date()) OR with associated date (for example, TLS data)
	if (processing_level == IOUtils::raw){
		iohandler.readDEM(grid2D);
	} else {
		if (buffer.get(grid2D, "/:DEM"))
			return;

		iohandler.readDEM(grid2D);
		buffer.push(grid2D, "/:DEM");
	}
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

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const std::string& name)
{
	iohandler.write2DGrid(grid2D, name);
}

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write2DGrid(grid2D, parameter, date);
}

void GridsManager::write3DGrid(const Grid3DObject& grid_out, const std::string& options)
{
	iohandler.write3DGrid(grid_out, options);
}

void GridsManager::write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write3DGrid(grid_out, parameter, date);
}

/** @brief Create a list of virtual stations from each grid point
 * @details Please not that the two options are mutually exclusive.
 * @param[in] dem the Digital Elevation Model to use (each grid point will be made into a virtual station)
 * @return a vector of virtual stations
 */
std::vector<StationData> GridsManager::initVirtualStationsAtAllGridPoints(const DEMObject& dem) const
{
	std::vector<StationData> v_stations;
	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
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
			const std::string id_str( static_cast<ostringstream*>( &(ostringstream() << stat_id) )->str() );
			StationData sd(curr_point, "VIR"+id_str, "Virtual_Station_"+id_str);
			sd.setSlope(dem.slope(ii,jj), dem.azi(ii,jj));
			v_stations.push_back( sd );
		}
	}
	
	return v_stations;
}

/** @brief Create a list of virtual stations from the user-provided input
 * @details Please not that the two options are mutually exclusive.
 * @param[in] dem the Digital Elevation Model to use
 * @param[in] adjust_coordinates should the coordinates be recomputed to match DEM cells?
 * @param[in] fourNeighbors pick the surrounding four nodes instead of only the exact one?
 * @return a vector of virtual stations
 */
std::vector<StationData> GridsManager::initVirtualStations(const DEMObject& dem, const bool& adjust_coordinates, const bool& fourNeighbors) const
{
	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	Coords llcorner( dem.llcorner );
	llcorner.setProj(coordin, coordinparam); //make sure the DEM and the VStations are in the same projection
	const double dem_easting = llcorner.getEasting();
	const double dem_northing = llcorner.getNorthing();

	//read the provided coordinates, remove duplicates and generate metadata
	std::vector<StationData> v_stations;
	const std::vector< std::pair<std::string, std::string> > vecStation( cfg.getValues("Vstation", "INPUT") );
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		if (vecStation[ii].first.find('_') != std::string::npos) continue; //so we skip the other vstations_xxx parameters
		
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		Coords curr_point(coordin, coordinparam, vecStation[ii].second);
		if (curr_point.isNodata()) continue;

		if (!dem.gridify(curr_point))
			throw NoDataException("Virtual station "+vecStation[ii].second+" is not contained in provided DEM "+dem.toString(DEMObject::SHORT), AT);

		size_t i = curr_point.getGridI(), j = curr_point.getGridJ();
		if (fourNeighbors) { //pick the surrounding four nodes
			if (i==dem.getNx()) i--;
			if (j==dem.getNy()) j--;
			
			const std::string id_num( vecStation[ii].first.substr(string("Vstation").length()) );
			for (size_t mm=i; mm<=i+1; mm++) {
				for (size_t nn=j; nn<=j+1; nn++) {
					const double easting = dem_easting + dem.cellsize*static_cast<double>(mm);
					const double northing = dem_northing + dem.cellsize*static_cast<double>(nn);
					curr_point.setXY(easting, northing, dem(mm,nn));
					curr_point.setGridIndex(static_cast<int>(mm), static_cast<int>(nn), IOUtils::inodata, true);
					const std::string sub_id( static_cast<ostringstream*>( &(ostringstream() << (mm-i)*2+(nn-j)+1) )->str() );
					const std::string grid_pos( static_cast<ostringstream*>( &(ostringstream() << mm << "-" << nn) )->str() );
					StationData sd(curr_point, "VIR"+id_num+"_"+sub_id, "Virtual_Station_"+grid_pos);
					sd.setSlope(dem.slope(mm,nn), dem.azi(mm,nn));
					v_stations.push_back( sd );
				}
			}
		} else { //pick the exact node
			if (adjust_coordinates) { //adjust coordinates to match the chosen cell
				const double easting = dem_easting + dem.cellsize*static_cast<double>(i);
				const double northing = dem_northing + dem.cellsize*static_cast<double>(j);
				curr_point.setXY(easting, northing, dem(i,j));
				curr_point.setGridIndex(static_cast<int>(i), static_cast<int>(j), IOUtils::inodata, true);
			} else {
				curr_point.setAltitude(dem(i,j), false);
			}
			
			//extract vstation number, build the station name and station ID
			const std::string id_num( vecStation[ii].first.substr(string("Vstation").length()) );
			StationData sd(curr_point, "VIR"+id_num, "Virtual_Station_"+id_num);
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
* @return a vector of meteodata for the configured virtual stations at the provided date, for the provided parameters
*/
METEO_SET GridsManager::getVirtualStationsFromGrid(const DEMObject& dem, const std::vector<size_t>& v_params, const std::vector<StationData>& v_stations, const Date& date)
{
	//HACK handle extra parameters when possible
	const size_t nrStations = v_stations.size();
	METEO_SET vecMeteo;
	
	//create stations without measurements
	for (size_t ii=0; ii<nrStations; ii++) {
		MeteoData md(date, v_stations[ii]);
		vecMeteo.push_back( md );
	}
	
	for (size_t param=0; param<v_params.size(); param++) { //loop over required parameters
		const MeteoGrids::Parameters grid_param = static_cast<MeteoGrids::Parameters>( v_params[param] );
		Grid2DObject grid;
		read2DGrid(grid, grid_param, date);
		
		if (!grid.isSameGeolocalization(dem))
			throw InvalidArgumentException("In GRID_EXTRACT, the DEM and the source grid don't match for '"+MeteoGrids::getParameterName(grid_param)+"' on "+date.toString(Date::ISO));
		
		for (size_t ii=0; ii<nrStations; ii++) { //loop over all virtual stations
			const size_t grid_i = v_stations[ii].position.getGridI(); //this should work since invalid stations have been removed in init
			const size_t grid_j = v_stations[ii].position.getGridJ();
			
			//check if this is a standard MeteoData parameter
			const size_t meteo_param = vecMeteo[ii].getParameterIndex( MeteoGrids::getParameterName(grid_param) ); //is this name also a meteoparameter?
			if (meteo_param!=IOUtils::npos)
				vecMeteo[ii]( static_cast<MeteoData::Parameters>(meteo_param) ) = grid(grid_i, grid_j);
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
* @return a vector of meteodata for the configured virtual stations at the provided date, for the provided parameters
*/
std::vector<METEO_SET> GridsManager::getVirtualStationsFromGrid(const DEMObject& dem, const std::vector<size_t>& v_params, const std::vector<StationData>& v_stations, const Date& dateStart, const Date& dateEnd)
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
		const METEO_SET vecMeteo( getVirtualStationsFromGrid(dem, v_params, v_stations, it->first) ); //the number of stations can not change
		for (size_t ii=0; ii<nrStations; ii++) 
			vecvecMeteo[ii].push_back( vecMeteo[ii] );
		
		if (it->first>=dateEnd) break;
	}
	
	return vecvecMeteo;
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
