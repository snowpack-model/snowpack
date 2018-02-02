/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/Timer.h>

using namespace std;

namespace mio {

/**
 * @page spatial_resampling Spatial resampling handling
 * It is possible to use spatially interpolated meteorological fields or time series of 2D grids to extract meteorological time series for a set of points.
 * This is handled as "spatial resampling" and the data **will seem to originate from these virtual stations points** where no station is present. This obviously
 * comes at the cost of much higher run times. Several strategies are available (with the *RESAMPLING_STRATEGY* keyword, after setting RESAMPLING to TRUE):
 *     + VSTATIONS: points measurements are spatially interpolated at the chosen locations;
 *     + GRID_EXTRACT: gridded values are extracted for the cells containing the given locations;
 *     + GRID_ALL: all grid points are extracted.
 * 
 * Currently, it is necessary to provide a hint on how often the data should be extrated versus temporally interpolated between extracted point. This is described
 * by providing a refresh rate and an offset (both in seconds, with the VSTATIONS_REFRESH_RATE and VSTATIONS_REFRESH_OFFSET keywords, respectively)
 * \image html vstations_sampling.png "Resampling workflow"
 * \image latex vstations_sampling.eps "Resampling workflow" width=0.9\textwidth
 *
 * @section vstations VSTATIONS
 * The data from real input stations (as read by the plugin defined with the METEO key in the [input] section) is filtered/processed, temporally interpolated and 
 * spatially interpolated as defined in the configuration file. Then time series are reconstructed from these grids at a set of defined points (which will receive
 * station IDs such as <i>VIR#</i> for each station). This behavior is configured by the following keys (in the [Input] section):
 *    + RESAMPLING set to *true*;
 *    + RESAMPLING_STRATEGY set to VSTATIONS;
 *    + VSTATION# : provide the lat, lon and altitude or easting, northing and altitude for a virtual station (see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax);
 *    + VIRTUAL_PARAMETERS: list of MeteoData::Parameters that have to be interpolated to populate the virtual stations;
 *    + VSTATIONS_REFRESH_RATE: how often to rebuild the spatial interpolations, in seconds;
 *    + VSTATIONS_REFRESH_OFFSET: time offset to the stations' refresh rate, in seconds;
 *    + INTERPOL_USE_FULL_DEM: should the spatial interpolations be performed on the whole DEM? (this is necessary for some algorithms, for example WINSTAL);
 * 
 * Currently, a DEM also has to be provided since this will be used to retrieve the elevation, slope and azimuth of the virtual stations.
 * 
 * In the example provided below, 4 stations provide the original data that will be spatially interpolated at 2 points (or virtual stations, VIR1 and VIR2) for
 * 7 meteorological parameters. Every 6 hours, with starting offset of on hour, the original data will be spatially interpolated (so at 01:00, 07:00, 13:00 and 19:00).
 * Any data requested at other time steps will be temporally resampled from the spatially interpolated data.
 * @code
 * DEM = ARC
 * DEMFILE = ./input/surface-grids/davos.asc
 *
 * #here, the real data as measured by some stations 
 * METEO		= IMIS
 * DBNAME		= sdbo
 * DBUSER		= xxx
 * DBPASS		= xxx
 * STATION1	= *WFJ
 * STATION2	= STB2
 * STATION3	= WFJ2
 * STATION4	= *DAV
 * 
 * #here the locations where the data will be generated. The caller will only see these stations!
 * RESAMPLING = true
 * RESAMPLING_STRATEGY = VSTATIONS
 * VSTATION1 = latlon 46.793029 9.821343 1987
 * VSTATION2 = latlon 46.793031 9.831572 2300
 * Virtual_parameters = TA RH PSUM ILWR P VW RSWR
 * VSTATIONS_REFRESH_RATE = 21600
 * VSTATIONS_REFRESH_OFFSET = 3600
 * @endcode
 * 
 * @section grids_extract From gridded data
 * The meteorological time series are extracted from time series of user-provided grids. therefore a plugin for 2D grids must have been defined (with the GRID2D key in
 * the [Input] section). The following keys control this downscaling process:
 *    + RESAMPLING set to *true*;
 *    + RESAMPLING_STRATEGY set to either *GRID_EXTRACT* or *GRID_ALL*;
 *    + VSTATION# : provide the lat, lon and (optionally) the epsg code for a virtual station;
 *    + VIRTUAL_PARAMETERS: list of MeteoData::Parameters that have to be interpolated to populate the virtual stations;
 *    + VSTATIONS_REFRESH_RATE: how often to rebuild the spatial interpolations, in seconds;
 *    + VSTATIONS_REFRESH_OFFSET: time offset to the stations' refresh rate, in seconds;
 *
 * Currently, a DEM has to be provided in order to check the position of the stations and the consistency of the grids.
 * @code
 * DEM = ARC
 * DEMFILE = ./input/surface-grids/davos.asc
 *
 * #here, the real data as measured by some stations 
 * METEO		= IMIS
 * DBNAME		= sdbo
 * DBUSER		= xxx
 * DBPASS		= xxx
 * STATION1	= *WFJ
 * STATION2	= STB2
 * STATION3	= WFJ2
 * STATION4	= *DAV
 * 
 * #here the locations where the data will be generated. The caller will only see these stations!
 * RESAMPLING = true
 * RESAMPLING_STRATEGY = GRID_ALL
 * Virtual_parameters = TA RH PSUM ILWR P VW ISWR
 * VSTATIONS_REFRESH_RATE = 21600
 * VSTATIONS_REFRESH_OFFSET = 3600
 * @endcode
 * 
 * @section resampling_behind_the_scene Behind the scene
 * Behind the scene, this is a two stages setup: the IOManager that normally requests its data to its internal TimeSeriesManager will first ask the 
 * Meteo2DInterpolator to perform the spatial interpolations. Once this is done, the Meteo2DInterpolator (that used its own, internal TimeSeriesManager 
 * to get the real, measured data) pushes the computed timeseries into the raw data buffer of the IOManager's internal TimeSeriesManager. 
 * The IOManager then request the data as usual, from its TimeSeriesManager.
 * \image html vstations_workflow.png "Resampling internal workflow"
 * \image latex vstations_workflow.eps "Resampling internal workflow" width=0.9\textwidth
 * 
 */

Meteo2DInterpolator::Meteo2DInterpolator(const Config& i_cfg, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                    : cfg(i_cfg), tsmanager(&i_tsmanager), gridsmanager(&i_gridsmanager),
                      grid_buffer(0), internal_dem(), mapAlgorithms(),
                      v_params(), v_coords(), v_stations(), vstations_refresh_rate(1), vstations_refresh_offset(0.), resampling_strategy(NONE),
                      algorithms_ready(false), use_full_dem(false)
{
	std::string resampling_strategy_str( "NONE" );
	cfg.getValue("Resampling_strategy", "Input", resampling_strategy_str, IOUtils::nothrow);
	IOUtils::toUpper( resampling_strategy_str );
	if (resampling_strategy_str=="VSTATIONS") {
		resampling_strategy = VSTATIONS;
	} else if (resampling_strategy_str=="GRID_EXTRACT") {
		resampling_strategy = GRID_EXTRACT;
	} else if (resampling_strategy_str=="GRID_ALL") {
		resampling_strategy = GRID_ALL;
	} else if (resampling_strategy_str=="GRID_SMART") {
		resampling_strategy = GRID_SMART;
	}
	
	if (resampling_strategy!=NONE) {
		tsmanager = new TimeSeriesManager(i_tsmanager.getIOHandler(), i_cfg);
		gridsmanager = new GridsManager(i_gridsmanager.getIOHandler(), i_cfg);
		cfg.getValue("VSTATIONS_REFRESH_RATE", "Input", vstations_refresh_rate, IOUtils::nothrow);
		cfg.getValue("VSTATIONS_REFRESH_OFFSET", "Input", vstations_refresh_offset, IOUtils::nothrow);
	}

	size_t max_grids = 10; //default number of grids to keep in buffer
	cfg.getValue("BUFF_GRIDS", "Interpolations2D", max_grids, IOUtils::nothrow);
	grid_buffer.setMaxGrids(max_grids);
	
	setAlgorithms();
}

Meteo2DInterpolator::Meteo2DInterpolator(const Meteo2DInterpolator& source)
           : cfg(source.cfg), tsmanager(source.tsmanager), gridsmanager(source.gridsmanager),
                      grid_buffer(source.grid_buffer), internal_dem(source.internal_dem), mapAlgorithms(source.mapAlgorithms),
                      v_params(source.v_params), v_coords(source.v_coords), v_stations(source.v_stations), 
                      vstations_refresh_rate(1), vstations_refresh_offset(0.), resampling_strategy(source.resampling_strategy),
                      algorithms_ready(source.algorithms_ready), use_full_dem(source.use_full_dem)
{}

Meteo2DInterpolator& Meteo2DInterpolator::operator=(const Meteo2DInterpolator& source) {
	if (this != &source) {
		//cfg = source.cfg;
		tsmanager = source.tsmanager;
		gridsmanager = source.gridsmanager;
		grid_buffer = source.grid_buffer;
		internal_dem = source.internal_dem;
		mapAlgorithms = source.mapAlgorithms;
		v_params = source.v_params;
		v_coords = source.v_coords;
		v_stations = source.v_stations;
		vstations_refresh_rate = source.vstations_refresh_rate;
		vstations_refresh_offset = source.vstations_refresh_offset;
		resampling_strategy = source.resampling_strategy;
		algorithms_ready = source.algorithms_ready;
		use_full_dem = source.use_full_dem;
	}
	return *this;
}

Meteo2DInterpolator::~Meteo2DInterpolator()
{
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		const vector<InterpolationAlgorithm*>& vecAlgs( iter->second );
		for (size_t ii=0; ii<vecAlgs.size(); ++ii)
			delete vecAlgs[ii];
	}
	if (resampling_strategy!=NONE) {
		delete tsmanager;
		delete gridsmanager;
	}
}

/* By reading the Config object build up a list of user configured algorithms
* for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, psum, ...)
* Concept of this constructor: loop over all MeteoData::Parameters and then look
* for configuration of interpolation algorithms within the Config object.
*/
void Meteo2DInterpolator::setAlgorithms()
{
	const std::set<std::string> set_of_used_parameters( getParameters(cfg) );

	std::set<std::string>::const_iterator it  = set_of_used_parameters.begin();
	for (; it != set_of_used_parameters.end(); ++it) {
		const std::string parname( *it );
		const std::vector<std::string> tmpAlgorithms( getAlgorithmsForParameter(cfg, parname) );
		const size_t nrOfAlgorithms = tmpAlgorithms.size();

		std::vector<InterpolationAlgorithm*> vecAlgorithms( nrOfAlgorithms );
		for (size_t jj=0; jj<nrOfAlgorithms; jj++) {
			const std::vector< std::pair<std::string, std::string> > vecArgs( getArgumentsForAlgorithm(parname, tmpAlgorithms[jj], "Interpolations2D") );
			vecAlgorithms[jj] = AlgorithmFactory::getAlgorithm( tmpAlgorithms[jj], *this, vecArgs, *tsmanager, *gridsmanager, parname);
		}

		if (nrOfAlgorithms>0) {
			mapAlgorithms[parname] = vecAlgorithms;
		}
	}
	algorithms_ready = true;
}

//get a list of all meteoparameters referenced in the Interpolations2D section
std::set<std::string> Meteo2DInterpolator::getParameters(const Config& cfg)
{
	const std::vector<std::string> vec_keys( cfg.getKeys("::algorithms", "Interpolations2D", true) );

	std::set<std::string> set_parameters;
	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			if (vec_keys[ii].length()<=(found+2))
				throw InvalidFormatException("Invalid syntax: \""+vec_keys[ii]+"\"", AT);
			if (vec_keys[ii][found+1]!=':')
				throw InvalidFormatException("Missing ':' in \""+vec_keys[ii]+"\"", AT);
			std::string tmp( vec_keys[ii].substr(0,found) );
			IOUtils::toUpper( tmp );
			set_parameters.insert( tmp );
		}
	}

	return set_parameters;
}


std::vector<std::string> Meteo2DInterpolator::getAlgorithmsForParameter(const Config& cfg, const std::string& parname)
{
	// This function retrieves the user defined interpolation algorithms for
	// parameter 'parname' by querying the Config object
	std::vector<std::string> vecAlgorithms;
	const std::vector<std::string> vecKeys( cfg.getKeys(parname+"::algorithms", "Interpolations2D") );

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.empty()) return vecAlgorithms;

	cfg.getValue(vecKeys[0], "Interpolations2D", vecAlgorithms, IOUtils::nothrow);
	return vecAlgorithms;
}

std::vector< std::pair<std::string, std::string> > Meteo2DInterpolator::getArgumentsForAlgorithm(const std::string& parname,
                                                     const std::string& algorithm, const std::string& section) const
{
	const std::string key_prefix( parname+"::"+algorithm+"::" );
	std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getValues(key_prefix, section) );

	//clean the arguments up (ie remove the {Param}::{algo}:: in front of the argument key itself)
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		const size_t beg_arg_name = vecArgs[ii].first.find_first_not_of(":", key_prefix.length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+vecArgs[ii].first+"'", AT);
		vecArgs[ii].first = vecArgs[ii].first.substr(beg_arg_name);
	}

	return vecArgs;
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result)
{
	std::string InfoString;
	interpolate(date, dem, meteoparam, result, InfoString);
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result, std::string& InfoString)
{
	if (!algorithms_ready)
		setAlgorithms();

	const std::string param_name( MeteoData::getParameterName(meteoparam) );
	
	//Get grid from buffer if it exists
	std::ostringstream grid_hash;
	grid_hash << dem.llcorner.toString(Coords::LATLON) << " " << dem.getNx() << "x" << dem.getNy() << " @" << dem.cellsize << " " << date.toString(Date::ISO) << " " << param_name;
	if (grid_buffer.get(result, grid_hash.str(), InfoString))
		return;

	//Show algorithms to be used for this parameter
	const std::map<string, vector<InterpolationAlgorithm*> >::iterator it = mapAlgorithms.find(param_name);
	if (it==mapAlgorithms.end())
		throw IOException("No interpolation algorithms configured for parameter "+param_name, AT);

	//look for algorithm with the highest quality rating
	const std::vector<InterpolationAlgorithm*>& vecAlgs( it->second );
	double maxQualityRating = -1.;
	size_t bestalgorithm = 0;
	for (size_t ii=0; ii < vecAlgs.size(); ++ii){
		const double rating = vecAlgs[ii]->getQualityRating(date);
		if ((rating != 0.0) && (rating > maxQualityRating)) { //we use ">" so that in case of equality, the first choice will be kept
			bestalgorithm = ii;
			maxQualityRating = rating;
		}
	}

	//finally execute the algorithm with the best quality rating or throw an exception
	if (maxQualityRating<=0.0)
		throw IOException("No interpolation algorithm with quality rating >0 found for parameter "+param_name+" on "+date.toString(Date::ISO_TZ), AT);
	vecAlgs[bestalgorithm]->calculate(dem, result);
	InfoString = vecAlgs[bestalgorithm]->getInfo();

	//Run soft min/max filter for RH, PSUM and HS
	if (meteoparam == MeteoData::RH){
		Meteo2DInterpolator::checkMinMax(0.0, 1.0, result);
	} else if (meteoparam == MeteoData::PSUM){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::HS){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::VW){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	}

	//save grid in buffer
	grid_buffer.push(result, grid_hash.str(), InfoString);
}

//NOTE make sure that skip_virtual_stations = true before calling this method when using virtual stations!
void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	result.clear();
	std::vector<Coords> vec_coords( in_coords );

	if (use_full_dem) {
		Grid2DObject result_grid;
		interpolate(date, dem, meteoparam, result_grid, info_string);
		const bool gridify_success = dem.gridify(vec_coords);
		if (!gridify_success)
			throw InvalidArgumentException("Coordinate given to interpolate is outside of dem", AT);

		for (size_t ii=0; ii<vec_coords.size(); ii++) {
			//we know the i,j are positive because of gridify_success
			const size_t pt_i = static_cast<size_t>( vec_coords[ii].getGridI() );
			const size_t pt_j = static_cast<size_t>( vec_coords[ii].getGridJ() );
			result.push_back( result_grid(pt_i,pt_j) );
		}
	} else {
		for (size_t ii=0; ii<vec_coords.size(); ii++) {
			const bool gridify_success = dem.gridify(vec_coords[ii]);
			if (!gridify_success)
				throw InvalidArgumentException("Coordinate given to interpolate is outside of dem", AT);

			//we know the i,j are positive because of gridify_success
			const size_t pt_i = static_cast<size_t>( vec_coords[ii].getGridI() );
			const size_t pt_j = static_cast<size_t>( vec_coords[ii].getGridJ() );

			//Make new DEM with just one point, namely the one specified by vec_coord[ii]
			//Copy all other properties of the big DEM into the new one
			DEMObject one_point_dem(dem, pt_i, pt_j, 1, 1, false);

			one_point_dem.min_altitude = dem.min_altitude;
			one_point_dem.max_altitude = dem.max_altitude;
			one_point_dem.min_slope = dem.min_slope;
			one_point_dem.max_slope = dem.max_slope;
			one_point_dem.min_curvature = dem.min_curvature;
			one_point_dem.max_curvature = dem.max_curvature;

			Grid2DObject result_grid;
			interpolate(date, one_point_dem, meteoparam, result_grid, info_string);
			result.push_back(result_grid(0,0));
		}
	}
}

void Meteo2DInterpolator::checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj)
{
	for (size_t ii=0; ii<gridobj.size(); ii++){
		double& value = gridobj(ii);
		if (value == IOUtils::nodata) continue;

		if (value < minval) {
			value = minval;
		} else if (value > maxval) {
			value = maxval;
		}
	}
}

void Meteo2DInterpolator::check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo)
{
	//check that the stations are using the same projection as the dem
	for (size_t ii=0; ii<vec_meteo.size(); ii++) {
		const StationData& meta( vec_meteo[ii].meta );
		if (!meta.position.isSameProj(dem.llcorner)) {
			std::string type, args;
			meta.position.getProj(type, args);
			const std::string station_str( "Station "+meta.stationID+" is using projection ("+type+" "+args+") " );
			dem.llcorner.getProj(type, args);
			const std::string dem_str( "while DEM is using projection ("+type+" "+args+") " );
			throw IOException(station_str+dem_str, AT);
		}
	}
}

//////////////////////////// Starting from here, virtual stations / resampling stuff
void Meteo2DInterpolator::pushVirtualMeteoData(const Date& i_date, TimeSeriesManager &user_tsmanager)
{
	//find the nearest sampling points (vstations_refresh_rate apart) around the requested point
	const Date i_date_down( Date::rnd(i_date-vstations_refresh_offset, vstations_refresh_rate, Date::DOWN) + vstations_refresh_offset );
	const Date i_date_up( Date::rnd(i_date-vstations_refresh_offset, vstations_refresh_rate, Date::UP) + vstations_refresh_offset );
	const Date buff_start( user_tsmanager.getRawBufferStart() );
	const Date buff_end( user_tsmanager.getRawBufferEnd() );
	const double half_range = (vstations_refresh_rate)/(3600.*24.*2.);

	if (buff_start.isUndef() || i_date_down<buff_start || i_date_down>buff_end) {
		METEO_SET vecMeteo;
		getVirtualMeteoData(i_date_down, vecMeteo);
		user_tsmanager.push_meteo_data(IOUtils::raw, i_date_down - half_range, i_date_down + half_range, vecMeteo);
	}

	if (i_date_down!=i_date_up && (buff_start.isUndef() || i_date_up<buff_start || i_date_up>buff_end)) { //sometimes, up and down are rounded up to the same value...
		METEO_SET vecMeteo;
		getVirtualMeteoData(i_date_up, vecMeteo);
		user_tsmanager.push_meteo_data(IOUtils::raw, i_date_up - half_range, i_date_up + half_range, vecMeteo);
	}
}

size_t Meteo2DInterpolator::getVirtualMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	if (v_stations.empty()) {
		if (resampling_strategy==GRID_ALL) {
			initVirtualStationsAtAllGridPoints();
		} else {
			const bool adjust_coordinates = (resampling_strategy==GRID_EXTRACT);
			initVirtualStations(adjust_coordinates);
		}
	}
	
	if (resampling_strategy==VSTATIONS) {
		//this reads station data, interpolates the stations and extract points from the interpolated grids
		return getVirtualStationsData(i_date, vecMeteo);
	} else if (resampling_strategy==GRID_EXTRACT || resampling_strategy==GRID_ALL) {
		//This reads already gridded data and extract points from the grids at the provided locations
		//(the virtual stations must have been initialized before, see above)
		return getVirtualStationsFromGrid(i_date, vecMeteo);
	} else if (resampling_strategy==GRID_SMART) {
		
	}

	throw UnknownValueException("Unknown virtual station strategy", AT);
}

void Meteo2DInterpolator::initVirtualStationsAtAllGridPoints()
{
	if (!cfg.keyExists("DEM", "Input"))
		throw NoDataException("In order to use virtual stations, please provide a DEM!", AT);
	if (internal_dem.empty()) gridsmanager->readDEM(internal_dem);

	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	internal_dem.llcorner.setProj(coordin, coordinparam); //make sure the DEM and the VStations are in the same projection
	const double dem_easting = internal_dem.llcorner.getEasting();
	const double dem_northing = internal_dem.llcorner.getNorthing();
	const double cellsize = internal_dem.cellsize;
	
	Coords curr_point(coordin, coordinparam);
	size_t stat_id=0;
	for (size_t jj=0; jj<internal_dem.getNy(); jj++) {
		for (size_t ii=0; ii<internal_dem.getNx(); ii++) {
			const double easting = dem_easting + cellsize*static_cast<double>(ii);
			const double northing = dem_northing + cellsize*static_cast<double>(jj);
			curr_point.setXY(easting, northing, internal_dem(ii,jj));
			curr_point.setGridIndex(static_cast<int>(ii), static_cast<int>(jj), IOUtils::inodata, true);
			
			//extract vstation number, build the station name and station ID
			stat_id++;
			const std::string id_str( static_cast<ostringstream*>( &(ostringstream() << stat_id) )->str() );
			StationData sd(curr_point, "VIR"+id_str, "Virtual_Station_"+id_str);
			sd.setSlope(internal_dem.slope(ii,jj), internal_dem.azi(ii,jj));
			v_stations.push_back( sd );
		}
	}
}

/** @brief read the list of virtual stations
 * @param adjust_coordinates should the coordinates be recomputed to match DEM cells?
 */
void Meteo2DInterpolator::initVirtualStations(const bool& adjust_coordinates)
{
	if (!cfg.keyExists("DEM", "Input"))
		throw NoDataException("In order to use virtual stations, please provide a DEM!", AT);
	if (internal_dem.empty()) gridsmanager->readDEM(internal_dem);

	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	internal_dem.llcorner.setProj(coordin, coordinparam); //make sure the DEM and the VStations are in the same projection
	const double dem_easting = internal_dem.llcorner.getEasting();
	const double dem_northing = internal_dem.llcorner.getNorthing();
	const double cellsize = internal_dem.cellsize;

	//read the provided coordinates, remove duplicates and generate metadata
	const std::vector< std::pair<std::string, std::string> > vecStation( cfg.getValues("Vstation", "INPUT") );
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		Coords curr_point(coordin, coordinparam, vecStation[ii].second);

		if (!curr_point.isNodata()) {
			v_coords.push_back( curr_point ); //so we can check for duplicates

			if (!internal_dem.gridify(curr_point)) {
				ostringstream ss;
				ss << "Virtual station \"" << vecStation[ii].second << "\" is not contained in provided DEM " << internal_dem.toString(DEMObject::SHORT);
				throw NoDataException(ss.str(), AT);
			}

			const size_t i = curr_point.getGridI(), j = curr_point.getGridJ();
			if (adjust_coordinates) { //adjust coordinates to match the chosen cell
				const double easting = dem_easting + cellsize*static_cast<double>(i);
				const double northing = dem_northing + cellsize*static_cast<double>(j);
				curr_point.setXY(easting, northing, internal_dem(i,j));
				curr_point.setGridIndex(static_cast<int>(i), static_cast<int>(j), IOUtils::inodata, true);
			} else {
				curr_point.setAltitude(internal_dem(i,j), false);
			}

			//remove duplicate stations, ie stations that have same easting,northing,altitude
			bool is_duplicate = false;
			for (size_t jj=0; jj<ii; jj++) {
				if (curr_point==v_coords[jj]) {
					std::cout << "[W] removing VSTATION" << ii+1 << " as a duplicate of VSTATION" << jj+1 << "\n";
					is_duplicate = true;
					break;
				}
			}
			if (!is_duplicate) {
				//extract vstation number, build the station name and station ID
				const std::string id_num( vecStation[ii].first.substr(string("Vstation").length()) );
				StationData sd(curr_point, "VIR"+id_num, "Virtual_Station_"+id_num);
				sd.setSlope(internal_dem.slope(i,j), internal_dem.azi(i,j));
				v_stations.push_back( sd );
			}
		}
	}

	if (v_stations.empty())
		throw NoDataException("No virtual stations provided", AT);
	
	//cleaning up v_coords so it is still usable for the VSTATIONS (ie no more duplicates)
	v_coords.resize(v_stations.size());
	for (size_t ii=0; ii<v_stations.size(); ii++) {
		v_coords[ii] = v_stations[ii].position;
	}
}

size_t Meteo2DInterpolator::getVirtualStationsMeta(const Date& /*date*/, STATIONS_SET& vecStation)
{
	vecStation = v_stations;
	return vecStation.size();
}

size_t Meteo2DInterpolator::getVirtualStationsData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	//get data from real input stations
	METEO_SET vecTrueMeteo;
	tsmanager->getMeteoData(i_date, vecTrueMeteo);
	if (vecTrueMeteo.empty()) return 0;

	if (v_params.empty()) {
		//get parameters to interpolate if not already done
		//we need valid data in order to handle extra parameters
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for (size_t ii=0; ii<vecStr.size(); ii++) {
			v_params.push_back( vecTrueMeteo[0].getParameterIndex(vecStr[ii]) );
		}
	}

	//create stations without measurements
	for (size_t ii=0; ii<v_stations.size(); ii++) {
		MeteoData md(i_date, v_stations[ii]);
		vecMeteo.push_back( md );
	}

	//fill meteo parameters
	if (internal_dem.empty()) gridsmanager->readDEM(internal_dem); //this is not a big deal since it will be in the buffer
	
	std::string info_string;
	for (size_t param=0; param<v_params.size(); param++) {
		std::vector<double> result;
		interpolate(i_date, internal_dem, static_cast<MeteoData::Parameters>(v_params[param]), v_coords, result, info_string);
		for (size_t ii=0; ii<v_coords.size(); ii++) {
			vecMeteo[ii](v_params[param]) = result[ii];
		}
	}

	return vecMeteo.size();
}

/**
* @brief Extract time series from grids at the specified points (virtual stations).
* @param i_date when to extract the virtual stations
* @param vecMeteo a vector of meteodata for the configured virtual stations
*/
size_t Meteo2DInterpolator::getVirtualStationsFromGrid(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	if (v_params.empty()) { //get parameters to interpolate if not already done
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for (size_t ii=0; ii<vecStr.size(); ii++) {
			const size_t param_idx = MeteoGrids::getParameterIndex( vecStr[ii] );
			if (param_idx==IOUtils::npos)
				throw InvalidArgumentException("Invalid parameter '" + vecStr[ii] + "', only standard parameters can be extracted from grids for virtual stations! ", AT);
			v_params.push_back( param_idx );
		}
	}

	//create stations without measurements
	for (size_t ii=0; ii<v_stations.size(); ii++) {
		MeteoData md(i_date, v_stations[ii]);
		vecMeteo.push_back( md );
	}
	
	if (internal_dem.empty()) gridsmanager->readDEM(internal_dem); //this is not a big deal since it will be in the buffer
	
	for (size_t param=0; param<v_params.size(); param++) { //loop over required parameters
		const MeteoGrids::Parameters grid_param = static_cast<MeteoGrids::Parameters>(v_params[param]);
		Grid2DObject grid;
		gridsmanager->read2DGrid(grid, grid_param, i_date);
		
		if (!grid.isSameGeolocalization(internal_dem))
			throw InvalidArgumentException("In GRID_EXTRACT, the DEM and the source grid don't match for '"+MeteoGrids::getParameterName(grid_param)+"' on "+i_date.toString(Date::ISO));
		
		for (size_t ii=0; ii<v_stations.size(); ii++) { //loop over all virtual stations
			const size_t grid_i = v_stations[ii].position.getGridI(); //this should work since invalid stations have been removed in init
			const size_t grid_j = v_stations[ii].position.getGridJ();
			
			//check if this is a standard MeteoData parameter
			const size_t meteo_param = vecMeteo[ii].getParameterIndex( MeteoGrids::getParameterName(grid_param) ); //is this name also a meteoparameter?
			if (meteo_param!=IOUtils::npos)
				vecMeteo[ii]( static_cast<MeteoData::Parameters>(meteo_param) ) = grid(grid_i, grid_j);
		}
	}
	
	return vecMeteo.size();
}

const std::string Meteo2DInterpolator::toString() const
{
	ostringstream os;
	os << "<Meteo2DInterpolator>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "TimeSeriesManager& tsmanager = "  << hex << &tsmanager << dec << "\n";
	os << "GridsManager& gridsmanager = "  << hex << &gridsmanager << dec << "\n";

	os << "Spatial resampling algorithms:\n";
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::const_iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		os << setw(10) << iter->first << "::";
		for (size_t jj=0; jj<iter->second.size(); jj++) {
			os << iter->second[jj]->algo << " ";
		}
		os << "\n";
	}

	//cache content
	os << grid_buffer.toString();
	os << "</Meteo2DInterpolator>\n";
	return os.str();
}

} //namespace
