/***********************************************************************************/
/*  Copyright 2018-2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <alpine3d/TechSnowA3D.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/AlpineMain.h>

#include <errno.h>

#include <algorithm> // for std::find
#include <iterator>

using namespace std;
using namespace mio;

/**
 * @brief Reading the slope file and slope conditions.
 * @details It reads in the slope file and the defined slope conditions for three different priorities
 * @param[in] cfg is used to get 'TechSnow' conditions for snow production and reading the slope conditions file, defined in a separate file
 * @param[in] dem file to check whether the geolocalization of the slope file is the same
 */
TechSnowA3D::TechSnowA3D(const mio::Config& cfg, const mio::DEMObject& dem) 
              : skiRunsMap(), psum_tech(dem, IOUtils::nodata), slope_conditions(),
              start_season(), end_season(), earliest_production(), priorities(),
              max_snowgun_water(), mean_elevation(dem.grid2D.getMean()), slope_open(), slope_closed()
{
	const double TZ = cfg.get("TIME_ZONE", "Input");
	
	const std::string mapFile = cfg.get("SKIRUNS_FILE", "TechSnow");
	const bool isMaster = MPIControl::instance().master();
	if (isMaster) {
		mio::IOManager io( cfg );
		io.read2DGrid(skiRunsMap, "/" + mapFile);
		if (!skiRunsMap.isSameGeolocalization( dem )) 
			throw InvalidArgumentException("The ski runs map does not has the same geolocalization as the DEM", AT);
	}
	MPIControl::instance().broadcast( skiRunsMap );
	
	// Read slope conditions
	setSlopeConditions(cfg.get("SLOPE_CONDITIONS", "TechSnow"));
	cfg.getValue("SEASON_OPENING", "TechSnow", start_season, TZ);
	cfg.getValue("SEASON_CLOSING", "TechSnow", end_season, TZ);
	cfg.getValue("SLOPE_OPEN", "TechSnow", slope_open);
	cfg.getValue("SLOPE_CLOSED", "TechSnow", slope_closed);
	cfg.getValue("MAX_SNOWGUN_WATER", "TechSnow", max_snowgun_water);

	for (unsigned int ii=1; ii<10000; ii++) {
		const std::string prio_key( "PRIO"+mio::IOUtils::toString(ii)+"::START_PROD" );
		if (!cfg.keyExists(prio_key, "TechSnow") ) break; //stop looking for other priorities
		
		priorities.push_back( setSnowStrategy(cfg, TZ, ii) );
		if (earliest_production.isUndef()) earliest_production = priorities.back().startProd;
		if (earliest_production > priorities.back().startProd) earliest_production = priorities.back().startProd;
	}
}

TechSnowA3D::snowStrategy TechSnowA3D::setSnowStrategy(const mio::Config& cfg, const double& TZ, const unsigned int& nr)
{
	const std::string root_key( "PRIO"+mio::IOUtils::toString(nr)+"::" );
	TechSnowA3D::snowStrategy ppt;
	
	cfg.getValue(root_key+"start_prod", "TechSnow", ppt.startProd, TZ);
	cfg.getValue(root_key+"end_prod", "TechSnow", ppt.endProd, TZ);
	cfg.getValue(root_key+"start_aim", "TechSnow", ppt.startAim);
	cfg.getValue(root_key+"end_aim", "TechSnow", ppt.endAim);
	cfg.getValue(root_key+"gun_operation", "TechSnow", ppt.gunOperation);
	return ppt;
}

/**
 * @brief Creating a vector to store all conditions for the defined slope sections.
 * @details Creating a vector to store all conditions (slope number, slope area, number of snowguns, priority, minimum snow height, wetbulb treshold)
 * for each slope section.
 * @param[in] filename of the file containing the slope conditions for the defined slope sections.
*/
void TechSnowA3D::setSlopeConditions(const std::string& filename)
{
	//first, compute the surface of each slope
	for (size_t ii=0; ii<skiRunsMap.size(); ii++) {
		if (skiRunsMap(ii)==mio::IOUtils::nodata) continue;
		slope_conditions[ getSlopeNumber(skiRunsMap(ii)) ].slope_area++;
	}
	
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "error opening slope conditions file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	
	const char eoln = mio::FileUtils::getEoln(fin); //get the end of line character for the file
	size_t lcount=0;
	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			lcount++;
			mio::IOUtils::stripComments(line);
			mio::IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);

			size_t slope_number;
			unsigned int number_snowguns, priority;
			double min_height, wet_bulb_thresh;
			iss >> std::skipws >> slope_number;
			iss >> std::skipws >> number_snowguns;
			iss >> std::skipws >> priority;
			iss >> std::skipws >> min_height;
			iss >> std::skipws >> wet_bulb_thresh;
			
			if ( !iss || !iss.eof()) {
				std::ostringstream ss;
				ss << "TechSnow: invalid line in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			
			slope_conditions[ slope_number ].setUserProperties(number_snowguns, priority, min_height, wet_bulb_thresh);
		} while (!fin.eof());
		
		fin.close();
	} catch (const std::exception&){
		//close fin if open
		if (fin.is_open()) fin.close();
		throw;
	}
	
	const double pixel_area = mio::Optim::pow2( skiRunsMap.cellsize );
	for (std::map<size_t, condition>::iterator it = slope_conditions.begin(); it != slope_conditions.end(); ++it) {
		if (it->second.slope_area==0.) {
			std::cerr << "[W] ski slope " << it->first << " declared in the SLOPE_CONDITIONS file but not found in the SKIRUNS_FILE grid\n";
			slope_conditions.erase( it );
			continue;
		}
		if (it->second.priority==0.)
			throw mio::InvalidArgumentException("Could not find slope " + mio::IOUtils::toString(it->first) + " in the SLOPE_CONDITIONS file", AT);
		
		it->second.slope_area *= pixel_area;
		
		 std::cout << "Found slope " << it->first << " -> " << it->second.toString() << "\n";
	}
}

/**
 * @brief Get the slope number for each pixel.
 * @details Getting the slope number for each pixel. The pixel code is 1xxxx. This function gets rid of the 1 before xxxx.
 * @param[in] dbl_code of the slope map.
 * @return value of the slope pixel defined in the SLOPE_CONDITIONS file.
*/
size_t TechSnowA3D::getSlopeNumber(const double& dbl_code)
{
	static const double epsilon = 0.0001;
	return static_cast<size_t>( (dbl_code - 10000. + epsilon) );
}

/**
 * @brief Setting priority of the technical snow production
 * @details The priority of the technical snow production is set of the defined slope conditions in the [TechSnow] in the io.ini file
 * @param[in] date actual date
 * @param[in] ppt definition of the snow production strategy
 * @param[in] snow_height actual snow height
 * @param[in] slope properties of the slope section
 * @param[in] date_hour time
 * @return psum_technical_snow sum of the technical snow production for each pixel.
*/
double TechSnowA3D::setPriority(const mio::Date& date, 
							 const TechSnowA3D::snowStrategy &ppt, const double& snow_height, 
							 const TechSnowA3D::condition& slope,
							 const int date_hour) const
{
	double psum_technical_snow = IOUtils::nodata;

	// Start snow production of winter season; ski resort still closed
	// Technical snow production the whole day
	// Check min snow height is reached before opening
	if (date >= ppt.startProd && date < start_season && date <= ppt.endProd) {
		if (snow_height < slope.min_height * ppt.startAim) 
			psum_technical_snow = max_snowgun_water * slope.number_snowguns / (slope.slope_area) * dt_main * (ppt.gunOperation);
	} else if (date >= ppt.startProd && date >= start_season && date <= ppt.endProd) {
		// Continue snow production after season opening; ski resort is open
		// Technical snow production only when ski resort is closed
		if (date_hour > slope_closed || date_hour < slope_open) {
			if (snow_height < slope.min_height * ppt.endAim)
				psum_technical_snow = max_snowgun_water * slope.number_snowguns / (slope.slope_area) * dt_main * (ppt.gunOperation);
		}
	}
	
	// Technical snow production parameters [mm]
	return psum_technical_snow;
}

/**
 * @brief Get the grooming and amount of technical snow production map
 * @details Setting psum_tech and grooming for each pixel based on the defined slope conditions
 * @param[in] ta temperature map
 * @param[in] rh relative humidity map
 * @param[in] hs snow height map
 * @param[in] date actual date
*/
void TechSnowA3D::setMeteo(const mio::Grid2DObject& ta, const mio::Grid2DObject& rh, const mio::Grid2DObject& hs, const mio::Date& date)
{
	// Cleaning data from last timestep
	psum_tech = IOUtils::nodata;
	
	// Get actual time (hour and minutes)
	int date_hour, date_min;
	date.getTime(date_hour, date_min);

	//Period of technical snow production?
	if (date < earliest_production || date > end_season) return;
	
	std::cout << "[I] Production of technical snow has started\n";
	for (size_t ii=0; ii<skiRunsMap.size(); ii++) {
		// Check if slope or not -> if not, do nothing
		if (skiRunsMap(ii)==IOUtils::nodata) continue;
		// Check slope and section
		const size_t slope_nr = getSlopeNumber(skiRunsMap(ii));
		// Wet bulb temperature, in Celsius
		const double T_wet = IOUtils::K_TO_C( mio::Atmosphere::wetBulbTemperature(ta(ii), rh(ii), mean_elevation) );
		
		// Threshold wet bulb temperature
		double psum_technical = IOUtils::nodata;				// [mm] Technical snow height in water equivalent
		
		if (T_wet < slope_conditions[slope_nr].wet_bulb_thresh) {
			// Check priority. Since we have cleaned up slope_conditions in setSlopeConditions(), it should always be valid
			const unsigned int priority = slope_conditions[slope_nr].priority - 1; //since vectors start at 0 and not 1!
			psum_technical = setPriority(date, priorities[ priority ], hs(ii), slope_conditions[slope_nr], date_hour);
		}
		
		// Technical snow production
		psum_tech(ii) = psum_technical;
	}
	std::cout << "[I] Production of technical snow has finished\n";
}


mio::Grid2DObject TechSnowA3D::getGrid(const SnGrids::Parameters& param) const
{
	switch (param) {
		case SnGrids::GROOMING:
			return skiRunsMap; //the caller will decide if it is the right time for grooming or not
		case SnGrids::PSUM_TECH:
			return psum_tech;
		default:
			throw mio::InvalidArgumentException("The requested grid can not be provided by TechSnow", AT);
	}
}
