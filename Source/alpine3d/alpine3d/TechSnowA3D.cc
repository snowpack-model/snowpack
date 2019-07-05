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

TechSnowA3D::TechSnowA3D(const mio::Config& cfg, const mio::DEMObject& dem) 
              : skiRunsMap(), grooming(dem, IOUtils::nodata), psum_tech(dem, IOUtils::nodata), slope_condition(),
              start_season(), end_season(), start_prod1(), end_prod1(), start_prod2(), end_prod2(), start_prod3(), end_prod3(),
              start_aim1(), end_aim1(), gun_operation1(), start_aim2(), end_aim2(), gun_operation2(), start_aim3(), end_aim3(), gun_operation3(),
              wet_bulb_thresh(), max_snowgun_water(), mean_elevation(dem.grid2D.getMean()), slope_open(), slope_closed(), number_of_slopes()
{
	const std::string mapFile = cfg.get("SKIRUNS", "Input");
	const double TZ = cfg.get("TIME_ZONE", "Input");
	const bool isMaster = MPIControl::instance().master();
	if (isMaster) {
		mio::IOManager io( cfg );
		io.read2DGrid(skiRunsMap, "/" + mapFile);
		if (!skiRunsMap.isSameGeolocalization( dem )) 
			throw InvalidArgumentException("The ski runs map does not has the same geolocalization as the DEM", AT);
	}
	MPIControl::instance().broadcast( skiRunsMap );
	
	const mio::Config tech_cfg("./io_technical_snow_production.ini");
	
	// Read slope conditions
	tech_cfg.getValue("NUMBER_OF_SLOPES", "Input", number_of_slopes);
	slope_condition = readSlopeConditions(number_of_slopes, tech_cfg.get("SLOPE_CONDITIONS", "Input"));
	
	tech_cfg.getValue("WET_BULB_THRESH", "Slope_Priority", wet_bulb_thresh);
	tech_cfg.getValue("SEASON_OPENING", "Slope_Priority", start_season, TZ);
	tech_cfg.getValue("SEASON_CLOSING", "Slope_Priority", end_season, TZ);
	tech_cfg.getValue("SLOPE_OPEN", "Slope_Priority", slope_open);
	tech_cfg.getValue("SLOPE_CLOSED", "Slope_Priority", slope_closed);
	tech_cfg.getValue("MAX_SNOWGUN_WATER", "Slope_Priority", max_snowgun_water);

	tech_cfg.getValue("PRIO1::start_prod", "Slope_Priority", start_prod1, TZ);
	tech_cfg.getValue("PRIO1::end_prod", "Slope_Priority", end_prod1, TZ);
	tech_cfg.getValue("PRIO1::start_aim", "Slope_Priority", start_aim1);
	tech_cfg.getValue("PRIO1::end_aim", "Slope_Priority", end_aim1);
	tech_cfg.getValue("PRIO1::gun_operation", "Slope_Priority", gun_operation1);
	
	tech_cfg.getValue("PRIO2::start_prod", "Slope_Priority", start_prod2, TZ);
	tech_cfg.getValue("PRIO2::end_prod", "Slope_Priority", end_prod2, TZ);
	tech_cfg.getValue("PRIO2::start_aim", "Slope_Priority", start_aim2);
	tech_cfg.getValue("PRIO2::end_aim", "Slope_Priority", end_aim2);
	tech_cfg.getValue("PRIO2::gun_operation", "Slope_Priority", gun_operation2);
	
	tech_cfg.getValue("PRIO3::start_prod", "Slope_Priority", start_prod3, TZ);
	tech_cfg.getValue("PRIO3::end_prod", "Slope_Priority", end_prod3, TZ);
	tech_cfg.getValue("PRIO3::start_aim", "Slope_Priority", start_aim3);
	tech_cfg.getValue("PRIO3::end_aim", "Slope_Priority", end_aim3);
	tech_cfg.getValue("PRIO3::gun_operation", "Slope_Priority", gun_operation3);
}

std::vector<TechSnowA3D::condition> TechSnowA3D::readSlopeConditions(const int& numbers_of_slopes, const std::string& filename)
{
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	
	const char eoln = mio::FileUtils::getEoln(fin); //get the end of line character for the file
	std::vector<condition> slope_conditions;
	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			mio::IOUtils::stripComments(line);
			mio::IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);

			condition tmp_cond;
			iss >> std::skipws >> tmp_cond.slope_number;
			iss >> std::skipws >> tmp_cond.slope_area;
			iss >> std::skipws >> tmp_cond.number_snowguns;
			iss >> std::skipws >> tmp_cond.priority;
			iss >> std::skipws >> tmp_cond.min_height;
			
			slope_conditions.push_back( tmp_cond );
		} while (!fin.eof());
		// Check wether total number of slopes are correct
		if ((signed)slope_conditions.size() == numbers_of_slopes) {
			fin.close();
		} else {
			throw mio::InvalidArgumentException("Inconsistent number of slopes", AT);
		}
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
	
	return slope_conditions;
}


short int TechSnowA3D::getSlopeNumber(const double& dbl_code)
{
	static const double epsilon = 0.0001;
	return static_cast<short int>( (dbl_code - 10000. + epsilon) );
}


short int TechSnowA3D::findSlope(const int& numbers_of_slopes, const std::vector<condition>& slope_condition, 
                              const int& findSlope)
{
	static const double epsilon = 0.001;
	int index = -1;
	for (int i=0; i<numbers_of_slopes; i++) {
		const double slopeFile = slope_condition[i].slope_number*100. + epsilon;		
		if ((int)slopeFile == findSlope) {
			index = i;
			break;
		}
	}
	if (index == -1)
		throw mio::InvalidArgumentException("Could not find slope "+mio::IOUtils::toString(findSlope), AT);
	
	return index;
}	


double TechSnowA3D::setPriority(const mio::Date& date, 
							 const mio::Date& start_prod, 
							 const mio::Date& end_prod, const double& start_aim, const double& end_aim,
							 const double& gun_operation, const double& snow_height, 
							 const double& slope_area, const double& nr_snowguns, const double& min_height,
							 const int date_hour) const
{
	double psum_technical_snow = IOUtils::nodata;

	// Start snow production of winter season; ski resort still closed
	// Technical snow production the whole day
	// Check min snow height is reached before opening
	if (date >= start_prod && date < start_season && date <= end_prod) {
		if (snow_height >= min_height * start_aim) psum_technical_snow = IOUtils::nodata;
		else psum_technical_snow = max_snowgun_water * nr_snowguns / (slope_area) * dt_main * (gun_operation);
	} else if (date >= start_prod && date >= start_season && date <= end_prod) {
		// Continue snow production after season opening; ski resort is open
		// Technical snow production only when ski resort is closed
		if (date_hour <= slope_closed && date_hour >= slope_open) {
			if (snow_height >= min_height * end_aim) psum_technical_snow = IOUtils::nodata;
			else psum_technical_snow = max_snowgun_water * nr_snowguns / (slope_area) * dt_main * (gun_operation);
		} else {
			psum_technical_snow = IOUtils::nodata; // No production during opening hours of the ski resort
		}
	}
	
	// Technical snow production parameters [mm]
	return psum_technical_snow;
}


void TechSnowA3D::setMeteo(const mio::Grid2DObject& ta, const mio::Grid2DObject& rh, const mio::Grid2DObject& hs, const mio::Date& date)
{
	// Cleaning data from last timestep
	grooming = IOUtils::nodata;
	psum_tech = IOUtils::nodata;
	
	// Get actual time (hour and minutes)
	int date_hour, date_min;
	date.getTime(date_hour, date_min);

	//Period of technical snow production and grooming?
	if (date < start_prod1 || date > end_season) return;
	
	std::cout << "[I] Production of technical snow has started\n";
	for (size_t ii=0; ii<skiRunsMap.size(); ii++) {
		// Check if slope or not -> if not, do nothing
		if (skiRunsMap(ii)==IOUtils::nodata) continue;
		// Check slope and section
		const int slope_nr = getSlopeNumber(skiRunsMap(ii));
		// Get slope section properties
		const int slope_index = findSlope(number_of_slopes, slope_condition, slope_nr);
		// Wet bulb temperature, in Celsius
		const double T_wet = IOUtils::K_TO_C( mio::Atmosphere::wetBulbTemperature(ta(ii), rh(ii), mean_elevation) );
		
		// Threshold wet bulb temperature
		double psum_technical = IOUtils::nodata;				// [mm] Technical snow height in water equivalent
		
		if (T_wet < wet_bulb_thresh) {
			// Check priority
			if (slope_condition[slope_index].priority == 1) {
				psum_technical = setPriority(date, 
								start_prod1, end_prod1, 
								start_aim1, end_aim1,
								gun_operation1, hs(ii),
								slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
								slope_condition[slope_index].min_height,
								date_hour);
			} else if (slope_condition[slope_index].priority == 2) {
				psum_technical = setPriority(date, 
								start_prod2, end_prod2, 
								start_aim2, end_aim2,
								gun_operation2, hs(ii),
								slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
								slope_condition[slope_index].min_height,
								date_hour);

			} else if (slope_condition[slope_index].priority == 3) {
				psum_technical = setPriority(date, 
								start_prod3, end_prod3, 
								start_aim3, end_aim3,
								gun_operation3, hs(ii),
								slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
								slope_condition[slope_index].min_height,
								date_hour);

			}
		}
		
		// Grooming and technical snow production
		if (hs(ii) >= 0.4) grooming(ii) = 1;
		else grooming(ii) = IOUtils::nodata;
		psum_tech(ii) = psum_technical;
	}
	std::cout << "[I] Production of technical snow has finished\n";
}


mio::Grid2DObject TechSnowA3D::getGrid(const SnGrids::Parameters& param) const
{
	switch (param) {
		case SnGrids::GROOMING:
			return grooming;
		case SnGrids::PSUM_TECH:
			return psum_tech;
		default:
			throw mio::InvalidArgumentException("The requested grid can not be provided by TechSnow", AT);
	}
}
