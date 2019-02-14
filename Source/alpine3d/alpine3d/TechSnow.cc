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
#include <alpine3d/TechSnow.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/AlpineMain.h>

#include <errno.h>

#include <algorithm> // for std::find
#include <iterator>

using namespace std;
using namespace mio;

TechSnow::TechSnow(const mio::Config& cfg, const mio::DEMObject& dem) 
              : skiRunsMap(), grooming(dem, IOUtils::nodata), psum_tech(dem, IOUtils::nodata), isMaster( MPIControl::instance().master() )
{
	const std::string mapFile = cfg.get("SKIRUNS", "Input");
	if (isMaster) {
		mio::IOManager io( cfg );
		if (istringstream(mapFile)) {	// Check IO.ini if path to a SKIRUNS-File is given
			cout << "Read SKIRUNS" << endl;
			io.read2DGrid(skiRunsMap, "/" + mapFile);
			if (!skiRunsMap.isSameGeolocalization( dem )) 
				throw InvalidArgumentException("The ski runs map does not has the same geolocalization as the DEM", AT);
		}
	}
		
	MPIControl::instance().broadcast( skiRunsMap );
}

void TechSnow::readSlopeConditions(const int& numbers_of_slopes, std::vector<condition>& slope_condition, const std::string& filename)
{	
	int jj=0;
	
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	
	const char eoln = mio::FileUtils::getEoln(fin); //get the end of line character for the file

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

			iss >> std::skipws >> slope_condition[jj].slope_number;
			iss >> std::skipws >> slope_condition[jj].slope_area;
			iss >> std::skipws >> slope_condition[jj].number_snowguns;
			iss >> std::skipws >> slope_condition[jj].priority;
			iss >> std::skipws >> slope_condition[jj].min_height;
		
			jj++;
		} while (!fin.eof());
		// Check wether total number of slopes are correct
		if (jj == numbers_of_slopes) fin.close();
		else {
			std::cout << "Error: Please check if the total number of slopes are correct";
			std::terminate();
		}
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}


short int TechSnow::getSlopeNumber(const double& dbl_code)
{
	static const double epsilon = 0.0001;
	return static_cast<short int>( (dbl_code - 10000. + epsilon) );
}


short int TechSnow::findSlope(const int& numbers_of_slopes, const std::vector<condition>& slope_condition, 
                              const int& findSlope)
{
	static const double epsilon = 0.001;
	int index = -1;
	for (int i=0; i<=numbers_of_slopes; i++)
	{
		double slopeFile = slope_condition[i].slope_number*100. + epsilon;		
		if ((int)slopeFile == findSlope)
		{
			index = i;
			break;
		}
	}
	if (index == -1)
	{
		std::cout << "Slope and section number could not be found. Please check slope file." << endl;
		std::cout << "Slope number to find:   " << findSlope << endl;
		std::terminate();
	}
	return index;
}	


double TechSnow::setPriority(const std::string& date, const std::string& season_opening, 
							 const std::string& start_prod, 
							 const std::string& end_prod, const double& start_aim, const double& end_aim,
							 const double& gun_operation, const double& snow_height, const double& V_water, 
							 const double& slope_area, const double& nr_snowguns, const double& min_height,
							 const int date_hour, const int slope_open, const int slope_closed)
{
	double psum_technical_snow = IOUtils::nodata;

	// Start snow production of winter season; ski resort still closed
	// Technical snow production the whole day
	// Check min snow height is reached before opening
	if (date >= start_prod && date < season_opening && date <= end_prod) 
	{
		if (snow_height >= min_height * start_aim) psum_technical_snow = IOUtils::nodata;
		else psum_technical_snow = V_water * nr_snowguns / (slope_area) * dt_main * (gun_operation);
	}
	// Continue snow production after season opening; ski resort is open
	// Technical snow production only when ski resort is closed
	else if (date >= start_prod && date >= season_opening && date <= end_prod) 
	{	
		if (date_hour <= slope_closed && date_hour >= slope_open)
		{
			if (snow_height >= min_height * end_aim) psum_technical_snow = IOUtils::nodata;
			else psum_technical_snow = V_water * nr_snowguns / (slope_area) * dt_main * (gun_operation);

		}
		else psum_technical_snow = IOUtils::nodata; // No production during opening hours of the ski resort
	}
	
	// Technical snow production parameters [mm]
	return psum_technical_snow;
}


void TechSnow::setMeteo(const mio::Grid2DObject& ta, const mio::Grid2DObject& rh, const mio::Grid2DObject& hs, const mio::Date& timestamp)
{
	// Load input-file for technical snow production
	Config cfg("./io_technical_snow_production.ini");

	// Cleaning data from last timestep
	grooming = IOUtils::nodata;
	psum_tech = IOUtils::nodata;

	// Initiate parameters
	int slope_nr;
	int slope_index;
	double psum_technical = IOUtils::nodata;				// [mm] Technical snow height in water equivalent
	static double T_wet;			// [C] Wet bulb temperature

	const double wet_bulb_thresh = cfg.get("WET_BULB_THRESH", "Slope_Priority");		// [C] wet bulb temperature
	const std::string start_prod = cfg.get("PRIO1::start_prod", "Slope_Priority");	// [-] earliest snow production date
	const std::string end_season = cfg.get("SEASON_CLOSING", "Slope_Priority");	// [-] latest production date
	const std::string date = timestamp.toString(Date::ISO);		// [-] actual simulation date
	int slope_open = cfg.get("SLOPE_OPEN", "Slope_Priority");
	int slope_closed = cfg.get("SLOPE_CLOSED", "Slope_Priority");
	const int number_of_slopes = cfg.get("NUMBER_OF_SLOPES", "Input");
	//condition slope_condition[number_of_slopes];
	std::vector<condition> slope_condition(number_of_slopes);
	
	// Get actual time (hour and minutes)
	int date_hour, date_min;
	timestamp.getTime(date_hour, date_min);

	// Read slope conditions
	readSlopeConditions(number_of_slopes, slope_condition, cfg.get("SLOPE_CONDITIONS", "Input"));
	
	//Period of technical snow production and grooming
	if (date >= start_prod && date <= end_season){
		std::cout << "Production of technical snow has started" << endl;
		for (size_t ii=0; ii<skiRunsMap.size(); ii++) {
			// Check if slope or not -> if not, do nothing
			if (skiRunsMap(ii)==IOUtils::nodata) continue;
			// Check slope and section
			slope_nr = getSlopeNumber(skiRunsMap(ii));
			// Get slope section properties
			slope_index = findSlope(number_of_slopes, slope_condition, slope_nr);
			// Wet bulb temperature				
			T_wet = (IOUtils::K_TO_C(ta(ii))) * atan(0.151977 * pow(rh(ii)*100. + 8.313659, 0.5)) + atan(IOUtils::K_TO_C(ta(ii)) + rh(ii)*100.) - atan(rh(ii)*100. - 1.676331) + 0.00391838 * pow(rh(ii)*100, 1.5) * atan(0.023101 * rh(ii)*100) - 4.686035;
			// Threshold wet bulb temperature
			if (T_wet >= wet_bulb_thresh) psum_technical = IOUtils::nodata;
			// Check priority
			else
			{
				if (slope_condition[slope_index].priority == 1){
					psum_technical = setPriority(date, cfg.get("SEASON_OPENING", "Slope_Priority"), 
								  cfg.get("PRIO1::start_prod", "Slope_Priority"), cfg.get("PRIO1::end_prod", "Slope_Priority"), 
								  cfg.get("PRIO1::start_aim", "Slope_Priority"), cfg.get("PRIO1::end_aim", "Slope_Priority"),
						    	  cfg.get("PRIO1::gun_operation", "Slope_Priority"), hs(ii),
						    	  cfg.get("MAX_SNOWGUN_WATER", "Slope_Priority"), 
						    	  slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
						    	  slope_condition[slope_index].min_height,
						    	  date_hour, slope_open, slope_closed);
				}
				else if (slope_condition[slope_index].priority == 2){
					psum_technical = setPriority(date, cfg.get("SEASON_OPENING", "Slope_Priority"), 
								  cfg.get("PRIO2::start_prod", "Slope_Priority"), cfg.get("PRIO2::end_prod", "Slope_Priority"), 
								  cfg.get("PRIO2::start_aim", "Slope_Priority"), cfg.get("PRIO2::end_aim", "Slope_Priority"),
						    	  cfg.get("PRIO2::gun_operation", "Slope_Priority"), hs(ii),
						    	  cfg.get("MAX_SNOWGUN_WATER", "Slope_Priority"), 
						    	  slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
						    	  slope_condition[slope_index].min_height,
						    	  date_hour, slope_open, slope_closed);

				}
				else if (slope_condition[slope_index].priority == 3){
					psum_technical = setPriority(date, cfg.get("SEASON_OPENING", "Slope_Priority"), 
								  cfg.get("PRIO3::start_prod", "Slope_Priority"), cfg.get("PRIO3::end_prod", "Slope_Priority"), 
								  cfg.get("PRIO3::start_aim", "Slope_Priority"), cfg.get("PRIO3::end_aim", "Slope_Priority"),
						    	  cfg.get("PRIO3::gun_operation", "Slope_Priority"), hs(ii),
						    	  cfg.get("MAX_SNOWGUN_WATER", "Slope_Priority"), 
						    	  slope_condition[slope_index].slope_area, slope_condition[slope_index].number_snowguns, 
						    	  slope_condition[slope_index].min_height,
						    	  date_hour, slope_open, slope_closed);

				}
				else psum_technical = IOUtils::nodata;
			}
			// Grooming and technical snow production
			if (hs(ii) >= 0.4) grooming(ii) = 1;
			else grooming(ii) = IOUtils::nodata;
			psum_tech(ii) = psum_technical;
		}
		std::cout << "Production of technical snow has finished" << endl;
	}
}


mio::Grid2DObject TechSnow::getGrid(const SnGrids::Parameters& param) const
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

