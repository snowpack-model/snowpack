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
#ifndef TECHSNOW_H
#define TECHSNOW_H

class SnowpackInterfaceWorker;

#include <meteoio/MeteoIO.h>
#include <alpine3d/MeteoObj.h>

/**
 * @page techsnow Technical Snow
 */
class TechSnow
{
	public:
		TechSnow(const mio::Config& cfg, const mio::DEMObject& dem);

		void setMeteo(const mio::Grid2DObject& ta,
		              const mio::Grid2DObject& rh,
		              const mio::Grid2DObject& hs,
		              const mio::Date& timestamp);
		
		mio::Grid2DObject getGrid(const SnGrids::Parameters& param) const;
		
		static std::string getGridsRequirements() { return "TA RH"; }

	private:
		typedef struct CONDITION {
			double slope_number;
			double slope_area;
			double number_snowguns;
			int priority;
			double min_height;
		} condition;
		
		static void readSlopeConditions(const int& column, std::vector<condition>& slope_condition, const std::string& filename);
		static short int getSlopeNumber(const double& dbl_code);
		static short int findSlope(const int& numbers_of_slopes, const std::vector<condition>& slope_condition, const int& findSlope);
		static double setPriority(const std::string& date, const std::string& season_opening, const std::string& start_prod, const std::string& end_prod, const double& start_aim, const double& end_aim, const double& gun_operation, const double& snow_height, const double& V_water, const double& slope_area, const double& nr_snowguns, const double& min_height, const int date_hour, const int slope_open, const int slope_closed);
							 	  
		mio::Grid2DObject skiRunsMap; ///< All ski runs pixels are tagged by ski run and section number
		mio::Grid2DObject grooming, psum_tech;
		
		const bool isMaster;
};

#endif
