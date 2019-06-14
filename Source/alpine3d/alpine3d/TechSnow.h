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
#ifndef TECHSNOWA3D_H
#define TECHSNOWA3D_H

class SnowpackInterfaceWorker;

#include <meteoio/MeteoIO.h>
#include <alpine3d/MeteoObj.h>

/**
 * @page techsnow Technical Snow
 */
class TechSnowA3D
{
	public:
		TechSnowA3D(const mio::Config& cfg, const mio::DEMObject& dem);

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
		
		static std::vector<TechSnowA3D::condition> readSlopeConditions(const int& column, const std::string& filename);
		static short int getSlopeNumber(const double& dbl_code);
		static short int findSlope(const int& numbers_of_slopes, const std::vector<condition>& slope_condition, const int& findSlope);
		double setPriority( const mio::Date& date, const mio::Date& start_prod, const mio::Date& end_prod, const double& start_aim, const double& end_aim, const double& gun_operation, const double& snow_height, const double& slope_area, const double& nr_snowguns, const double& min_height, const int date_hour) const;
		
		mio::Grid2DObject skiRunsMap; ///< All ski runs pixels are tagged by ski run and section number
		mio::Grid2DObject grooming, psum_tech;
		
		std::vector<condition> slope_condition;
		mio::Date start_season, end_season;		///< [-] latest production date
		mio::Date start_prod1, end_prod1;
		mio::Date start_prod2, end_prod2;
		mio::Date start_prod3, end_prod3;
		double start_aim1, end_aim1, gun_operation1;
		double start_aim2, end_aim2, gun_operation2;
		double start_aim3, end_aim3, gun_operation3;
		double wet_bulb_thresh;		//< [C] wet bulb temperature
		double max_snowgun_water;	///< max l/s that a snogun can provide
		const double mean_elevation;
		int slope_open, slope_closed; ///< local hour for slope opening / closing
		int number_of_slopes;
};

#endif
