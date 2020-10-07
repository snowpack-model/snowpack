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
#include <snowpack/libsnowpack.h>
#include <alpine3d/MeteoObj.h>

/**
 * @page techsnowA3D Technical Snow
 * @details
 * Implementation of technical snow production and grooming. Slopes are tagged in a provided 2D grid
 * (same geolocalization as the DEM) where each pixel either has its slope number or nodata.
 * Then different snow production priorities (1, 2, 3 etc) are available and can be attributed to
 * the individual ski slopes. The number of priorities is unlimited, their relative order has no
 * special meaning but it is not possible to have gaps in the numbering when defining the properties
 * of the priorities (see below for the configuration keys).Please also keep in mind that during 
 * the ski season, snow is only produced when the slopes are closed to the public.
 * 
 * This module works hand in hand with the Snowpack TechSnow module, so it also uses the configuration keys defined in 
 * <a href="https://models.slf.ch/docserver/snowpack/html/classTechSnow.html">Snowpack::TechSnow</a> to control the grooming.
 *
 * Then, it relies on the following configuration keys, all in the [TechSnow] section:
 *  - SNOW_PRODUCTION: if set to true, enables this module (default: false);
 *  - SKIRUNS_FILE: the grid where all the pixels are tagged either with their slope number or nodata;
 *  - SLOPE_CONDITIONS: CSV file that gives the properties of each slope
 *  - SEASON_OPENING: Ski resort opening dates of the skiing season
 *  - SEASON_CLOSING: Ski resort ending dates of the skiing season
 *  - SLOPE_OPEN: at what local time do the slopes open to the public?
 *  - SLOPE_CLOSED: at what local time do the slopes close to the public?
 *  - MAX_SNOWGUN_WATER: max l/s that a snogun can provide
 * 
 * For each snow production priority, the following keys are defined (unlimited number of priorities 
 * but please, no gaps in the numbering!):
 *  - PRIO\#\::start_prod: snow production starting date
 *  - PRIO\#\::end_prod: date of the end of snow production
 *  - PRIO\#\::start_aim: Snow production aim in snow height [m] until ski resort opening (factor x snow_prod_min)
 *  - PRIO\#\::end_aim: Snow production aim in snow height [m] until ski resort ending (factor x snow_prod_min)
 *  - PRIO\#\::gun_operation: Percentage of snow gun operation [%]
 *
 * The *slope conditions* file is a CSV file (the '#' or ';' characters are used to comment a line or part of it) that
 * contains the following: the slope number, the slope area, the number of snow guns, its priority, the minimum snow
 * height and wetbulb threshold temperature that should be present on the slope. For example:
 * @code
 * #slope_number	nr_snowguns	priority	minimum_hs	wetbulb_threshold
 * #[-] 		[-]		[-]		[m]		[°C]
 * 1		10		1		0.6		-2.0
 * @endcode
 * 
 *@author Mathias Bavay, Pirmin Ebner and others
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
		
		static std::string getGridsRequirements() { return "TA RH HS"; }

	private:
		typedef struct CONDITION {
			CONDITION() 
			   : slope_area(0.), number_snowguns(0), priority(mio::IOUtils::unodata), min_height(0.)  {}
			
			void setUserProperties(const unsigned int& nr_snowguns, const unsigned int& i_priority, const double& i_min_height, const double& wb_thresh) {
				number_snowguns = nr_snowguns;
				priority = i_priority;
				min_height = i_min_height;
				wet_bulb_thresh = wb_thresh;
			}
			
			const std::string toString() {std::ostringstream os; os << "[" << slope_area << "m², " << number_snowguns << " snowguns, " << priority << " priority, " << min_height << "m min, min wbt: " << wet_bulb_thresh << "°C]"; return os.str();}
			
			double slope_area;		//area of the slope section
			unsigned int number_snowguns;	//number of snow guns per slope section
			unsigned int priority;	//snow production priority for this slope
			double min_height;		//minimum snow height for technical snow production
			double wet_bulb_thresh;	//wet bulb temperature threshold (°C) below which snow can be produced
		} condition;
		
		typedef struct SNOWSTRATEGY {
			SNOWSTRATEGY() : startProd(), endProd(), startAim(0.), endAim(0.), gunOperation(0.) {}
			const std::string toString() {std::ostringstream os; os << "[" << startProd.toString(mio::Date::ISO) << "-" << endProd.toString(mio::Date::ISO) << " aim: " << startAim << " -> " << endAim << " @ " << gunOperation << "]"; return os.str();}
			
			mio::Date startProd;	//start date of snow production
			mio::Date endProd;		//end date of snow production
			double startAim;		//snow production aim in snow height [m] until ski resort opening (factor x snow_prod_min)
			double endAim;			//snow production aim in snow height [m] between ski resort opening and end date of snow production (factor x snow_prod_min)
			double gunOperation;	//percentage of snow gun operation
		} snowStrategy;
		
		static TechSnowA3D::snowStrategy setSnowStrategy(const mio::Config& cfg, const double& TZ, const unsigned int& nr);
		void setSlopeConditions(const std::string& filename);
		static size_t getSlopeNumber(const double& dbl_code);
		double setPriority(const mio::Date& date, const TechSnowA3D::snowStrategy &ppt, const double& snow_height, const TechSnowA3D::condition& slope, const int date_hour) const;
		
		mio::Grid2DObject skiRunsMap; ///< All ski runs pixels are tagged by ski run and section number
		mio::Grid2DObject psum_tech;
		
		std::map<size_t, condition> slope_conditions;
		mio::Date start_season, end_season;		///< start and end of ski season (snow is only produced at off hours during ski season)
		mio::Date earliest_production;			///< this is only used to optimize the processing
		std::vector<snowStrategy> priorities;
		
		double max_snowgun_water;	///< max l/s that a snogun can provide
		const double mean_elevation;
		int slope_open, slope_closed; ///< local hour for slope opening / closing
};

#endif
