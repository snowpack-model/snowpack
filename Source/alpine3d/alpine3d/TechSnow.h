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
		              const mio::Date& timestamp);
		
		mio::Grid2DObject getGrid(const SnGrids::Parameters& param) const;
		
		static std::string getGridsRequirements() { return "TA RH"; }

	private:
		static short int getRunNumber(const double& dbl_code);
		
		mio::Grid2DObject skiRunsMap; ///< All ski runs pixels are tagged by ski run and section number
		mio::Grid2DObject grooming, psum_tech;
		
		const bool isMaster;
};

#endif
