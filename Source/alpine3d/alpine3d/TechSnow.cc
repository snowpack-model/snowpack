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

using namespace std;
using namespace mio;

TechSnow::TechSnow(const mio::Config& cfg, const mio::DEMObject& dem) 
              : skiRunsMap(), grooming(), psum_tech(), isMaster( MPIControl::instance().master() )
{
	const std::string inpath = cfg.get("GRID2DPATH", "Input");
	const std::string mapFile = cfg.get("SKIRUNS", "Input");
	
	if (isMaster) {
		mio::IOManager io( cfg );
		io.read2DGrid(skiRunsMap, inpath + "/" + mapFile);
		
		if (!skiRunsMap.isSameGeolocalization( dem )) 
			throw InvalidArgumentException("The ski runs map does not has the same geolocalization as the DEM", AT);
	}
	
	MPIControl::instance().broadcast( skiRunsMap );
}

short int TechSnow::getRunNumber(const double& dbl_code)
{
	static const double epsilon = 0.0001;
	return static_cast<short int>( (dbl_code + epsilon) / 100. );
}

void TechSnow::setMeteo(const mio::Grid2DObject& ta, const mio::Grid2DObject& rh, const mio::Date& timestamp)
{
	grooming.clear();
	psum_tech.clear();
	
	const unsigned short weekNr = timestamp.getISOWeekNr();
	
	//according to ta and rh, fill grooming and psum_tech grids
	for (size_t ii=0; ii<grooming.size(); ii++) {
		grooming(ii) = IOUtils::nodata;
		psum_tech(ii) = IOUtils::nodata;
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

