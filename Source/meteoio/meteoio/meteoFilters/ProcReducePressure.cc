// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcReducePressure.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <cmath>
#include <stdio.h>

using namespace std;

namespace mio {

ProcReducePressure::ProcReducePressure(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
                       : ProcessingBlock(vecArgs, name, cfg)
{
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcReducePressure::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;

	for (size_t ii=0; ii<ovec.size(); ii++) {
		double& tmp = ovec[ii](param);
		if (tmp==IOUtils::nodata) continue;
		
		const double lat = ovec[ii].meta.position.getLat();
		const double alt = ovec[ii].meta.position.getAltitude();
		
		if (lat!=IOUtils::nodata && alt!=IOUtils::nodata) {
			tmp = Atmosphere::reducedAirPressure(tmp, lat, alt);
		} else {
			tmp = IOUtils::nodata;
		}
	}
}

} //end namespace
