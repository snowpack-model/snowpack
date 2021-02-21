/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/TerrainRadiationSimple.h>
#include <alpine3d/ebalance/TerrainRadiationHelbig.h>
#include <alpine3d/ebalance/TerrainRadiationComplex.h>

using namespace std;
using namespace mio;

TerrainRadiationAlgorithm::~TerrainRadiationAlgorithm() {}
// FELIX: const RadiationField* radfield
//please document these keywords in EnergyBalance.h
TerrainRadiationAlgorithm* TerrainRadiationFactory::getAlgorithm(const Config& cfg, const DEMObject &dem,
                                                                 const int& nbworkers)
{
	string method = "SIMPLE";
	cfg.getValue("Terrain_Radiation_Method", "EBalance", method, IOUtils::nothrow);
	IOUtils::toUpper(method);

	if (method == "SIMPLE") {
		return new TerrainRadiationSimple(cfg, dem, method);
	} else if (method == "HELBIG") {
		return new TerrainRadiationHelbig(cfg, dem, nbworkers, method);
	} else if (method == "COMPLEX") {
		return new TerrainRadiationComplex(cfg, dem, method);
	} else {
		throw IOException("The terrain radiation method '"+method+"' is not implemented/activated" , AT);
	}
}
