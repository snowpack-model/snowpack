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
#ifndef TERRAINRADIATIONALGORITHM_H
#define TERRAINRADIATIONALGORITHM_H

#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/RadiationField.h>
#include <alpine3d/ebalance/SolarPanel.h>

//Class for shooting cells
class CellsList {
	public:
		double radiation;
		int x;
		int y;
};

inline bool operator_greater(const CellsList& a, const CellsList& b) {
	return a.radiation > b.radiation;
}

class TerrainRadiationAlgorithm {
	public:
		TerrainRadiationAlgorithm(const std::string& i_algo) : algo(i_algo), _hasSP(false) {}
		virtual ~TerrainRadiationAlgorithm();
		const bool hasSP(){return _hasSP;}
		virtual void getRadiation(mio::Array2D<double>& direct, mio::Array2D<double>& diffuse,
                              mio::Array2D<double>& terrain, const mio::Array2D<double>& direct_unshaded_horizontal,
                              const mio::Array2D<double>& total_ilwr, mio::Array2D<double>& sky_ilwr,
                              mio::Array2D<double>& terrain_ilwr, double solarAzimuth, double solarElevation) = 0;
		virtual void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta) = 0;

		const std::string algo;
		virtual void setSP(const mio::Date timestamp, const double solarAzimuth, const double solarElevation){};
		virtual void writeSP(const unsigned int max_steps){};
		virtual void getSkyViewFactor(mio::Array2D<double> &o_sky_vf) = 0;

	protected:
		bool _hasSP;
};

class TerrainRadiationFactory {
	public:
		// FELIX: const RadiationField* radfield
		static TerrainRadiationAlgorithm* getAlgorithm(const mio::Config& cfg, const mio::DEMObject &dem, const int& nbworkers);
};

#endif
