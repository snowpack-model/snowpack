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
#ifndef RADIATIONFIELD_H
#define RADIATIONFIELD_H

#include <meteoio/MeteoIO.h>

class RadiationField {
	public:
		RadiationField();
		RadiationField(const mio::DEMObject& in_dem, const size_t& in_startx, const size_t& in_nx);

		void setDEM(const mio::DEMObject& in_dem);
		void setDEM(const mio::DEMObject& in_dem, const size_t& in_startx, const size_t& in_nx);
		void setStations(const std::vector<mio::MeteoData>& vecMeteo, const mio::Grid2DObject& albedo);

		void setMeteo(const mio::Grid2DObject& in_ta, const mio::Grid2DObject& in_rh, const mio::Grid2DObject& in_p, const mio::Grid2DObject& in_albedo);

		void getPositionSun(double& o_solarAzimuth, double& o_solarElevation) const;
		void getRadiation(mio::Array2D<double>& o_direct, mio::Array2D<double>& o_diffuse) const;
		void getBandOffsets(size_t& o_startx, size_t& o_stopx) const;

	private:
		mio::Date date;
		mio::DEMObject dem;
		mio::Grid2DObject dem_band, direct, diffuse;
		mio::SunObject Sun;
		std::vector<mio::StationData> vecMeta;
		std::vector<double> vecMd, vecCorr;
		double dem_mean_altitude;
		double cellsize;
		size_t dem_dimx, band_dimx, dimy;
		size_t startx;
		bool day, night; ///>is it day or night? It can be neither! (ie: at dawn or dusk)
};

#endif
