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
#ifndef METEOOBJ_H
#define METEOOBJ_H

#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

#include <alpine3d/MPIControl.h>
#include <alpine3d/Glaciers.h>

class MeteoObj
{
	public:
		MeteoObj(const mio::Config& config, const mio::DEMObject& in_dem);
		~MeteoObj();
		
		void setSkipWind(const bool& i_skipWind);
		void prepare(const mio::Date& in_date);
		void get(const mio::Date& in_date,
		         mio::Grid2DObject& ta,
		         mio::Grid2DObject& rh,
		         mio::Grid2DObject& psum,
		         mio::Grid2DObject& psum_ph,
		         mio::Grid2DObject& vw,
		         mio::Grid2DObject& p,
		         mio::Grid2DObject& ilwr);
		void get(const mio::Date& in_date, std::vector<mio::MeteoData>& o_vecMeteo);
		void checkMeteoForcing(const mio::Date& calcDate);
		void setGlacierMask(const mio::Grid2DObject& glacierMask);
		double getTiming() const;

	private:
		static void checkLapseRate(const std::vector<mio::MeteoData>& i_vecMeteo, const mio::MeteoData::Parameters& param);
		static void checkGridRange(const mio::Date& calcDate, const mio::Grid2DObject& grid, const mio::MeteoData::Parameters& param);
		static void checkInputsRequirements(std::vector<mio::MeteoData>& vecData);
		void fillMeteoGrids(const mio::Date& calcDate);
		void getMeteo(const mio::Date& calcDate);

		mio::Timer timer;
		const mio::Config &config;
		mio::IOManager io;
		const mio::DEMObject &dem;
		mio::Grid2DObject ta, rh, psum, psum_ph, vw, p, ilwr;
		mio::Grid2DObject sum_ta, sum_rh, sum_rh_psum, sum_psum, sum_psum_ph, sum_vw, sum_ilwr;
		std::vector<mio::MeteoData> vecMeteo;
		mio::Date date;
		Glaciers *glaciers;
		unsigned int count_sums, count_precip;
		bool skipWind; ///<should the grids be filled or only the data vectors returned?
};

#endif
