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

#include <iostream>

class SnGrids {
	public:
		/// \anchor SnGrids this enum provides names for possible Snowpack grids
		enum Parameters {firstparam=0,
                                          TA=firstparam, ///< Air temperature
                                          RH, ///< Relative humidity
                                          VW, ///< Wind velocity
                                          VW_DRIFT, ///< Wind velocity
                                          DW, ///< Wind direction (deg)
                                          ISWR, ///< Incoming short wave radiation
                                          ISWR_DIFF, ///< Incoming short wave, diffuse
                                          ISWR_DIR, ///< Incoming short wave, direct
                                          ILWR, ///< Incoming long wave radiation
                                          HS, ///< Height of snow
                                          PSUM, ///< Water equivalent of precipitations, either solid or liquid
                                          PSUM_PH, ///<  Precipitation phase, between 0 (fully solid) and 1 (fully liquid)
                                          PSUM_TECH, ///< Water equivalent precipitation from artificial snow production
                                          GROOMING, ///< Used as a boolean flag to decide whever a pixel is scheduled for grooming or not
                                          TSG, ///< Temperature ground surface
                                          TSS, ///< Temperature snow surface
                                          TS0, ///< Temperature soil surface
                                          TSNOW, ///< Snow temperature at depth xxx m
                                          TSNOW_AVG, ///< Average snow temperature in the top xxx m
                                          RHOSNOW_AVG, ///< Average snow density in the top xxx m
                                          SWE, ///< Snow Water Equivalent
                                          RSNO, ///< Snow mean density
                                          TOP_ALB, ///< Albedo from the top (ie above canopy)
                                          SURF_ALB, ///< Albedo of the surface (ie below canopy)
                                          SP, ///< sphericity
                                          RB, ///< bond radius
                                          RG, ///< grain radius
                                          N3, ///< grain Coordination number
                                          MS_SNOWPACK_RUNOFF, ///< runoff on the surface of the soil (vitual lysimeter)
                                          MS_SOIL_RUNOFF, ///< runoff at the bottom of the snow/soil column	
                                          MS_WATER, ///< The total amount of water in the snowpack at the present time
                                          SFC_SUBL, ///< The mass loss or gain of the top element due to snow (ice) sublimating
                                          STORE, ///< internal usage (precipitation events that are delayed because they are too small)
                                          ERODEDMASS, ///< wind eroded mass (kg/m2)
                                          WINDEROSIONDEPOSITION, ///< wind erosion and deposition (kg/m2)
                                          GLACIER, ///< mask showing the glaciated pixels
                                          GLACIER_EXPOSED, ///< mask showing the exposed glaciated pixels (ie not snow covered)
                                          TSOIL1, TSOIL2, TSOIL3, TSOIL4, TSOIL5, ///< Temperature within the soil, at a given depth
                                          lastparam=TSOIL5};


		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const size_t& parindex);
		static size_t getParameterIndex(const std::string& parname);

	private:
		static std::vector<std::string> paramname;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

class MeteoObj
{
	public:
		MeteoObj(const mio::Config& config, const mio::DEMObject& in_dem);
		~MeteoObj();

		void setSkipWind(const bool& i_skipWind);
		void prepare(const mio::Date& in_date);
		void get(const mio::Date& in_date,
		         mio::Grid2DObject& ta,
		         mio::Grid2DObject& tsg,
		         mio::Grid2DObject& rh,
		         mio::Grid2DObject& psum,
		         mio::Grid2DObject& psum_ph,
		         mio::Grid2DObject& vw,
		         mio::Grid2DObject& vw_drift,
		         mio::Grid2DObject& dw,
		         mio::Grid2DObject& p,
		         mio::Grid2DObject& ilwr);
		void get(const mio::Date& in_date, std::vector<mio::MeteoData>& o_vecMeteo);
		void checkMeteoForcing(const mio::Date& calcDate);
		void setGlacierMask(const mio::Grid2DObject& glacierMask);
		double getTiming() const;

	private:
		static void checkLapseRate(const std::vector<mio::MeteoData>& i_vecMeteo, const mio::MeteoData::Parameters& param);
		static void checkGridRange(const mio::Date& calcDate, const mio::Grid2DObject& grid, const mio::MeteoData::Parameters& param);
		static void checkInputsRequirements(std::vector<mio::MeteoData>& vecData, const bool& soil_flux);
		void fillMeteoGrids(const mio::Date& calcDate);
		void getMeteo(const mio::Date& calcDate);

		mio::Timer timer;
		const mio::Config &config;
		mio::IOManager io;
		const mio::DEMObject &dem;
		mio::Grid2DObject ta, tsg, rh, psum, psum_ph, vw, vw_drift, dw, p, ilwr;
		mio::Grid2DObject sum_ta, sum_rh, sum_rh_psum, sum_psum, sum_psum_ph, sum_vw, sum_ilwr;
		std::vector<mio::MeteoData> vecMeteo;
		mio::Date date;
		Glaciers *glaciers;
		unsigned int count_sums, count_precip;
		bool skipWind; ///<should the grids be filled or only the data vectors returned?
		bool soil_flux;
		bool enable_simple_snow_drift;
};

#endif
