/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef WINSTRAL_LISTON_DRIFT_ALGORITHM_H
#define WINSTRAL_LISTON_DRIFT_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

class WinstralListonDriftAlgorithm : public InterpolationAlgorithm {
	public:
		WinstralListonDriftAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
		                               GridsManager& i_gdm, Meteo2DInterpolator& i_mi);
		virtual double getQualityRating(const Date& i_date);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		static bool windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static void getSynopticWind(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, double& VW, double& DW);

		Meteo2DInterpolator& mi;
		GridsManager& gdm;
		std::string base_algo_user, ref_station;
		bool inputIsAllZeroes;
		double dmax;
};

} //end namespace mio

#endif
