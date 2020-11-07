/***********************************************************************************/
/*  Copyright 2020 ALPsolut.eu                                                     */
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

#ifndef SNOWLINE_ALGORITHM_H
#define SNOWLINE_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

#include <meteoio/thirdParty/tinyexpr.h>
#include <string>
#include <utility>
#include <vector>

namespace mio {

class SnowlineAlgorithm : public InterpolationAlgorithm {
	public:
		SnowlineAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs,
		    const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
		    GridsManager& i_gdm, Meteo2DInterpolator& i_mi);
		virtual double getQualityRating(const Date& i_date);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);

	private:
		struct aspect {
			double deg_beg;
			double deg_end;
			double sl_elevation;
			aspect() : deg_beg(0.), deg_end(0.), sl_elevation(IOUtils::nodata) {}
			//to check if a vector of aspects is sorted:
			bool operator< (const aspect &other) const { return this->deg_beg < other.deg_beg; }
			bool operator<= (const aspect &other) const { return this->deg_beg <= other.deg_beg; }
		};

		void baseInterpol(const double& snowline, const DEMObject& dem, Grid2DObject& grid);
		void assimilateCutoff(const double& snowline, const DEMObject& dem, Grid2DObject& grid);
		void assimilateBands(const double& snowline, const DEMObject& dem, Grid2DObject& grid);
		void assimilateFormula(const double& snowline, const DEMObject& dem, Grid2DObject& grid);
		Grid2DObject mergeSlopes(const DEMObject& dem, const std::vector<Grid2DObject>& azi_grids);
		void initExpressionVars(const std::vector< std::pair<std::string, double> >& substitutions, te_variable* vars) const;
		te_expr* compileExpression(const std::string& expression, const te_variable* te_vars, const size_t& sz) const;
		std::vector<aspect> readSnowlineFile();
		void getSnowlines();
		double probeTrend();
		double calculateRate(const double& snowline);
		std::vector< std::pair<std::string, std::string> > prepareBaseArgs(const double& snowline, const std::string& base_alg_fallback);
		void msg(const std::string& message);

		template <class T> static bool isSorted(std::vector<T> vec, const bool& strict = false)
		{ //check if a vector is sorted (with 'strict', duplicates are not allowed either)
			if (vec.empty())
				return true;
			for (typename std::vector<T>::iterator it = vec.begin() + 1; it != vec.end(); ++it) {
				if (strict) {
					if (*it <= *(it - 1))
						return false;
				} else {
					if (*it < *(it - 1))
						return false;
				}
			}
			return true;
		}

		GridsManager& gdm_;
		Meteo2DInterpolator& mi_;
		std::string base_alg_;
		std::vector< std::pair<std::string, std::string> > input_args_; //store ini input args from constructor
		std::string where_;
		void (SnowlineAlgorithm::*assimilateFunction_)(const double& snowline, const DEMObject& dem, Grid2DObject& grid); //points to enabled assimilation function
		bool smoothing_;
		std::vector<aspect> snowlines_;
		std::string snowlines_file_;
		bool enforce_positive_rate_; //deduce lapse rate from snowline if data has "reversed" rate
		bool calc_base_rate_; //deduce lapse rate from snowline in any case
		double fallback_rate_;
		double cutoff_val_;
		double band_height_;
		unsigned int band_no_;
		std::string formula_;
		bool has_warned_deduced_trend_;
		bool has_warned_fixed_trend_;
		bool verbose_; //suppress warnings?
};

} //end namespace mio

#endif //SNOWLINE_ALGORITHM_H
