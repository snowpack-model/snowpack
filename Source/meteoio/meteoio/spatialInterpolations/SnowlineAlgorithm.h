// SPDX-License-Identifier: LGPL-3.0-or-later
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

/**
 * @class SnowlineAlgorithm
 * @ingroup spatialization
 * @brief Assimilation of snowline elevation information into snow maps. 
 * @details
 * This algorithm applies a base interpolation method and reads in snowline elevation information to be able to<br>
 * a) set snow heights below the snowline elevation to zero<br>
 * b) check the lapse rate of the base algorithm, and if this is "reversed" inject one derived from the snowline elevation.
 *
 * A problem that may turn up for practitioners when creating snow maps is that a station higher up is more eposed to wind,
 * or may measure less snow than a station further down in the valley for any other reason. Then, if applying an algorithm
 * that computes a lapse rate, snow is transported downwards resulting in maps that show little snow on the mountain tops
 * and a lot of snow in the valleys.
 *
 * To combat this, external knowledge is required. In case of this SNOWLINE algorithm the user provides a text file containing
 * information about the (estimated) snowline elevation, i. e. the elevation above which there is snow. This can for example
 * come from satellite images, webcams, or manual observations. If it is deemed unphysical in the calculated domain to have
 * less snow further up, then the trend can be forced to be the other way around.
 *
 * The algorithm takes the following parameters:
 *  - BASE: Name of the base interpolation algorithm. Can be any 2d interpolation available. Of course, injecting a new rate only works if the base algorithm uses a lapse rate (default: `IDW_LAPSE`).
 *  - SNOWLINEFILE: Text file containing the snowline elevation(s), potentially split up into multiple slope aspects (see below).
 *  - SNOWLINE: The snowline elevation can be supplied in the ini file directly which is most useful for standalone experiments and testing (single aspect only).
 *  - ENFORCE_POSITIVE_RATE: If the trend is "reversed" as explained above, calculate it from the snowline elevation and the highest available station (default: false).
 *  - CALC_BASE_RATE: Use a lapse rate deduced from the snowline in any case, i. e. not only when it is "reversed" in the data (default: false).
 *  - FALLBACK_RATE: If the trend is reversed but no snowline information is available, use this lapse rate.
 *  - METHOD: Decide how to cut off snow below the snowline elevation. Can be `CUTOFF` (hard cut to zero below snowline, default), `BANDS` (simple linear gradient from zero
 *    to the measured data), or `FORMULA` (evaluate all pixels according to a custom formula).
 *  - BAND_NO: For `BANDS` smoothing, use this many steps to go from zero snow to measured snow (default: 10).
 *  - BAND_HEIGHT: For `BANDS` smoothing, use bands of this height (default: 10 m).
 *  - FORMULA: Expression for `FORMULA` smoothing. The same mathematics as in FilterMaths are available. In addition, you can use the following substitutions:
 *    `snowline` (the snowline elevation), `altitude` (the pixel's altitude) and `param` (the original value as calculated by the base algorithm).
 *  - VERBOSE: Print some info messages (default: true)? 
 *  - SET: Below the snow line, do not set snow to 0 but to this value.
 *
 * In addition, you can set the base algorithm's parameters as usual.
 *
 * The minimal settings to control the default `IDW_LAPSE` interpolation by setting `HS` to zero below the snowline elevation are the following:
 * @code
 * HS::ALGORITHMS             = SNOWLINE
 * HS::SNOWLINE::SNOWLINEFILE = /home/snow/sat/snowlines.txt
 * @endcode
 *
 * To calculate the lapse rate from the snowline and inject it in an `IDW_SLOPES` interpolation if necessary (using a fixed rate if the snowline is n/a) use the following setup:
 * @code
 * HS::ALGORITHMS                      = SNOWLINE
 * HS::SNOWLINE::BASE                  = IDW_SLOPES
 * HS::SNOWLINE::METHOD                = BANDS
 * HS::SNOWLINE::SNOWLINEFILE          = /home/snow/sat/snowlines.txt
 * HS::SNOWLINE::ENFORCE_POSITIVE_RATE = TRUE
 * HS::SNOWLINE::FALLBACK_RATE         = 0.0005
 * @endcode
 *
 * The format specification of the text file containing snowline elevations per slope aspect is as follows: `aspect_#,beg_azimuth end_azimuth snowline`
 * with one line per slope aspect. For example:
 * @code
 * #Format: aspect_#,beg_azimuth end_azimuth snowline
 * #Snowline in meters, azimuth in degrees.
 * aspect_1,45 135 2200
 * aspect_2,135 225 2400
 * aspect_3,225 315 2200
 * aspect_4,315 45 2100
 * @endcode
 *
 * Set elevations to the `nodata` value (usually -999) to mark them as unknown.
 *
 * A minimal file for a single slope aspect (complete circle from 0 to 360/0 degrees) might look like this for a snowline elevation of 2500 m:
 * @code
 * aspect_1,0 0 2500
 * @endcode
 *
 * @author Michael Reisecker
 * @date 2020-09
 */
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
