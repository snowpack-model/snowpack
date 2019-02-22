/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef FILTERHUNHEATEDPSUM_H
#define FILTERHUNHEATEDPSUM_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterUnheatedPSUM
 * @ingroup processing
 * @brief Filters out snow melting in an unheated rain gauge.
 * @details
 * This filter can ONLY be applied to precipitation. Non-zero measurements are accepted only if they take place
 * when the relative humidity is greater than 0.5 and (TA-TSS) < 3, otherwise they get reset to 0.
 * It can take two optional arguments overwriting these thresholds (THRESH_RH and THRESH_DT). If none of these conditions could be tested
 * (for lack of data), then the precipitation is reset to nodata. On the contrary, if the "soft" option is given,
 * the lack of validation data keeps the precipitation as it is.
 *
 * @code
 * PSUM::filter2 = unheated_raingauge
 * PSUM::arg2::soft = TRUE
 * PSUM::arg2::thresh_rh = 0.5
 * PSUM::arg2::thresh_dt = 3.
 * @endcode
 */

class FilterUnheatedPSUM : public ProcessingBlock {
	public:
		FilterUnheatedPSUM(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		double thresh_rh, thresh_Dt;
		bool is_soft;
};

} //end namespace

#endif
