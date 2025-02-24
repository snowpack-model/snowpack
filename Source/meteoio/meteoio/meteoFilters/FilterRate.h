// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef FILTERRATE_H
#define FILTERRATE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterRate
 * @ingroup processing
 * @brief Rate of change filter.
 * @details
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point. It takes the following arguments:
 *  - MIN: minimum permissible rate of change (per seconds, optional);
 *  - MAX: either the absolute value of the maximum permissible rate of change (per seconds) if no other argument is provided, or maximum
 * permissible rate of change (per seconds) if a MIN was provided.
 *  - METHOD: which points are taken into account when computing the rate of change (optional, default is LEFT). It can be either one of:
 *      - LEFT: looking for the rate of change on the left of the current point, ie between the last valid point and the current one (default);
 *      - RIGHT: looking for the rate of change on the right of the current point, ie between the current point and the next valid one;
 *      - LEFT_AND_RIGHT: excluding a point if the rate of change exceeds MAX or MIN both on the left **and** on the right of the current point;
 *      - LEFT_OR_RIGHT: excluding a point if the rate of change exceeds MAX or MIN either on the left **or** on the right of the current point;
 *
 * So depending if MIN and MAX were provided or only MAX, every point where the local rate of change is outside <em>[ MIN , MAX]</em> or
 * every point outside  <em>[ -MAX , MAX]</em> is rejected.
 *
 * @code
 * TA::filter1   = rate
 * TA::arg1::MIN = -0.01
 * TA::arg1::MAX = 0.015
 * @endcode
 */
class FilterRate : public ProcessingBlock {
	public:
		FilterRate(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec) override;

	private:
         typedef enum IMPLEMENTATION_TYPE {
			LEFT,
			RIGHT,
            LEFT_AND_RIGHT,
            LEFT_OR_RIGHT
		} implementation_type;

		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
        static size_t findNextPoint(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& start_idx);
        static double getRate(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& idx, const size_t& cmp_idx);
        bool filterOut(const std::vector<MeteoData>& vecM, const unsigned int& param, const size_t& idx, const size_t& last_good, const size_t& next_good) const;

		double min_rate_of_change, max_rate_of_change;
        implementation_type methodParam; ///< controls which implementation of the filter is used
};

} //end namespace

#endif
