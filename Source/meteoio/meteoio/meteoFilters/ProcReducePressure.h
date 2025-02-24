// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCREDUCEPRESSURE_H
#define PROCREDUCEPRESSURE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcReducePressure
 * @ingroup processing
 * @brief Convert local pressure to sea level pressure.
 * @details
 * This implements a correction for the local pressure to be converted to a sea level pressure, 
 * according to the equation given in Atmosphere::reducedAirPressure().
 *
 * @code
 * P::filter1    = REDUCE_PRESSURE
 * @endcode
 */

class ProcReducePressure : public ProcessingBlock {
	public:
		ProcReducePressure(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec) override;

	private:
};

} //end namespace

#endif
