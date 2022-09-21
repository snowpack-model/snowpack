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
#ifndef SYNTHETICIO_H
#define SYNTHETICIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {

/**
 * @class SynthIO
 * @brief This plugin generate synthetic data
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2022-09-19
 */
class SynthIO : public IOInterface {
	public:
		SynthIO(const std::string& configfile);
		SynthIO(const SynthIO&);
		SynthIO(const Config& cfgreader);

		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void init();
		
		const Config cfg;
		std::vector<StationData> vecStations;
		Date dt_start, dt_end;
		double dt_step;
};

} //namespace
#endif
