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
#ifndef BUFRIO_H
#define BUFRIO_H

#include "BUFRFile.h"

#include <string>

namespace mio {

/**
 * @class BUFRIO
 * @brief This (empty) class is to be used as a bufrio for developing new plugins
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2010-06-14
 */
class BUFRIO : public IOInterface {
	public:
		BUFRIO(const std::string& configfile);
		BUFRIO(const Config& cfgreader);

		virtual void readStationData(const Date &date, std::vector<StationData> &vecStation);
        virtual void readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecMeteo);

        virtual void writeMeteoData(const std::vector<std::vector<MeteoData>> &vecMeteo, const std::string &name = "");

	private:
		const Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters

		std::vector<BUFRFile> station_files;
        std::vector<std::string> additional_params;
        
        std::string outpath;
        bool separate_stations;

		static const std::string template_filename;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999

        void parseInputSection();
        void parseOutputSection();


};

} //namespace
#endif
