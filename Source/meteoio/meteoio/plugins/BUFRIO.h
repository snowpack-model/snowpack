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
 *
 * @ingroup plugins
 * @author Patrick Leibersperger
 * @date   2024-08-06
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
        bool verbose_out;

        // Cryo Station specific
        bool write_cryo;
        long wigos_id_series, wigos_issuer, wigos_issue_no, station_type, surface_type, snow_depth_method;
        std::string wigos_local_id;
        
		static const std::string template_filename;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999

        void parseInputSection();
        void parseOutputSection();
        void setStationData(CodesHandlePtr &message, const StationData &station, const Coords &position, const std::string &subset_prefix);
        void setWIGOSId(CodesHandlePtr &message, const std::string &subset_prefix);


};

} //namespace
#endif
