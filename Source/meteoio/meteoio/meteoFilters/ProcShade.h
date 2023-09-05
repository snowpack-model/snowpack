// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCSHADE_H
#define PROCSHADE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcShade
 * @ingroup processing
 * @brief Apply a shading mask to the Incoming or Reflected Short Wave Radiation
 * @details
 * A shading mask that is either computed from the DEM or read from a separate file will be applied to the radiation
 * and combined with the radiation splitting model in order to properly compute the shading effects on the measurement point. 
 * This mask will be linearly interpolated between the provided points in order to be applied to the true sun position.
 * 
 * When providing the shading mask in a separate file, the same mask will be applied to all stations. It simply need to
 * contain the horizon elevation (in deg.) as a function of azimuth (in deg.), prefixed by either a stationID or "*" as fallback 
 * wildcard for all stationIDs:
 * @code
 * DAV1	0	5
 * DAV1	15	25
 * DAV1	45	12
 * DAV1	180	30
 * DAV1	270	20
 * SLF2	0	7.5
 * SLF2	180	10
 * *	0	5
 * *	180	5
 * @endcode
 *
 * Therefore, the filter supports the following arguments:
 *  - INFILE: file (and path) to read the mask from (optional);
 *  - OUTFILE: file (and path) to write the mask to (optional, this is usefull if some masks have been computed from the DEM).
 * 
 * Then the filter is declared with the file name containing the horizon mask as argument:
 * @code
 * ISWR::filter1    = SHADE
 * ISWR::arg1::infile = ../input/iswr_mask.dat
 * @endcode
 *
 * If no arguments are provided, then it will compute the mask from the Digital Elevation Model. In such as case, 
 * a DEM must be declared in the [Input] section and must contain the stations of interest as a mask will be computed for
 * each station. Please make sure that the extent of the DEM is appropriate to correctly compute the shading effects!
 */

class ProcShade : public ProcessingBlock {
	public:
		ProcShade(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config &i_cfg);
		~ProcShade();

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		static std::vector< std::pair<double,double> > computeMask(const DEMObject& i_dem, const StationData& sd);

		const Config &cfg;
		DEMObject dem;
		std::map< std::string , std::vector< std::pair<double,double> > > masks;
		std::string horizons_outfile;
		bool has_dem;

		static const double diffuse_thresh;
};

} //end namespace

#endif
