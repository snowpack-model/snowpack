// SPDX-License-Identifier: LGPL-3.0-or-later
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
#ifndef USERINTERPOLATION_H
#define USERINTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class USERInterpolation
 * @ingroup spatialization
 * @brief Reads user provided gridded data on the disk.
 * @details
 * By default, the grids are all in the GRID2DPATH directory given in the [Input] section and the file extension is assumed to
 * be ".asc". But the following arguments allow overriding it:
 *  - SUBDIR: look for grids in the provided subdirectory of GRID2DPATH;
 *  - EXT: use another file extension.
 *  - NAMING: provide file name template using YYYY for year, MM for month, DD for day, hh for hour, mm for minute, ss for second (optional) and PARAM for parameter name.
 *		Examples:
 *		- param::user::naming = YYYYMMDDhhmmss_PARAM (default) searches for files named: <b>{numeric date with second resolution}_{capitalized meteo parameter}.{ext}</b>, for example 20081201150000_TA.asc.
 *		- TA::user::naming = YYYY-MM-DDThh.mm_PARAM searches for files named 2008-12-01T15.00_TA.asc.
 *  - TIME_CONSTANT: if true, use the same grid for all timesteps (default: false). If TIME_CONSTANT is set to true, files must be named according to:
 *    <b>{capitalized meteo parameter}.{ext}</b>, for example TA.asc.
 *
 *  If no grid exists for a given timestamp and parameter, the algorithm returns a zero rating so any other interpolation algorithm can pickup
 * and provide a fallback. Therefore, it is not necessary to provide grids for all time steps but one can focuss on only the relevant and interesting
 * time steps. The following argument changes this behavior:
 *  - LOWEST_PRIORITY: by default, when a user provided grid exists it will have priority over any other interpolation algorithm. Setting LOWEST_PRIORITY 
 * to TRUE gives user provided grids the lowest priority so such grids are only used when no other spatial interpolation algorithm can provide interpolations.
 *
 * The meteo parameters can be found in \ref meteoparam "MeteoData". Examples of use:
 * @code
 * TA::algorithms = USER	# read grids from GRID2DPATH using the GRID2D plugin
 *
 * VW::algorithms   = USER	# read grids from GRID2DPATH/wind
 * VW::user::subdir = wind
 *
 * PSUM::algorithms   = USER	# read grids from GRID2DPATH/precip with the ".dat" extension
 * PSUM::user::subdir = precip
 * PSUM::user::ext    = .dat
 * PSUM::user::naming = YYYY-MM-DDThh.mm.ss_PARAM
 *
 * TSG::algorithms    = USER     # read grids from GRID2DPATH using the GRID2D plugin
 * TSG::user::ext     = .asc
 * TSG::user::time_constant  = TRUE  # use the same grid for all timesteps.
 * @endcode
 */
class USERInterpolation : public InterpolationAlgorithm {
	public:
		USERInterpolation(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
		                                GridsManager& i_gdm);
		virtual double getQualityRating(const Date& i_date);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::string getGridFileName() const;
		GridsManager& gdm;
		std::string filename, grid2d_path;
		std::string subdir, file_ext, naming;
		bool time_constant, lowest_priority;
};

} //end namespace mio

#endif
