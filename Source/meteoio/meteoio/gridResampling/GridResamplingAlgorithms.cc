// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 MobyGIS Srl, Trento, Italy                                      */
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

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>
#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/gridResampling/GridTimeseriesResampling.h>

namespace mio {
/**
 * @page grid_resampling Grid resampling overview
 * In addition to spatial interpolations (from stations or grids with different resolutions) it is also possible to temporally resample grids.
 * Hence, if a set of grids is available at certain dates (from weather models for example), a number of resampling algorithms can be
 * configured to provide data at dates/times in between.
 *
 * @section grid_resampling_section Grid resampling section
 * First, grid resampling must be enabled by setting the REGRIDDING_STRATEGY to GRID_1DINTERPOLATE in the [InputEditing] section.
 * Note that this currently forces the path through the temporal resampling (regridding) algorithms, and any spatial resampling
 * configured (from point timeseries) will not take action.
 * Then, grid resampling is configured in the [GridInterpolations1D] section for each meteo parameter separately.
 *
 * For example:
 * @code
 * [InputEditing]
 * REGRIDDING_STRATEGY = GRID_1DINTERPOLATE
 *
 * [GridInterpolations1D]
 * TA::RESAMPLE = LINEAR
 * @endcode
 *
 * @section grid_algorithms_available Available Grid resampling algorithms
 * The following grid resampling algorithms are implemented:
 * - LINEAR: linear data resampling, see GridLinearResampling
 * - TIMESERIES: extract time series at all grid points and use Interpolations1D algorithm, see GridTimeseriesResampling
 *
 * @note By default, linear resampling will be used. It is possible to enable/disable all resampling with the ENABLE_GRID_RESAMPLING key
 * (in the [GridInterpolations1D] section). You can make MeteoIO write all grids that were resampled/regridded this way to the file system
 * by setting WRITE_RESAMPLED_GRIDS = TRUE (also in [GridInterpolations1D]).
 */

/**
 * @brief Facade constructor for a generic grid resampling algorithm.
 * @param[in] algorithm The current algorithm's semantic name.
 * @param[in] i_parname The current meteo parameter's identifier.
 * @param[in] dflt_window_size The default grid resampling window size.
 * @param[in] vecArgs vector of arguments (user settings) for this algorithm.
 */
GridResamplingAlgorithm::GridResamplingAlgorithm(const std::string& algorithm, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: algo(algorithm), parname(i_parname), grid_window_size(dflt_window_size)
{
	//do nothing
	(void)vecArgs;
}

/**
 * @brief Set this algorithm's window size to something other than the default value.
 * @param[in] window_size Desired window size in seconds.
 */
void GridResamplingAlgorithm::setWindowSize(const double& window_size)
{
	if (window_size <= 0.)
		throw InvalidArgumentException("Invalid WINDOW_SIZE for grid resampling algorithm", AT);
	grid_window_size = window_size / 86400.; //end user enters seconds, Julian internally
}

/**
 * @brief Object factory for temporal grid resampling algorithms.
 * @param[in] i_algorithm Semantic name of algorithm (as given in the INI file) to build.
 * @param[in] parname Meteo parameter to build the algorithm for.
 * @param[in] grid_window_size Standard window size for temporal grid resampling.
 * @param[in] vecArgs The algorithm's parameters as parsed from the user setings.
 * @param[in] cfg Config object in order to be able to read configuration keys.
 * @return A resampling algorithm object for the desired parameter.
 */
GridResamplingAlgorithm* GridResamplingAlgorithmsFactory::getAlgorithm(const std::string& i_algorithm, const std::string& parname,
	const double& grid_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	const std::string algorithm( IOUtils::strToUpper(i_algorithm) );
	if (algorithm == "LINEAR")
		return new GridLinearResampling(algorithm, parname, grid_window_size, vecArgs);
	else if (algorithm == "TIMESERIES")
		return new GridTimeseriesResampling(algorithm, parname, grid_window_size, vecArgs, cfg);
	else
		throw IOException("The temporal grid resampling algorithm '" + algorithm + "' is not implemented", AT);
}

} //namespace

