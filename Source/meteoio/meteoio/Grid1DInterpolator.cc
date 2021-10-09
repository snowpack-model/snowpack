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

#include <meteoio/Grid1DInterpolator.h>

#include <iostream>
#include <utility>

using namespace std;

namespace mio {

const std::string Grid1DInterpolator::section_name("GridInterpolations1D");

/**
 * @brief Constructor for a grid interpolator.
 * @details On initialization, a Grid1DInterpolator instance prepares a map of per-parameter
 * interpolation algorithms which are then ready to do the work when called.
 * @param[in] in_cfg The current simulation configuration.
 */
Grid1DInterpolator::Grid1DInterpolator(const Config& in_cfg) : algorithm_map(), cfg(in_cfg)
{
	cfg.getValue("WINDOW_SIZE", section_name, grid_window_size, IOUtils::nothrow);
	if (grid_window_size <= 1.)
		throw IOException("WINDOW_SIZE for grids not valid, it should be a duration in seconds at least greater than 1", AT);
	grid_window_size /= 86400.; //user inputs seconds, internally we use Julian
	cfg.getValue("ENABLE_GRID_RESAMPLING", section_name, enable_grid_resampling, IOUtils::nothrow);

	//create the grid resampling algorithms for each MeteoData::Parameters entry:
	for (size_t ii = 0; ii < MeteoData::nrOfParameters; ++ii) { //loop over all MeteoData member variables
		const std::string parname( MeteoData::getParameterName(ii) ); //current semantic parameter name
		const std::string algo_name( IOUtils::strToUpper(getGridAlgorithmForParameter(parname)) );

		//parse the algorithm's INI file settings and create an object for each algorithm:
		const std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getArgumentsForAlgorithm(
			parname, algo_name, section_name) );
		algorithm_map[parname] = GridResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname,
			grid_window_size, vecArgs);

		//set generic parameters available for all algorithms:
		const std::string where( "GridInterpolations1D::" + parname + "::" + algo_name );
		for (size_t jj = 0; jj < vecArgs.size(); ++jj) {
			if (vecArgs[jj].first == "WINDOW_SIZE") {
				double algo_window_size;
				IOUtils::parseArg(vecArgs[jj], where, algo_window_size);
				algorithm_map[parname]->setWindowSize(algo_window_size);
			}
		} //endfor jj
	} //endfor ii
}

/**
 * @brief Destructor for a Grid1DInterpolator instance.
 * @details Upon destruction all algorithms previously created through the object factory
 * are deleted here.
 */
Grid1DInterpolator::~Grid1DInterpolator()
{
	for (auto it = algorithm_map.begin(); it != algorithm_map.end(); ++it)
		delete it->second;
}

/**
 * @brief Wrapper function to call the desired grid resampling algorithm.
 * @details The grid processor calls this function which in turn redirects the
 * request to the proper algorithm object.
 * @param[in] date The date to resample to.
 * @param[in] parameter The meteo parameter to resample.
 * @param[in] available_grids A list of grids on the hard drive (or generated ones) which
 * should be used for temporal resampling inbetween them, as well as the corresponding dates.
 * @param[out] resampled_grid The grid that was temporally resampled.
 * @return True if resampling could be performed. This means that the setup is in order,
 * but nodata may still be returned.
 */
bool Grid1DInterpolator::resampleData(const Date& date, const MeteoGrids::Parameters& parameter,
	const std::map<Date, Grid2DObject>& available_grids, Grid2DObject& resampled_grid)
{
	//translate MeteoGrids parameter to MeteoData parameter:
	const MeteoData::Parameters mpar( MeteoData::findMeteoParam(parameter) );
	const std::string mparname( MeteoData::getParameterName(mpar) );
	if (algorithm_map.find(mparname) == algorithm_map.end())
		return false; //no algorithm configured for this parameter
	algorithm_map[mparname]->resample(date, available_grids, resampled_grid);
	return true; //was able to call interpolation routine
}

/**
 * @brief Parse the algorithm to use for a parameter from the INI file.
 * @param[in] parname Meteo parameter to find the interpolation algorithm for.
 * @return Name of the algorithm to use (default: LINEAR).
 */
std::string Grid1DInterpolator::getGridAlgorithmForParameter(const std::string& parname) const
{
	std::string algorithm( "linear" ); //default algorithm
	cfg.getValue(parname + "::resample", section_name, algorithm, IOUtils::nothrow);
	return algorithm;
}

} //namespace
