/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file ReSolver1d.h
 * @version 10.02
 */

#ifndef SALINITY_H
#define SALINITY_H

#include <meteoio/MeteoIO.h>

/**
 * @class Salinity
 * @version 18.07
 * @author Nander Wever
 * @brief This module contains the solver for the diffusion-advection equation for the transport of salinity
 */
class Salinity {

	public:
		bool SetDomainSize(size_t nE);

		std::vector<double> flux_up;	//Flux with element above (negative=upward, positive=downward)
		std::vector<double> flux_down;	//Flux with element below (negative=upward, positive=downward)

	private:
		size_t NumberOfElements;
};
#endif //End of Salinity.h
