/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OMPCONTROL_H
#define OMPCONTROL_H
#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

#include <cstdio>

namespace OMPControl
{
    /**
		* @brief Returns the parameters for splitting an array in several, balanced sub-arrays.
		* This is mostly usefull for parallel calculations, where an array will be split and sent to different
		* workers.
		* @param[in] dimx number of cells in the desired dimension
		* @param[in] nbworkers total number of slices
		* @param[in] idx_wk current slice index (starting at 0)
		* @param[out] startx_sub calculated start index for the current slice
		* @param[out] nx_sub calculated number of cells (in the desired dimension) of the current slice
		*/
  void getArraySliceParams(const size_t& dimx, const size_t& nbworkers, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub);
  void getArraySliceParamsOptim(const size_t& nbworkers, const std::vector<SnowStation*>&, const mio::DEMObject& mpi_sub_dem,
                                const mio::Grid2DObject& mpi_sub_landuse, std::vector<std::vector<size_t> >& omp_snow_stations_ind);

}
#endif
