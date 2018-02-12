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
#ifndef CHECKSUM_H
#define CHECKSUM_H

#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>
#include <alpine3d/snowdrift/SnowDrift.h>

double checksum(const CDoubleArray &x);
double checksum(const CDoubleArray &x, int start, int step);

double checksum(const CElementArray &x);
double checksum(const mio::Array2D<double> &x);
double checksum_rows(const mio::Array2D<double> &x, const size_t& from, size_t to);
double checksum_cols(const mio::Array2D<double> &x, const size_t& from, size_t to);
double checksum_c(const mio::Grid3DObject &grid);

//data structures from Snowpack
double checksum(const mio::Array1D<SnowStation> &x);
double checksum(const std::vector<ElementData>& x, const size_t n);
double checksum(const std::vector<NodeData>& x, const size_t n);
double checksum(const CanopyData &x);

#endif
