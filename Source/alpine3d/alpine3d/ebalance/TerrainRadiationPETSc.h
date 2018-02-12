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
#ifndef TERRAINRADIATIONPETSC_H
#define TERRAINRADIATIONPETSC_H

#include <petscksp.h>

#include <meteoio/MeteoIO.h>

#include <alpine3d/MPIControl.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/ViewFactors.h>

class TerrainRadiationPETSc : public TerrainRadiationAlgorithm {

	public:
		TerrainRadiationPETSc(const mio::Config& i_cfg, const mio::DEMObject& dem_in, const int& i_nbworkers, const std::string& method);
		~TerrainRadiationPETSc();

		void getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain);
		void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta,
		              const mio::Array2D<double>& rh, const mio::Array2D<double>& ilwr);

	private:
		const mio::DEMObject& dem;
		const unsigned int dimx, dimy;

		Mat M, VF;
		Vec ALBEDO, x, b;
		KSP ksp;

		PetscInt Mstart, Mend;
		PetscInt Istart, Iend, N, B;

		int multiIndexToRowIndex(const PetscInt& i, const PetscInt& N);
		int multiIndexToColIndex(const PetscInt& i, const PetscInt& N);
		PetscInt indexToMultiIndex(const unsigned int& i, const unsigned int& j, const PetscInt& N);
};

#endif
