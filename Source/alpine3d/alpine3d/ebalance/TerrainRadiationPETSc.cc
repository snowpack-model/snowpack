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
#include <alpine3d/ebalance/TerrainRadiationPETSc.h>

using namespace std;
using namespace mio;

TerrainRadiationPETSc::TerrainRadiationPETSc(const mio::Config& cfg, const mio::DEMObject& dem_in, const int&, const std::string& method)
                      : TerrainRadiationAlgorithm(method), dem(dem_in), dimx(dem_in.getNx()), dimy(dem_in.getNy())
{
	cout << "[i] Calculating radiation with TerrainRadiationPETSc method" << endl;
	bool isMaster = MPIControl::instance().master();
	N = dimx*dimy;
	B = dimy;

	ViewFactors viewFactors(cfg, dem_in);
	//viewFactors.setBoundarys(Istart, Iend, 0, dimy);
	//viewFactors.generate();

	MatCreate(MPI_COMM_WORLD, &VF);
	MatSetSizes(VF, PETSC_DECIDE, PETSC_DECIDE, N, N);
	MatSetFromOptions(VF);
	MatSetUp(VF);
	MatGetOwnershipRange(VF, &Istart, &Iend);

	PetscInt DIAGONALnonZeroValues[Iend-Istart];
	PetscInt OFFDIAGONALnonZeroValues[Iend-Istart];

	for (int Ii = Istart; Ii < Iend ; Ii++) {
		DIAGONALnonZeroValues[Ii-Istart] = 0;
		OFFDIAGONALnonZeroValues[Ii-Istart] = 0;
	}

	viewFactors.generate_visible_matrix();
	for (int Ii = Istart; Ii < Iend ; Ii++) {
		const unsigned int i = multiIndexToRowIndex(Ii,B);
		const unsigned int j = multiIndexToColIndex(Ii,B);
		for (int Ij = 0; Ij < N; Ij++) {
			const unsigned int a = multiIndexToRowIndex(Ij,B);
			const unsigned int b = multiIndexToColIndex(Ij,B);
			const bool entries = viewFactors.getVisibleMatrixEntry(a,b,i,j);

			if ((entries) || (Ii == Ij)) {
				if (Ij>=Istart && Ij < Iend) {
					DIAGONALnonZeroValues[Ii-Istart]++;
				} else {
					OFFDIAGONALnonZeroValues[Ii-Istart]++;
				}
			}
		}
	}

	MatMPIAIJSetPreallocation(VF, N, DIAGONALnonZeroValues, N, OFFDIAGONALnonZeroValues);

	PetscInt Ii, Ij;
	PetscScalar v;

	if (isMaster) cout << "[i] TerrainRadiation PETSc, copying view factors:" << endl;
	for (Ii=Istart; Ii<Iend; Ii++) {
		const unsigned int i_from = multiIndexToRowIndex(Ii,B);
		const unsigned int j_from = multiIndexToColIndex(Ii,B);
		if (isMaster && (Ii%(Iend/100)==0))
			cout << (Ii-Istart)*100./(Iend-Istart) << "%" << "\r";
		for (Ij=0; Ij<N; Ij++) {
			const unsigned int i_to = multiIndexToRowIndex(Ij,B);
			const unsigned int j_to = multiIndexToColIndex(Ij,B);

			v = viewFactors.GetViewfactor(i_from,j_from,i_to,j_to);

			if ((v != 0) || (Ij==Ii))
				MatSetValuesBlocked(VF,1,&Ij,1,&Ii,&v,INSERT_VALUES);
		}

	}


	MatAssemblyBegin(VF,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(VF,MAT_FINAL_ASSEMBLY);
	MatDuplicate(VF, MAT_COPY_VALUES, &M);
	VecCreate(MPI_COMM_WORLD, &b);

	VecSetSizes(b, PETSC_DECIDE, N);
	VecSetFromOptions(b);
	VecDuplicate(b, &x);
	VecDuplicate(b, &ALBEDO);
	VecSet(ALBEDO, 0.5);

	KSPCreate(MPI_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, M, M, SAME_NONZERO_PATTERN);
	KSPSetFromOptions(ksp);
}

TerrainRadiationPETSc::~TerrainRadiationPETSc()
{
	KSPDestroy(&ksp);
	VecDestroy(&b);
	VecDestroy(&x);
	VecDestroy(&ALBEDO);
	MatDestroy(&M);
}


void TerrainRadiationPETSc::getRadiation(const mio::Array2D<double>& direct, mio::Array2D<double>& diffuse, mio::Array2D<double>& terrain)
{
	Array2D<double> Qin_raw(dimx,dimy,0);
	Qin_raw += diffuse;
	Qin_raw += direct;

	PetscInt bi;
	PetscScalar s;

	for (bi=Istart; bi<Iend; bi++) {
		const unsigned int i = multiIndexToRowIndex(bi,B);
		const unsigned int j = multiIndexToColIndex(bi,B);
		s = -1.0 * Qin_raw(i,j);

		VecSetValues(b,1,&bi,&s,INSERT_VALUES);
	}

	MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(b);
	VecAssemblyEnd(b);

	KSPSetOperators(ksp, M, M, SAME_NONZERO_PATTERN);
	KSPSolve(ksp, b, x);

	terrain.resize(dimx,dimy,0);
	PetscInt xi;

	for (xi=Istart; xi<Iend; xi++) {
		const unsigned int i = multiIndexToRowIndex(xi,B);
		const unsigned int j = multiIndexToColIndex(xi,B);
		VecGetValues(x, 1, &xi, &s);
		const double value = s-Qin_raw(i,j);
		if (value >= 0)
			terrain(i,j) = value;
		else
			terrain(i,j) = 0.;
	}

	//MPIControl::instance().allreduce_sum(terrain);
	return terrain;
}

int TerrainRadiationPETSc::multiIndexToRowIndex(const PetscInt& i, const PetscInt& N)
{
	return floor(i/N);
}

int TerrainRadiationPETSc::multiIndexToColIndex(const PetscInt& i, const PetscInt& N)
{
	return i % N;
}

PetscInt TerrainRadiationPETSc::indexToMultiIndex(const unsigned int& i, const unsigned int& j, const PetscInt& N)
{
	return N*i + j;
}

void TerrainRadiationPETSc::setMeteo(const mio::Array2D<double> &albedo,const mio::Array2D<double> &/*ta*/,const mio::Array2D<double> &/*rh*/,const mio::Array2D<double> &/*ilwr*/)
{
	PetscInt ALBEDOi;
	PetscScalar s;
	for (ALBEDOi=Istart; ALBEDOi<Iend; ALBEDOi++) {
		const unsigned int i = multiIndexToRowIndex(ALBEDOi,B);
		const unsigned int j = multiIndexToColIndex(ALBEDOi,B);
		s = albedo(i,j);
		VecSetValues(ALBEDO,1,&ALBEDOi,&s,INSERT_VALUES);
	}

	MatCopy(VF, M, SAME_NONZERO_PATTERN);

	VecAssemblyBegin(ALBEDO);
	VecAssemblyEnd(ALBEDO);

	MatDiagonalScale(M, ALBEDO, NULL);
	MatShift(M,-1.0);
}
