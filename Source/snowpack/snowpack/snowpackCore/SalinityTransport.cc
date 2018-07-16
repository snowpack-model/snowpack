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

#include <snowpack/snowpackCore/SalinityTransport.h>
#include <snowpack/Utils.h>
#include <stdio.h>

#ifdef CLAPACK
	// Matching C data types with FORTRAN data types (taken from f2c.h):
	typedef long int integer;
	typedef double doublereal;

	// Declare the function interfaces with the LAPACK library (taken from clapack.h):
	extern "C" {
		/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
		doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
		ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *info);

		/* Subroutine */ int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
		a, integer *lda, doublereal *s, doublereal *u, integer *ldu,
		doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
		integer *iwork, integer *info);

		/* Subroutine */ int dgtsv_(integer *n, integer *nrhs, doublereal *dl,
		doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer
		*info);
	}
#endif


#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif


/**
 * @brief Class for solving diffusion-advection equation for salinity using the Crank-Nicolson implicit method\n
 * Solve Richards Equation \n
 * @author Nander Wever
 * @param nE Domain size
 */
SalinityTransport::SalinityTransport(const size_t nE)
           : flux_up(), flux_down(), dz_(), dz_up(), dz_down(), theta1(), theta2(), BrineSal(), D(), sb(), NumberOfElements(0)
{
	SetDomainSize(nE);
}


/**
 * @brief Resizing vectors to match given domain size \n
 * @author Nander Wever
 * @param nE Domain size
 */
void SalinityTransport::SetDomainSize(size_t nE) {
	NumberOfElements = nE;

	flux_up.resize(nE, 0.);
	flux_down.resize(nE, 0.);
	dz_.resize(nE, 0.);
	dz_up.resize(nE, 0.);
	dz_down.resize(nE, 0.);
	theta1.resize(nE, 0.);
	theta2.resize(nE, 0.);
	BrineSal.resize(nE, 0.);
	D.resize(nE, 0.);
	sb.resize(nE, 0.);
	return;
}


/**
 * @brief Solve diffusion-advection equation using the Crank-Nicolson implicit, or fully implicit method\n
 * @author Nander Wever
 * @par This function solves the following equation ($n$ and $i$ denoting time and spatial level, respectively), in Latex code:
\begin{multline}
\frac{ \left ( \theta^{n+1}_i S_{\mathrm{b}, i}^{n+1} - \theta^{n}_i S_{\mathrm{b}, i}^{n} \right ) } { \Delta t } \\
- f \left [ \left ( \frac{ 2 D_{i+1}^{n} \theta^{n+1}_{i+1} S_{\mathrm{b}, i+1}^{n+1} }{ \Delta z_{\mathrm{up}} \left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } - \frac{ 2 D_{i}^{n} \theta^{n+1}_{i} S_{\mathrm{b}, i}^{n+1} }{\left ( \Delta z_{\mathrm{up}} \Delta z_{\mathrm{down}} \right ) } + \frac{ D_{i-1}^{n} \theta^{n+1}_{i-1} S_{\mathrm{b}, i-1}^{n+1} }{ \Delta z_{\mathrm{down}} \left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } \right ) \right ] \\
- \left ( 1-f \right ) \left [ \left ( \frac{ 2 D_{i+1}^{n} \theta^{n}_{i+1} S_{\mathrm{b}, i+1}^{n} }{ \Delta z_{\mathrm{up}} \left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } - \frac{ 2 D_{i}^{n} \theta^{n}_{i} S_{\mathrm{b}, i}^{n} }{\left ( \Delta z_{\mathrm{up}} \Delta z_{\mathrm{down}} \right ) } + \frac{ D_{i-1}^{n} \theta^{n}_{i-1} S_{\mathrm{b}, i-1}^{n} }{ \Delta z_{\mathrm{down}} \left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } \right ) \right ] \\
- f \left [ \left ( \frac{q^{n}_{i+1} S_{\mathrm{b},i+1}^{n+1} - q^{n}_{i-1} S_{\mathrm{b},i-1}^{n+1}}{\left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } \right ) \right ] - \left ( 1-f \right ) \left [ \left ( \frac{q^{n}_{i+1} S_{\mathrm{b},i+1}^{n} - q^{n}_{i-1} S_{\mathrm{b},i-1}^{n}}{\left ( \Delta z_{\mathrm{up}} + \Delta z_{\mathrm{down}} \right ) } \right ) \right ] - s_{\mathrm{sb}} = 0
\end{multline}
Here, $f=1$ results in the fully implicit scheme, whereas $f=0.5$ corresponds to the Crank-Nicolson scheme. The implicit scheme is first order accurate, whereas the Crank-Nicolson scheme is second order accurate. Furthermore, both are unconditionally stable and suffer only minimal numerical diffusion for the advection part. As with many other common schemes, the advection part is not perfectly conserving sharp transitions. Futhermore, the reason to not use the fully implicit or the Crank Nicolson scheme is the occurrence of spurious oscillations in the solution, which negatively impact the accuracy of the simulations more than the negative effect on computational efficiency imposed by the CFL criterion required for the explicit method (see SalinityTransport::SolveSalinityTransportEquationExcplicit).
 * @param dt Time step (s)
 * @param DeltaSal Result vector (change in salinity over time step)
 */
bool SalinityTransport::SolveSalinityTransportEquationImplicit(const double dt, std::vector <double> &DeltaSal) {

	if(NumberOfElements==0) return false;	// Nothing to do

	const double f = .5;			// Factor: f = .5 is Crank-Nicolson, f = 1. is fully implicit.

	const bool WriteDebugOutput = false;
	if(WriteDebugOutput) setvbuf(stdout, NULL, _IONBF, 0);

	// Declare and initialize l.h.s. matrix and r.h.s. vector
	std::vector<double> ad(NumberOfElements, 0.);		// Matrix diagonal
	std::vector<double> adu(NumberOfElements-1, 0.);	// Matrix upper diagonal
	std::vector<double> adl(NumberOfElements-1, 0.);	// Matrix lower diagonal
	std::vector<double> b(NumberOfElements, 0.);		// Vector

	std::vector<double> theta12(NumberOfElements, 0.);	// Theta at iteration level for diffusion
	theta12=theta1;						// Assign iteration level theta

	// Fill matrix and r.h.s. vector
	for(size_t i = 0; i < NumberOfElements; i++) {
		// The matrix diagonal, the time derivative:
		ad[i] += theta2[i] / dt;
		// The matrix diagonal, the diffusion part:
		ad[i] += (2. * D[i] * theta12[i]) / (dz_up[i] * dz_down[i]);

		// The lower diagonal
		if(i!=0) {
			// the diffusion part:
			adl[i-1] += -f * 2. * D[i-1] * theta12[i-1] / (dz_down[i] * (dz_up[i] + dz_down[i]));

			// the advection part:
			adl[i-1] += f * flux_down[i] / (dz_up[i] + dz_down[i]);
		}

		// The upper diagonal
		if(i!=NumberOfElements-1) {
			// the diffusion part:
			adu[i] += -f * 2. * D[i+1] * theta12[i+1] / (dz_up[i] * (dz_up[i] + dz_down[i]));

			// the advection part:
			adu[i] += -f * flux_up[i] / (dz_up[i] + dz_down[i]);
		}


		// The r.h.s. vector time derivative:
		b[i] += (theta1[i] * BrineSal[i]) / dt;

		// The r.h.s. vector diffusion part:
		if(i==NumberOfElements-1) {
			b[i] += 0.;
		} else {
			b[i] += (1. - f) * (2. * D[i+1] * theta12[i+1] * BrineSal[i+1]) / (dz_up[i] * (dz_up[i] + dz_down[i]));
		}
		b[i] += -2. * D[i] * theta12[i] * BrineSal[i] / (dz_up[i] * dz_down[i]);
		if(i==0) {
			b[i] += (1. - f) * 2. * D[i] * SeaIce::OceanSalinity / (dz_down[i] * (dz_up[i] + dz_down[i]));
		} else {
			b[i] += (1. - f) * 2. * D[i-1] * theta12[i-1]  * BrineSal[i-1] / (dz_down[i] * (dz_up[i] + dz_down[i]));
		}

		//The r.h.s. vector advection part:
		if(i==0 && i==NumberOfElements-1) {
			// TODO: What to do in the case of only 1 element??
			throw;
		} else if (i==0) {
			b[i] += (1. - f) * (flux_up[i] * BrineSal[i+1] - flux_down[i] * SeaIce::OceanSalinity) / (dz_up[i] + dz_down[i]);
		} else if (i==NumberOfElements-1) {
			b[i] += (1. - f) * (flux_up[i] * 0. - flux_down[i] * BrineSal[i-1]) / (dz_up[i] + dz_down[i]);
		} else {
			b[i] += (1. - f) * (flux_up[i] * BrineSal[i+1] - flux_down[i] * BrineSal[i-1]) / (dz_up[i] + dz_down[i]);
		}

		// The r.h.s. vector source/sink term:
		b[i] += -sb[i];
	}


	// Deal with boundary conditions:

	// Add the terms from "out of boundary" diffusion
	b[0] += -f * (2. * D[0] * SeaIce::OceanSalinity) / (dz_down[0] * (dz_up[0] + dz_down[0]));
	b[NumberOfElements-1] += f * (2. * D[NumberOfElements-1] * 0.) / (dz_up[NumberOfElements-1] * (dz_up[NumberOfElements-1] + dz_down[NumberOfElements-1]));

	// Add the terms from "out of boundary" advection
	b[0] += -f * (flux_down[0] * SeaIce::OceanSalinity) / (dz_up[0] + dz_down[0]);
	b[NumberOfElements-1] += f * (flux_down[NumberOfElements-1] * 0.) / (dz_up[NumberOfElements-1] + dz_down[NumberOfElements-1]);


	if(WriteDebugOutput) {
		std::cout << "SalinityTransport.cc > Coefficients:\n";
		for(size_t i = 0; i < NumberOfElements; i++) {
			if(i==NumberOfElements-1) {
				std::cout << i << ": " << flux_up[i] << " " << flux_down[i] << " " << BrineSal[i] << " " << ad[i] << " " << "---" << " " << "---" << " " << b[i] << "\n";
			} else {
				std::cout << i << ": " << flux_up[i] << " " << flux_down[i] << " " << BrineSal[i] << " " << ad[i] << " " << adl[i] << " " << adu[i] << " " << b[i] << "\n";
			}
		}
	}


	// Call solver
	const int matrixdimensions=int(NumberOfElements);	// Cast from size_t to int is necessary, to interface correctly with LAPACK dgtsv_.
#ifdef CLAPACK
	// Call LAPACK DGTSV: Solver for tridiagonal matrices, with partial pivoting.
	int info=0;
	const int vectordimensions=1;
	dgtsv_( (integer*) &matrixdimensions, (integer*) &vectordimensions, &adl[0], &ad[0], &adu[0], &b[0], (integer*) &matrixdimensions, (integer*) &info );

	if(info!=0) {
		//= 0: successful exit
		//< 0: if INFO = -i, the i-th argument had an illegal value
		//> 0: if INFO = i, U(i,i) is exactly zero, and the solution
		//    has not been computed.  The factorization has not been
		//    completed unless i = N.
		std::cout << "[E] Error in SalinityTransport.cc: DGTSV failed [info = " << info << "].\n";
		return false;
	}
#else
	// Call TDMASolver: Thomas algorithm for tidiagonal matrices. Not the recommended choice, but useful when LAPACK is not available.
	std::vector<double> b_ = b;
	const int ret = ReSolver1d::TDMASolver(matrixdimensions, &adl[0], &ad[0], &adu[0], &b[0], &b_[0]);
	b=b_;
	if (ret != 0) {
		std::cout << "[E] Error in SalinityTransport.cc: TDMA failed.\n";
		std::cout << "    Using LAPACK (see compile options) may increase numerical stability in SalinityTransport.\n";
		return false;
	}
#endif


	// Apply solution
	if(WriteDebugOutput) std::cout << "SalinityTransport.cc > Solution vector:\n";
	for(size_t i=0; i<NumberOfElements; i++) {
		if(WriteDebugOutput) std::cout << i << ": " << b[i] << "\n";
		DeltaSal[i]=b[i]-BrineSal[i];
		BrineSal[i]=b[i];
	}

	return true;
}


/**
 * @brief Solve diffusion-advection equation using the upwind explicit method\n
 * @author Nander Wever
 * @param dt Time step (s)
 * @param DeltaSal Result vector (change in salinity over time step)
 */
bool SalinityTransport::SolveSalinityTransportEquationExplicit(const double dt, std::vector <double> &DeltaSal) {

	if(NumberOfElements==0) return false;	// Nothing to do

	const bool WriteDebugOutput = false;
	if(WriteDebugOutput) setvbuf(stdout, NULL, _IONBF, 0);

	// Declare vectors
	std::vector<double> b(NumberOfElements, 0.);		// Solution vector
	std::vector<double> theta12(NumberOfElements, 0.);	// Theta at iteration level for diffusion
	theta12=theta1;						// Assign iteration level theta

	// Fill matrix and r.h.s. vector
	for(size_t i = 0; i < NumberOfElements; i++) {
		// Explicit upwind scheme for advection:
		const double tmp_flux = (flux_up[i] * dz_up[i] + flux_down[i] * dz_down[i]) / (dz_up[i] + dz_down[i]);
		// We assume that the incoming flux at the top element consists of rain and has no salinity, and the incoming flux at the bottom element consists of ocean salinity.
		const double TopFluxSalinity = 0.;
		b[i] += (theta1[i] * BrineSal[i]) +  (  (tmp_flux > 0.)   ? (tmp_flux * dt * (((i==NumberOfElements-1) ? (TopFluxSalinity) : (BrineSal[i+1])) - BrineSal[i]) / dz_up[i])            : (tmp_flux * dt * ((BrineSal[i] - ((i==0) ? (SeaIce::OceanSalinity) : (BrineSal[i-1]))) / dz_down[i]))  );

		// Explicit scheme for diffusion
		b[i] += dt * ( ((i==NumberOfElements-1) ? (D[i] * 0.) : (theta12[i+1] * D[i+1] * BrineSal[i+1])) / (dz_up[i]*(dz_up[i]+dz_down[i])) - (2. * theta12[i] * D[i] * BrineSal[i]) / (dz_up[i]+dz_down[i]) + (((i==0) ? (D[i] * SeaIce::OceanSalinity) : (theta12[i-1] * D[i-1] * BrineSal[i-1]))) / (dz_down[i]*(dz_up[i]+dz_down[i])) );

		// Source/sink term
		b[i] += -sb[i];
	}

	// Apply solution
	for(size_t i=0; i<NumberOfElements; i++) {
		DeltaSal[i]=b[i] / theta2[i] - BrineSal[i];
		BrineSal[i]=b[i] / theta2[i];
	}

	return true;
}


/**
 * @brief Check for CFL criterion\n
 * @author Nander Wever
 * @param dt Time step (s)
 * @return true when provided time step dt satisfies CFL criterion, false otherwise.
 */
bool SalinityTransport::VerifyCFL(const double dt)
{
	const double CFL_limit = 0.98;
	for(size_t i = 0; i < NumberOfElements; i++) {
		if (std::max(fabs(flux_up[i]), fabs(flux_down[i])) * dt / std::min(dz_up[i], dz_down[i]) > CFL_limit) return false;
	}
	return true;
}
