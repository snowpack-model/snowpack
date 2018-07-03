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

#include <snowpack/snowpackCore/Salinity.h>
#include <snowpack/Utils.h>

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

using namespace std;
using namespace mio;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif


bool Salinity::SetDomainSize(size_t nE) {
	NumberOfElements=nE;
	flux_up.resize(nE);
	flux_down.resize(nE);
	return true;
}
