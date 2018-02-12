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
#include <alpine3d/snowdrift/checksum.h>

using namespace mio;
using namespace std;

double checksum(const CDoubleArray &x)
{
	double t=0;
	size_t n=x.getNx();
	for (size_t i=0;i<n;i++) {
		t+=x[i];
	}
	return t;
}

double checksum(const CDoubleArray &x, int start, int step)
{
	double t=0;
	size_t n=x.getNx();
	if (step<=0) {
		step=1;
	}
	for (size_t i=static_cast<size_t>(start);i<n;i+=static_cast<size_t>(step)) {
		t+=x[i];
	}
	return t;
}

double checksum(const CElementArray &x)
{
	double t=0;
	size_t n, tmp;
	x.size(n, tmp);
	for (size_t i=0;i<n;i++) {
		for (size_t j=0;j<8;j++) {
			t+=x(i,j);
		}
	}
	return t;
}

double checksum(const mio::Array2D<double> &x)
{
	size_t nx, ny;
	x.size(nx,ny);
	double t=0;
	for (size_t i=0;i<nx;i++) {
		for (size_t j=0;j<ny;j++) {
			t+=x(i,j);
		}
	}

	return t;
}

double checksum_rows(const mio::Array2D<double> &x, const size_t& from, size_t to)
{
	size_t nx, ny;
	x.size(nx,ny);
	double t=0;

	if (to>=nx) {
		to=nx-1;
	}

	for (size_t i=from;i<=to;i++) {
		for (size_t j=0;j<ny;j++) {
			t+=x(i,j);
		}
	}

	return t;
}

double checksum_cols(const mio::Array2D<double> &x, const size_t& from, size_t to)
{
	size_t nx, ny;
	x.size(nx,ny);
	double t=0;

	if (to>=ny) {
		to=ny-1;
	}

	for (size_t i=0;i<nx;i++) {
		for (size_t j=from;j<=to;j++) {
			t+=x(i,j);
		}
	}
	return t;
}

double checksum_c(const Grid3DObject &grid)
{
    const double t = grid.grid3D.getMean()*static_cast<double>(grid.grid3D.getCount());
    return t;

}

double checksum(const mio::Array1D<SnowStation>&x)
{
	const size_t n = x.getNx();
	double t1=0;
	bool first=true;

	for (size_t ii=0; ii<n; ii++)
	{
		double t=0;
		t+=x[ii].Albedo+x[ii].SoilAlb+x[ii].BareSoil_z0+static_cast<double>(x[ii].SoilNode)+x[ii].cH+x[ii].mH+x[ii].Ground;
		t+=x[ii].hn+x[ii].rho_hn+x[ii].windward+static_cast<double>(x[ii].ErosionLevel);
		t+=x[ii].S_class1+x[ii].S_class2+x[ii].S_d+x[ii].z_S_d+x[ii].S_s+x[ii].z_S_s+
			x[ii].S_n+x[ii].z_S_n+x[ii].S_4+x[ii].z_S_4+x[ii].S_5+x[ii].z_S_5+static_cast<double>(x[ii].getNumberOfNodes());

		t+=checksum(x[ii].Ndata, x[ii].getNumberOfNodes());
		t+=checksum(x[ii].Edata, x[ii].getNumberOfElements());
		if (x[ii].Kt!=NULL) {
			t+=0.1;
		}

		t+=checksum(x[ii].Cdata);

		if (t>1.0E20 && first) {
			first=false;
			cerr << "STATION VALUE TO BIG AT INDEX #" << ii << "\n";
		}

		t1+=t;
	}
	return t1;

}

double checksum(const std::vector<NodeData>& x, const size_t n)
{
	double t=0;
	for (size_t i=0;i<n;i++) {
		t+=x[i].z+x[i].u+x[i].f+x[i].udot+x[i].T+x[i].S_n+x[i].S_s+x[i].hoar;
	}
	return t;

}

double checksum(const std::vector<ElementData>& x, const size_t n)
{
	double t=0;
	for (size_t i=0; i<n; i++) {
		t+=x[i].depositionDate.getJulian()+x[i].L0+x[i].L+x[i].Te+x[i].gradT+x[i].Rho+x[i].M;
		for (int j=0;j<N_COMPONENTS;j++) {
			t+=x[i].theta[j];
			for (unsigned int k=0;k<SnowStation::number_of_solutes;k++) {
				t += (x[i]).conc(j,k);
			}
		}
		for (int j=0;j<N_SN_FIELDS;j++) {
			if (x[i].k[j] < Constants::big) {
				t+=x[i].k[j]+x[i].c[j]+x[i].soil[j];
			}
		}

		t+=x[i].sw_abs+x[i].rg+x[i].dd+x[i].sp+x[i].rb+x[i].N3+static_cast<double>(x[i].mk+x[i].type);

		t+=x[i].dEps+x[i].Eps+x[i].Eps_e+x[i].Eps_v+x[i].Eps_Dot+x[i].Eps_vDot+x[i].S+x[i].C+ x[i].S_dr+x[i].hard;
	}
	return t;
}

double checksum(const CanopyData &x)
{
	const size_t n = sizeof(CanopyData)/sizeof(double);
	double *xp=(double *)(&x);
	double t=0;
	for (size_t i=0;i<n;i++, xp++) {
		t+=*xp;
	}
	return t;
}

