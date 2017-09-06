#include <time.h>
#include <meteoio/MeteoIO.h>

using namespace std;
using namespace mio;

bool array1d(const unsigned int& n) {
	cout << "Testing Array1D<double>\n";
	bool status = true;
	srand((unsigned)time(0));
	const double range = 10.;

	Array1D<double> grid1(n);

	for(unsigned int ii=0; ii<grid1.getNx(); ii++) {
		grid1(ii) = (double)rand()/(double)RAND_MAX*range;
	}

	Array1D<double> grid2=grid1;
	grid2 -= 8.;
	grid2 /= 2.;
	grid2 += 4.;
	grid2 *= 2.;
	if(!grid2.checkEpsilonEquality(grid1, 1e-6)) {
		cout << "\terror: basic operations with constants fail!\n";
		status=false;
	}

	grid2.resize(n, 0.);
	grid2 += grid1;
	grid2 /= grid1;
	if(!IOUtils::checkEpsilonEquality(grid2.getMean(), 1., 1e-6)) {
		cout << "\terror: grids addition or division fails!\n";
		status=false;
	}
	grid2 *= grid1;
	if(!grid2.checkEpsilonEquality(grid1, 1e-6)) {
		cout << "\terror: grids multiplication fails!\n";
		status=false;
	}
	grid2 -= grid1;
	if(!IOUtils::checkEpsilonEquality(grid2.getMean(), 0., 1e-6)) {
		cout << "\terror: grids substraction fails!\n";
		status=false;
	}

	return status;
}

bool grid2d(const unsigned int& n) {
	cout << "Testing Grid2DObject\n";
	bool status = true;
	srand((unsigned)time(0));
	const double range = 10.;
	Coords llcorner("CH1903","");
	llcorner.setXY(785425. , 191124., 1400.);

	Grid2DObject grid1(n, n, 100, llcorner);

	for(unsigned int jj=0; jj<grid1.getNy(); jj++) {
		for(unsigned int ii=0; ii<grid1.getNx(); ii++) {
			grid1.grid2D(ii,jj) = (double)rand()/(double)RAND_MAX*range;
		}
	}

	Grid2DObject grid2(100, llcorner, grid1.grid2D);
	grid2.grid2D -= 8.;
	grid2.grid2D /= 2.;
	grid2.grid2D += 4.;
	grid2.grid2D *= 2.;
	if(!grid2.grid2D.checkEpsilonEquality(grid1.grid2D, 1e-6)) {
		cout << "\terror: basic operations with constants fail!\n";
		status=false;
	}

	grid2.set(n, n, 100, llcorner, 0.);
	grid2.grid2D += grid1.grid2D;
	grid2.grid2D /= grid1.grid2D;
	if(!IOUtils::checkEpsilonEquality(grid2.grid2D.getMean(), 1., 1e-6)) {
		cout << "\terror: grids addition or division fails!\n";
		status=false;
	}
	grid2.grid2D *= grid1.grid2D;
	if(!grid2.grid2D.checkEpsilonEquality(grid1.grid2D, 1e-6)) {
		cout << "\terror: grids multiplication fails!\n";
		status=false;
	}
	grid2.grid2D -= grid1.grid2D;
	if(!IOUtils::checkEpsilonEquality(grid2.grid2D.getMean(), 0., 1e-6)) {
		cout << "\terror: grids substraction fails!\n";
		status=false;
	}

	return status;
}

bool grid3d(const unsigned int& n) {
	cout << "Testing Grid3DObject\n";
	bool status = true;
	srand((unsigned)time(0));
	const double range = 10.;
	Coords llcorner("CH1903","");
	llcorner.setXY(785425. , 191124., 1400.);

	Grid3DObject grid1(n, n, n, 100, llcorner);

	for(unsigned int kk=0; kk<grid1.getNz(); kk++) {
		for(unsigned int jj=0; jj<grid1.getNy(); jj++) {
			for(unsigned int ii=0; ii<grid1.getNx(); ii++) {
				grid1.grid3D(ii,jj,kk) = (double)rand()/(double)RAND_MAX*range;
			}
		}
	}

	Grid3DObject grid2(100, llcorner, grid1.grid3D);
	grid2.grid3D -= 8.;
	grid2.grid3D /= 2.;
	grid2.grid3D += 4.;
	grid2.grid3D *= 2.;
	if(!grid2.grid3D.checkEpsilonEquality(grid1.grid3D, 1e-6)) {
		cout << "\terror: basic operations with constants fail!\n";
		status=false;
	}

	grid2.set(n, n, n, 100, llcorner, 0.);
	grid2.grid3D += grid1.grid3D;
	grid2.grid3D /= grid1.grid3D;
	if(!IOUtils::checkEpsilonEquality(grid2.grid3D.getMean(), 1., 1e-6)) {
		cout << "\terror: grids addition or division fails!\n";
		status=false;
	}
	grid2.grid3D *= grid1.grid3D;
	if(!grid2.grid3D.checkEpsilonEquality(grid1.grid3D, 1e-6)) {
		cout << "\terror: grids multiplication fails!\n";
		status=false;
	}
	grid2.grid3D -= grid1.grid3D;
	if(!IOUtils::checkEpsilonEquality(grid2.grid3D.getMean(), 0., 1e-6)) {
		cout << "\terror: grids substraction fails!\n";
		status=false;
	}

	return status;
}

bool matrix(const size_t& n) {
	cout << "Testing Matrix\n";
	bool status=true;
	Matrix I(n,1.); //build an n*n identity matrix
	Matrix m1(n,n);
	m1.random(10.);

	const double det = m1.det();
	const Matrix m1_trans = m1.getT();
	const double det_trans = m1_trans.det();

	if(!IOUtils::checkEpsilonEquality(det, det_trans, Matrix::epsilon_mtr*fabs(det))) {
		cout << "\terror: m1.det != m1T.det\n";
		status=false;
	}

	Matrix m2 = m1-8.;
	m2 /= 2.;
	m2 += 4.;
	m2 *= 2.;
	if(m2!=m1) {
		cout << "\terror: basic operations with constants fail!\n";
		status=false;
	}

	m2 = m1+m1;
	m2 -= m1;
	if(m2 != m1) {
		cout << "\terror when adding/substracting matrix\n";
		status=false;
	}


	m2 = m1.getInv();
	Matrix m3=m1*m2;
	if(m3.isIdentity()!=true) {
		cout << "\terror: m1*inv(m1) is NOT identity matrix\n";
		status=false;
	}

	m3 = m1*m1;
	if(m2*m3 != m1) {
		cout << "\terror when multiplying matrix\n";
		status=false;
	}

	Matrix m4=Matrix::solve(m1,I); //solve m1*X=I
	if(m1*m4 != I){
		cout << "\terror when solving A*X=B\n";
		status=false;
	}

	Matrix L,U;
	if(m1.LU(L,U)==false) {
		cout << "\terror: LU decomposition could NOT be computed\n";
		status=false;
	}

	return status;
}


int main() {
	const unsigned int n=50;

	const bool grid1d_status = array1d(n);
	const bool grid2d_status = grid2d(n);
	const bool grid3d_status = grid3d(n);
	const bool matrix_status = matrix(n);
	if(grid1d_status!=true || grid2d_status!=true || grid3d_status!=true || matrix_status!=true) throw IOException("Grid/Matrix error", AT);
	return 0;
}
