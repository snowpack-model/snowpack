#include <meteoio/MeteoIO.h>

using namespace std;
using namespace mio; //The MeteoIO namespace is called mio

//please give a size for the square matrix on the command line.
//(in this example, we only consider square matrices even if this is not mandatory)
//Obviously, a large number will fill your screen with the output.
int main(int /*argc*/, char** argv) {

	size_t n;
	IOUtils::convertString(n, argv[1]);

	//we create a suqare matrix of size n, filled with random numbers
	Matrix m1(n,n);
	m1.random(10.);

	std::cout << "Randomly filled (" << n << "," << n << ") matrix\n";
	std::cout << "\t\tm1=" << m1.toString(); //print the matrix on the screen
	std::cout << "=> det(m1)=" << m1.det() << "\n"; //determinant
	std::cout << "\t\tT(m1)=" << m1.getT().toString(); //transpose
	std::cout << "\t\tm1-2=" << (m1-2.).toString(); //arithmetic with scalars
	std::cout << "\t\tm1*m1=" << (m1*m1).toString(); //arithmetic with matrices
	Matrix m2= m1.getInv();
	std::cout << "\t\tinv(m1)=" << m2.toString(); //inverse
	Matrix m3=m1*m2;
	if(m3.isIdentity()==true) //test if this is the identity matrix
		std::cout << "=> m1*inv(m1) is identity matrix\n";
	else
		std::cout << "=> m1*inv(m1) is NOT identity matrix\n";

	Matrix I(n,1.); //build an n*n identity matrix
	Matrix m4=Matrix::solve(m1,I); //solve m1*X=I
	std::cout << "\tIn m1*X=I, X=" << m4.toString(); //print the matrix on the screen

	Matrix L,U;
	if(m1.LU(L,U)==false) { //LU factorization
		std::cout << "\tLU factorization of m1 can not be computed\n";
	} else {
		std::cout << "\tLU factorization of m1: L=" << L.toString() << "\t\tU=" << U.toString();
	}

	return 0;
}
