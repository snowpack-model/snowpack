#include <meteoio/MeteoIO.h>
#include <snowpack/snowpackCore/Solver.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mio;

// Forward declarations
void EL_INCID(const size_t &e, int Ie[]);

// PARAMETERS
const double tol_inf = 1e-6;  // Tolerance for the infinite-norm error
const double tol_2 = 1e-6;  // Tolerance for the 2-norm error

/********** Input file format **********/
/********** START **********/
// NAME OF DATA SET
// Symmetric Matrix: 0-th diagonal and 1-st diagonal
// 0-TH DIAGONAL AS LINE
// 1-ST DIAGONAL AS LINE
// Solution vector:
// SOLUTION VECTOR AS LINE
/********** END **********/

int main(int argc, char *argv[]) {

  /* Testing solver.h */

  if (argc < 1) {
    cout << "solverTest called without argument. Error.\n";
    cerr << "solverTest called without argument. Error.\n";
    exit(1);
  } else if (argc > 2) {
    cout << "solverTest called with too many arguments. Error.\n";
    cerr << "solverTest called with too many arguments. Error.\n";
    exit(1);
  }

  /* General memory allocations and declarations */
  // The FEM discretization uses 1D elements with two nodes
  void *Kt = NULL;
  int Ie[2];
  double Se[2][2];

  /* Read in test matrices */
  // Open input stream
  string str = argv[1];
  std::ifstream ifs;
  ifs.open(str.c_str());
  if (!ifs.is_open()) {
    cerr << "\nCould not open data file " << argv[1] << "\n";
    exit(1);
  }

  // Read in test name
  string testname, line;
  getline(ifs, testname);

  // Read in 0-th diagonal
  getline(ifs, line);
  getline(ifs, line);
  istringstream iss(line);
  std::vector<double> diag0, diag1, sol;
  double value;
  while (iss >> value) {
    diag0.push_back(value);
  }

  // Read in 1-th diagonal
  getline(ifs, line);
  iss.str(line);
  iss.clear();
  while (iss >> value) {
    diag1.push_back(value);
  }

  // Read in exact solution
  getline(ifs, line);
  getline(ifs, line);
  iss.str(line);
  iss.clear();
  while (iss >> value) {
    sol.push_back(value);
  }

  // Read in exact solution
  getline(ifs, line);
  getline(ifs, line);
  iss.str(line);
  iss.clear();
  bool rankDeficient;
  iss >> rankDeficient;

  // Close stream
  ifs.close();

  // Size of matrix and solution
  const size_t nN = diag0.size();
  const size_t nE = nN - 1;

  // Initialization of the solver base on 2 nodes 1-dim finite element with "pyramid" basis function
  if (Kt != NULL) {
    ds_Solve(ReleaseMatrixData, (SD_MATRIX_DATA*) Kt, 0);
  }
  ds_Initialize(nN, (SD_MATRIX_DATA**) &Kt);
  for (size_t e = 0; e < nE; e++) {
    int Nodes[2] = { (int) e, (int) e + 1 };
    ds_DefineConnectivity((SD_MATRIX_DATA*) Kt, 2, Nodes, 1, 0);
  }
  ds_Solve(SymbolicFactorize, (SD_MATRIX_DATA*) Kt, 0);

  /* Filling matrix A and compute RHS */
  Matrix A(nN, nN);
  Matrix x(nN, (size_t) 1);

  // Lower right block
  EL_INCID(nE - 1, Ie);
  Se[1][1] = A(nE + 1, nE + 1) = diag0.back();
  diag0.pop_back();
  Se[0][0] = A(nE, nE) = diag0.back();
  diag0.pop_back();
  Se[0][1] = Se[1][0] = A(nE, nE + 1) = A(nE + 1, nE) = diag1.back();
  diag1.pop_back();

  ds_AssembleMatrix((SD_MATRIX_DATA*) Kt, 2, Ie, 2, (double*) Se);

  // Last entry of solution vector
  x(nE + 1, 1) = sol.back();
  sol.pop_back();

  // Rest of A matrix and solution vector
  // Idea is to empty the diagonal vectors one entry after the other
  for (size_t e = nE - 1; !diag0.empty(); e--) {
    EL_INCID(e - 1, Ie);
    Se[0][0] = A(e, e) = diag0.back();
    diag0.pop_back();
    Se[0][1] = Se[1][0] = A(e, e + 1) = A(e + 1, e) = diag1.back();
    diag1.pop_back();
    Se[1][1] = 0.0;
    ds_AssembleMatrix((SD_MATRIX_DATA*) Kt, 2, Ie, 2, (double*) Se);

    x(e + 1, 1) = sol.back();
    sol.pop_back();
  }

  // Top entry of the solution vector
  x(1, 1) = sol.back();
  sol.pop_back();

  // Build RHS and test-solve with meteoIO matrix class
  Matrix b = A * x;
  vector<double> dU;
  for (size_t i = 1; i <= b.getNy(); i++) {
    dU.push_back(b(i, 1));
  }

  // Solve with solver.h
  if (ds_Solve(ComputeSolution, (SD_MATRIX_DATA*) Kt, dU.data())) {
    if (!rankDeficient) {  // Testing solver behavior if matrix is rank deficient!
      cerr << "Matrix from " << testname << " could not be inverted" << endl;
      exit(1);
    } else {  // if A is rank deficient and solver is unsuccessful then test is passed!
      return 0;
    }
  }

  // Retrieve result and for nan value
  Matrix xStar(dU.size(), (size_t) 1);
  size_t j = 1;
  for (vector<double>::iterator it = dU.begin(); it != dU.end(); it++) {
    xStar(j++, 1) = *(it);
  }


  /* Result comparison */
  // 2-norm and maximum/infinite norm
  Matrix error = x - xStar;
  double tmp = 0.0;
  double norm_error_inf = 0.0;
  for (size_t i = 1; i <= error.getNy(); i++) {
    tmp += error(i, 1) * error(i, 1);  // adding square terms
    if (abs(error(i, 1)) > norm_error_inf) {  // largest error so far
      norm_error_inf = error(i, 1);
    }
  }
  double norm_error_2 = sqrt(tmp);  // square root of the squares' sum

  // Check against given tolerances
  if (norm_error_2 > tol_2 || norm_error_inf > tol_inf) {
    cerr << setprecision(12) << "Error for Test " << testname
         << ": \n- inf-norm: " << norm_error_inf << " > " << tol_inf
         << "\n- 2-norm: " << norm_error_2 << " > " << tol_2 << "\n";
    exit(1);
  }

  return 0;
}

void EL_INCID(const size_t &e, int Ie[]) {
  Ie[0] = static_cast<int>(e);
  Ie[1] = static_cast<int>(e + 1);

}
