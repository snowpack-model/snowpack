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
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <ctime>

#include <meteoio/MeteoIO.h>

#include <alpine3d/snowdrift/SnowDrift.h>

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Function    : setQuadraturePoints
// Authors     : Henning Loewe
//--------------------------------------------------------------------------------
// Description : sets the SnowDrift variable qValues which holds the
// quadrature points for a one point gaussian quadrature
//
//        IP_8-----IP_7
//        /|       /|              /\  z
//       / |      / |               |
//      /  |     /  |               |
//    IP_5-----IP_6 |               |------>  y
//     |   |    |   |              /
//     |  IP_4--|--IP_3           /
//     |  /     |  /             /
//     | /      | /              V  x
//     |/       |/
//    IP_1------IP_2
//--------------------------------------------------------------------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SnowDriftA3D::setQuadraturePoints( void )
{
  double IP[3][8]={ { 1, 1,-1,-1, 1, 1,-1,-1},
		    {-1, 1, 1,-1,-1, 1, 1,-1},
		    {-1,-1,-1,-1, 1, 1, 1, 1} };

  double qValue = 0.577350269189626;
  for ( int i = 0; i < 3; i++)	{
      for ( int j = 0; j < 8; j++)	{
	  		qPoint(i,j) = qValue * IP[i][j];
			}
  }
}

/**
 * @brief Compute Jacobian
 * This function calculates the Jacobian at the given
 * point P and and stores the value of its determinant at point P in
 * the global variable DETERMINANTJ[k] the input value is the
 * respective element, the position where to store the
 * DETERMINANTJ[k], and the function is of type void, the Jacobian
 * matrix J and the integration point P and the DETERMINANTJ array
 * @param DETERMINANTJ Determinant of Jacobian matrix
 * @param J Jacobian matrix
 * @param element element number
 * @param P point where Jacobian has to be calculated
 * @param k for storing Jacobian at DETERMINANTJ(k)
 */
void SnowDriftA3D::Jacobian(double *DETERMINANTJ,
			 double J[][3],
			 const int element,
			 const double *P,
			 const int k,
			 const int ix,
			 const int iy,
			 const int iz)
{
  double G[3][8];  //dummy array

  // the partial derivatives of phi evaluated at the integration point P
  TTfun(G,P);
  // set Jacobian J to zero
  for (int i = 0; i < 3; i++) {
      for ( int j = 0; j < 3; j++) {
	  		J[i][j] = 0;
			}
  }

  int ii = 0, jj = 0, kk = 0;

  for (int j=0; j<3; j++){
	int count=0;
	while (count<8){
	    if (count==0 || count==3 || count==4 || count==7){
		ii=ix;
	    } else{
		ii=ix+1;
	    }
	    if (count==0 || count==1 || count==4 || count==5){
		jj=iy;
	    } else{
		jj=iy+1;
	    }
	    if (count<=3){
		kk=iz;
	    } else{
		kk=iz+1;
	    }

	    J[0][j] = J[0][j] + nodes_x.grid3D(ii,jj,kk) * G[j][count];
	    J[1][j] = J[1][j] + nodes_y.grid3D(ii,jj,kk) * G[j][count];
	    J[2][j] = J[2][j] + nodes_z.grid3D(ii,jj,kk) * G[j][count];
	    count++;
	}
  }


  // now calculate the determinant
  DETERMINANTJ[k] = J[0][0] * (J[1][1]*J[2][2]-J[2][1]*J[1][2])
    -J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])
    +J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]);

  if (DETERMINANTJ[k]<0){
    cout<<"the determinant of the element "<<element<<" is negative, watch out"<<endl;
  }

}

/**
 * @brief J0 function
 * computes the matrix of cofactors of the input matrix J and stores it in J0
 * @param J0 matrix of cofactors (transpose of the adjoint of J)
 * @param J Jacobian matrix
 */
void SnowDriftA3D::J0fun( double J0[3][3], const double J[3][3])
{
  J0[0][0] = J[1][1]*J[2][2]-J[2][1]*J[1][2];
  J0[0][1] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2]);
  J0[0][2] = J[1][0]*J[2][1]-J[1][1]*J[2][0];

  J0[1][0] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2]);
  J0[1][1] = J[0][0]*J[2][2]-J[2][0]*J[0][2];
  J0[1][2] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1]);

  J0[2][0] = J[0][1]*J[1][2]-J[1][1]*J[0][2];
  J0[2][1] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2]);
  J0[2][2] = J[0][0]*J[1][1]-J[1][0]*J[0][1];
}

/**
 * @brief GQIntB
 * returns the integral of the product of basis
 * functions i and j over a hexahedral element using gaussian
 * quadrature, DETERMINANTK[0],...DETERMINANTK[7] contains the
 * determinant of the jacobian of the mapping onto the unit cube
 * evaluated at the corners
 * @param DETERMINANTJ determinant of Jacobian
 * @param i nr shape function
 * @param j nr shape function
 */
double SnowDriftA3D::GQIntB( double *DETERMINANTJ,const int i, const int j)
{
  double q_loc=0;
  double PHI[8];
  int k;
  double cad[3];

  for ( k = 0; k < 8; k++) {
      cad[0]=qPoint(0,k);
      cad[1]=qPoint(1,k);
      cad[2]=qPoint(2,k);

      phi(PHI,cad);

      q_loc += PHI[i] * PHI[j] * fabs(DETERMINANTJ[k]);  // i and j have already been adjustet!
  }
  return q_loc;
}

/**
 * @brief GQIntC
 * returns the integral of the product of the gradients
 * of basis functions i and j over a hexahedral element using gaussian
 * quadrature, DETERMINANTK[0],...DETERMINANTK[7] contains the
 * determinant of the jacobian of the mapping onto the unit cube
 * evaluated at the corners
 * @param DETERMINANTJ determinant of Jacobian
 * @param J0M
 * @param i pointer/index
 * @param j pointer/index
 * @param b drift vector
 * @param K Diffusion tensor
 */
double SnowDriftA3D::GQIntC(double * DETERMINANTJ,
			 const double J0M[3][3][8],
			 const int i, const int j,
			 const double b[3],
			 const double K[3][3])
{
  double q_loc=0;
  double PHI[8], TT[3][8], Hi[3], Hj[3];
  double A[3], B[3], C[3];   //pure dummy vectors
  int k,m;
  double cad[3];


  for ( k = 0;k<8;k++){

    cad[0]=qPoint(0,k);		//the integration point
    cad[1]=qPoint(1,k);
    cad[2]=qPoint(2,k);

    phi(PHI,cad);
    TTfun(TT,cad);

    Hi[0]=TT[0][i];		//dummy vector containing the gradient of the ith shape function
    Hi[1]=TT[1][i];
    Hi[2]=TT[2][i];

    Hj[0]=TT[0][j];		// dummy vector containing the gradient of the jth shape function
    Hj[1]=TT[1][j];
    Hj[2]=TT[2][j];

    /* the first part is due to diffusion */

    /* first calculate the intermediate vectors*/

    for (m=0;m<3;m++){   /* set to zero*/
      A[m]=0;
      B[m]=0;
      C[m]=0;
    }

    for (m=0;m<3;m++){
      A[0]=A[0]+J0M[0][m][k]*Hj[m];   /* this is J0M*Hj */
      A[1]=A[1]+J0M[1][m][k]*Hj[m];
      A[2]=A[2]+J0M[2][m][k]*Hj[m];
    }

    for (m=0;m<3;m++){
      B[0]=B[0]+K[0][m]*A[m];  /* this is K*A  */
      B[1]=B[1]+K[1][m]*A[m];
      B[2]=B[2]+K[2][m]*A[m];

    }

    for (m=0;m<3;m++){
      C[0]=C[0]+J0M[0][m][k]*Hi[m];	/* this is J0M*Hi  */
      C[1]=C[1]+J0M[1][m][k]*Hi[m];
      C[2]=C[2]+J0M[2][m][k]*Hi[m];
    }


    for (m=0;m<3;m++) {
      q_loc=q_loc+B[m]*C[m]*(1/fabs(DETERMINANTJ[k]));    /* add up the diffusion part */
		}

    /* now we add the part that is due to convection*/

    // A is the matrix product J0M(k)*Hj and has been calculated previously, no need to repeat this calculation

    if (fabs(DETERMINANTJ[k])<1e-10)
      printf("\n the determinant of J at point k is zero or almost zero, watch out\n");

    for (m=0;m<3;m++){
      q_loc=q_loc+b[m]*A[m]*PHI[i]*(DETERMINANTJ[k]/fabs(DETERMINANTJ[k]));
		}
  }	//end loop in k

  return q_loc;

}

/**
 * @brief GQIntApdx
 * returns the integral of the product of a basis
 * function i and the gradient of basis function j over a hexahedral
 * element using gaussian quadrature,
 * DETERMINANTK[0],...DETERMINANTK[7] contains the determinant of the
 * jacobian of the mapping onto the unit cube evaluated at the corners
 * @param DETERMINANTJ determinant of Jacobian
 * @param J0M
 * @param i nr basis function
 * @param j nr basis function
 * @param b drift vector
 * @param deltak parameter SUPG method
 */
double SnowDriftA3D::GQIntApdx(double *DETERMINANTJ,
			    const double J0M[3][3][8],
			    const int i, const int j,
			    double b[],
			    const double deltak)
{
  double q_loc=0;
  double A[3];
  double PHI[8], TT[3][8], Hi[3];
  double cad[3];
  int m,k;

  for (k=0;k<8;k++){

    cad[0]=qPoint(0,k);
    cad[1]=qPoint(1,k);
    cad[2]=qPoint(2,k);

    phi(PHI,cad);
    TTfun(TT,cad);

    Hi[0]=TT[0][i];		//dummy vector containing the gradient of the ith shape function
    Hi[1]=TT[1][i];
    Hi[2]=TT[2][i];

    for (m=0;m<3;m++){
      A[m]=0;
		}

    for (m=0;m<3;m++){
      A[0]+=J0M[0][m][k]*Hi[m];
      A[1]+=J0M[1][m][k]*Hi[m];
      A[2]+=J0M[2][m][k]*Hi[m];
    }


    for (m=0;m<3;m++){
      q_loc+=PHI[j]*b[m]*A[m]*(DETERMINANTJ[k]/fabs(DETERMINANTJ[k]));
    }
	}
  q_loc=q_loc*deltak;
  return q_loc;
}

/**
 * @brief GQIntAdxdx
 * ...
 * @param DETERMINANTJ determinant of Jacobian
 * @param J0M
 * @param i pointer
 * @param j pointer
 * @param *b drift  vector
 * @param deltak parameter SUPG method
 */
double SnowDriftA3D::GQIntAdxdx( double *DETERMINANTJ,
			      double J0M[3][3][8],
			      const int i, const int j,
			      double *b,
			      const double deltak)
{
  double q_loc=0; double z_loc,y_loc;
  double TT[3][8], Hi[3], Hj[3];
  double A[3],C[3];
  double cad[3];

  for (int k = 0; k < 8; k++)	{

      cad[0] = qPoint(0,k);   //the kth integration point
      cad[1] = qPoint(1,k);
      cad[2] = qPoint(2,k);

      TTfun(TT,cad);

      Hi[0] = TT[0][i];		//dummy vector containing the gradient of the ith shape function
      Hi[1] = TT[1][i];
      Hi[2] = TT[2][i];

      Hj[0] = TT[0][j];		// dummy vector containing the gradient of the jth shape function
      Hj[1] = TT[1][j];
      Hj[2] = TT[2][j];

    /* the first part is due to diffusion */

    /* first calculate the intermediate vectors*/
      for (int m = 0; m < 3; m++)	{
	  			A[m] = 0;
	  			C[m] = 0;
			}
      z_loc = 0;
      y_loc = 0;

      for (int m = 0; m < 3; m++)	{
	 			 A[0] = A[0] + J0M[0][m][k] * Hj[m];   /* this is J0M*Hj */
	 			 A[1] = A[1] + J0M[1][m][k] * Hj[m];
				 A[2] = A[2] + J0M[2][m][k] * Hj[m];
			}

      for (int m = 0;m < 3; m++)	{
			  C[0]=C[0]+J0M[0][m][k]*Hi[m];	/* this is J0M*Hi  */
			  C[1]=C[1]+J0M[1][m][k]*Hi[m];
	 			C[2]=C[2]+J0M[2][m][k]*Hi[m];
			}

      for (int m = 0;m < 3; m++)	{
			  z_loc+=b[m]*A[m];
			  y_loc+=b[m]*C[m];
			}

      q_loc+=z_loc*y_loc*(1/fabs(DETERMINANTJ[k]));    /* add up the diffusion part */
  }

  q_loc = q_loc*deltak;
  return q_loc;

}

/**
 * @brief TT function
 * TTfun is a function calculating the values of the
 * partial derivative of phi at the point P and putting them in a
 * (3,8) array
 * @param TT partial derivative of phi at point P (result of this function)
 * @param P pointer
 */
void SnowDriftA3D::TTfun(double TT[3][8],const double *P)	// P is a pointer pointing at an adress containing a double
{
  const double p_loc = P[0];
  const double q_loc = P[1];
  const double r = P[2];		//the three coordiantes of the point P

  TT[0][0] = 1/8.*(1-q_loc)*(1-r);
  TT[0][1] = 1/8.*(1+q_loc)*(1-r);
  TT[0][2] = -1/8.*(1+q_loc)*(1-r);
  TT[0][3] = -1/8.*(1-q_loc)*(1-r);
  TT[0][4] = 1/8.*(1-q_loc)*(1+r);
  TT[0][5] = 1/8.*(1+q_loc)*(1+r);
  TT[0][6] = -1/8.*(1+q_loc)*(1+r);
  TT[0][7] = -1/8.*(1-q_loc)*(1+r);

  TT[1][0] = -1/8.*(1+p_loc)*(1-r);
  TT[1][1] = 1/8.*(1+p_loc)*(1-r);
  TT[1][2] = 1/8.*(1-p_loc)*(1-r);
  TT[1][3] = -1/8.*(1-p_loc)*(1-r);
  TT[1][4] = -1/8.*(1+p_loc)*(1+r);
  TT[1][5] = 1/8.*(1+p_loc)*(1+r);
  TT[1][6] = 1/8.*(1-p_loc)*(1+r);
  TT[1][7] = -1/8.*(1-p_loc)*(1+r);

  TT[2][0] = -1/8.*(1+p_loc)*(1-q_loc);
  TT[2][1] = -1/8.*(1+p_loc)*(1+q_loc);
  TT[2][2] = -1/8.*(1-p_loc)*(1+q_loc);
  TT[2][3] = -1/8.*(1-p_loc)*(1-q_loc);
  TT[2][4] = 1/8.*(1+p_loc)*(1-q_loc);
  TT[2][5] = 1/8.*(1+p_loc)*(1+q_loc);
  TT[2][6] = 1/8.*(1-p_loc)*(1+q_loc);
  TT[2][7] = 1/8.*(1-p_loc)*(1-q_loc);

}

/**
 * @brief Phi - shape functions
 * Calculates the value of the 8 shape functions at point P
 * @param *PHI ...................shape functions
 * @param *P pointer
 */
void SnowDriftA3D::phi(double *PHI,double *P)
{
  double p_loc=P[0], q_loc=P[1], r=P[2];

  PHI[0] = 1/8.*(1+p_loc)*(1-q_loc)*(1-r);
  PHI[1] = 1/8.*(1+p_loc)*(1+q_loc)*(1-r);
  PHI[2] = 1/8.*(1-p_loc)*(1+q_loc)*(1-r);
  PHI[3] = 1/8.*(1-p_loc)*(1-q_loc)*(1-r);
  PHI[4] = 1/8.*(1+p_loc)*(1-q_loc)*(1+r);
  PHI[5] = 1/8.*(1+p_loc)*(1+q_loc)*(1+r);
  PHI[6] = 1/8.*(1-p_loc)*(1+q_loc)*(1+r);
  PHI[7] = 1/8.*(1-p_loc)*(1-q_loc)*(1+r);
}

/**
 * @brief matmult
 * computes a matrix vector product for a sparse matrix
 * @param res result (matrix vector product)
 * @param x
 * @param sm sparse matrix
 * @param ijm
 */
void SnowDriftA3D::matmult(CDoubleArray& res, const CDoubleArray& x_loc, double* sm, int* ijm)
{
  int n= ijm[1]-2;

  for (int i=1;i<=n;i++){
    res[i]=0;
	}

  for (int i = 1; i <= n; i++)	{
      res[i]=sm[i]*x_loc[i];
      for (int k = ijm[i]; k <= (ijm[i+1]-1); k++) {
				  //res[i]=ijm[5];//fake line
				  //res[i]=x[ijm[k]];//fake line
				  res[i] = res[i] + sm[k]*x_loc[ijm[k]];  //real line
			}
  }
}

/**
 * @brief matmult
 * computes a matrix vector product for a sparse matrix of CSR format
 * @param y result (matrix vector product)
 * @param x
 * @param sA
 * @param colind column index
 * @param rowPtr row index
 */
void SnowDriftA3D::matmult(CDoubleArray& y_loc,
			const CDoubleArray& x_loc,
			const CDoubleArray& sA_loc,
			const CIntArray& colInd,
			CIntArray& rowPtr )
{
	const size_t dim = rowPtr.getNx() - 1;

	for (size_t i = 0; i < dim; i++) {
		y_loc[i] = 0;
		for (int j = rowPtr[i]; j < rowPtr[i+1] ; j++) {
			y_loc[i] += sA_loc[ j ] * x_loc[ colInd[j] ];
		}
	}
}

/**
 * @brief transmult
 * this is the transmult function that multiplies the vector x by the transpose of the matrix M in its
 * sparse form sm, ijm
 * @param res result
 * @param x
 * @param sm
 * @param ijm index
 */
void SnowDriftA3D::transmult(CDoubleArray& res, const CDoubleArray& x_loc, double* sm, int* ijm)
{
  int n=ijm[1]-2;

  int j;
  for (int i = 1; i <= n; i++) {
      res[i] = sm[i]*x_loc[i];
  }

  for (int i = 1; i<= n;i++) {
      for (int k = ijm[i]; k <= (ijm[i+1]-1); k++) {
	  		j = ijm[k];
	  		res[j] = res[j] + sm[k] * x_loc[i];
			}
  }
}

/**
 * @brief bicgStab  iterative equation solver
 * iterative equation solver
 * Tests : Tested by the followin procedure: given a sparse matrix A
 * (CRS-format) characterized by colA and rowA and an arbitrary, or
 * rather: a few nontrivial examples of a vector x. For each x
 * matmult(y,x,sA,rowA,colA) and bicgStab(result,y,sA,colA,rowA,...)
 * have been computed and then verified that result=x
 * @param res result
 * @param rhs
 * @param sA
 * @param colA
 * @param rowA
 * @param nmax max number of iterations
 * @param tol tolerance
 */
void SnowDriftA3D::bicgStab(CDoubleArray& result,
			 CDoubleArray& rhs_loc,
			 const CDoubleArray& sA_loc,
			 const CIntArray& colA_loc,
			 CIntArray& rowA_loc,
			 const int nmax,
			 const double tol,
			 double& testres)
{

  //dimension of the system
  const size_t  n = rowA_loc.getNx() - 1;

  CDoubleArray res1;
  res1.resize( n);

  CDoubleArray r_0;
  r_0.resize( n );

  CDoubleArray r;
  r.resize( n );

  CDoubleArray p_loc;
  p_loc.resize( n );

  CDoubleArray v;
  v.resize( n );

  CDoubleArray aux1;
  aux1.resize( n );

  CDoubleArray auxhat;
  auxhat.resize( n );

  CDoubleArray phat;
  phat.resize( n );

  CDoubleArray aux2;
  aux2.resize( n );

  double rho_old = 1;
  double rho_new = 1;

  double omega = 1;
  double alpha = 1;
  double beta = 0;

  double residual,norm1,norm2;
  double res4,res5;
  double tmp_res=1;
  int iterations=0;

  CDoubleArray res_ok_state;


  srand(static_cast<unsigned>(time(NULL)));
  //intitialization
  for (size_t i=0;i<n;i++) {
      result[i]=0;//(rand()%RAND_MAX)/(1.0*RAND_MAX);
      v[i] = 0;
      p_loc[i] = 0;
  }

  matmult(res1,result,sA_loc,colA_loc,rowA_loc);  // multiply Bx and store it into the dummy res1

  norm1 = 0;
  norm2 = 0;
  for ( size_t i = 0; i < n; i++ )  {
      r_0[i] = rhs_loc[i]-res1[i];
	//cout<<"r_0="<<r_0[i]<<" rhs= "<<rhs[i]<<" res1="<<res1[i]<<"at i"<<i<<endl;

      r[i] = rhs_loc[i]-res1[i];
      norm1 += rhs_loc[i] * rhs_loc[i];
      norm2 += r_0[i] * r_0[i];
  }

  if ( rho_new == 0 ){
    printf("Iteration fails\n");  //in this case, can't take x=0 as initial guess
	}

  int k = 0;
  double mark = 0;  //as soon as mark==1 you can stop the iteration, good approximatin is attained

  residual = sqrt(norm2);

  if ( residual < tol) {		//stopping criterion
      mark=1;
      printf("-------> Lucky failure within %d steps with residual: %f\n", k, residual);

  }

  //main loop
  while ( (k<=nmax) && (mark==0) ) {

      rho_new = 0;
      for (size_t i=0;i<n;i++)	{
	  		rho_new += r_0[i]*r[i];
			}

      beta = rho_new / rho_old * alpha / omega;

      for (size_t i=0;i<n;i++) {
	 			p_loc[i] *= beta;
	  		p_loc[i] += (r[i] - beta * omega *v[i] );

	  		//precond
	  		phat[i] = p_loc[i]/precond[i];
			}

      matmult(v,phat,sA_loc,colA_loc,rowA_loc);		// put B*p into v

      res4 = 0;
      for ( size_t i = 0; i < n ; i++ )	{
	  		res4 += v[i] * r_0[i];
			}

      //HERE test print if v or res4 don't make sense anymore
      if (!(fabs(v[n-2])<1e20)){
				printf("-------> LH: v too large?\n");
			}
      if (!(fabs(res4)<1e40)){
				printf("-------> LH: res4 too large? res4=%f \n", res4);
			}
      alpha = rho_new / res4;

      for ( size_t i = 0; i < n; i++ ) {
	  		aux1[i] = r[i] - alpha * v[i];
	  		//precond
	  		auxhat[i] = aux1[i]/precond[i];
      }

      matmult(aux2,auxhat,sA_loc,colA_loc,rowA_loc);

      res4 = 0;
      res5 = 0;
      for ( size_t i = 0; i < n; i++ ) {
	 			 res4 += ( aux2[i] * aux1[i] );
	 			 res5 += ( aux2[i] * aux2[i] );
      }

      omega = res4 / res5;

      for ( size_t i = 0; i < n; i++ ) {
	  		result[i] += ( alpha * phat[i] + omega * auxhat[i] );
	 			r[i] = aux1[i] - omega * aux2[i];
      }

      matmult(res1,result,sA_loc,colA_loc,rowA_loc);

      norm1=0;
      norm2=0;
      for ( size_t i = 0; i < n; i++ ) {
	  		norm1 += ( rhs_loc[i] * rhs_loc[i] );
			norm2 +=((rhs_loc[i]-res1[i])*(rhs_loc[i]-res1[i])) ; //norm2 += pow( rhs[i] - res1[i], 2 );
      }

      residual = sqrt(norm2) / sqrt(norm1);

      if ( residual <= tol)	{		//stopping criteria!
	 			 mark=1;
	 			 printf("-------> Convergence within %d steps ", k);
   			 printf("with residual: %f \n", residual);
      }

      if ( residual >= 1.e4 && tmp_res>=1.) {		//stopping criteria!
	  			mark=1;
	  			printf("-------> Hopeless after %d steps ", k);
	  printf("with residual: %f, try again! \n ", residual);
      }

      if ( residual < tmp_res){
	  //copy this state temporarily
	  res_ok_state=result;
	  tmp_res=residual;
	  iterations=k;
      }

      rho_old = rho_new;
//       //LH-CHANGES BEGIN
// //       1) comment out the following line
// 	  cout<<" this is k: "<<k<< " and the residual: "<< residual<<endl;
// //       2) redirect output to logfile
//         FILE *logfile;
//  	string fname=string("../results/numerics-tests/residual.dat");
//         logfile = fopen(fname.c_str(), "a");
//         fprintf(logfile, "%d %.15lf\n", k, residual);
//         fclose (logfile);
//        //LH-CHANGES END
      k++;
  }
  if (mark==0) {
	printf("-------> No convergence within maximal number of iterations (residual = %f)\n",residual);
	if (residual > 1e7*tol && tmp_res<1.){
	    printf("Use a previous step with residual = %f, #iterations= %d)\n", tmp_res, iterations);
	    //copy previous result and residual
	    result=res_ok_state;
	    residual=tmp_res;
      }
  }
    testres=residual;

}//end of function
