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

#include <alpine3d/snowdrift/Cell.h>
#include <meteoio/MeteoIO.h>

#include <alpine3d/snowdrift/SnowDrift.h>

using namespace mio;

/**
 *@brief Assemble System
 * loop over all elements, updates the system matrices A, B and the 'bare' right hand side (rhs)
 * of the system, in other words: prepares the system prior to the inclusion of dirichlet
 * boundary conditions
 * @param colA column index to locate within sparse matrix
 * @param rowA row index
 * @param sA "system matrix"
 * @param sB "system matrix"
 * @param Psi vector for incorporating inhomogeneous Dirichlet BC
 * @param f source in diffusion equation
 * @param dt
*/
void SnowDriftA3D::assembleSystem( CIntArray& colA_loc,
				CIntArray& rowA_loc,
				CDoubleArray& sA_loc,
				CDoubleArray& sB_loc,
				CDoubleArray& Psi_loc,
				CDoubleArray& f_loc,
				const double dt)
{
  //-------------------------------------------------------------
  //variables for integrals over elements
  //-------------------------------------------------------------
  double DETERMINANTJ[8];   // to store the determinant of the
			    // isoparametric point transf. at the 8
			    // int. points
  double J0M[3][3][8];      // to store the eight J0 matrices at the 8
			    // int. points
  double J[3][3];	    // the Jacobian matrix
  double J0[3][3];	    // the J0 matrix
  //element variables
  double b[3];		// the wind field
  double K[3][3];	// the diffusion matrix
  //the element matrix
  double Ael[9][9];
  double Del[9][9];

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ELEMENT MATRICES AND ASSEMBLING
  // now run over all elements: first interior of domain, then faces, corners
  // and bars
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  double deltak; //the SUPG parameter
  int element;  //used as element counter

  //the array dofNode contains the local node indices the first
  //nDofNode entries are the degrees of freedom of that element and
  //the following nBoundaryNodes = 8-nDofNodes entries are the
  //boundary nodes.  the number of boundary and dof-nodes depends on
  //the location of the element (face,bar,corner) and is set by the Cell class
  int dofNode[8];
  int nDofNodes;
  int nBoundaryNodes;

  Cell cell;
  std::string cellType;

  //LH_BC
  //loop over all elements
	for ( unsigned int iz = 0; iz < nz-1;iz++)	{
		for ( unsigned int iy = 0; iy < ny-1;iy++)	{
			for ( unsigned int ix = 0; ix < nx-1;ix++)	{

				element = iz*(nx-1)*(ny-1)+iy*(nx-1) + ix;
	   			cellType = "interior";
	   			cell.classifyNodes(dofNode, &nDofNodes, &nBoundaryNodes, cellType,0);

					computeDriftVector(b,ix,iy,iz);

					computeDiffusionTensor(K,ix,iy,iz);

					computeElementParameters(element, DETERMINANTJ,J0M,J0,J,b,K,deltak,qualla,ix,iy,iz);

					computeElementSystem(element,nDofNodes,dofNode,Ael,Del,STATIONARY, DETERMINANTJ,J0M,b,K,deltak,dt,f_loc,Psi_loc);
   				addElementMatrix(sA_loc, colA_loc, rowA_loc, Ael,element,dofNode,nDofNodes);//LH
					for ( int i = 0; i < 8; i++)	{
		  				precond[ nodeMap[element][i] ] += Ael[i+1][i+1];
					}
  	      			addElementMatrix(sB_loc, colA_loc, rowA_loc, Del,element,dofNode,nDofNodes);//LH

	      			//boundary conditions
				int boundaryFaceNodes[6][4]=
					{
						{0, 1, 4, 5},//neg y direction in global system
						{2, 3, 6, 7},//pos y direction in global system
						{1, 2, 5, 6},//pos x direction in global system
						{0, 3, 4, 7},//neg x direction in global system
						{4, 5, 6, 7},//pos z direction in global system
						{0, 1, 2, 3} //neg z direction in global system
					};
				int isBFace[6];
				isBFace[0] = (iy == 0)   ?   1 : 0;
				isBFace[1] = (iy == ny-2)?  -1 : 0;
				isBFace[2] = (ix == nx-2)?   1 : 0;
				isBFace[3] = (ix == 0)   ?  -1 : 0;
				isBFace[4] = (iz == nz-2)?   1 : 0;
				isBFace[5] = (iz == 0)   ?  -1 : 0;

				double BCel[9][9];
				double surfaceMetric;
				double qp[3];
				double PHI[8];
				double rhsel[8];

				for ( int i = 0; i < 8; i++)	{
					rhsel[i] = 0;
					for ( int j = 0; j < 8; j++)	{
						BCel[i+1][j+1] = 0;
					}
				}

				//loop over element faces
				for ( unsigned int bf = 0; bf < 6; bf++) {
					//if b-face
					if ( isBFace[bf] != 0 )	{
						//loop over quadrature points
						//compute face averages
						double gc = 0;
						double gN = 0;
						double gD = 0;
						for ( unsigned int k = 0; k < 4; k++ ) {
							gc += 0.25*gamma[ nodeMap[element][ boundaryFaceNodes[bf][k] ] ];
							gN += 0.25*gNeumann[ nodeMap[element][ boundaryFaceNodes[bf][k] ] ];
							gD += 0.25*gDirichlet[ nodeMap[element][ boundaryFaceNodes[bf][k] ] ];
						}
		 				//compute surface integral
		 				//loop over quadrature points
						for ( unsigned int k = 0; k < 4; k++ ) {
							qp[0]= qPoint(0, boundaryFaceNodes[bf][k]);
							qp[1]= qPoint(1, boundaryFaceNodes[bf][k]);
							qp[2]= qPoint(2, boundaryFaceNodes[bf][k]);

							//coordinate direction
							int cDir = bf/2;
							qp[ cDir ] = isBFace[bf];
							phi(PHI,qp);
							Jacobian(DETERMINANTJ,J,element,qp,k,ix,iy,iz);
							J0fun(J0,J);

							surfaceMetric = 0;
			  				for ( unsigned int l = 0; l < 3; l++) {
							    surfaceMetric += J0[l][cDir]*J0[l][cDir];
			    				}
			  				surfaceMetric = sqrt(surfaceMetric);

			 				for ( unsigned int i = 0; i < 8; i++){

								rhsel[i] += surfaceMetric*(gc*gD-gN)*PHI[i];
			      					for ( unsigned int j = 0; j < 8; j++)	{
				  					BCel[i+1][j+1] += ( surfaceMetric * gc* PHI[i]*PHI[j]);
								}
			    				}
						}
		    			} //if b-face
				}//loop element faces

				for ( unsigned int i = 0; i < 8; i++)	{
					rhs[ nodeMap[element][i] ] += rhsel[i];
					precond[ nodeMap[element][i] ] += BCel[i+1][i+1];
				}

				addElementMatrix(sA_loc,colA_loc,rowA_loc,BCel,element,dofNode,nDofNodes);//LH
			}//end of ix
		}// end of iy
	}//end of iz
}

/**
 *@brief applyBoundaryValues
 * Apply boundary values to all elements
 * Author: Marc Ryser
 * @param var00 initial variable
 * @param Psi vector for incorporating inhomogeneous Dirichlet BC
*/
void SnowDriftA3D::applyBoundaryValues(CDoubleArray& var00,
				    CDoubleArray& Psi_loc)
{
  //variables for integrals over elements
  double DETERMINANTJ[8];   // to store the determinant of the
			    // isoparametric point transf. at the 8
			    // int. points
  double J0M[3][3][8];      // to store the eight J0 matrices at the 8
			    // int. points
  double J[3][3];	    // the Jacobian matrix
  double J0[3][3];	    // the J0 matrix
  //element variables
  double b[3];		// the wind field
  double K[3][3];	// the diffusion matrix
  double deltak; //the SUPG parameter
  int element;  //used as element counter

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Loop over all boundary elements
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //the array dofNode contains the local node indices the first
  //nDofNode entries are the degrees of freedom of that element and
  //the following nBoundaryNodes = 8-nDofNodes entries are the
  //boundary nodes.  the number of boundary and dof-nodes depends on
  //the location of the element (face,bar,corner) and is set by the Cell class
  int dofNode[8];
  int nDofNodes;
  int nBoundaryNodes;

  Cell cell;
  std::string cellType;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Step 1: ELEMENTS AND ELEMENT MATRICES ON FACES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  size_t nz_face_size, tmp;
  nz_face.size(nz_face_size, tmp);
  //for (int face = 0; face < nz_face.size(); face++)
  for (size_t face = 0; face < nz_face_size; face++) {
	for (int layer = nz_face[face][0]; layer<=nz_face[face][1]; layer++ ){
	    for (int depth = ny_face[face][0]; depth<=ny_face[face][1]; depth++ )	{
		for (int width = nx_face[face][0]; width<=nx_face[face][1]; width++ )	{
			element = layer*(nx-1)*(ny-1)+depth*(nx-1)+width;
			computeDriftVector(b,width,depth,layer);
			computeDiffusionTensor(K,width,depth,layer);
			cellType = "face";
			cell.classifyNodes(dofNode,&nDofNodes,&nBoundaryNodes,cellType,(int)face);
			computeElementParameters(element, DETERMINANTJ,J0M,J0,J,b,K,deltak,qualla,width,depth,layer);

			//DirichletBoundaryValues hasn't been implemented
			computeDirichletBoundaryValues( element,DETERMINANTJ,J0M,J0,b,K,deltak,dofNode,nDofNodes,nBoundaryNodes,var00,Psi_loc);
		} // end of width loop
	    }// end of depth loop
	}// end of layer loop
  }// end of face loop


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Step 3: ELEMENTS AND ELEMENT MATRICES ON BOUNDARY CORNERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for (unsigned int corner = 0; corner < n_corner.getNx() ;corner++ ) {
      element = n_corner[corner];
      cellType = "corner";
      cell.classifyNodes(dofNode,&nDofNodes,&nBoundaryNodes,cellType,corner);
      unsigned int ix=0;
      unsigned int iy=0;
      unsigned int iz=0;

      if (corner<=3){
	    iz=0;
	    if (corner==0 || corner==3){
		ix=0;
	    }else{
		ix=nx-2;
	    }
	    if (corner==2 || corner==3){
		iy=ny-2;
	    } else{
		iy=0;
	    }
      }else if (corner<=7){
	    iz=nz-1;
	    if (corner==4 || corner==7){
		ix=0;
	    }else{
		ix=nx-2;
	    }
	    if (corner==6 || corner==7){
		iy=ny-2;
	    } else{
		iy=0;
	    }
      }

      computeDriftVector(b, ix, iy, iz);

      computeDiffusionTensor(K, ix, iy, iz);

      computeElementParameters(element,DETERMINANTJ,J0M,J0,J,b,K,deltak,qualla,ix,iy,iz);

      //DirichletBoundaryValues hasn't been implemented
      computeDirichletBoundaryValues( element,DETERMINANTJ,J0M,J0,b,K,deltak,dofNode,nDofNodes,nBoundaryNodes,var00,Psi_loc);
  }//end corner loop

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Step 4: ELEMENTS AND ELEMENT MATRICES ON BARS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  size_t nz_bar_size, tmp2;
  nz_bar.size(nz_bar_size, tmp2);
  //for (int bar = 0; bar < nz_bar.size(); bar++ )
    for (size_t bar = 0; bar < nz_bar_size; bar++ ){
	for (int layer=nz_bar[bar][0];layer<=nz_bar[bar][1];layer++){
	    for (int depth=ny_bar[bar][0];depth<=ny_bar[bar][1];depth++) {
		for (int width=nx_bar[bar][0];width<=nx_bar[bar][1];width++){

			element = layer*(nx-1)*(ny-1) + depth*(nx-1) + width;
			computeDriftVector(b, width, depth, layer);
			computeDiffusionTensor(K,width, depth, layer);
			cellType = "bar";
			cell.classifyNodes(dofNode,&nDofNodes,&nBoundaryNodes,cellType,(int)bar);

			computeElementParameters(element, DETERMINANTJ,J0M,J0,J,b,K,deltak,qualla,width,depth,layer);
			//DirichletBoundaryValues hasn't been implemented >>remove?
			computeDirichletBoundaryValues( element,DETERMINANTJ,J0M,J0,b,K,deltak,dofNode,nDofNodes,nBoundaryNodes,var00,Psi_loc);
		}// end width loop
	    }// end depth loop
	}// end layer loop
    } // end bar loop
}

/**
 * @brief Compute diffusion tensor
 * Computes the diffusion tensor of an element, accesses the nodes array and the nodeMap array
 * @param element element number
 * @param K diffusion coefficient
 * @param layer layer of the element
*/
void SnowDriftA3D::computeDiffusionTensor(double K[3][3], const unsigned int ix, const unsigned int iy, const unsigned int iz)
{
 double Length_Scale;
 double K_model,z_loc,K_min,du;
 double K_max=1.6;

 int version=1; //1: old version with mixing length, 2: take K from ARPS, 3: version where K depends on standard deviation of wind speed within a layer, 4: K is calculated after the windfield is read (see SnowDrift.cc),calulation was moved because it is time consuming and is now calculated only once.
 if (READK){
    version=2;
 }

 for (int i = 0; i < 3; i++ ) {

	for (int j = 0; j < 3; j++ )	{
	    K[i][j] = 0;
	}
 }

 K_model = 0.; K_min = 0.;

    if (version==1){
	//loop over 4 nodes
	unsigned int ii, jj;
	for ( ii=ix; ii<=ix+1; ii++){
	    for (jj = iy; jj<=iy+1; jj++) {

		Length_Scale = fabs(nodes_z.grid3D(ii,jj,iz)-nodes_z.grid3D(ii,jj,iz+1));

		//vertical height
		z_loc = (nodes_z.grid3D(ii,jj,iz) - nodes_z.grid3D(ii,jj,0)) + 0.5*Length_Scale;

		/* Approximation of ABL neutral: 1/8*dudz*SQR(1/((1/k(z+z0))+(1/LengthScale))) */
		du = sqrt( Optim::pow2( nodes_u.grid3D(ii,jj,iz+1) - nodes_u.grid3D(ii,jj,iz))
			    + Optim::pow2( nodes_v.grid3D(ii,jj,iz+1) - nodes_v.grid3D(ii,jj,iz) )
			    + Optim::pow2( nodes_w.grid3D(ii,jj,iz+1) - nodes_w.grid3D(ii,jj,iz) )
			    );

		K_model += 0.25*du/Length_Scale*(Optim::pow2(1./((1./(0.4*z_loc))+(1./Length_Scale))));

		/* Approximation of k*ustar*z for the surface layer */
		K_min += 0.25*0.4*z_loc*0.1;
	    }
	}
	    K_min = std::max(0.1,K_min);
    }else if (version==2){
	    //use the diffusivity as given in ARPS >> will only work if you use wind fields of ARPS with a somewhat longer integration step (thus incl. turbulence). Some tests with K of mean //flow fields > solver couldn't find steady state of c
	    double temp_K=0.;
	    unsigned int ii,jj;
	    for ( ii=ix; ii<=ix+1; ii++){
		for (jj = iy; jj<=iy+1; jj++) {
		    temp_K=temp_K+nodes_K.grid3D(ii,jj,iz)+nodes_K.grid3D(ii,jj,iz+1);

		}
	    }
	    //now we have added all possible K of the 8 nodes around element, take the average:
	    K_model=temp_K/8.;
	    K_min=0.1;


    } 
    //make sure K_model is in range min-max
    K_model=std::max(K_min,K_model);
    K_model=std::min(K_max,K_model);

    K[2][2] = K_model;
    K[1][1] = K_model;
    K[0][0] = K_model;
    //extra test
    K[0][0]=std::min(std::max(K[0][0],K_min),50.*K_model);
    K[1][1]=std::min(std::max(K[1][1],K_min),50.*K_model);
    K[2][2]=std::min(std::max(K[2][2],K_min),10.*K_model);//ori 10
}

/**
 * @brief Compute the drift vector of an element
 * @param element element number
 * @param b drift vector
 */
void SnowDriftA3D::computeDriftVector(double b[3], const unsigned int ix, const unsigned int iy, const unsigned int iz)
{
  //reset
  for (int i = 0; i < 3; i++) {
      b[i]=0.;
  }

  unsigned int ii, jj, kk;
    for ( ii=ix; ii<=ix+1; ii++){
	for (jj = iy; jj<=iy+1; jj++) {
	    for (kk = iz; kk<=iz+1; kk++) {
		b[0]+=1/8.*nodes_u.grid3D(ii,jj,kk);
		b[1]+=1/8.*nodes_v.grid3D(ii,jj,kk);
		b[2]+=1/8.*(nodes_w.grid3D(ii,jj,kk)-nodes_wstar.grid3D(ii,jj,kk));
	    }
	}
    }
}

/**
 * @brief Compute parameters of the SUPG method
 * Compute the parameters of the SUPG method for each element
 * @param element Element number
 * @param DETERMINANTJ Determinant of Jacobian matrix
 * @param J0M
 * @param J0 JO matrix
 * @param J Jacobian matrix
 * @param b drift vector
 * @param K Diffusion tensor
 * @param &deltak parameter SUPG method
 * @param &qualla parameter SUPG method (to vary delta_k)
 */
void SnowDriftA3D::computeElementParameters(const int& element,
					 double DETERMINANTJ[8],
					 double J0M[3][3][8],
					 double J0[3][3],			// the J0 matrix
					 double J[3][3],			// the Jacobian matrix
					 double b[3],
					 double K[3][3],
					 double &deltak,
					 double &qualla_loc,
					 const int ix,
					 const int iy,
					 const int iz)
{
  //cad is a dummy vector the integration point for the Gaussian
  //quadrature is assigned to
  double  cad[3];
  double epsilon;
  double hk;      // the diameter (=longest side) of the respective
		  // element
  double Pe;      // the local Peclet number

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //element infos
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  epsilon = 1/3.* (K[0][0]+K[1][1]+K[2][2]);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //hk
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hk=fabs(nodes_x.grid3D(ix,iy,iz+1)-nodes_x.grid3D(ix+1,iy,iz+1));

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // deltak
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Pe = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]) * hk/12/epsilon; /*calculate the local Peclet number*/
  if (Pe<1){
      deltak=hk*hk/12/epsilon*qualla_loc; 	//qualla is parameter to vary delta_k, standard:1
  }else{
      deltak=hk/2/sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])*qualla_loc;
  }

  for (int i = 0 ; i < 8; i++){  //loop over all 8 points
      cad[0] = qPoint(0, i);
      cad[1] = qPoint(1, i);
      cad[2] = qPoint(2, i);

      // this function will calculate the Jacobian of the
      // isoparametric point transformation on the
      // respective element and stores its determinant in
      // DETERMINANTJ[i] for each of the integration
      // points and store the Jacobian J at the respective
      // point in J; //note that you'll only use J once
      // more for J0
      Jacobian(DETERMINANTJ,J,element,cad,i,ix,iy,iz);


      // this function will calculate the J0 matrix at the
      // integration point and store it in the variable
      // J0*/
      J0fun(J0,J);

      for (int k = 0; k < 3; k++ ){
	  for (int j = 0; j < 3; j++){
	      J0M[k][j][i] = J0[k][j];
	  }
       }
  }
}


/**
 * @brief compute element system
 * ....
 * @param element element number
 * @param nDofNodes degrees of freedom of element
 * @param Ael
 * @param Del
 * @param stationary
 * @param DETERMINANTJ
 * @param J0M
 * @param b wind field
 * @param K diffusion tensor
 * @param deltak parameter
 * @param dt
 * @param f source
 * @param Psi
*/
void SnowDriftA3D::computeElementSystem(int &element,
				     int &nDofNodes,
				     int* dofNode,
				     double Ael[9][9],
				     double Del[9][9],
				     bool stationary,
				     double DETERMINANTJ[8],
				     double J0M[3][3][8],
				     double b[3],
				     double K[3][3],
				     double &deltak,
				     const double &dt,
				     CDoubleArray& f_loc,
				     CDoubleArray& Psi_loc)
{
  //element matrices
  double Bel[9][9];
  double Cel[9][9];
  double Apdxel[9][9];
  double Adxdxel[9][9];

  //reset element matrices
  for (int i = 0; i<9; i++ ){
      for (int j = 0; j<9; j++ ){
	Ael[i][j] = 0; //LH
	Bel[i][j] = 0;
	Cel[i][j] = 0;
	Apdxel[i][j] = 0;
	Adxdxel[i][j] = 0;
      }
  }

  //------------------------------------
  // calculation of the element matrices
  //------------------------------------
  for (int i = 0; i < nDofNodes ; i++ )	{
    double auxSum = 0;
    for (int j = 0; j < nDofNodes ; j++)	{
	//element matrices are indexed from 1..8
	Bel[ 1+dofNode[i] ][ 1+dofNode[i] ] = GQIntB(DETERMINANTJ,dofNode[i],dofNode[j]);
	Cel[ 1+dofNode[i] ][ 1+dofNode[j] ] = GQIntC(DETERMINANTJ,J0M,dofNode[i],dofNode[j],b,K);
	Adxdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] = GQIntAdxdx(DETERMINANTJ,J0M,dofNode[i],dofNode[j],b,deltak);
	Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] = GQIntApdx(DETERMINANTJ,J0M,dofNode[i],dofNode[j],b,deltak);

	//assemble matrix according to crank nicolson
	if ( stationary ) {			//LH
		Ael[ 1+dofNode[i] ][ 1+dofNode[j] ] =
		( Adxdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Cel[ 1+dofNode[i] ][ 1+dofNode[j] ] );
	}else{
		Ael[ 1+dofNode[i] ][ 1+dofNode[j] ] =
		( Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Bel[ 1+dofNode[i] ][ 1+dofNode[j] ] ) / dt
		+ ( Adxdxel[ 1+dofNode[i] ][ 1+dofNode[j]] + Cel[ 1+dofNode[i] ][ 1+dofNode[j] ] ) * 0.5;
	}

	//assemble right hand matrix according to crank nicolson
	if ( stationary ){//LH
		Del[ 1+dofNode[i] ][ 1+dofNode[j] ] = 0;
	} else {
	      	Del[ 1+dofNode[i] ][ 1+dofNode[j] ] =
		( Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Bel[ 1+dofNode[i] ][ 1+dofNode[j] ] ) / dt
		- ( Adxdxel[ 1+dofNode[i] ][ 1+dofNode[j]] + Cel[ 1+dofNode[i] ][ 1+dofNode[j] ] ) * 0.5;
	}

//			assemble right hand side
// 	  	if ( stationary )//LH
// 	    {
	 auxSum += ( Bel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] )
	  	* f_loc[ nodeMap[element][ dofNode[j] ] ] ;

// 	    }
// 	  	else
// 	    {
// 	      auxSum += ( ( ( Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Bel[ 1+dofNode[i] ][ 1+dofNode[j] ] )/dt
// 			    - ( Cel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Adxdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] )*0.5
// 			    ) * c[ nodeMap[element][ dofNode[j] ] ]
// 			  + ( Bel[ 1+dofNode[i] ][ 1+dofNode[j] ] + Apdxel[ 1+dofNode[i] ][ 1+dofNode[j] ] )
// 			  * f[ nodeMap[element][ dofNode[j] ] ]
// 			  );
// 		  }
    }
    Psi_loc[ nodeMap[element][ dofNode[i] ] ] += auxSum;
 }
}

/**
 * @brief Add element matrix
 * given a sparse matrix A in CSR format ( specified by
 * sA, colInd, rowPtr) and the element matrix Bel, this function adds
 * the element contributions to the global ones ; the first
 * length_spec entries of the vector spec contains the indices of the
 * matrix which are inserted
 * Comments : Beware of the enumeration of the matrix sA (size is one
 * bigger than required, zero index is unused => dangerous
 * accesses nodeMap, nodeMap
 * @param sA
 * @param colInd
 * @param rowPtr
 * @param Bel
 * @param element
 * @param spec
 * @param length_spec
*/
void SnowDriftA3D::addElementMatrix( CDoubleArray& sA_loc,
				  const CIntArray& colInd,
				  const CIntArray& rowPtr,
				  const double Bel[9][9],
				  const int element,
				  const int* spec,
				  const int length_spec)
{
  //global indices
  int I;
  int J;
  for ( int i = 0; i < length_spec; i++)	{
      I = nodeMap[element][ spec[i] ];

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //insert entry (i,i)
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //binary search of I's column index
      int first = rowPtr[I];
      int last = rowPtr[I+1] - 1;
      int mid;
      while ( first <= last ){
	    mid = (first + last) / 2;  // compute mid point.
	    if ( I > colInd[mid] ){
		first = mid + 1;  // repeat search in top half.
	    }else if ( I < colInd[mid] ) {
		last = mid - 1; // repeat search in bottom half.
	    }else {
		sA_loc[ mid ] += Bel[ spec[i]+1 ][ spec[i]+1 ]; // got it
	    break;
	    }
      }

      for ( int j = i + 1; j < length_spec; j++)	{
	    J = nodeMap[element][ spec[j] ];
	    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    //insert entry (i,j)
	    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	    //binary search of J's column index
	    first = rowPtr[I];
	    last = rowPtr[I+1] - 1;
	    while (first <= last){
		mid = (first + last) / 2;  // compute mid point.
		if ( J > colInd[mid] ){
			first = mid + 1;  // repeat search in top half.
		}else if ( J < colInd[mid] ) {
			last = mid - 1; // repeat search in bottom half.
		}else{
			sA_loc[ mid ] += Bel[ spec[i]+1 ][ spec[j]+1 ]; // got it
			break;
		}
	    }

	    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    //insert entry (j,i)
	    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    //binary search of I's column index
	    first = rowPtr[J];
	    last = rowPtr[J+1] - 1;
	    while ( first <= last ) {
		    mid = (first + last) / 2;  // compute mid point.
		    if ( I > colInd[mid] ) {
			    first = mid + 1;  // repeat search in top half.
		    }else if ( I < colInd[mid] ) {
			    last = mid - 1; // repeat search in bottom half.
		    }else {
			    sA_loc[ mid ] += Bel[ spec[j]+1 ][ spec[i]+1 ]; // got it
			    break;
		    }
	    }
	}//i
    }//j
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function    : computeDirichletBoundaryValues
// Authors     : Marc Ryser
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SnowDriftA3D::computeDirichletBoundaryValues(int element,
					       double DETERMINANTJ[8],
					       double J0M[3][3][8],
					       double J0[3][3],			// the J0 matrix
					       double b[3],
					       double K[3][3],
					       double deltak,
					       int node[8],
					       int nDofNodes,
					       int nBoundaryNodes,
					       CDoubleArray& c00_loc,
					       CDoubleArray& Psi_loc )
{
	//since this method is currently empty, we just remove the warnings...
	//TODO: implement this method or get rid of it
	(void)element;
	(void)DETERMINANTJ;
	(void)J0M;
	(void)J0;
	(void)b;
	(void)K;
	(void)deltak;
	(void)node;
	(void)nDofNodes;
	(void)nBoundaryNodes;
	(void)c00_loc;
	(void)Psi_loc;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // incorporation of inhomogeneous boundary conditions
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//   for (int i = 0; i < nDofNodes; i++ )
//     {
//       for (int j = nDofNodes; j < nDofNodes + nBoundaryNodes; j++)
// 	{
// 	  Psi[ nodeMap[ element][ node[i] ] ] +=
// 	    c00[ nodeMap[ element][ node[j] ] ]
// 	    * GQIntC( DETERMINANTJ,J0M, node[i], node[j],b,K);

// 	  Psi[ nodeMap[element][ node[i] ] ] +=
// 	    c00[ nodeMap[element][ node[j] ] ]
// 	    * GQIntAdxdx(DETERMINANTJ,J0M, node[i], node[j],b,deltak);

// 	  rhs[ nodeMap[element][ node[j] ] ] = c00[ nodeMap[element][ node[j] ] ];
// // 	  printf("c00[%d]=%lf\n", nodeMap[element][ node[j]], c00[ nodeMap[element][ node[j] ] ] );
// 	  zeroRow( nodeMap[element][ node[j] ] );
// 	}
//     }

  //if bottom
//   if ( element <= (nx-1)*(ny-1) )
//     {
      //enumeration from 1..8
//       double BCel[9][9];
//       double surfaceMetric;
//       double qp[3];
//       double PHI[8];

//       for ( int i = 0; i < 8; i++)
// 	{
// 	  for ( int j = 0; j < 8; j++)
// 	    {
// 	      BCel[i+1][j+1] = 0;

// 	      //determine element normal to bottom
// 	      for ( int k = 0; k < 4; k++ )
// 		{
// 		}

// 	      //loop over quadrature points
// 	      for ( int k = 0; k < 4; k++ )
// 		{
// 		  qp[0]= qPoint[0][k];
// 		  qp[1]= qPoint[1][k];
// 		  qp[2]= -1;

// 		  phi(PHI,qp);
// 		  Jacobian(DETERMINANTJ,J,element,qp,k);
// 		  J0fun(J0,J);
// 		  surfaceMetric = sqrt( J0[1][3]*J0[1][3]
// 					+J0[2][3]*J0[2][3]
// 					+J0[3][3]*J0[3][3] );





// 		  BCel[i+1][j+1] += ( surfaceMetric
// 				      * gamma
// 				      * b[2]
// 				      * (PHI[i]-gDirichlet[ nodeMap[element][k]] )
// 				      + gNeumann[ nodeMap[element][k] ]
// 				      ) * PHI[j];

// 		}
// 	    }
// 	}
//       addElementMatrix(sA,colA,rowA,BCel,element,dofNode,nDofNodes+nBoundaryNodes);//LH
// //     }



}



/**
 * @brief compute deposition flux
 * Authors     : Marc Ryser
 * Description : calculates the deposition flux i.e. (b c + K grad c) at the
 * artificial layer of nodes and write it to the flux_x,y,z variable
 * Comments : Beware of the enumeration of elements and nodes, Must be
 * consistent with SnowDriftA3D::SnowMassChange(...)  Note: this solution
 * has temporary character since flux_x has size (nx-1)*(ny-1) and
 * therefore only internal fluxes are written
 * @param CDoubleArray&c Concentration of snow
 * @param theta height above bottom of element
*/
void SnowDriftA3D::computeDepositionFlux(const CDoubleArray& concentration,
				      const double theta_l)
{
  double DETERMINANTJ[8];
  double J0[3][3];
  double J[3][3];
  double TT[3][8];
//   double b[3];
  double K[3][3];

  //cad is a dummy vector the integration point for the Gaussian
  //quadrature is assigned to
  double  cad[3];

  double hivec[3];
  double gradc[3];

  int element;

  // loop over all the second layer
  // elements having the bottom node 0
  // in the interior of the domain
  for (unsigned int k = 1; k < ny-1; k++ ){
      for (unsigned int l = 1; l < nx-1; l++){
	  // the respective element of the second layer
	  element = (nx-1)*(ny-1) + k*(nx-1) + l;
// 	  int element0 = k*(nx-1) + l;
	  int node = nodeMap[element][0]-nx*ny;
	  //reset fluxes
	  flux_x[node] = 0;
	  flux_y[node] = 0;
	  flux_z[node] = 0;

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // advective contribution
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	flux_x[node] += concentration[nodeMap[element][0]] * nodes_u.grid3D(l,k,1);
	flux_y[node] += concentration[nodeMap[element][0]] * nodes_v.grid3D(l,k,1);
	flux_z[node] += concentration[nodeMap[element][0]] * ( nodes_w.grid3D(l,k,1)-nodes_wstar.grid3D(l,k,1));

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // diffusive contribution
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  computeDiffusionTensor(K, l,k,1); //take layer consistent with element

	  // This Is The Vector Containing The Point Where We Want To
	  // Calculate K gradc dot n theta is the height above bottom of
	  // the element, very small if possible to get as close as
	  // possible to the point
	  cad[0] = 1 - theta_l;
	  cad[1] = -1 + theta_l;
	  cad[2] = -1 + theta_l;

	  // calculation of J at the point cad and storing the determinant in the
	  // 0th (=kth component) component of DETERMINANTJ
	  Jacobian( DETERMINANTJ, J, element, cad, 0 ,l,k,1);

	  // calculation of J0 at the point cad
	  J0fun(J0,J);

	  //calculation of the gradient of all 8 shape functions at the point cad in the reference element
	  TTfun(TT,cad);

	  for (int j = 0; j < 3; j++){ // the 3 components of gradc
	    gradc[j] = 0;
	    for (int m = 0; m < 8; m++){ //add over all 8 shape functions
		hivec[j] = 0;// need to clear hivec[j] for the matrix product of J0 with gradphi(m)
		for (int n = 0; n < 3; n++){ // dummy loop for matrix product
		    hivec[j] += J0[j][n] * TT[n][m];//store  J0 * gradphi(m)
		}
		gradc[j] += concentration[nodeMap[element][m]] * (1/DETERMINANTJ[0]) * hivec[j];//gradc=sum_m
		//(c(m)*1/detJ
		//*
		//J0*gradphi(m))
	    }
	  }
	  // k=1, l=2, and element given above implies
	  // nodeMap[element][0]=1, ie flux_x is filled up from the index 1
	  // k=ny-2, l=nx-1, and element given above implies
	  // nodeMap[element][0]=3721=(nx-2)*(ny-2)
	  // note: flux_x,y,z have size (nx-1)*(ny-1) such that the
	  // remaining memory is unaltered (better: resize flux_x)
	  flux_x[ node ] -= K[0][0] * gradc[0];
	  flux_y[ node ] -= K[1][1] * gradc[1];
	  flux_z[ node ] -= K[2][2] * gradc[2];
	} //end of kicking all elements on 1st interior element layer in l!!!
    } //dito in k!!
}

/**
 * @brief compute deposition flux
 * Authors     : Marc Ryser
 * Description : calculates the deposition flux i.e. (b c + K grad c) at the
 * artificial layer of nodes and write it to the flux_x,y,z variable
 * Comments : Beware of the enumeration of elements and nodes, Must be
 * consistent with SnowDriftA3D::SnowMassChange(...)  Note: this solution
 * has temporary character since flux_x has size (nx-1)*(ny-1) and
 * therefore only internal fluxes are written
 * @param CDoubleArray&c Concentration of snow
 * @param theta height above bottom of element
*/
void SnowDriftA3D::computeDepositionFluxSublimation(const CDoubleArray& concentration,
				      const double theta_l)
{
  double DETERMINANTJ[8];
  double J0[3][3];
  double J[3][3];
  double TT[3][8];
//   double b[3];
  double K[3][3];

  //cad is a dummy vector the integration point for the Gaussian
  //quadrature is assigned to
  double  cad[3];

  double hivec[3];
  double gradc[3];

  int element;

  // loop over all the second layer
  // elements having the bottom node 0
  // in the interior of the domain
  for (unsigned int k = 1; k < ny-1; k++ ){
      for (unsigned int l = 1; l < nx-1; l++){
	  // the respective element of the second layer
	  element = (nx-1)*(ny-1) + k*(nx-1) + l;
// 	  int element0 = k*(nx-1) + l;
	  int node = nodeMap[element][0]-nx*ny;
	  //reset fluxes
	  flux_x_subl[node] = 0;
	  flux_y_subl[node] = 0;
	  flux_z_subl[node] = 0;

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // advective contribution
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  flux_x_subl[node] += concentration[nodeMap[element][0]] * nodes_u.grid3D(l,k,1);
	  flux_y_subl[node] += concentration[nodeMap[element][0]] * nodes_v.grid3D(l,k,1);
	  flux_z_subl[node] += concentration[nodeMap[element][0]] * ( nodes_w.grid3D(l,k,1)-nodes_wstar.grid3D(l,k,1));

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  // diffusive contribution
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  computeDiffusionTensor(K, l,k,1); //take layer consistent with element

	  // This Is The Vector Containing The Point Where We Want To
	  // Calculate K gradc dot n theta is the height above bottom of
	  // the element, very small if possible to get as close as
	  // possible to the point
	  cad[0] = 1 - theta_l;
	  cad[1] = -1 + theta_l;
	  cad[2] = -1 + theta_l;

	  // calculation of J at the point cad and storing the determinant in the
	  // 0th (=kth component) component of DETERMINANTJ
	  Jacobian( DETERMINANTJ, J, element, cad, 0 ,l,k,1);

	  // calculation of J0 at the point cad
	  J0fun(J0,J);

	  //calculation of the gradient of all 8 shape functions at the point cad in the reference element
	  TTfun(TT,cad);

	  for (int j = 0; j < 3; j++){ // the 3 components of gradc
	    gradc[j] = 0;
	    for (int m = 0; m < 8; m++){ //add over all 8 shape functions
		hivec[j] = 0;// need to clear hivec[j] for the matrix product of J0 with gradphi(m)
		for (int n = 0; n < 3; n++){ // dummy loop for matrix product
		    hivec[j] += J0[j][n] * TT[n][m];//store  J0 * gradphi(m)
		}
		gradc[j] += concentration[nodeMap[element][m]] * (1/DETERMINANTJ[0]) * hivec[j];//gradc=sum_m
		//(c(m)*1/detJ
		//*
		//J0*gradphi(m))
	    }
	  }
	  // k=1, l=2, and element given above implies
	  // nodeMap[element][0]=1, ie flux_x is filled up from the index 1
	  // k=ny-2, l=nx-1, and element given above implies
	  // nodeMap[element][0]=3721=(nx-2)*(ny-2)
	  // note: flux_x,y,z have size (nx-1)*(ny-1) such that the
	  // remaining memory is unaltered (better: resize flux_x)
	  flux_x_subl[ node ] -= K[0][0] * gradc[0];
	  flux_y_subl[ node ] -= K[1][1] * gradc[1];
	  flux_z_subl[ node ] -= K[2][2] * gradc[2];
	} //end of kicking all elements on 1st interior element layer in l!!!
    } //dito in k!!
}

void SnowDriftA3D::zeroRow(int node)
{
  //version 3
  int I = node;
  int nnz = rowA[I+1]-rowA[I];

  //loop over all indices in row I
  for (int k = 0; k < nnz; k++){
      int J = colA[ rowA[I] + k ];
      //set diagonal element to 1
      if ( I == J ){
	  sA[rowA[I] + k] = 1;
      }else{
	  //zero off diagonal elements
	  //IJ
	  sA[rowA[I] + k] = 0;

	  //JI
	  //binary search of I's column index
	  int first = rowA[J];
	  int last = rowA[J+1] - 1;
	  int mid;
	  while (first <= last) {
	      mid = (first + last) / 2;  // compute mid point.
	      if ( I > colA[mid] ) {
		  first = mid + 1;  // repeat search in top half.
	      } else if ( I < colA[mid] ) {
		  last = mid - 1; // repeat search in bottom half.
	      }else{
		  sA[ mid ] = 0; // got it
		  break;
	      }
	  }
      }
 }
}
