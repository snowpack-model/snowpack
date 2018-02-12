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

#include <meteoio/MeteoIO.h>

#include <alpine3d/snowdrift/SnowDrift.h>

using namespace mio;

/**
 * @brief Initializes system size and all data members which are required for the finite element solution
 */
void SnowDriftA3D::InitializeFEData( void )
{
	nNodes = nx * ny * nz;
	nElements = (nx-1) * (ny-1) * (nz-1);
	nDOF = (nx-2) * (ny-2) * (nz-2);

	//LH_BC
	nDOF=nNodes;

	nNZ = numberOfNonzeros();

	nz_face.resize(6,2);
	ny_face.resize(6,2);
	nx_face.resize(6,2);

	nz_bar.resize(12,2);
	ny_bar.resize(12,2);
	nx_bar.resize(12,2);

	nz_interior.resize(1,2);
	ny_interior.resize(1,2);
	nx_interior.resize(1,2);

	n_corner.resize(8);

	qPoint.resize(3,8);

	classifySubdomain();
	setQuadraturePoints();

	//set some parameter
	//parameter to play with
	qualla = 1.0;   // paramater to vary the deltak coefficient of
			// the SUPG approach

	theta = 1.0; // here you can vary the heigth above the point used in
		// the diffusional part of the computeDepositionFlux function

	sA.resize( nNZ );
	sB.resize( nNZ );
	colA.resize( nNZ );
	rowA.resize( nDOF + 1 );
	Psi.resize( nDOF );
	c.resize( nDOF );
	q.resize( nDOF );
	T.resize( nDOF );
	T00.resize( nDOF );

	rhs.resize( nDOF );
	precond.resize( nDOF );
	f.resize( nDOF );
	c00.resize( nNodes );

	q00.resize( nNodes );
	nnzA.resize( nDOF );
	adjA.resize( nNZ );

	//LH_BC
	gDirichlet.resize( nDOF );
	gNeumann.resize( nDOF );
	gamma.resize( nDOF );

}


/**
 * @brief Description...
 * this functions sets summation limits for the loops
 * over boundary elements (faces,bars,corners), presently the assemble
 * routines depend on these limits => fragile
 * the values of the limits is also affected by the
 * enumeration of the elements presently the the first
 * element has the number 1 (instead of 0) => even more
 * fragile

 * This setting will become important for the domain
 * decomposition, since then an element on the domain
 * boundary is not necessarily an element of the
 * physical boundary, and these array account for that.
 * Note: these summation limits are related to the
 * number of nodes, the number of elements, the number
 * of dofs and accordingly the number of nonzero matrix
 * elemens in the system matrix. It would be better to
 * have one class which holds all these information
 * consistently.
 */
void SnowDriftA3D::classifySubdomain( void )
{
	// determine the summation limits for layer, depth and width
	nz_face(0,0)=1;      nz_face(0,1)=(nz-3);
	nz_face(1,0)=1;      nz_face(1,1)=(nz-3);
	nz_face(2,0)=1;      nz_face(2,1)=(nz-3);
	nz_face(3,0)=1;      nz_face(3,1)=(nz-3);
	nz_face(4,0)=0;      nz_face(4,1)=0;
	nz_face(5,0)=(nz-2); nz_face(5,1)=(nz-2);

	ny_face(0,0)=0;      ny_face(0,1)=0;
	ny_face(1,0)=1;      ny_face(1,1)=(ny-3);
	ny_face(2,0)=(ny-2); ny_face(2,1)=(ny-2);
	ny_face(3,0)=1;      ny_face(3,1)=(ny-3);
	ny_face(4,0)=1;      ny_face(4,1)=(ny-3);
	ny_face(5,0)=1;      ny_face(5,1)=(ny-3);

	nx_face(0,0)=1;      nx_face(0,1)=(nx-3);
	nx_face(1,0)=(nx-2); nx_face(1,1)=(nx-2);
	nx_face(2,0)=1;      nx_face(2,1)=(nx-3);
	nx_face(3,0)=0;      nx_face(3,1)=0;
	nx_face(4,0)=1;      nx_face(4,1)=(nx-3);
	nx_face(5,0)=1;      nx_face(5,1)=(nx-3);


	nz_bar(0,0)=0;      nz_bar(0,1)=0;
	nz_bar(1,0)=0;      nz_bar(1,1)=0;
	nz_bar(2,0)=0;      nz_bar(2,1)=0;
	nz_bar(3,0)=0;      nz_bar(3,1)=0;
	nz_bar(4,0)=(nz-2); nz_bar(4,1)=(nz-2);
	nz_bar(5,0)=(nz-2); nz_bar(5,1)=(nz-2);
	nz_bar(6,0)=(nz-2); nz_bar(6,1)=(nz-2);
	nz_bar(7,0)=(nz-2); nz_bar(7,1)=(nz-2);
	nz_bar(8,0)=1;      nz_bar(8,1)=(nz-3);
	nz_bar(9,0)=1;      nz_bar(9,1)=(nz-3);
	nz_bar(10,0)=1;     nz_bar(10,1)=(nz-3);
	nz_bar(11,0)=1;     nz_bar(11,1)=(nz-3);

	ny_bar(0,0)=0;       ny_bar(0,1)=0;
	ny_bar(1,0)=1;       ny_bar(1,1)=(ny-3);
	ny_bar(2,0)=(ny-2);  ny_bar(2,1)=(ny-2);
	ny_bar(3,0)=1;       ny_bar(3,1)=(ny-3);
	ny_bar(4,0)=0;       ny_bar(4,1)=0;
	ny_bar(5,0)=1;       ny_bar(5,1)=(ny-3);
	ny_bar(6,0)=(ny-2);  ny_bar(6,1)=(ny-2);
	ny_bar(7,0)=1;       ny_bar(7,1)=(ny-3);
	ny_bar(8,0)=0;       ny_bar(8,1)=0;
	ny_bar(9,0)=0     ;  ny_bar(9,1)=0;
	ny_bar(10,0)=(ny-2); ny_bar(10,1)=(ny-2);
	ny_bar(11,0)=(ny-2); ny_bar(11,1)=(ny-2);

	nx_bar(0,0)=1;       nx_bar(0,1)=(nx-3);
	nx_bar(1,0)=(nx-2);  nx_bar(1,1)=(nx-2);
	nx_bar(2,0)=1;       nx_bar(2,1)=(nx-3);
	nx_bar(3,0)=0;       nx_bar(3,1)=0;
	nx_bar(4,0)=1;  		 nx_bar(4,1)=(nx-3);
	nx_bar(5,0)=(nx-2);  nx_bar(5,1)=(nx-2);
	nx_bar(6,0)=1; 			 nx_bar(6,1)=(nx-3);
	nx_bar(7,0)=0;  		 nx_bar(7,1)=0;
	nx_bar(8,0)=0;    	 nx_bar(8,1)=0;
	nx_bar(9,0)=(nx-2);  nx_bar(9,1)=(nx-2);
	nx_bar(10,0)=(nx-2); nx_bar(10,1)=(nx-2);
	nx_bar(11,0)=0; 		 nx_bar(11,1)=0;


	n_corner[0] =        0*(nx-1)*(ny-1) +      0*(nx-1) + 0;
	n_corner[1] =	     0*(nx-1)*(ny-1) +      0*(nx-1) + (nx-2);
	n_corner[2] =	     0*(nx-1)*(ny-1) + (ny-2)*(nx-1) + (nx-2);
	n_corner[3] =	     0*(nx-1)*(ny-1) + (ny-2)*(nx-1) + 0;
	n_corner[4] =	(nz-2)*(nx-1)*(ny-1) +      0*(nx-1) + 0;
	n_corner[5] =	(nz-2)*(nx-1)*(ny-1) +      0*(nx-1) + (nx-2);
	n_corner[6] =	(nz-2)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + (nx-2);
	n_corner[7] =   (nz-2)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + 0;

	nz_interior(0,0) = 1; nz_interior(0,1) = nz-3;
	ny_interior(0,0) = 1; ny_interior(0,1) = ny-3;
	nx_interior(0,0) = 1; nx_interior(0,1) = nx-3;

}

/**
 * @brief Number of nonzero elements of the finite element system matrix
 * Description : computes the number of nonzero elements of the finite
 * element system matrix, assumes a rectangular domain with simple
 * cubic connectivity of size nx ny nz and linear elements. Then the
 * support of the basis function at each of the (nx-2)*(ny-2)*(nz-2)
 * degrees of freedom overlaps with the support of the basis functions
 * of a certain number of nodes. The number of these nodes depends on the
 * location of the degree of freedom. This defines four classes of
 * dofs: interior (= adjacent to 0 boundary nodes) face (= adjacent to
 * one boundary node) bar (= adjacent to two boundary nodes) and
 * corner (= adjacent to three boundary nodes) where here adjacency refers
 * to the simple cubic lattice graph.
 *
 * accesses the SnowDrift variable nNZ
 * Comments : For parallelization issues this must be in accordance
 * with the summation limits set in SnowDriftA3D::classifySubdomain()
 */
int SnowDriftA3D::numberOfNonzeros()
{
	int nnz = 27 * (nx-2) * (ny-2) * (nz-2) +
		18 * 2 * ( (nx-2) * (ny-2) + (nz-2) * (ny-2) + (nz-2) * (nx-2) ) +
		12 * 4 * ( (nx-2) + (ny-2) + (nz-2) ) +
		8 * 8;

	return nnz;
}

/**
 * @brief prepare sparse matrix
 * Description : generates the auxiliary arrays (i.e. the arrays
 * row-pointer and column-index) for a sparse matrix of CSR
 * format. The sparsity pattern is assumed to be given by a finite
 * element problem on a nodal mesh which is topologically equivalent
 * to a simple cubic lattice with N=nx*ny*nz number of nodes. For
 * these hexahedral elements each node belongs to 8 elements. Hence,
 * the support of the basis function of that given node overlaps with
 * the support of 27 nodes' basis functions (including its own) which
 * implies 27 nonzero entries per node in the interior of the
 * domain. On the surface one has 18, 12 and 8 nonzero entries
 * depending on the location (interior, edge, corner) From that the
 * total number of nonzero elements can be computed.  Usually, the
 * dimensions of colA and rowA are nNZ and nDOF + 1,respectively where
 * the last entry of rowA holds the number of nonzeros. Note: Here the
 * sizes are nNZ + 1 for colA nDOF + 2 for rowA. This stems from the
 * fact that originally Numerical recipes routines are used for the
 * sparse matrix arrays in numerical recipes start with index one. In
 * order to adopt them an the size was simply enlarged by one and the
 * zero-index entry is unused. => fragile
*/
void SnowDriftA3D::prepareSparseMatrix( CIntArray& colA_loc,
				     CIntArray& rowA_loc,
				     CDoubleArray& adjA_loc)
{

	// init aux variables for the loop
	int nnz = 0;

	rowA_loc[0]=0;
	rowA_loc[nDOF]=colA_loc.getNx();

	//loop over all degrees of freedom to compute the number of non-zero
	//matrix elements per dof
	for ( unsigned int iz = 0; iz < nz; iz++) {
		for ( unsigned int iy = 0; iy < ny; iy++) {
			for ( unsigned int ix = 0; ix < nx; ix++) {
				//node takes values from 1...nDOF
				const unsigned int node = iz * nx*ny + iy * nx + ix;
				nnz = 0;

				unsigned int diagonal = rowA_loc[node] + nnz; //HACK: to initialize with something...

				// loop over adjacent nodes ( adjacency defined by graph of the
				// sparsity pattern
				for ( int kz = (signed)iz-1; kz <= (signed)iz+1; kz++) {
					for ( int ky = (signed)iy-1; ky <= (signed)iy+1; ky++) {
						for ( int kx = (signed)ix-1; kx <= (signed)ix+1; kx++) {
							const unsigned int neighborNode = kz * (nx)*(ny) + ky * (nx) + kx;

							// if adjacent node is not a degree of freedom
							if ( ( kx < 0 || kx >= (signed)nx ) || ( ky < 0 || ky >= (signed)ny )
								|| ( kz < 0 || kz >= (signed)nz ) ) {
							//skip the node
							} else {
								//insert column index
								colA_loc[ rowA_loc[node] + nnz ] = neighborNode;

								//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
								//LH_DEBUG: Algorithm test BEGIN
								adjA_loc[ rowA_loc[node] + nnz ] = -1;

								if ( (neighborNode == node) ) {
									diagonal = rowA_loc[node] + nnz;
								}
								//LH_DEBUG: Algorithm test END
								//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

								nnz++;
							}
						}
					}
				}
				rowA_loc[ node + 1 ] = rowA_loc[ node ] + nnz;

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				//LH_DEBUG: Algorithm test BEGIN
			// 	      nnzA[ node ] = nnz;
				adjA_loc[ diagonal ] = nnz+node/(nx*ny);
			// 	      precond[ node ] = nnz+node/(nx*ny);//nnz;
				//LH_DEBUG: Algorithm test END
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			}
		}
	}

}

/**
 * @brief Set Bottom Boundary Condition
 * Set the boundary condition of bottom layer +intialize concentration in rest of domain.
 * @param c00 initial concentration/humidity at bottom
 */
void SnowDriftA3D::setBC_BottomLayer(CDoubleArray& var00,  const param_type param)
{

	if (param==CON){
		//set BC for the bottom layer ( ext. node layer number one)
		for (unsigned int iy = 0; iy < ny; iy++) {
			for (unsigned int ix = 0; ix < nx; ix++) {
				if (saltation(ix,iy) > 0) {
					var00[ iy*nx + ix ] = S_TO_H(psum.grid2D(ix,iy))/WS0 + c_salt(ix , iy);
				} else {
					var00[ iy*nx + ix ] = S_TO_H(psum.grid2D(ix,iy))/WS0;
				}
				for (unsigned int iz = 1; iz < nz; iz++) {
					var00[ nx*ny*iz + iy*nx + ix ] = S_TO_H(psum.grid2D(ix,iy))/WS0;
				}
			}
		}
		printf("BC value (0,0): %f\n", S_TO_H(psum.grid2D(0,0))/WS0);
	}else if (param==HUM){
		values_nodes_to_elements(nodes_q_ini, var00 );
	}else if (param==TEM){
		Grid3DObject tempT;
		tempT.set(nx,ny,nz,ta.cellsize,ta.llcorner);
		for (unsigned int iz = 0; iz<nz; iz++){
			for (unsigned int iy = 0; iy < ny; iy++) {
				for (unsigned int ix = 0; ix < nx; ix++) {
					tempT.grid3D(ix,iy,iz)=(nodes_Tair_ini.grid3D(ix,iy,iz)+(Cst::gravity/Constants::specific_heat_air)*nodes_z.grid3D(ix,iy,iz));
				}
			}
		}
		values_nodes_to_elements(tempT,var00);
	}
}

/**
* @brief Set Robin Boundary Condition
* @param aspect Location of boundary (NSWE or top/bottom)
* @param gamma_val (double) constant to approach Dirichlet or Neumann condition
* @param ix, iy, iz location
* @param var00 Initial 3D field (T, q or c) that needs to be set
* @param param Type of field (T, q or c)
*/
void SnowDriftA3D::setRobinBoundaryCondition(const aspect_type aspect, const double gamma_val, const int ix, const int iy, const int iz, CDoubleArray& var00, const param_type param)
{
	const int node=nx*ny*iz+nx*iy+ix;

	gamma[ node ] = gamma_val;
	if (aspect==BOTTOM && (param==HUM || param==TEM)){
		gamma[node]=0.;
	}else if (aspect==BOTTOM && param==CON){
		gamma[ node ] = gamma_val;
	}
	gDirichlet[ node ] = var00[node];
	gNeumann[ node ] = 0.;
}

/**
 * @brief Initialize system
 * initializes vectors for the algorithm required for each time step
 * @param colA array determining the sparsity pattern
 * @param rowA array determining the sparsity pattern
 * @param sA system matrix in compressed sparse row (CSR) storage format
 * @param sB system matrix in compressed sparse row (CSR) storage format
 * @param rhs vector which contains the right hand side of the linear system
 * @param f source term in the advection/diffusion equation (sublimation amount)
 * @param Psi auxiliary vector for incorporating inhomogeneous Dirichlet boundary conditions
 * @param var variable for which equation has to be solved (Humidity, Concentration or Temperature)
 * @param var00 boundary condition at bottom + initial condition interior domain of variable var
 */

void SnowDriftA3D::initializeSystem( CIntArray& colA_loc, CIntArray& rowA_loc, CDoubleArray& sA_loc, CDoubleArray& sB_loc, CDoubleArray& rhs_loc, CDoubleArray& f_loc, CDoubleArray& Psi_loc,CDoubleArray& var, CDoubleArray& var00,const param_type param)
{
	for (unsigned int i = 0; i < nDOF; i++) {
		// Initialization of source/sink terms on the interior of the domain

		precond[i] = 0;
		Psi_loc[i] = 0;
		colA_loc[i] = 0;
		rowA_loc[i] = 0;
		rhs_loc[i] = 0;
		sA_loc[i] = 0;
		sB_loc[i] = 0;
		var[i]=0;
		var00[i]=0;
		f_loc[i] = 0;
	}

	rowA[nDOF] = 0;
	for ( unsigned int i = nDOF; i < nNZ; i++) {
		sA_loc[i] = 0;
		colA_loc[i] = 0;
		adjA[i] = 0;
	}

	setBC_BottomLayer(var00, param);

	for (unsigned int iz = 0; iz < nz; iz++) {
		for (unsigned int iy = 0; iy < ny; iy++) {
			for (unsigned int ix = 0; ix < nx; ix++) {
				//east
				if (	( ix == 0 )
					|| ( ix == nx-1 )
					|| ( iy == ny-1 && !(ix==0 || ix==nx-1) )
					|| ( iy == 0 && !(ix==0 || ix==nx-1))
					|| ( iz == nz-1 && !(iy==0 || iy==ny-1 || ix==0 || ix==nx-1))) {

						//if East, West, North South or Top
						setRobinBoundaryCondition(OTHER, 1.e5, ix,iy,iz, var00, param);
				}
				//if bottom
				if ( iz == 0 && !(iy==0 || iy==ny-1 || ix==0 || ix==nx-1) ) {
					double K[3][3];
					computeDiffusionTensor(K,ix,iy, iz);
					double K_perp = K[2][2];
					setRobinBoundaryCondition(BOTTOM, K_perp,ix,iy,iz, var00, param);
				} else {
			//    		  gamma[ node ] = 1e5;
			//    		  gDirichlet[ node ] = 0.001;//c00[node];
			//    		  gNeumann[ node ] = 0;//0.5 * c_salt[ ix ][ iy ];//u_perp *c00[node];
			//    		  gDirichlet[ node ] = ( S_TO_H(sn_MdataT.nswc/WS0)
			//    					 *( 1 + sin( 2*M_PI * ix/(1.0*(nx-1)) ) )
			//   					 *( 1 + sin( 2*M_PI * iy/(1.0*(ny-1)) ) ) );
				}
			}//iz
		}//iy
	}//ix
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function    : resetArray
// Authors     : Henning Loewe
// Description : set Array
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SnowDriftA3D::resetArray( CIntArray& cda )
{
	for (unsigned int i = 0; i < cda.getNx(); i++) {
		cda[i] = 0;
	}
}

void SnowDriftA3D::resetArray( CDoubleArray& cda )
{
	for (unsigned int i = 0; i < cda.getNx(); i++) {
		cda[i] = 0;
	}
}

/**
 * @brief nodes_to_elements
 * @param nodesGrid 3Dgrid of nodes
 * @param elementsArray array with elements
 */
void SnowDriftA3D::values_nodes_to_elements(const Grid3DObject& nodesGrid, CDoubleArray& elementsArray )
{
	const int nx_grid=nodesGrid.getNx();
	const int ny_grid=nodesGrid.getNy();
	const int nxy=nx_grid*ny_grid;
	const int nz_grid=nodesGrid.getNz();
	const unsigned int Nelems=elementsArray.getNx();

	for (int i=0; i<(signed)Nelems; i++){
		//find the nodes for this element
		int iz = (int)floor(i/nxy);
		int iy = (int)floor(i-iz*nxy) / nx_grid;
		int ix = i-iz*nxy-iy*nx_grid;
		//elements should always be surrounded by 8 nodes, but there are some extra elements
		iz = std::min(iz, nz_grid-2);
		iy = std::min(iy, ny_grid-2);
		ix = std::min(ix, nx_grid-2);

		double value=0;
		unsigned int count=0;
		for (int ii=ix; ii<=ix+1; ii++){
			for (int jj=iy; jj<=iy+1; jj++){
				for (int kk=iz; kk<=iz+1; kk++){
				value += nodesGrid.grid3D(ii,jj,kk);
				count++;
				}
			}
		}
		elementsArray(i) = value/count;
	}
}

/**
 * @brief elements_to_nodes
 * Extract a 3D-grid of nodes from an array of elements
 * @param nodesGrid grid of nodes that should be filled
 * @param elementsArray array with elements
 */
void SnowDriftA3D::values_elements_to_nodes(Grid3DObject& nodesGrid, const CDoubleArray& elementsArray )
{
    const int ncols=nodesGrid.getNx();
    const int nrows=nodesGrid.getNy();
    const int ndepths=nodesGrid.getNz();
    int ixmin=0;
    int ixmax=0;
    int iymin=0;
    int iymax=0;
    int izmin=0;
    int izmax=0;

    for (int ii=0; ii<=ncols-1; ii++){
	for (int jj=0; jj<=nrows-1; jj++){
	    for (int kk=0; kk<=ndepths-1; kk++){

		if (ii==0){
		    ixmin=ii;
		    ixmax=ii;
		}else if (ii==ncols-1){
		    ixmin=ii-1;
		    ixmax=ii-1;
		}else{
		    ixmin=ii-1;
		    ixmax=ii;
		}

		if (jj==0){
		    iymin=jj;
		    iymax=jj;
		}else if (jj==nrows-1){
			iymin=jj-1;
			iymax=jj-1;//iymax=jj;
		}else{
		    iymin=jj-1;
		    iymax=jj;
		}

		if (kk==0){
		    izmin=kk;
		    izmax=kk;
		}else if (kk==ndepths-1){
		    izmin=kk-1;
		    izmax=kk-1;
		}else{
		    izmin=kk-1;
		    izmax=kk;
		}

		double value = 0.;
		int count = 0;
		for ( int ix=ixmin; ix<=ixmax; ix++){
		    for (int iy = iymin; iy<=iymax; iy++) {
			for (int iz = izmin; iz<=izmax; iz++) {
			    int element = iz*(ncols)*(nrows)+iy*(ncols)+ix;
			    value += elementsArray(element);
			    count+=1;
			}
		    }
		}
		value=value/count;
		nodesGrid.grid3D(ii,jj,kk)=value;
	    }
	}
    }

}
