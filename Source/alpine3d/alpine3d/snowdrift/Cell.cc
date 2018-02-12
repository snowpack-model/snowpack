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
#include <alpine3d/snowdrift/Cell.h>

Cell::Cell()
{
	initialize();
}

/**
 * @brief arrays to select the interior nodes of a boundary element located either on
 * faces, corners, or bars. Note: the corners are enumerated from the
 * bottom layer to the top layer in the usual order of enumeration
 * used for the local node indices in an element. The faces are
 * ordered according to: front, right, back, lef,t bottom, top.  the
 * bars are in order: 1) horizontal bars in the bottom: front right
 * back left 2) horizontal bars on the top: same order 3) vertical
 * bars: starting from front-left
 */
void Cell::initialize()
{
	const int aux1[6][4] = { {2, 3, 6, 7},
			{0, 3, 4, 7},
			{0, 1, 4, 5},
			{1, 2, 5, 6},
			{4, 5, 6, 7},
			{0, 1, 2, 3} };

	for ( int i = 0; i < 6; i++) {
		for ( int j = 0; j < 4; j++) {
			FACE[i][j] = aux1[i][j];
		}
	}


	const int aux2[12][2] = { {6,7}, {4,7}, {4,5},
			{5,6}, {2,3}, {0,3},
			{0,1}, {1,2}, {2,6},
			{3,7}, {0,4}, {1,5} };

	for ( int i = 0; i < 12; i++) {
		for ( int j = 0; j < 2; j++) {
			BAR[i][j] = aux2[i][j];
		}
	}


	const int aux3[8] = {6, 7, 4, 5, 2, 3, 0, 1};
	
	for ( int i = 0; i < 8; i++) {
		CORNER[i] = aux3[i];
	}

	// array to pick out the boundary nodes of an element located on a face,
	// hence: NODE[i][] = {0,1,2..7} \ FACE[i][]
	const int aux4[6][4] = { {0, 1, 4, 5},
			{1, 2, 5, 6},
			{2, 3, 6, 7},
			{0, 3, 4, 7},
			{0, 1, 2, 3},
			{4, 5, 6, 7} };

	for ( int i = 0; i < 6; i++) {
		for ( int j = 0; j < 4; j++) {
			NODE[i][j] = aux4[i][j];
		}
	}
}

void Cell::classifyNodes(int* dofNode,
			 int* nDofNodes,
			 int* nBoundaryNodes,
			 const std::string& type,
			 int number) const
{
	if ( type == "interior" ) {
		for (int i = 0; i < 8; i++) {
			dofNode[i]=i;
		}
		*nDofNodes = 8;
		*nBoundaryNodes = 0;
	} else if ( type == "face" ) {
		//specify the nodes of the element which are interior nodes
		for (int i = 0; i < 4; i++) {
			dofNode[i] = FACE[number][i];
		}
		*nDofNodes=4;

		//specify the nodes of the element which are boundary nodes
		for (int i = 0; i < 4; i++) {
			dofNode[i+*nDofNodes] = NODE[number][i];
		}
		*nBoundaryNodes = 4;
	} else if ( type == "bar" ) {
		//specify the nodes of the element which are interior nodes
		dofNode[0] = BAR[number][0];
		dofNode[1] = BAR[number][1];

		*nDofNodes=2;

		//specify the nodes of the element which are boundary nodes
		//works only iff BAR[bar][0] < BAR[bar][1] for all bar
		for (int i = 0; i < BAR[number][0]; i++) {
			dofNode[i+*nDofNodes]=i;
		}
		for (int i = BAR[number][0]; i < BAR[number][1] - 1; i++) {
			dofNode[i+*nDofNodes]=i+1;
		}
		for (int i = BAR[number][1]-1; i < 6; i++) {
			dofNode[i+*nDofNodes]=i+2;
		}
		*nBoundaryNodes = 6;
	} else if ( type == "corner" ) {
		//specify the nodes of the element which are interior nodes
		dofNode[0] = CORNER[number];
		*nDofNodes = 1;

		//specify the nodes of the element which are boundary nodes
		for (int i = 0; i < CORNER[number]; i++) {
			dofNode[i+*nDofNodes]=i;
		}
		for (int i = CORNER[number]; i < 8; i++) {
			dofNode[i+*nDofNodes]=i+1;
		}
		*nBoundaryNodes = 7;
	} else {
		exit(1);
	}
}
