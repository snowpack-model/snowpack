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
#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>

class Cell
{
	public:
		Cell();
		
		void classifyNodes(int* node,
				int* nDofNodes,
				int* nBoundaryNodes,
				const std::string& type,
				int number) const;
		
	private:
		void initialize();
		
		int FACE[6][4];
		int NODE[6][4];
		int CORNER[8];
		int BAR[12][2];
};

#endif
