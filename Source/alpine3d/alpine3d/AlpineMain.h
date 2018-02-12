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
#ifndef ALPINEMAIN_H
#define ALPINEMAIN_H

// _VERSION is given as a compilation flag to tell us what is the version number
// Please only use A3D_VERSION in the code
#ifndef A3D_VERSION
	//here below, the double-expansion stringification macro trick...
	#define STR1(x) #x
	#define STR2(x) STR1(x)
	#define A3D_VERSION STR2( _VERSION )
#endif

//Only for ALPINE3d the time steps are defined here
#define dt_main 3600. /* Large Calculation step length for EB, Saltation, Snowpack etc.  */

#endif
