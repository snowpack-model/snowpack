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
#ifndef EBSTRUCT_H
#define EBSTRUCT_H

// Ssolar contains the information for sun position
 typedef struct
 {

	double ecc_corr;       // correction due to the eccentricity of the earth's orbit
	double decl;           // solar declination: angle between the vector earth/sun
			       // and the equatorial plane
	double eq_time;        // equation of time
	double azi;            // solar azimuth
	double elev;           // solar elevation
	double toa_h;          // incoming irradiance on horizontal area on top of the atmosphere
	double solar_time;     // true solar time
	double hr_angle;

} Ssolar;

// Smeteo contains the meteo input information as well as time information and atmosphere correction factors
 typedef struct
 {
	int hh;                // actual hour
	double gr;             // incoming global radiation
	double ea;             // atmospheric emissivity (longwave radiation, cloud cover)
	double tdirhor;        // at station or biggest theoretical incident direct solar radiation on a flat field
	double tdiff;          // at station or biggest theoretical clear sky diffuse radiation
	double Md;             // splitting coefficient for direct and diffuse solar radiation

 } Smeteo;

 #endif
