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
#ifndef TERRAINRadiationComplex_H
#define TERRAINRadiationComplex_H

#include <meteoio/MeteoIO.h>

#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
#include <alpine3d/ebalance/RadiationField.h>
#include <alpine3d/ebalance/SolarPanel.h>
#include <alpine3d/ebalance/SnowBRDF.h>
#include <array>

/**
 * @page TerrainRadiationComplex

 * This module calculates the radiative transfer of SW radiation in snow-covered terrain. It can take into account anisotropic reflection
 * of light on snow (forward scattering) and multiple reflections in the terrain. The price for this complexity is high:
 * The initialization time is extremely long (for a DEM with 30'000 pixels it can be 12 hours. However, if a ViewList-file is generated
 * during the initialization, it can be bypassed in later simulations. ) In the main computation loop, the effort is considerable too
 * and may dominate the simulation time in many cases. Therefore, this module is recommended for explicit radiation investigations based
 * on relatively small, high-resolution DEM's. It can also be coupled to the SolarPanel-module for precise simulation of radiation on
 * photovoltaic panels.
 * The algorithm is described in the master thesis of Felix von Rütte, entitled "Radiative Transfer Model for Snowy Mountains". In the
 * comments is often referred to by the abbreviation "MT". If you don't find a version online, ask Michael Lehning, he was the supervisor.
 *
 * @keys in io-file
 * COMPLEX_ANISOTROPY [true or false] 		: Whether a snow BRDF should be included. False means isotropic scattering.
 * COMPLEX_MULTIPLE [true or false]			: Whether multiple scattering in the terrain should be taken into account.
 * COMPLEX_WRITE_VIEWLIST [true or false]	: Whether the initialization stuff should be written to file. (Make sure you have folder "output")
 * COMPLEX_READ_VIEWLIST [true or false]	: Whether an existing initialization file should be read in; bypassing the initialization.
 * COMPLEX_VIEWLISTFILE [<path>/<filename>]	: Path to the ViewList file if existing. (e.g ../input/surface-grids/ViewList_Totalp_30x30.rad)
 *
 *
 *
 * @code
 * [EBalance]
 * TERRAIN_RADIATION = TRUE
 * TERRAIN_RADIATION_METHOD = COMPLEX
 *
 * COMPLEX_ANISOTROPY	= TRUE
 * COMPLEX_MULTIPLE	= TRUE
 * COMPLEX_WRITE_VIEWLIST	= FALSE
 * COMPLEX_READ_VIEWLIST	= TRUE
 * COMPLEX_VIEWLISTFILE = ../input/surface-grids/ViewList_Totalp_30x30.rad
 * @endcode
 *
 */

class TerrainRadiationComplex : public TerrainRadiationAlgorithm
{

public:
	TerrainRadiationComplex(const mio::Config &cfg, const mio::DEMObject &dem_in, const std::string &method);
	~TerrainRadiationComplex();

	virtual void getRadiation(mio::Array2D<double>& direct, mio::Array2D<double>& diffuse,
                            mio::Array2D<double>& terrain, const mio::Array2D<double>& direct_unshaded_horizontal,
                            const mio::Array2D<double>& total_ilwr, mio::Array2D<double>& sky_ilwr,
                            mio::Array2D<double>& terrain_ilwr, double solarAzimuth, double solarElevation);
	virtual void setMeteo(const mio::Array2D<double>& albedo, const mio::Array2D<double>& ta);

	void getSkyViewFactor(mio::Array2D<double> &o_sky_vf);
	void setSP(const mio::Date timestamp, const double solarAzimuth, const double solarElevation);
	void writeSP(const unsigned int max_steps);

private:
	typedef std::array<double, 3> Vec3D;

	// Initialisation Functions
	void initBasicSetHorizontal();
	void initBasicSetRotated();
	void initViewList();
	void initRList();
	void initSortList();
	void WriteViewList();
	bool ReadViewList();

	// auxiliary functions
	void TriangleNormal(size_t ii_dem, size_t jj_dem, int which_triangle, Vec3D &v_out);
	double IntersectionRayTriangle(const Vec3D &v_view, size_t ii_0, size_t jj_0, size_t ii_dem, size_t jj_dem, size_t which_triangle);
	size_t vectorToSPixel(const Vec3D &vec_in, size_t ii_dem, size_t jj_dem, size_t which_triangle);
	void initSkyViewFactor();
	double computeSkyViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle);
	double getLandViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle);
	double getSkyViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle);
	void getVectorSun(double solarAzimuth, double solarElevation, Vec3D &v_out);
	double TerrainBiggestDifference(const mio::Array3D<double> &terrain_old, const mio::Array3D<double> &terrain_new);

	// Standard Vector operations
	double NormOfVector(const Vec3D &vec1);
	void normalizeVector(const Vec3D &vec1, Vec3D &v_out);
	double VectorScalarProduct(const Vec3D &vec1, const Vec3D &vec2);
	void VectorCrossProduct(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out);
	void VectorSum(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out);
	void VectorDifference(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out);
	void VectorStretch(const Vec3D &vec1, double factor, Vec3D &v_out);
	void RotN(const Vec3D &axis, const Vec3D &vec_in, double rad, Vec3D &v_out);
	void ProjectVectorToPlane(const Vec3D &vec1, const Vec3D &plane_normal, Vec3D &v_out);
	double AngleBetween2Vectors(const Vec3D &vec1, const Vec3D &vec2);

	// PVP functions
	void readSP();


	// Output functions
	void PrintProgress(double percentage);

	// Variables
	const size_t dimx, dimy;
	size_t startx, endx;
	mio::DEMObject dem;
	const mio::Config &cfg;

	SnowBRDF BRDFobject;
	SolarPanel SP;
	std::vector<std::vector<double> > pv_points;

	mio::Array3D<std::vector<double>> SortList;			  // Used for speedup in Terrain Iterations
	std::vector<Vec3D> BasicSet_Horizontal; // Horizontal Basic Set [MT 2.1.1 Basic Set]
	mio::Array4D<Vec3D> BasicSet_rotated;	  // Basic Set rotated in Triangular Pixel Plane [MT 2.1.3 View List, eq. 2.38]
	mio::Array4D<std::vector<double>> ViewList;			  // Stores all information of network between pixels [MT 2.1.3 View List, eq. 2.47]
	mio::Array2D<double> RList;							  // List pre-storage of BRDF values
	mio::Array2D<double> albedo_grid;					  // Albedo value for each square Pixel

	mio::Array2D<double> sky_vf_mean;		  // Grid to store view factor
	std::vector<mio::Array2D<double>> sky_vf; // Array to stor grid of view factor for both values of which_triangle
	unsigned int M_epsilon;					  // Number of small circles in Basic Set [MT fig. 2.1]
	unsigned int M_phi;						  // Number of vectors per small circle of Basic Set [MT fig. 2.1]
	unsigned int S;							  // Number of vectors per Basic Set [MT fig. 2.1]
	double delta_F_max;						  // Stopping Treshold for Iteration in W/m2 [MT eq. 2.100]

	// Keys from io-file
	bool if_anisotropy = false;		// Anisotropic or Isotropic Snow Model ?
	bool if_multiple = false;		// Do Multiple Scattering in Terrain or not ?
	bool if_write_view_list = true; // Write ViewList to file ?
	bool if_read_view_list = false; // Read existing View-list file? -> Speeds up Initialisation by factor ~200
};

#endif
