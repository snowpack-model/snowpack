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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOB.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <alpine3d/ebalance/TerrainRadiationComplex.h>
#include <alpine3d/MPIControl.h>

//file operations
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace mio;

TerrainRadiationComplex::TerrainRadiationComplex(const mio::Config &cfg_in, const mio::DEMObject &dem_in,
                                                 const std::string &method)
  : TerrainRadiationAlgorithm(method), dimx(dem_in.getNx()), dimy(dem_in.getNy()), startx(0), endx(dimx),
    dem(dem_in), cfg(cfg_in), BRDFobject(cfg_in), pv_points(),
    albedo_grid(dem_in.getNx(), dem_in.getNy(), IOUtils::nodata),
    sky_vf(2, mio::Array2D<double>(dimx, dimy, IOUtils::nodata)), sky_vf_mean(dimx, dimy, IOUtils::nodata)
{
	// PRECISION PARAMETERS
	// ####################
	M_epsilon = 30;	 // Number of small circles in Basic Set
	M_phi = 30;		 // Number of vectors per small circle of Basic Set
	delta_F_max = 1; // Stopping Treshold for Iteration in W/m2 [MT eq. 2.100]
	// ####################

	size_t nx = dimx;
	MPIControl::instance().getArraySliceParams(dimx, startx, nx);
	endx = startx + nx;
	if (startx == 0)
	{
		startx = 1;
	}

	if (endx >= dimx - 1)
	{
		endx = dimx - 1;
	}

	std::cout << "[i] In TerrainRadiationComplex: Process " << MPIControl::instance().rank() << " working on slice (" << startx << ", " << endx << ")" << std::endl;

	S = M_epsilon * M_phi; // Number of vectors per Basic Set

	if (cfg.keyExists("Complex_Anisotropy", "Ebalance"))
		cfg.getValue("Complex_Anisotropy", "Ebalance", if_anisotropy);
	else
		std::cout << "[i] In TerrainRadiationComplex: No flag <<Complex_Anisotropy>> set in [Ebalance]. Use default "
		<< "Complex_Anisotropy = " << if_anisotropy << " .\n";

	if (cfg.keyExists("Complex_Multiple", "Ebalance"))
		cfg.getValue("Complex_Multiple", "Ebalance", if_multiple);
	else
		std::cout << "[i] In TerrainRadiationComplex: No flag <<Complex_Multiple>> set in [Ebalance]. Use default "
		<< "Complex_Multiple = " << if_multiple << " .\n";

	if (cfg.keyExists("Complex_Write_Viewlist", "Ebalance"))
		cfg.getValue("Complex_Write_Viewlist", "Ebalance", if_write_view_list);
	else
		std::cout << "[i] In TerrainRadiationComplex: No flag <<Complex_Write_Viewlist>> set in [Ebalance]. Use default "
		<< "Complex_Write_Viewlist = " << if_write_view_list << " .\n";

	if (cfg.keyExists("Complex_Read_Viewlist", "Ebalance"))
	{
		if (cfg.get("Complex_Read_Viewlist", "EBalance"))
			if_read_view_list = true;
	}
	else
		std::cout << "[i] In TerrainRadiationComplex: No flag <<Complex_Read_Viewlist>> set in [Ebalance]. Use default "
		<< "Complex_Read_Viewlist = " << if_read_view_list << " .\n";

	if (cfg.keyExists("PVPFILE", "EBalance"))
	{
		//load PVP data
		readSP();
		//performs some validation on loaded data
		std::vector<Coords> co_vec;
		for (size_t ii = 0; ii < pv_points.size(); ii++)
		{
			Coords point;
			point.setXY(pv_points[ii][0], pv_points[ii][1], pv_points[ii][2]);
			co_vec.push_back(point);
		}
		bool master = MPIControl::instance().master();
		if (!dem.gridify(co_vec, true))
		{ //keep invalid points
			if (master)
				std::cerr << "[E] Some PVP are invalid or outside the DEM:\n";
			for (size_t ii = 0; ii < co_vec.size(); ii++)
				if (!co_vec[ii].indexIsValid() && master)
					std::cout << "[E] Point " << ii << "\t" << co_vec[ii].toString(Coords::CARTESIAN) << "\n";
			throw InvalidArgumentException("Invalid PVP, please check in the logs", AT);
		}
		else if (master)
			std::cout << "[i] Using " << pv_points.size() << " PVP\n";

		SP = SolarPanel(cfg, dem, pv_points);
		_hasSP=true;
	}

	initBasicSetHorizontal();
	std::cout << "[i] Initialized BasicSetHorizontal" << std::endl;
	initBasicSetRotated();
	std::cout << "[i] Initialized BasicSetRotated" << std::endl;

	// Initialise based on ViewList
	bool succesful_read = false;
	if (if_read_view_list)
	{
		succesful_read = ReadViewList();
	}
	// Initialise from scratch
	if (!succesful_read)
	{
		initViewList();
	}
	std::cout << "[i] Initialized ViewList" << std::endl;

	// Initialise Speed-up
	initRList();
	std::cout << "[i] Initialized RList" << std::endl;
	initSortList();
	initSkyViewFactor();
	if (if_write_view_list)
		WriteViewList(); // Write ViewList to file

	if (_hasSP)
		SP.initTerrain(M_epsilon, M_phi); // Link SolarPanel-object to ViewList

}

TerrainRadiationComplex::~TerrainRadiationComplex() {}

//########################################################################################################################
//                                             INITIALISATION FUNCTIONS
//########################################################################################################################

/**
* @brief Initializes set of Vectors that point to equal solid angels for horizontal hemisphere (BasicSet) [MT 2.1.1 Basic Set]
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::initBasicSetHorizontal()
{

	double psi_0 = 0, psi_1 = 0, epsilon = 0;
	double delta = 0;
	double phi = 0;
	int i = 0;

	BasicSet_Horizontal.clear();

	for (size_t n = 0; n < M_epsilon; ++n)
	{
		psi_0 = psi_1;
		if (n == (M_epsilon - 1))
		{
			psi_1 = 3.1415926 / 2;
		}
		else
		{
			delta = acos(-2 / float(M_epsilon) + cos(2 * psi_0)) / 2 - psi_0; // [~ MT eq. 2.5]
			psi_1 = psi_0 + delta;
		}
		epsilon = (psi_1 + psi_0) / 2; // [MT eq. 2.4]

		for (size_t m = 0; m < M_phi; ++m)
		{
			phi = 2 * 3.1415926 * (m + 0.5) / M_phi; // [MT eq. 2.3]
			i += 1;
			BasicSet_Horizontal.push_back({cos(epsilon) * sin(phi), cos(epsilon) * cos(phi), sin(epsilon)}); // [MT eq. 2.1]
		}
	}
}

/**
* @brief Rotates BasicSet in Plane of triangular Surface Object		[MT 2.1.3 View-List]
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::initBasicSetRotated()
{

	BasicSet_rotated.resize(dimx, dimy, 2, S, {0, 0, 0});

#pragma omp parallel for
	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				Vec3D triangle_normal, axis, z_axis = {0, 0, 1};
				TriangleNormal(ii, jj, which_triangle, triangle_normal);
				VectorCrossProduct(z_axis, triangle_normal, axis);						 // [MT eq. 2.39]
				double inclination = acos(VectorScalarProduct(triangle_normal, z_axis)); // [MT eq. 2.40]

				if (NormOfVector(axis) == 0)
				{ // [MT see text after eq. 2.41]
					axis[0] = 1;
					axis[1] = 0;
					axis[2] = 0;
				}

				for (size_t solidangle = 0; solidangle < S; ++solidangle)
				{
					Vec3D rotated;
					const auto &solid_angle = BasicSet_Horizontal[solidangle];
					RotN(axis, solid_angle, inclination, rotated); // [MT eq. 2.38]
					auto &to_rotate = BasicSet_rotated(ii, jj, which_triangle, solidangle);
					for (auto i = 0; i < rotated.size(); ++i)
					{
						to_rotate[i] = rotated[i];
					}
				}
			}
		}
	}
}

/**
* @brief Assigns a pixel (or sky) to each space vector of each ProjectVectorToPlane		[MT 2.1.2 Surface Generation and 2.1.3 View-List]
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::initViewList()
{

	std::cout << "[i] Initialize Terrain Radiation Complex\n";

	ViewList.resize(dimx, dimy, 2, S, {0, 0, 0, 0, 0});

	int counter = 0; // For output bar
//loop over all triangles of surface
#pragma omp parallel for
	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
#pragma omp parallel for
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				for (size_t solidangle = 0; solidangle < S; ++solidangle) // loop over all vectors of basic set
				{
					// Initialize with Impossible Numbers
					size_t ii_temp = 999999;
					size_t jj_temp = 999999;
					size_t solidangle_temp = 99999;
					size_t which_triangle_temp = 9;
					double distance, minimal_distance = dem.cellsize * (dem.getNx() + dem.getNy());

					// Search for Intersection Candidates
					size_t ii_dem = ii, jj_dem = jj, nb_cells = 0;
					const auto &ray = BasicSet_rotated(ii, jj, which_triangle, solidangle);
					if (ray.size() != 3)
					{
						throw IndexOutOfBoundsException("Expected array of size 3");
					}
					Vec3D projected_ray, z_axis = {0, 0, 1};
					ProjectVectorToPlane(ray, z_axis, projected_ray); // [in MT 2.1.3:  projected_ray ~ v_view,yx]
					if (NormOfVector(projected_ray) != 0)
						normalizeVector(projected_ray, projected_ray);
					else
					{
						projected_ray[0] = 1;
						projected_ray[1] = 0;
						projected_ray[2] = 0;
					}

					// [MT Figure 2.4] Smart loop over intersection Candidates [Amanatides et al. (1987)]
					while (!(ii_dem < 1 || ii_dem > dem.getNx() - 2 || jj_dem < 1 || jj_dem > dem.getNy() - 2))
					{

						for (int k = -1; k < 2; ++k)
						{
							for (int kk = -1; kk < 2; ++kk)
							{

								ii_dem = ii + (int)round(((double)nb_cells) * projected_ray[0]) + k;  // [~ MT eq. 2.32]
								jj_dem = jj + (int)round(((double)nb_cells) * projected_ray[1]) + kk; // [~ MT eq. 2.32]
								if ((ii_dem < 1 || ii_dem > dem.getNx() - 2 || jj_dem < 1 || jj_dem > dem.getNy() - 2))
									continue;

								distance = IntersectionRayTriangle(ray, ii, jj, ii_dem, jj_dem, 0); // Triangles B: [MT Figure 2.2]
								if (distance != -999 && distance < minimal_distance)
								{
									minimal_distance = distance; // If intersection with sevaral trangles, take the closest intersection
									ii_temp = ii_dem;
									jj_temp = jj_dem;
									which_triangle_temp = 0;
								}

								distance = IntersectionRayTriangle(ray, ii, jj, ii_dem, jj_dem, 1); // Triangles A: [MT Figure 2.2]
								if (distance != -999 && distance < minimal_distance)
								{
									minimal_distance = distance; // If intersection with sevaral trangles, take the closest intersection
									ii_temp = ii_dem;
									jj_temp = jj_dem;
									which_triangle_temp = 1;
								}
							}
						}
						ii_dem = ii + (size_t)round(((double)nb_cells) * projected_ray[0]);
						jj_dem = jj + (size_t)round(((double)nb_cells) * projected_ray[1]);
						nb_cells += 2.7; // [MT eq. 2.37]
					}
					if (minimal_distance == dem.cellsize * (dem.getNx() + dem.getNy()))
						minimal_distance = -999; // No intersection found    =>    -999 == "SKY"
					if (ii_temp > 0 && ii_temp < dimx)
					{
						Vec3D ray_stretched;
						VectorStretch(ray, -1, ray_stretched);
						solidangle_temp = vectorToSPixel(ray_stretched, ii_temp, jj_temp, which_triangle_temp);
					}
					ViewList(ii, jj, which_triangle, solidangle) = {(double)ii_temp, (double)jj_temp, (double)which_triangle_temp, minimal_distance, (double)solidangle_temp}; // [MT eq. 2.47]
				}
				counter++;
				if (counter % 10 == 0)
					PrintProgress((double)counter / (double)(dimx - 2) / (double)(dimy - 2) / 2.);
			}
		}
	}
	PrintProgress(1);
	std::cout << "  Done.\n";
}

/**
* @brief +SPEEDUP+ Projects BRDF on Basic Set. Instead of Calculating/Intrapolating all the time
* 					new BRDF factors, store discrete number covering whole hemisphere
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::initRList()
{

	RList.resize(S, S, 0);

// Every in/out Vector-Pair of Basic Set is assigned BRDF Value (more precise: Anisotropy Factor)
#pragma omp parallel for
	for (size_t number_solid_in = 0; number_solid_in < S; ++number_solid_in)
	{
		for (size_t number_solid_out = 0; number_solid_out < S; ++number_solid_out)
		{
			const auto &solid_in = BasicSet_Horizontal[number_solid_in];
			const auto &solid_out = BasicSet_Horizontal[number_solid_out];

			// For Syntax See SnowBRDF-Class
			Vec3D z_axis = {0, 0, 1}, projected_in, projected_out;
			double cth_i = VectorScalarProduct(solid_in, z_axis);
			double cth_v = VectorScalarProduct(solid_out, z_axis);
			ProjectVectorToPlane(solid_in, z_axis, projected_in);
			normalizeVector(projected_in, projected_in);
			ProjectVectorToPlane(solid_out, z_axis, projected_out);
			normalizeVector(projected_out, projected_out);
			double cphi = VectorScalarProduct(projected_in, projected_out);

			RList(number_solid_in, number_solid_out) = BRDFobject.get_RF(cth_i, cphi, cth_v);
		}
	}
}

/**
* @brief +SPEEDUP+: For most terrain a large part of the basicSet points in the sky. SortList Stores land-pointing vectors only.
* Only these need a iterative radiation computation. [not discussed in MT and somewhat confusing syntax in MT eq. 2.95, sorry. Better Look @ Paper ????]
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::initSortList()
{

	SortList.resize(dimx, dimy, 2);

	// Store for every triangle, the actually used Vectors (by othe triangles). Different triangles may use same vector, therefore redundancy produced
	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				for (size_t solidangle = 0; solidangle < S; ++solidangle)
				{

					double distance = ViewList(ii, jj, which_triangle, solidangle)[3];
					if (distance == -999)
						continue;

					size_t ii_source = ViewList(ii, jj, which_triangle, solidangle)[0];
					size_t jj_source = ViewList(ii, jj, which_triangle, solidangle)[1];
					size_t which_triangle_source = ViewList(ii, jj, which_triangle, solidangle)[2];
					size_t solidangle_source = ViewList(ii, jj, which_triangle, solidangle)[4];

					SortList(ii_source, jj_source, which_triangle_source).push_back(solidangle_source);
				}
			}
		}
	}

	// get rid of redundancy / Delete doubles
	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				sort(SortList(ii, jj, which_triangle).begin(), SortList(ii, jj, which_triangle).end());
				SortList(ii, jj, which_triangle).erase(unique(SortList(ii, jj, which_triangle).begin(), SortList(ii, jj, which_triangle).end()), SortList(ii, jj, which_triangle).end());
			}
		}
	}
}

/**
* @brief Writes Viewlist to file, use SMET format only
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::WriteViewList()
{
	std::ofstream GL_file;
	const std::string filename = cfg.get("Complex_ViewListFile", "Ebalance");
	GL_file.open(filename);
	GL_file.precision(std::numeric_limits<double>::digits10);
	GL_file << "SMET 1.1 ASCII\n[HEADER]\n";
	GL_file << "ncols =\t" << dimx << "\n";
	GL_file << "nrows =\t" << dimy << "\n";
	GL_file << "cellsize =\t" << dem.cellsize << "\n";
	GL_file << "easting =\t" << dem.llcorner.getEasting() << "\n";
	GL_file << "northing =\t" << dem.llcorner.getNorthing() << "\n";
	GL_file << "Elevation Points =\t" << M_epsilon << "\n";
	GL_file << "Azimuth Points =\t" << M_phi << "\n";

	GL_file << "epsg\t= 0\nnodata\t= -999\naltitude\t=0\nstation_id\t= ViewList\nstation_name\t= TerrainList\n";
	GL_file << "fields =\t ii\t jj\t which_triangle\t solidangle\t ii_seen\t jj_seen\t which_triangle_seen\t distance\t soldiangle_seen\n";

	GL_file << "[DATA]\n";

	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				for (size_t solidangle = 0; solidangle < S; ++solidangle)
				{
					size_t ii_source = ViewList(ii, jj, which_triangle, solidangle)[0];
					size_t jj_source = ViewList(ii, jj, which_triangle, solidangle)[1];
					size_t which_triangle_source = ViewList(ii, jj, which_triangle, solidangle)[2];
					double distance = ViewList(ii, jj, which_triangle, solidangle)[3];
					size_t solidangle_source = ViewList(ii, jj, which_triangle, solidangle)[4];

					GL_file << ii << "\t" << jj << "\t" << which_triangle << "\t" << solidangle << "\t" << ii_source << "\t" << jj_source << "\t" << which_triangle_source << "\t" << distance << "\t" << solidangle_source << "\n";
				}
			}
		}
	}

	GL_file.close();
	std::cout << "[i] ViewList written to file.\n";
}

/**
* @brief Reads ViewList from file, use SMET format only
* @param[in] -
* @param[out] -
*
*/
bool TerrainRadiationComplex::ReadViewList()
{
	const std::string filename = cfg.get("Complex_ViewListFile", "Ebalance");

	if (!FileUtils::fileExists(filename))
	{
		throw NotFoundException(filename, AT);
		return false;
	}

	smet::SMETReader myreader(filename);
	std::vector<double> vec_data;
	myreader.read(vec_data);
	const size_t nr_fields = myreader.get_nr_of_fields();
	size_t dimx_file = static_cast<size_t>(myreader.get_header_doublevalue("ncols"));
	size_t dimy_file = static_cast<size_t>(myreader.get_header_doublevalue("nrows"));
	double cellsize_file = myreader.get_header_doublevalue("cellsize");
	double llx = myreader.get_header_doublevalue("easting");
	double lly = myreader.get_header_doublevalue("northing");
	size_t azimuth_points = static_cast<size_t>(myreader.get_header_doublevalue("Azimuth Points"));
	size_t elevation_points = static_cast<size_t>(myreader.get_header_doublevalue("Elevation Points"));

	M_epsilon = elevation_points;
	M_phi = azimuth_points;
	S = M_epsilon * M_phi;

	if (dimx_file != dimx)
	{
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: dimx\n";
		return false;
	}
	if (dimy_file != dimy)
	{
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: dimy\n";
		return false;
	}
	if (std::fabs(cellsize_file - dem.cellsize) >= std::numeric_limits<double>::epsilon())
	{
		std::cout.precision(std::numeric_limits<double>::digits10);
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: cellsize (got: " << cellsize_file << ", expected: " << dem.cellsize << ")\n";
		return false;
	}
	if (dimx_file != dimx)
	{
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: dimx\n";
		return false;
	}
	if (std::fabs(llx - dem.llcorner.getEasting()) >= std::numeric_limits<double>::epsilon())
	{
		std::cout.precision(std::numeric_limits<double>::digits10);
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: llx (got: " << llx << ", expected: " << dem.llcorner.getEasting() << ")\n";
		return false;
	}
	if (std::fabs(lly - dem.llcorner.getNorthing()) >= std::numeric_limits<double>::epsilon())
	{
		std::cout.precision(std::numeric_limits<double>::digits10);
		std::cout << "[E] in TerrainRadiationComplex::ReadViewList: DEM does not agree with TerrainList for field: lly (got: " << lly << ", expected: " << dem.llcorner.getNorthing() << ")\n";
		return false;
	}

	ViewList.resize(dimx, dimy, 2, S, {0, 0, 0, 0, 0});
	size_t ii_fd, jj_fd, which_triangle_fd, solidangle_fd, ii_seen_fd, jj_seen_fd, which_triangle_seen_fd, distance_fd, soldiangle_seen_fd;

	for (size_t kk = 0; kk < nr_fields; kk++)
	{
		const std::string tmp(myreader.get_field_name(kk));
		if (tmp == "ii")
			ii_fd = kk;
		if (tmp == "jj")
			jj_fd = kk;
		if (tmp == "which_triangle")
			which_triangle_fd = kk;
		if (tmp == "solidangle")
			solidangle_fd = kk;
		if (tmp == "ii_seen")
			ii_seen_fd = kk;
		if (tmp == "jj_seen")
			jj_seen_fd = kk;
		if (tmp == "which_triangle_seen")
			which_triangle_seen_fd = kk;
		if (tmp == "distance")
			distance_fd = kk;
		if (tmp == "soldiangle_seen")
			soldiangle_seen_fd = kk;
	}
	for (size_t ii = 0; ii < vec_data.size(); ii += nr_fields)
	{

		ViewList(vec_data[ii + ii_fd], vec_data[ii + jj_fd], vec_data[ii + which_triangle_fd], vec_data[ii + solidangle_fd]) = {(double)vec_data[ii + ii_seen_fd], (double)vec_data[ii + jj_seen_fd], (double)vec_data[ii + which_triangle_seen_fd], vec_data[ii + distance_fd], (double)vec_data[ii + soldiangle_seen_fd]};
	}

	std::cout << "[i] TerrainRadiationComplex: Initialized " << M_epsilon << "x" << M_phi << " ViewList from file\n";

	return true;
}

//########################################################################################################################
//                                             PUBLIC FUNCTIONS
//########################################################################################################################

/**
* @brief Computes direct, diffuse and terrain radiation for each gridpoint. Terrain radiation
* @param[in] -
* @param[out] -
*
*/
void TerrainRadiationComplex::getRadiation(mio::Array2D<double>& direct, mio::Array2D<double> &diffuse,
                                           mio::Array2D<double> &terrain, const mio::Array2D<double>
                                           &direct_unshaded_horizontal, const mio::Array2D<double>& total_ilwr,
                                           mio::Array2D<double>& sky_ilwr, mio::Array2D<double>& terrain_ilwr,
                                           double solarAzimuth, double solarElevation)
{
	MPIControl &mpicontrol = MPIControl::instance();

	if (_hasSP)
	{
		SP.setGridRadiation(albedo_grid, direct, diffuse, direct_unshaded_horizontal, solarAzimuth, solarElevation);
	}

	// Special T_Lists for Radation analysis
	mio::Array4D<double> TList_ms_old(dimx, dimy, 2, S, 0), TList_ms_new(dimx, dimy, 2, S, 0); // Total reflected radiance (W/m2/sr)  for all Vectors of basic set and all triangles of DEM
	mio::Array4D<double> TList_sky_aniso(dimx, dimy, 2, S, 0);								   // Anisotropic single-scattered radiance from sun&sky (W/m2/sr)
	mio::Array4D<double> TList_sky_iso(dimx, dimy, 2, S, 0);								   // Isotropic single-scattered radiance from sun&sky (W/m2/sr)
	mio::Array4D<double> TList_direct(dimx, dimy, 2, S, 0);									   // Anisotropic single-scattered radiance only from direct solar (W/m2/sr), (used in SolarPanel-module for shadow)

	// Direct, Diffuse, and Iterative Terrain Fux densities (triangular grid)
	mio::Array2D<double> direct_A(dimx, dimy, 0.), direct_B(dimx, dimy, 0.);					 // Direct Solar Flux density (W/m2) for triangles type A and B
	mio::Array2D<double> skydiffuse_A(dimx, dimy, 0.), skydiffuse_B(dimx, dimy, 0.);			 // Sky-diffuse Flux density (W/m2) for triangles type A and B
	mio::Array3D<double> terrain_flux_old(dimx, dimy, 2, 0), terrain_flux_new(dimx, dimy, 2, 0); // Total incident Flux density (W/m2) for all triangles of DEM (use new and old for iteration), [MT eq. 2.98]

	// Direct, Diffuse, and Terrain Fux densities (averaged to square grid)
	mio::Array2D<double> direct_temp(dimx, dimy, 0.), diffuse_temp(dimx, dimy, 0.), terrain_temp(dimx, dimy, 0.);

	Vec3D a_sun;
	getVectorSun(solarAzimuth, solarElevation, a_sun);
	Vec3D z_axis = {0, 0, 1};

// Calculate [MT eq. 2.96] and all Elements thereof
// --> Initialize direct_temp, diffuse_temp, TList_sky_aniso, TList_sky_iso, TList_direct
#pragma omp parallel for
	for (size_t ii = startx; ii < endx; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			//////////////////////////////////
			///// DIRECT [in MT eq. 2.96: F_direct,t]
			size_t solidangle_sun;
			double distance_closest_triangle;

			// Triangle-A, Geometric projection: horizontal radiation -> beam radiation -> triangle radiation
			Vec3D triangle;
			TriangleNormal(ii, jj, 1, triangle);
			if (a_sun[2] > 0 && VectorScalarProduct(triangle, a_sun) > 0)
			{
				double proj_to_ray, proj_to_triangle;

				proj_to_ray = 1. / VectorScalarProduct(a_sun, z_axis);
				proj_to_triangle = VectorScalarProduct(a_sun, triangle);

				direct_A(ii, jj) = direct_unshaded_horizontal(ii, jj) * proj_to_ray * proj_to_triangle;

				solidangle_sun = vectorToSPixel(a_sun, ii, jj, 1);
				distance_closest_triangle = ViewList(ii, jj, 1, solidangle_sun)[3];

				if (distance_closest_triangle != -999)
					direct_A(ii, jj) = 0;
			}
			else
				direct_A(ii, jj) = 0;

			// Triangle-B, Geometric projection: horizontal radiation -> beam radiation -> triangle radiation
			TriangleNormal(ii, jj, 0, triangle);
			if (a_sun[2] > 0 && VectorScalarProduct(triangle, a_sun) > 0)
			{
				double proj_to_ray, proj_to_triangle;

				proj_to_ray = 1. / VectorScalarProduct(a_sun, z_axis);
				proj_to_triangle = VectorScalarProduct(a_sun, triangle);

				direct_B(ii, jj) = direct_unshaded_horizontal(ii, jj) * proj_to_ray * proj_to_triangle;

				solidangle_sun = vectorToSPixel(a_sun, ii, jj, 0);
				distance_closest_triangle = ViewList(ii, jj, 0, solidangle_sun)[3];

				if (distance_closest_triangle != -999)
					direct_B(ii, jj) = 0;
			}
			else
				direct_B(ii, jj) = 0;

			direct_temp(ii, jj) = (direct_A(ii, jj) + direct_B(ii, jj)) / 2; // Average both triangles to one DEM-Gridpoint-Value for further Alpine3D use

			///////////////////////////////////////
			///// DIFFUSE [in MT eq. 2.96: F_diffuse,t]

			skydiffuse_A(ii, jj) = diffuse(ii, jj) * getSkyViewFactor(ii, jj, 1);
			skydiffuse_B(ii, jj) = diffuse(ii, jj) * getSkyViewFactor(ii, jj, 0);
			diffuse_temp(ii, jj) = (skydiffuse_A(ii, jj) + skydiffuse_B(ii, jj)) / 2; // Average both triangles to one DEM-Gridpoint-Value for further Alpine3D use

			/////////////////////////////////////////////
			///// PREPARE LAND [MT 2.96: Put it together]
			double direct_t, diffuse_t;
			double albedo_temp = albedo_grid(ii, jj);

			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				if (a_sun[2] < 0)
					continue; // might be to strict for Dawn..

				if (which_triangle == 1)
				{
					diffuse_t = skydiffuse_A(ii, jj);
					direct_t = direct_A(ii, jj);
				}
				else
				{
					diffuse_t = skydiffuse_B(ii, jj);
					direct_t = direct_B(ii, jj);
				}
				Vec3D triangle_normal;
				TriangleNormal(ii, jj, which_triangle, triangle_normal);
				if (VectorScalarProduct(triangle_normal, a_sun) > 0)
					solidangle_sun = vectorToSPixel(a_sun, ii, jj, which_triangle);
				else
					solidangle_sun = 0;

				// [MT eq. 2.96]: Original is TList_sky_aniso, the others are similar (used for Radation analysis in SolarPanel-Class)
				if (albedo_temp < 0.5 || !if_anisotropy)
				{ // Only Snow-Anisotropy if albedo > 0.5
					for (size_t solidangle_out = 0; solidangle_out < S; ++solidangle_out)
					{
						TList_sky_aniso(ii, jj, which_triangle, solidangle_out) = (direct_t + diffuse_t) * albedo_temp / Cst::PI;
						TList_sky_iso(ii, jj, which_triangle, solidangle_out) = TList_sky_aniso(ii, jj, which_triangle, solidangle_out);
						TList_direct(ii, jj, which_triangle, solidangle_out) = (direct_t)*albedo_temp / Cst::PI;
					}
				}
				else
				{
					for (size_t solidangle_out = 0; solidangle_out < S; ++solidangle_out)
					{
						TList_sky_aniso(ii, jj, which_triangle, solidangle_out) = (RList(solidangle_sun, solidangle_out) * direct_t + diffuse_t) * albedo_temp / Cst::PI;
						TList_sky_iso(ii, jj, which_triangle, solidangle_out) = (direct_t + diffuse_t) * albedo_temp / Cst::PI;
						TList_direct(ii, jj, which_triangle, solidangle_out) = (RList(solidangle_sun, solidangle_out) * direct_t) * albedo_temp / Cst::PI;
					}
				}
			}
		}
	}

	mpicontrol.allreduce_sum(TList_sky_aniso);

	///////////////////////////////
	///// LAND MAIN LOOP [MT eq. 2.97]

	bool another_round = 1;
	if (if_multiple == false)
		another_round = 0; // if multiple scattering is turned false, do just one iteration
	size_t number_rounds = 0;

	while (another_round && a_sun[2] > 0)
	{
		TList_ms_old = TList_ms_new + TList_sky_aniso;
		terrain_flux_old = terrain_flux_new;
		TList_ms_new = 0;
		terrain_flux_new = 0;

#pragma omp parallel for
		for (size_t ii = startx; ii < endx; ++ii)
		{
			for (size_t jj = 1; jj < dimy - 1; ++jj)
			{
				double albedo_temp = albedo_grid(ii, jj);

				for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
				{
					for (size_t solidangle_in = 0; solidangle_in < S; ++solidangle_in)
					{
						double Rad_solidangle;
						double distance_closest_triangle = ViewList(ii, jj, which_triangle, solidangle_in)[3];
						if (distance_closest_triangle == -999)
							continue;

						size_t ii_source = ViewList(ii, jj, which_triangle, solidangle_in)[0];
						size_t jj_source = ViewList(ii, jj, which_triangle, solidangle_in)[1];
						size_t which_triangle_source = ViewList(ii, jj, which_triangle, solidangle_in)[2];
						size_t solidangle_source = ViewList(ii, jj, which_triangle, solidangle_in)[4];

						Rad_solidangle = TList_ms_old(ii_source, jj_source, which_triangle_source, solidangle_source);
						terrain_flux_new(ii, jj, which_triangle) += Rad_solidangle / S * Cst::PI;

						Rad_solidangle = Rad_solidangle * albedo_temp / S;

						size_t solidangle_out = 0;
						// These are the most expensive loops... core of [MT eq. 2.97]
						if (albedo_temp < 0.5 || !if_anisotropy)
						{
							for (size_t kk = 0; kk < SortList(ii, jj, which_triangle).size(); ++kk)
							{
								solidangle_out = SortList(ii, jj, which_triangle)[kk];
								TList_ms_new(ii, jj, which_triangle, solidangle_out) += Rad_solidangle;
							}
						}
						else
						{
							for (size_t kk = 0; kk < SortList(ii, jj, which_triangle).size(); ++kk)
							{
								solidangle_out = SortList(ii, jj, which_triangle)[kk];
								TList_ms_new(ii, jj, which_triangle, solidangle_out) += Rad_solidangle * RList(solidangle_in, solidangle_out);
							}
						}
					}
				}
			}
		}

		mpicontrol.allreduce_sum(terrain_flux_new);
		mpicontrol.allreduce_sum(TList_ms_new);

		// Test if convergence is below treshold [MT eq. 2.100]
		if (number_rounds != 0 && TerrainBiggestDifference(terrain_flux_new, terrain_flux_old) < delta_F_max)
			another_round = 0;
		++number_rounds;
		std::cout << "[i] TerrainRadiationComplex: Done " << number_rounds << " Full Iteration(s). Domain Average is now " << terrain_flux_new.getMean() << " W/m2.\n";
	}
	// last (if multiple scattering) or only (if no multiple scattering) Iteration.
	// Very similar to Land-Main-Loop, but no speed up in core loop as data is prepared for SolarPanel-class
	if (a_sun[2] > 0)
	{
		TList_ms_old = TList_ms_new + TList_sky_aniso;
		terrain_flux_old = terrain_flux_new;
		TList_ms_new = 0;
		terrain_flux_new = 0;

#pragma omp parallel for
		for (size_t ii = startx; ii < endx; ++ii)
		{
			for (size_t jj = 1; jj < dimy - 1; ++jj)
			{
				double albedo_temp = albedo_grid(ii, jj);

				for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
				{
					for (size_t solidangle_in = 0; solidangle_in < S; ++solidangle_in)
					{
						double Rad_solidangle;
						double distance_closest_triangle = ViewList(ii, jj, which_triangle, solidangle_in)[3];
						if (distance_closest_triangle == -999)
							continue;

						size_t ii_source = ViewList(ii, jj, which_triangle, solidangle_in)[0];
						size_t jj_source = ViewList(ii, jj, which_triangle, solidangle_in)[1];
						size_t which_triangle_source = ViewList(ii, jj, which_triangle, solidangle_in)[2];
						size_t solidangle_source = ViewList(ii, jj, which_triangle, solidangle_in)[4];

						Rad_solidangle = TList_ms_old(ii_source, jj_source, which_triangle_source, solidangle_source);
						terrain_flux_new(ii, jj, which_triangle) += Rad_solidangle / S * Cst::PI;

						Rad_solidangle = Rad_solidangle * albedo_temp / S;

						if (_hasSP)
						{
							if (albedo_temp < 0.5 || !if_anisotropy)
							{
								for (size_t solidangle_out = 0; solidangle_out < S; ++solidangle_out)
								{
									TList_ms_new(ii, jj, which_triangle, solidangle_out) += Rad_solidangle;
								}
							}
							else
							{
								for (size_t solidangle_out = 0; solidangle_out < S; ++solidangle_out)
								{
									TList_ms_new(ii, jj, which_triangle, solidangle_out) += Rad_solidangle * RList(solidangle_in, solidangle_out);
								}
							}
						}
					}
				}
			}
		}
	}

	// If SolarPanel-module is used, send all needed data
	if (_hasSP){
		mpicontrol.allreduce_sum(TList_sky_iso);
		mpicontrol.allreduce_sum(TList_direct);
		mpicontrol.allreduce_sum(TList_ms_new);
		SP.setTLists(TList_direct, TList_sky_iso, TList_sky_aniso, TList_ms_new + TList_sky_aniso);
	}

// Average both triangles to one DEM-Gridpoint-Value for further Alpine3D use
#pragma omp parallel for
	for (size_t ii = startx; ii < endx; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			terrain_temp(ii, jj) = (terrain_flux_new(ii, jj, 1) + terrain_flux_new(ii, jj, 0)) / 2.;
		}
	}
	mpicontrol.allreduce_sum(diffuse_temp);
	mpicontrol.allreduce_sum(terrain_temp);
	diffuse = diffuse_temp;
	terrain = terrain_temp;

	std::cout << "[i] TerrainRadiationComplex converged in " << number_rounds + 1 << " Full Iteration(s).\n";
}

void TerrainRadiationComplex::setMeteo(const mio::Array2D<double> &albedo, const mio::Array2D<double> & /*ta*/)
{
	albedo_grid = albedo;
}

//########################################################################################################################
//                                               AUXILLARY FUNCTIONS
//########################################################################################################################

/**
* @brief returns normalized surface-normal vector [MT 2.1.2 Surface Generation]
* @param[in] ii_dem: easting DEM-Grid-Point
* @param[in] jj_dem: northing DEM-Grid-Point
* @param[in] which_triangle: 1 or 2
* @param[out] {n_x,n_y,n_z}
*/
void TerrainRadiationComplex::TriangleNormal(size_t ii_dem, size_t jj_dem, int which_triangle, Vec3D &v_out)
{
	double cellsize = dem.cellsize;
	Vec3D n;
	// Triangles: A=1, B=0  [MT Figure 2.2]
	if (which_triangle == 1)
	{
		Vec3D e_x = {-cellsize, 0, dem(ii_dem - 1, jj_dem) - dem(ii_dem, jj_dem)}; //[MT eq. 2.17]
		Vec3D e_y = {0, cellsize, dem(ii_dem, jj_dem + 1) - dem(ii_dem, jj_dem)};  //[MT eq. 2.18]
		VectorCrossProduct(e_y, e_x, n);										   //[MT eq. 2.21]
	}
	else
	{
		Vec3D e_x = {cellsize, 0, dem(ii_dem + 1, jj_dem) - dem(ii_dem, jj_dem)};  //[MT eq. 2.19]
		Vec3D e_y = {0, -cellsize, dem(ii_dem, jj_dem - 1) - dem(ii_dem, jj_dem)}; //[MT eq. 2.20]
		VectorCrossProduct(e_y, e_x, n);										   //[MT eq. 2.21]
	}

	return normalizeVector(n, v_out);
}

/**
* @brief returns for 2 triangles and viewing direction eighter -999 (no intersection) or the distance to intersection [MT 2.1.2 Surface Generation]
* @param[in] v_view: see [MT fig. 2.3]
* @param[in] nn, mm: origin of v_view. Why not specification of which_triangle in nn,mm? Because have identical base point (DEM-Grid_point) [MT Figure 2.2]
* @param[in] ii, jj, which_triangle: to check intersection with
* @param[out] distance
*/
double TerrainRadiationComplex::IntersectionRayTriangle(const Vec3D &v_view, size_t mm, size_t nn, size_t ii, size_t jj, size_t which_triangle)
{
	Vec3D intersection, e_x, e_y, e_xT, e_yT, n;

	double cellsize = dem.cellsize;
	double distance, P_3, r_3;
	double intersection_x, intersection_y;

	// Triangles: A=1, B=0  [MT Figure 2.2]
	if (which_triangle == 1)
	{
		//[MT eq. 2.17]
		e_x[0] = -cellsize;
		e_x[1] = 0;
		e_x[2] = dem(ii - 1, jj) - dem(ii, jj);
		//[MT eq. 2.18]
		e_y[0] = 0;
		e_y[1] = cellsize;
		e_y[2] = dem(ii, jj + 1) - dem(ii, jj);
	}
	else
	{
		//[MT eq. 2.19]
		e_x[0] = cellsize;
		e_x[1] = 0;
		e_x[2] = dem(ii + 1, jj) - dem(ii, jj);
		//[MT eq. 2.20]
		e_y[0] = 0;
		e_y[1] = -cellsize;
		e_y[2] = dem(ii, jj - 1) - dem(ii, jj);
	}
	TriangleNormal(ii, jj, which_triangle, n);

	if (VectorScalarProduct(n, v_view) > 0)
		return -999; // Must hit frontal

	Vec3D aufpunkt_ray = {mm * cellsize, nn * cellsize, dem(mm, nn) + 0.01}; // 0.01 : Place Basic Set slightly above surface to prevent trivial intersections with itself
	Vec3D balance_point_par = {ii * cellsize, jj * cellsize, dem(ii, jj)};	 // [ ~MT eq. 2.22]
	Vec3D v_diff;
	VectorDifference(aufpunkt_ray, balance_point_par, v_diff);
	P_3 = VectorScalarProduct(v_diff, n); // nominator in rhs of [MT eq. 2.25]
	r_3 = VectorScalarProduct(v_view, n); // denominator in rhs of [MT eq. 2.25]

	if (fabs(r_3) < 0.0000001)
		return -999; // Prevent overflow

	distance = -P_3 / r_3; // [MT eq. 2.25]

	if (distance < 0)
		return -999; // light cannot cross the ground

	Vec3D view_stretch, ray_view_sum;
	VectorStretch(v_view, distance, view_stretch);
	VectorSum(aufpunkt_ray, view_stretch, ray_view_sum);
	VectorDifference(ray_view_sum, balance_point_par, intersection); // [MT eq. 2.26] or better see [MT fig. 2.3]

	// Now slightl different (less formal) procedure than [MT eq. 2.27-2.30]:
	// Want elementary vector, that is normal to e_x => crossproduct
	// scalarproduct with [MT eq. 2.27] projects out intersection_x
	// same for e_y ..
	VectorCrossProduct(n, e_x, e_xT);
	VectorCrossProduct(n, e_y, e_yT);

	intersection_x = VectorScalarProduct(intersection, e_yT) / VectorScalarProduct(e_x, e_yT);
	intersection_y = VectorScalarProduct(intersection, e_xT) / VectorScalarProduct(e_y, e_xT);

	if (intersection_x >= 0 && intersection_y >= 0 && (intersection_x + intersection_y) <= 1)
		return distance; // [MT eq. 2.31]
	else
		return -999;
}

/**
* @brief returns for given triangle and vector the closest vector of the corresponding rotated basic set.
* 	Described after [MT fig. 2.5] in [MT 2.1.3 View-List]
* @param[in] vec_in: vector that is searching for closest rotated basic set-vector
* @param[in] ii, jj, which_triangle: triangle specification
* @param[out] list_index: Index of vector in Basic Set
*/
size_t TerrainRadiationComplex::vectorToSPixel(const Vec3D &vec_in, size_t ii_dem, size_t jj_dem, size_t which_triangle)
{

	double azimuth_flat;
	double delta, phi_temp = 0, phi;
	int m, n = -1;
	size_t list_index;

	Vec3D vec_horizontal;

	Vec3D triangle_normal, z_axis = {0, 0, 1};
	TriangleNormal(ii_dem, jj_dem, which_triangle, triangle_normal);
	double inclination = AngleBetween2Vectors(triangle_normal, z_axis);

	Vec3D axis;
	VectorCrossProduct(z_axis, triangle_normal, axis);
	if (NormOfVector(axis) == 0)
	{
		axis[0] = 1;
		axis[1] = 0;
		axis[2] = 0;
	}

	RotN(axis, vec_in, -inclination, vec_horizontal); // Backrotation in flat Basic Set [MT eq. 2.42]

	// If Vector is below trangle-plane
	if (vec_horizontal[2] < 0)
	{
		std::cout << "Vector lower than horizon TerrainRadiationComplex::vectorToSPixel\n";
		return -999;
	}

	Vec3D vec_projected_xy = {vec_horizontal[0], vec_horizontal[1], 0};
	azimuth_flat = AngleBetween2Vectors(vec_projected_xy, {0, 1, 0}); // [MT eq. 2.43]
	if (vec_projected_xy[0] < 0)
		azimuth_flat = 2 * Cst::PI - azimuth_flat;				  // AngleBetween2Vectors is not uniquely defined;
	phi = AngleBetween2Vectors(vec_horizontal, vec_projected_xy); // [MT eq. 2.44]

	m = (int)M_phi * (azimuth_flat) / 2 / Cst::PI; // [MT eq. 2.46]
	size_t i = 0;
	while (n == -1) // [MT eq. 2.45]
	{

		if (i != (M_epsilon - 1))
		{
			delta = acos(-2 / float(M_epsilon) + cos(2 * phi_temp)) / 2 - phi_temp;
			phi_temp += delta;
		}
		else
			phi_temp = Cst::PI / 2;
		if (phi <= phi_temp)
			n = i;
		i += 1;
		if (i > M_epsilon)
			std::cout << "[E] in TerrainRadiationComplex::vectorToSPixel\n";
	}

	list_index = n * M_phi + m; // [MT eq. 2.2] (ii starts @ 0)
	return list_index;
}

/**
* @brief returns land-view factor [MT 2.1.5 View Factors]
* @param[in] ii, jj, which_triangle: triangle specification
* @param[out] LandViewFactor
*/
double TerrainRadiationComplex::getLandViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle)
{
	return 1. - getSkyViewFactor(ii_dem, jj_dem, which_triangle);
}

void TerrainRadiationComplex::initSkyViewFactor()
{
	for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
	{
		for (size_t ii = 1; ii < dimx - 1; ++ii)
		{
			for (size_t jj = 1; jj < dimy - 1; ++jj)
			{
				sky_vf[which_triangle][ii][jj] = computeSkyViewFactor(ii, jj, which_triangle);
			}
		}
	}
	sky_vf_mean = (sky_vf[0] + sky_vf[1]) / 2;
}

/**
* @brief returns land-view factor [MT 2.1.5 View Factors]
* @param[in] ii, jj, which_triangle: triangle specification
* @param[out] LandViewFactor
*/
double TerrainRadiationComplex::computeSkyViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle)
{
	double sum = 0;

	for (size_t l = 0; l < S; ++l) // [MT eq. 2.60]
	{
		if (ViewList(ii_dem, jj_dem, which_triangle, l)[3] == -999) // [MT eq. 2.61]
		{
			sum++;
		}
	}

	return sum / S;
}

void TerrainRadiationComplex::getSkyViewFactor(mio::Array2D<double> &o_sky_vf)
{
	o_sky_vf = sky_vf_mean;
	MPIControl::instance().allreduce_sum(o_sky_vf);
}

double TerrainRadiationComplex::getSkyViewFactor(size_t ii_dem, size_t jj_dem, size_t which_triangle)
{
	return sky_vf[which_triangle][ii_dem][jj_dem];
}

/**
* @brief returns normalized vector pointing towards the sun
* @param[in] solarAzimuth: from north in degrees
* @param[in] solarElevation: from horizontal in degrees
* @param[out] v_out: normalized vector pointing towards the sun
*/
void TerrainRadiationComplex::getVectorSun(double solarAzimuth, double solarElevation, Vec3D &v_out)
{
	v_out[0] = cos(solarElevation * Cst::to_rad) * sin(solarAzimuth * Cst::to_rad);
	v_out[1] = cos(solarElevation * Cst::to_rad) * cos(solarAzimuth * Cst::to_rad);
	v_out[2] = sin(solarElevation * Cst::to_rad);
}

/**
* @brief returns the biggest point-like difference between two [MT eq. 2.100]
* @param[in] terrain_1, terrain_2: usually two Triangle-Arrays with Ground-reflected flux density
* @param[out] max_diff: biggest point-like difference
*/
double TerrainRadiationComplex::TerrainBiggestDifference(const mio::Array3D<double> &terrain_1, const mio::Array3D<double> &terrain_2)
{

	double max_diff = 0, diff;

	for (size_t ii = 1; ii < dimx - 1; ++ii)
	{
		for (size_t jj = 1; jj < dimy - 1; ++jj)
		{
			for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle) // Triangles: A=1, B=0  [MT Figure 2.2]
			{
				diff = fabs(terrain_1(ii, jj, which_triangle) - terrain_2(ii, jj, which_triangle));
				if (diff > max_diff)
					max_diff = diff;
			}
		}
	}
	return max_diff;
}

//########################################################################################################################
//                                               VECTOR OPERATIONS
//########################################################################################################################

/**
* @brief Standard vector operation: Calculates Norm
*/
double TerrainRadiationComplex::NormOfVector(const Vec3D &vec1)
{
	return sqrt(VectorScalarProduct(vec1, vec1));
}

/**
* @brief Standard vector operation: Normalization
*/
void TerrainRadiationComplex::normalizeVector(const Vec3D &vec1, Vec3D &v_out)
{
	double norm = NormOfVector(vec1);
	if (norm == 0)
	{
		std::cout << "Exception in TerrainRadiationComplex::normalizeVector : norm of (" << vec1[0] << ", " << vec1[1] << ", " << vec1[2] << ") is zero\n";
		return;
	}

	for (size_t i = 0; i < v_out.size(); ++i)
	{
		v_out[i] = vec1[i] / norm;
	}
}

/**
* @brief Standard vector operation: ScalarProduct
*/
double TerrainRadiationComplex::VectorScalarProduct(const Vec3D &vec1, const Vec3D &vec2)
{
	double sum = 0;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		sum += vec1[i] * vec2[i];
	}

	return sum;
}

/**
* @brief Standard vector operation: CrossProduct
*/
void TerrainRadiationComplex::VectorCrossProduct(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out)
{

	v_out[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	v_out[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	v_out[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

/**
* @brief Standard vector operation: Sum. Expects arrays of len 3.
*/
void TerrainRadiationComplex::VectorSum(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out)
{
	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out[i] = vec1[i] + vec2[i];
	}
}

/**
* @brief Standard vector operation: Difference. Expects arrays of len 3.
*/
void TerrainRadiationComplex::VectorDifference(const Vec3D &vec1, const Vec3D &vec2, Vec3D &v_out)
{
	Vec3D v_stretched;
	VectorStretch(vec2, -1, v_stretched);
	VectorSum(vec1, v_stretched, v_out);
}

/**
* @brief Standard vector operation: Multiplication with real number. Expects arrays of len 3.
*/
void TerrainRadiationComplex::VectorStretch(const Vec3D &vec1, double factor, Vec3D &v_out)
{
	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out[i] = vec1[i] * factor;
	}
}

/**
* @brief Standard vector operation: Rotation along given axis [MT eq. 2.41]
*/
void TerrainRadiationComplex::RotN(const Vec3D &axis, const Vec3D &vec_in, double rad, Vec3D &v_out)
{
	double c = cos(rad);
	double s = sin(rad);

	double n1 = axis[0] / NormOfVector(axis);
	double n2 = axis[1] / NormOfVector(axis);
	double n3 = axis[2] / NormOfVector(axis);

	v_out[0] = (n1 * n1 * (1 - c) + c) * vec_in[0] + (n1 * n2 * (1 - c) - n3 * s) * vec_in[1] + (n1 * n3 * (1 - c) + n2 * s) * vec_in[2];
	v_out[1] = (n2 * n1 * (1 - c) + n3 * s) * vec_in[0] + (n2 * n2 * (1 - c) + c) * vec_in[1] + (n2 * n3 * (1 - c) - n1 * s) * vec_in[2];
	v_out[2] = (n3 * n1 * (1 - c) - n2 * s) * vec_in[0] + (n3 * n2 * (1 - c) + n1 * s) * vec_in[1] + (n3 * n3 * (1 - c) + c) * vec_in[2];
}

/**
* @brief Standard vector operation: Project vector "vec_1" to plane defined by normal vector "plane_normal"
*/
void TerrainRadiationComplex::ProjectVectorToPlane(const Vec3D &vec1, const Vec3D &plane_normal, Vec3D &v_out)
{
	//std::vector<double> normal = plane_normal, v_out = vec1;
	Vec3D normal;
	normalizeVector(plane_normal, normal);
	double proj = VectorScalarProduct(normal, vec1);
	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out[i] = vec1[i] - proj * normal[i];
	}
}

/**
* @brief Standard vector operation: Angle between 2 Vectors [rad]
*/
double TerrainRadiationComplex::AngleBetween2Vectors(const Vec3D &vec1, const Vec3D &vec2)
{
	double sum = 0;
	double angle = 0;
	double norm1, norm2;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		sum += vec1[i] * vec2[i];
	}
	norm1 = NormOfVector(vec1);
	norm2 = NormOfVector(vec2);

	sum = sum / norm1 / norm2;
	if (sum > 1)
		sum = 1; // Had some issues with floating point
	if (sum < -1)
		sum = -1;

	angle = acos(sum);
	return angle;
}


void TerrainRadiationComplex::readSP()
{
	const std::string filename = cfg.get("PVPFILE", "EBalance");
	if (!FileUtils::fileExists(filename))
	{
		throw NotFoundException(filename, AT);
	}

	smet::SMETReader myreader(filename);
	std::vector<double> vec_data;
	myreader.read(vec_data);
	const size_t nr_fields = myreader.get_nr_of_fields();
	const int epsg = myreader.get_header_intvalue("epsg");
	const double smet_nodata = myreader.get_header_doublevalue("nodata");

	if (myreader.location_in_data(smet::WGS84) == true)
	{
		size_t lat_fd = IOUtils::unodata, lon_fd = IOUtils::unodata, wid_fd = IOUtils::unodata, hig_fd = IOUtils::unodata;
		size_t alt_fd = IOUtils::unodata, inc_fd = IOUtils::unodata, az_fd = IOUtils::unodata;
		for (size_t ii = 0; ii < nr_fields; ii++)
		{
			const std::string tmp(myreader.get_field_name(ii));
			if (tmp == "latitude")
				lat_fd = ii;
			if (tmp == "longitude")
				lon_fd = ii;
			if (tmp == "altitude")
				alt_fd = ii;
			if (tmp == "inclination")
				inc_fd = ii;
			if (tmp == "azimuth")
				az_fd = ii;
			if (tmp == "height")
				hig_fd = ii;
			if (tmp == "width")
				wid_fd = ii;
		}
		for (size_t ii = 0; ii < vec_data.size(); ii += nr_fields)
		{
			std::vector<double> point;
			point = {vec_data[ii + lat_fd], vec_data[ii + lon_fd], vec_data[ii + alt_fd], vec_data[ii + inc_fd], vec_data[ii + az_fd], vec_data[ii + hig_fd], vec_data[ii + wid_fd]};
			pv_points.push_back(point);
		}
	}
	else if (myreader.location_in_data(smet::EPSG) == true)
	{
		if (epsg == (int)floor(smet_nodata + 0.1))
			throw InvalidFormatException("In file \"" + filename + "\", missing EPSG code in header!", AT);

		size_t east_fd = IOUtils::unodata, north_fd = IOUtils::unodata, alt_fd = IOUtils::unodata;
		size_t inc_fd = IOUtils::unodata, az_fd = IOUtils::unodata, wid_fd = IOUtils::unodata, hig_fd = IOUtils::unodata;
		for (size_t ii = 0; ii < nr_fields; ii++)
		{
			const std::string tmp(myreader.get_field_name(ii));
			if (tmp == "easting")
				east_fd = ii;
			if (tmp == "northing")
				north_fd = ii;
			if (tmp == "altitude")
				alt_fd = ii;
			if (tmp == "inclination")
				inc_fd = ii;
			if (tmp == "azimuth")
				az_fd = ii;
			if (tmp == "height")
				hig_fd = ii;
			if (tmp == "width")
				wid_fd = ii;
		}
		if ((east_fd == IOUtils::unodata) || (north_fd == IOUtils::unodata) || (alt_fd == IOUtils::unodata))
			throw InvalidFormatException("File \"" + filename + "\" does not contain all data fields necessary for EPSG coordinates", AT);

		for (size_t ii = 0; ii < vec_data.size(); ii += nr_fields)
		{
			Coords coord_temp;
			coord_temp.setEPSG(epsg);
			coord_temp.setXY(vec_data[ii + east_fd], vec_data[ii + north_fd], vec_data[ii + alt_fd]);

			std::vector<double> point;
			//point={coord_temp.getLat(), coord_temp.getLon(), vec_data[ii+alt_fd], vec_data[ii+inc_fd], vec_data[ii+az_fd]};
			point = {coord_temp.getEasting(), coord_temp.getNorthing(), vec_data[ii + alt_fd], vec_data[ii + inc_fd], vec_data[ii + az_fd], vec_data[ii + hig_fd], vec_data[ii + wid_fd]};
			pv_points.push_back(point);
		}
	}
	else
	{
		throw InvalidFormatException("File \"" + filename + "\" does not contain expected location information in DATA section!", AT);
	}
}

void TerrainRadiationComplex::setSP(const mio::Date timestamp, const double solarAzimuth, const double solarElevation){
	SP.setSP(timestamp, solarAzimuth, solarElevation);
}

void TerrainRadiationComplex::writeSP(const unsigned int max_steps){
	SP.writeSP(max_steps);
}
//########################################################################################################################
//                                                     TEST OUTPUT
//########################################################################################################################

void TerrainRadiationComplex::PrintProgress(double percentage)
{
	char PBSTR[61] = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
	int PBWIDTH = 60;
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}
