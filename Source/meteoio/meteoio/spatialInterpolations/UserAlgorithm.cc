// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <meteoio/spatialInterpolations/UserAlgorithm.h>
#include <meteoio/FileUtils.h>

namespace mio {

USERInterpolation::USERInterpolation(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
                                                                 GridsManager& i_gdm)
                                : InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), gdm(i_gdm), filename(), grid2d_path(), subdir(), file_ext(), naming("YYYYMMDDhhmmss_PARAM"), allow_regridding(false), time_constant(false), lowest_priority(false)
{
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="SUBDIR") {
			subdir = vecArgs[ii].second;
		} else if (vecArgs[ii].first=="EXT") {
			file_ext = vecArgs[ii].second;
		} else if (vecArgs[ii].first=="NAMING") {
			naming = vecArgs[ii].second;
		} else if (vecArgs[ii].first=="REGRID") {
			const std::string where( "Interpolations2D::"+i_param+"::"+i_algo );
			IOUtils::parseArg(vecArgs[ii], where, allow_regridding);
		} else if (vecArgs[ii].first=="TIME_CONSTANT") {
			const std::string where( "Interpolations2D::"+i_param+"::"+i_algo );
			IOUtils::parseArg(vecArgs[ii], where, time_constant);
		} else if (vecArgs[ii].first=="LOWEST_PRIORITY") {
			const std::string where( "Interpolations2D::"+i_param+"::"+i_algo );
			IOUtils::parseArg(vecArgs[ii], where, lowest_priority);
		}
	}

	if (!subdir.empty()) subdir += "/";
	if (file_ext.empty()) file_ext = ".asc";

	gdm.getConfig().getValue("GRID2DPATH", "Input", grid2d_path);
}

double USERInterpolation::getQualityRating(const Date& i_date)
{
	if (time_constant) {
		filename = subdir + param + file_ext;
	} else {
		date = i_date;
		const std::string timenum = date.toString(Date::NUM);
		std::string timestr = naming;
		bool invalidnaming = false;
		size_t pos = 0;
		// Insert year
		pos = naming.find("YYYY");
		if (pos!=std::string::npos)
			timestr.replace(pos, 4, timenum.substr(0, 4));
		else
			invalidnaming = true;
		// Insert month
		pos = naming.find("MM");
		if (pos!=std::string::npos)
			timestr.replace(pos, 2, timenum.substr(4, 2));
		else
			invalidnaming = true;
		// Insert day
		pos = naming.find("DD");
		if (pos!=std::string::npos)
			timestr.replace(pos, 2, timenum.substr(6, 2));
		else
			invalidnaming = true;
		// Insert hours
		pos = naming.find("hh");
		if (pos!=std::string::npos)
			timestr.replace(pos, 2, timenum.substr(8, 2));
		else
			invalidnaming = true;
		// Insert minutes
		pos = naming.find("mm");
		if (pos!=std::string::npos)
			timestr.replace(pos, 2, timenum.substr(10, 2));
		else
			invalidnaming = true;
		// Insert seconds (optional)
		pos = naming.find("ss");
		if (pos!=std::string::npos)
			timestr.replace(pos, 2, timenum.substr(12, 2));
		// Insert param
		pos = naming.find("PARAM");
		if (pos!=std::string::npos)
			timestr.replace(pos, 5, param);
		else
			invalidnaming = true;

		if (invalidnaming) throw InvalidArgumentException("[E] Invalid NAMING key (" + naming + ") for "+algo+" interpolation algorithm.\n", AT);

		filename = subdir + timestr + file_ext;
	}

	if (!FileUtils::validFileAndPath(grid2d_path+"/"+filename)) {
		std::cerr << "[E] Invalid grid filename for "+algo+" interpolation algorithm: " << grid2d_path+"/"+filename << "\n";
		return 0.0;
	}

	const bool has_data = FileUtils::fileExists(grid2d_path+"/"+filename);

	if (!lowest_priority)
		return (has_data)? 1. : 0.;
	else
		return (has_data)? 1e-6 : 0.;
}

double USERInterpolation::bilinear_pixel(const Array2D<double> &i_grid, const size_t &org_ii, const size_t &org_jj, const size_t &org_nx, const size_t &org_ny, const double &x, const double &y)
{
	// From: LibResampling2D
	if (org_jj>=(org_ny-1) || org_ii>=(org_nx-1)) return i_grid(org_ii, org_jj);

	const double f_0_0 = i_grid(org_ii, org_jj);
	const double f_1_0 = i_grid(org_ii+1, org_jj);
	const double f_0_1 = i_grid(org_ii, org_jj+1);
	const double f_1_1 = i_grid(org_ii+1, org_jj+1);

	double avg_value = 0.;
	unsigned int avg_count = 0;
	if (f_0_0!=IOUtils::nodata) {
		avg_value += f_0_0;
		avg_count++;
	}
	if (f_1_0!=IOUtils::nodata) {
		avg_value += f_1_0;
		avg_count++;
	}
	if (f_0_1!=IOUtils::nodata) {
		avg_value += f_0_1;
		avg_count++;
	}
	if (f_1_1!=IOUtils::nodata) {
		avg_value += f_1_1;
		avg_count++;
	}

	if (avg_count==4) return f_0_0 * (1.-x)*(1.-y) + f_1_0 * x*(1.-y) + f_0_1 * (1.-x)*y + f_1_1 *x*y;

	//special cases: less than two neighbours or three neighbours
	if (avg_count<=2) return IOUtils::nodata;

	double value = 0.;
	const double avg = avg_value/(double)avg_count;
	if (f_0_0!=IOUtils::nodata) value += f_0_0 * (1.-x)*(1.-y);
	else value += avg * (1.-x)*(1.-y);
	if (f_1_0!=IOUtils::nodata) value += f_1_0 * x*(1.-y);
	else value += avg * x*(1.-y);
	if (f_0_1!=IOUtils::nodata) value += f_0_1 * (1.-x)*y;
	else value += avg * (1.-x)*y;
	if (f_1_1!=IOUtils::nodata) value += f_1_1 *x*y;
	else value += avg *x*y;

	return value;
}


Grid2DObject USERInterpolation::resample(const DEMObject& dem, const Grid2DObject& i_grid)
{
	// Initialize output grid
	Grid2DObject o_grid(dem, 0.);
	// Get relevant georeferencing info
	const double o_grid_ll_east = o_grid.llcorner.getEasting();
	const double o_grid_ll_north = o_grid.llcorner.getNorthing();
	const double o_grid_cz = o_grid.cellsize;
	const double i_grid_ll_east = i_grid.llcorner.getEasting();
	const double i_grid_ll_north = i_grid.llcorner.getNorthing();
	const double i_grid_cz = i_grid.cellsize;
	// Loop over coordinates in output grid
	for (size_t jj=0; jj<o_grid.getNy(); jj++) {
		// map y-coordinate from o_grid in i_grid
		const double j2 = ((jj*o_grid_cz + o_grid_ll_north) - i_grid_ll_north) / i_grid_cz;
		const int j2i = (int)j2;
		const double y = j2 - (double)j2i; //normalized y, between 0 and 1
		for (size_t ii=0; ii<o_grid.getNx(); ii++) {
			// map x-coordinate from o_grid in i_grid
			const double i2 = ((ii*o_grid_cz + o_grid_ll_east) - i_grid_ll_east) / i_grid_cz;
			const int i2i = (int)i2;
			if (i2i < 0 || j2i < 0 || i2i > i_grid.getNx()-1 || j2i > i_grid.getNy()-1) {
				o_grid(ii,jj) = IOUtils::nodata;
			} else {
				const double x = i2 - (double)i2i; //normalized x, between 0 and 1
				o_grid(ii,jj) = bilinear_pixel(i_grid.grid2D, i2, j2, i_grid.getNx(), i_grid.getNy(), x, y);
			}
		}
	}
	info << FileUtils::getFilename(filename) << ", regridded from " << i_grid.getNx() << "x" << i_grid.getNy() << " to " << o_grid.getNx() << "x" << o_grid.getNy() << ", cellsize " << i_grid_cz << " to " << o_grid_cz << ", ll (" << i_grid_ll_east << "," << i_grid_ll_north << ") to (" << o_grid_ll_east << "," << o_grid_ll_north << ") to match DEM!";
	return o_grid;
}

void USERInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	gdm.read2DGrid(grid, filename);
	if (!grid.isSameGeolocalization(dem)) {
		if (allow_regridding) {
			if (param == "DW") {
				Grid2DObject gridU(grid, IOUtils::nodata);
				Grid2DObject gridV(grid, IOUtils::nodata);
				for (size_t ii=0; ii<grid.size(); ii++) {
					if (grid(ii) != IOUtils::nodata) {
						gridU(ii) = IOUtils::VWDW_TO_U(1., grid(ii));
						gridV(ii) = IOUtils::VWDW_TO_V(1., grid(ii));
					}
				}
				info << "U: ";
				gridU = resample(dem, gridU);
				info << " V: ";
				gridV = resample(dem, gridV);
				grid.set(dem, 0.);
				for (size_t ii=0; ii<grid.size(); ii++) {
					if (gridU(ii) != IOUtils::nodata && gridV(ii) != IOUtils::nodata) {
						grid(ii) = IOUtils::UV_TO_DW(gridU(ii), gridV(ii));
					}
				}
			} else {
				grid = resample(dem, grid);
			}
		} else {
			throw InvalidArgumentException("[E] trying to load a grid(" + filename + ") that does not have the same georeferencing as the DEM!", AT);
		}
	} else {
		info << FileUtils::getFilename(filename);
	}
}

} //namespace
