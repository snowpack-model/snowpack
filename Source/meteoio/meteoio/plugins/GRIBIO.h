// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef GRIBIO_H
#define GRIBIO_H

#include <meteoio/IOInterface.h>
#include <meteoio/plugins/GRIBFile.h>

#include <string>

namespace mio {

	using namespace codes;
/**
 * @class GRIBIO
 * @brief This class is used to read grib files
 *
 * @ingroup plugins
 * @author Patrick Leibersperger
 * @date   2024-04-17
 * 
 */
class GRIBIO : public IOInterface {
	public:
		GRIBIO(const std::string& configfile);
		GRIBIO(const Config& cfgreader);

		virtual bool list2DGrids(const Date& /*start*/, const Date& /*end*/, std::map<Date, std::set<size_t> >& /*list*/) {return false;}
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);

		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		const Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters

		std::string meteopath_in, grid2dpath_in, table_path;
		std::string meteo_ext, meteo_pattern,grid2d_ext, grid_2d_pattern;
		bool recursive_search;
		bool verbose;

		bool update_dem;
		double bearing_offset, latitudeOfNorthernPole, longitudeOfNorthernPole;
		bool llcorner_initialized;
		Coords llcorner;
		double cellsize, factor_x, factor_y;


		bool grid_initialized, meteo_initialized;
		GRIBTable parameter_table;
		std::vector<GRIBFile> cache_meteo, cache_grid2d;

		std::vector<Coords> vecPts; //points to use for virtual stations if METEO=GRIB

		void setOptions();
		void initTable();
		void scanPath(const std::string &in_path, const std::string &in_ext, const std::string &in_pattern, std::vector<GRIBFile> &cache);
		void initialize();
		void readStations(std::vector<Coords> &vecPoints);
		Coords getGeolocalization(double &cellsize_x, double &cellsize_y, const std::map<std::string,double> &grid_params);

		// GRIDDED DATA
		void initializeLLCorner(const std::map<std::string, double>& grid_params);
		void read2Dlevel(CodesHandlePtr &h, Grid2DObject &grid_out, const std::map<std::string, double> &grid_params );
		void processSingleMessage(Grid2DObject& dem_out ,GRIBFile& dem_file, const GRIBTable& dem_table, const MeteoGrids::Parameters& parameter);

		// METEO DATA
		bool readMeteoMeta(GRIBFile& file ,std::vector<Coords>& vecPoints, std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons);
    	bool removeDuplicatePoints(std::vector<Coords> &vecPoints, std::vector<double> &lats, std::vector<double> &lons); // remove potential duplicates. Returns true if some have been removed

};

} //namespace
#endif
