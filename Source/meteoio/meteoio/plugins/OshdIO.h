// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2016 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef OSHDIO_H
#define OSHDIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {

/**
 * @class OshdIO
 * @brief This plugin reads Matlab binary files, relying on the MatIO library
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2016-03-24
 */
class OshdIO : public IOInterface {
	public:
		OshdIO(const std::string& configfile);
		OshdIO(const OshdIO&);
		OshdIO(const Config& cfgreader);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

		virtual bool list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list);
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& filename="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);

	private:
		struct file_index {
			file_index(const Date& i_date, const std::string& i_path, const std::string& i_filename, const std::string& i_run_date)
			                : date(i_date), run_date(i_run_date), path(i_path), filename(i_filename) {}
			bool operator<(const file_index& a) const {
				return date < a.date;
			}
			bool operator>(const file_index& a) const {
				return date > a.date;
			}
			std::string toString() const {
				return "<"+date.toString(Date::ISO)+" - "+run_date+" "+path+" "+filename+">";
			}
			Date date;
			std::string run_date;
			std::string path;
			std::string filename;
		};
		
		struct station_index {
			station_index(const std::string& i_id) 
			             : id(i_id), oshd_acro(), idx(IOUtils::npos) {}
			station_index(const std::string& i_id, const std::string& acro, const size_t& i_idx) 
			             : id(i_id), oshd_acro(acro), idx(i_idx) {}
			std::string id, oshd_acro;
			size_t idx;
		};
		
		void parseInputOutputSection();
		void readFromFile(const std::string& file_and_path, const Date& in_timestep, std::vector< std::vector<MeteoData> >& vecMeteo) const;
		void buildVecIdx(std::vector<std::string> vecAcro);
		void fillStationMeta();
		
		static size_t getFileIdx(const std::vector< struct file_index >& cache, const Date& start_date);
		static std::vector< struct file_index > scanMeteoPath(const std::string& meteopath_in, const bool& is_recursive);
		static double convertUnits(const double& val, const std::string& units, const size_t& param, const std::string& filename);
		
		const Config cfg;
		std::vector< struct file_index > cache_meteo_files; //cache of meteo files in METEOPATH
		std::vector< struct file_index > cache_grid_files; //cache of meteo files in METEOPATH
		std::vector<StationData> vecMeta;
		std::vector< struct station_index > statIDs;
		std::string coordin, coordinparam; //projection parameters
		std::string grid2dpath_in, in_meteopath, in_metafile;
		size_t nrMetadata; ///< How many stations have been provided in the metadata file
		bool debug; ///< write out extra information to help understand what is being read
		
		static const char* meteo_ext; //for the file naming scheme
		static std::map< std::string, MeteoGrids::Parameters > params_map; ///< parameters to extract from the files
		static std::map< MeteoGrids::Parameters, std::string > grids_map; ///< parameters to extract from the gridded data
		static const double in_dflt_TZ;     //default time zone, should be visible to matio for debugging
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

} //namespace
#endif
