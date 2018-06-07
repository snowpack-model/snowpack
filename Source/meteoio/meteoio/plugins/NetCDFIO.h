/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef NetCDFIO_H
#define NetCDFIO_H

#include <meteoio/IOInterface.h>
#include <meteoio/plugins/libncpp.h>

#include <string>

namespace mio {

class ncParameters {
	public:
		enum Mode {READ, WRITE};
		
		ncParameters(const std::string& filename, const Mode& mode, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug=false);
		
		std::pair<Date, Date> getDateRange() const;
		std::set<size_t> getParams() const;
		std::vector<Date> getTimestamps() const {return vecTime;}
		Grid2DObject read2DGrid(const size_t& param, const Date& date) const;
		Grid2DObject read2DGrid(const std::string& varname) const;
		
		void write2DGrid(const Grid2DObject& grid_in, ncpp::nc_variable& var, const Date& date);
		void write2DGrid(const Grid2DObject& grid_in, size_t param, std::string param_name, const Date& date);
		
		void writeMeteo(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx=IOUtils::npos);
		
	private:
		static std::map< std::string, std::vector<ncpp::var_attr> > initSchemasVars();
		static std::map< std::string, std::vector<ncpp::nc_dimension> > initSchemasDims();
		static std::vector<ncpp::var_attr> initUserSchemas(const Config& i_cfg);
		static std::vector<ncpp::nc_dimension> initUserDimensions(const Config& i_cfg);
		
		void initSchemaCst(const std::string& schema);
		void initFromFile(const std::string& filename, const std::string& schema);
		void initVariablesFromFile(const int& ncid, const std::string& schema_name);
		void initDimensionsFromFile(const int& ncid, const std::string& schema_name);
		void initFromSchema(const std::string& schema);
		
		Grid2DObject read2DGrid(const ncpp::nc_variable& var, const size_t& time_pos, const bool& m2mm=false, const bool& reZero=false) const;
		std::vector<Date> read_1Dvariable(const int& ncid) const;
		std::vector<double> read_1Dvariable(const int& ncid, const size_t& param) const;
		std::vector<std::string> read_1Dstringvariable(const int& ncid, const size_t& param) const;
		size_t read_1DvariableLength(const ncpp::nc_variable& var) const;
		size_t readDimension(const int& dimid) const;
		bool hasDimension(const size_t& dim) const;
		const ncpp::var_attr getSchemaAttributes(const std::string& var, const std::string& schema_name) const;
		const ncpp::var_attr getSchemaAttributes(const size_t& param, const std::string& schema_name) const;
		const ncpp::nc_dimension getSchemaDimension(const std::string& dimname, const std::string& schema_name) const;
		
		void writeGridMetadataHeader(const int& ncid, const Grid2DObject& grid_in) const;
		void writeMeteoMetadataHeader(const int& ncid, const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx) const;
		static Date getRefDate(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx);
		static void pushVar(std::vector<size_t> &nc_variables, const size_t& param);
		void addToVars(const size_t& param);
		void appendVariablesList(std::vector<size_t> &nc_variables, const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx);
		bool setAssociatedVariable(const int& ncid, const size_t& param, const Date& ref_date);
		size_t addTimestamp(const int& ncid, const Date& date);
		const std::vector<double> fillBufferForVar(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx, ncpp::nc_variable& var) const;
		static const std::vector<double> fillBufferForVar(const Grid2DObject& grid, ncpp::nc_variable& var);
		void applyUnits(Grid2DObject& grid, const std::string& units, const size_t& time_pos, const bool& m2mm) const;
		
		static std::map< std::string, std::vector<ncpp::var_attr> > schemas_vars; ///< all the variables' attributes for all schemas
		static std::map< std::string, std::vector<ncpp::nc_dimension> > schemas_dims; ///< all the dimensions' attributes for all schemas
		
		std::vector<ncpp::var_attr> user_schemas; ///< all the variables' attributes for the user defined schema
		std::vector<ncpp::nc_dimension> user_dimensions; ///< all the variables' attributes for the user defined schema
		std::map<size_t, ncpp::nc_variable> vars; ///< all the recognized variables for the selected schema_name and current file
		std::map<std::string, ncpp::nc_variable> unknown_vars; ///< all the unrecognized variables for the current file, as map< name, nc_variable>
		std::vector<Date> vecTime;
		std::vector<double> vecX, vecY;
		std::map<size_t, ncpp::nc_dimension> dimensions_map; ///< all the dimensions for the current schema, as found in the current file
		std::string file_and_path, current_schema;
		std::string coord_sys, coord_param;
		double TZ;
		double dflt_zref, dflt_uref; ///< default reference height for all data or wind data (respectively)
		double dflt_slope, dflt_azi; ///< default slope and azimuth
		double schema_nodata; ///< nodata value as defined in the schema
		int schema_dflt_type; ///< default data type as defined in the schema
		bool debug, isLatLon, force_station_dimension;
};

/**
 * @class NetCDFIO
 * @brief This plug-in allows reading and writing of NetCDF files for gridded data.
 *
 * @ingroup plugins
 */
class NetCDFIO : public IOInterface {
	public:
		NetCDFIO(const std::string& configfile);
		NetCDFIO(const NetCDFIO&);
		NetCDFIO(const Config& cfgreader);

		virtual bool list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list);
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);
		
		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name="");

	private:
		void parseInputOutputSection();
		void scanMeteoPath(const std::string& in_path, const std::string& nc_ext, std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files);
		void cleanMeteoCache(std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files);
		
		const Config cfg;
		std::vector< std::pair<std::pair<Date,Date>, ncParameters> > cache_grid_files; //cache of grid files in GRID2DPATH
		std::vector<MeteoGrids::Parameters> available_params;
		std::string in_schema, out_schema, in_grid2d_path, in_nc_ext, out_grid2d_path, grid2d_out_file;
		std::string out_meteo_path, in_meteo_path, out_meteo_file;
		double in_dflt_TZ, out_dflt_TZ;
		bool debug, out_single_file;
};

} //namespace
#endif
