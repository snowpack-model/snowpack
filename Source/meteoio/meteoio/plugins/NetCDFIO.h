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

#include <string>

namespace mio {
//https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_set_components.html

class ncParameters {
	public:
		enum Mode {READ, WRITE};
		enum Dimensions {firstdimension=MeteoGrids::AZI+10, NONE=firstdimension, TIME, LATITUDE, LONGITUDE, NORTHING, EASTING, STATION, lastdimension=STATION};
		
		typedef struct VAR_ATTR {
			VAR_ATTR() : name(), standard_name(), long_name(), units(), height(IOUtils::nodata), param(IOUtils::npos) {};
			VAR_ATTR(const std::string& i_name) : name(i_name), standard_name(), long_name(), units(), height(IOUtils::nodata), param(IOUtils::npos) {};
			VAR_ATTR(const size_t& prm, const std::string& str1, const double& hgt)
			                     : name(str1), standard_name(), long_name(), units(), height(hgt), param(prm) {};
			VAR_ATTR(const size_t& prm, const std::string& str1, const std::string& str2, const std::string& str3, const std::string& str4, const double& hgt)
			                     : name(str1), standard_name(str2), long_name(str3), units(str4), height(hgt), param(prm) {};
			std::string toString() const {std::ostringstream os; os << "["  << getParameterName(param) << " - " << name << " / " << standard_name << " / " << long_name << " , in " << units << " @ " << height << "]"; return os.str();};
  
			std::string name;
			std::string standard_name;
			std::string long_name;
			std::string units;
			double height;
			size_t param; //mapping to our MeteoGrids::Parameters or Dimensions
		} var_attr;

		typedef struct NC_VARIABLE {
			NC_VARIABLE() : attributes(), dimids(), scale(1.), offset(0.), nodata(IOUtils::nodata), varid(-1) {};
			NC_VARIABLE(const var_attr& attr)
			                   : attributes(attr), dimids(), scale(1.), offset(0.), nodata(IOUtils::nodata), varid(-1) {};
			NC_VARIABLE(const var_attr& attr, const double& i_scale, const double& i_offset, const double& i_nodata, const int& i_varid)
			                   : attributes(attr), dimids(), scale(i_scale), offset(i_offset), nodata(i_nodata), varid(i_varid) {};
			std::string toString() const {std::ostringstream os; os << "[" << varid << " - " << "\"" << attributes.name << "\" - packing( *" << scale << ", +" << offset << "), nodata=" << nodata << " - depends on ("; for(size_t ii=0; ii<dimids.size(); ii++) os << " " << dimids[ii]; os << ") ]"; return os.str();};
			
			var_attr attributes;
			std::vector<int> dimids;  //dimensions this variable depends on
			double scale, offset, nodata;
			int varid;
		} nc_variable;
		
		ncParameters(const std::string& filename, const Mode& mode, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug=false);
		
		std::pair<Date, Date> getDateRange() const;
		std::set<size_t> getParams() const;
		std::vector<Date> getTimestamps() const {return vecTime;}
		Grid2DObject read2DGrid(const size_t& param, const Date& date) const;
		Grid2DObject read2DGrid(const std::string& varname) const;
		
		void write2DGrid(const Grid2DObject& grid_in, nc_variable& var, const Date& date);
		void write2DGrid(Grid2DObject grid_in, const size_t& param, const Date& date);
		
		void writeMeteo(const std::vector< std::vector<MeteoData> >& vecMeteo);
		
	private:
		typedef struct NC_DIMENSION {
			NC_DIMENSION() : name(), length(0), dimid(-1), type(NONE), isUnlimited(false) {};
			NC_DIMENSION(const size_t& i_type, const std::string& i_name)
			                     : name(i_name), length(0), dimid(-1), type(i_type), isUnlimited(false) {};
			NC_DIMENSION(const size_t& i_type, const std::string& i_name, const size_t& len, const int& i_dimid, const bool& unlimited)
			                     : name(i_name), length(len), dimid(i_dimid), type(i_type), isUnlimited(unlimited) {};
			std::string toString() const {std::ostringstream os; os << getParameterName(type) << " -> [ " << dimid << " - " << name << ", length " << length; if (isUnlimited) os << ", unlimited"; os << "]"; return os.str();};
			
			std::string name;
			size_t length;
			int dimid;
			size_t type;
			bool isUnlimited;
		} nc_dimension;
		
		static std::vector<std::string> initDimensionNames();
		static std::map< std::string, std::vector<ncParameters::var_attr> > initSchemasVars();
		static std::map< std::string, std::vector<ncParameters::nc_dimension> > initSchemasDims();
		static std::vector<ncParameters::var_attr> initUserSchemas(const Config& i_cfg);
		static std::vector<ncParameters::nc_dimension> initUserDimensions(const Config& i_cfg);
		static void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, std::string& attr_value);
		static void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value);
		static void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier);
		static std::string getParameterName(const size_t& param);
		static size_t getParameterIndex(const std::string& param);
		
		void initFromFile(const std::string& filename, const std::string& schema);
		void initVariableFromFile(const int& ncid, nc_variable& var) const;
		void initVariablesFromFile(const int& ncid, const std::string& schema_name);
		void initDimensionsFromFile(const int& ncid, const std::string& schema_name);
		void initFromSchema(const std::string& schema);
		
		Grid2DObject read2DGrid(const nc_variable& var, const size_t& time_pos, const bool& m2mm=false, const bool& reZero=false) const;
		std::vector<Date> read_1Dvariable(const int& ncid) const;
		std::vector<double> read_1Dvariable(const int& ncid, const size_t& param) const;
		size_t read_1DvariableLength(const nc_variable& var) const;
		bool hasDimension(const size_t& dim) const;
		const ncParameters::var_attr getSchemaAttributes(const std::string& var, const std::string& schema_name) const;
		const ncParameters::nc_dimension getSchemaDimension(const std::string& dimname, const std::string& schema_name) const;
		double calculate_cellsize(double& factor_x, double& factor_y) const;
		double calculate_XYcellsize(double& factor_x, double& factor_y) const;
		void fill2DGrid(Grid2DObject& grid, const double data[], const double& nodata) const;
		
		size_t addTimestamp(const int& ncid, const Date& date);
		void fill_SpatialDimensions(const int& ncid, const Grid2DObject& grid_in);
		void fill_SpatialDimensions(const int& ncid, const std::vector< std::vector<MeteoData> >& vecMeteo);
		bool create_Dimension(const int& ncid, const size_t& param, const size_t& length);
		bool create_TimeDimension(const int& ncid, const Date& date, const size_t& length);
		static void create_variable(const int& ncid, nc_variable& var);
		
		static std::vector<std::string> dimnames;
		static std::map< std::string, std::vector<ncParameters::var_attr> > schemas_vars; ///< all the variables' attributes for all schemas
		static std::map< std::string, std::vector<ncParameters::nc_dimension> > schemas_dims; ///< all the dimensions' attributes for all schemas
		
		std::vector<ncParameters::var_attr> user_schemas; ///< all the variables' attributes for the user defined schema
		std::vector<ncParameters::nc_dimension> user_dimensions; ///< all the variables' attributes for the user defined schema
		std::map<size_t, nc_variable> vars; ///< all the recognized variables for the selected schema_name and current file
		std::map<std::string, nc_variable> unknown_vars; ///< all the unrecognized variables for the current file, as map< name, nc_variable>
		std::vector<Date> vecTime;
		std::vector<double> vecX, vecY;
		std::map<size_t, nc_dimension> dimensions_map; ///< all the dimensions for the current schema, as found in the current file
		std::string file_and_path;
		std::string coord_sys, coord_param;
		double TZ;
		bool wrf_hacks, debug, isLatLon;
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
		bool dem_altimeter, debug;
};

} //namespace
#endif
