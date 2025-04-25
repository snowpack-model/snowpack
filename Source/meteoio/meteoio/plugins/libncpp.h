// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef LIBNCPP_H
#define LIBNCPP_H

#include <meteoio/plugins/libacdd.h>
#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/MeteoData.h>

#include <string>
#include <vector>

namespace ncpp {
	/// This enum expands the parameters given in mio::MeteoGrids::Parameters and adds parameters used as dimensions in NetCDF files
	enum Dimensions {firstdimension=mio::MeteoGrids::lastparam+10, NONE=firstdimension, TIME, LATITUDE, LONGITUDE, NORTHING, EASTING, STATION, STATSTRLEN, DATESTRLEN, ZREF, UREF, lastdimension=UREF};
	
	//These methods are needed by the structures defined below
	std::string getParameterName(const size_t& param);
	std::string getParameterDescription(const size_t& param);
	std::string getParameterUnits(const size_t& param);
	size_t getParameterIndex(const std::string& param);
	
	/** This structure contains the metadata associated with a NetCDF variable that are schema specific (for example, CF-1) */
	typedef struct VAR_ATTR {
		VAR_ATTR() : name(), standard_name(), long_name(), units(), height(mio::IOUtils::nodata), param(mio::IOUtils::npos), type(-1) {} //please do NOT use this constructor!
		VAR_ATTR(const int& i_type) : name(), standard_name(), long_name(), units(), height(mio::IOUtils::nodata), param(mio::IOUtils::npos), type(i_type) {}
		VAR_ATTR(const std::string& i_name, const int& i_type) : name(i_name), standard_name(), long_name(), units(), height(mio::IOUtils::nodata), param(mio::IOUtils::npos), type(i_type) {}
		VAR_ATTR(const size_t& prm, const std::string& i_name, const double& hgt, const int& i_type)
								: name(i_name), standard_name(), long_name(), units(), height(hgt), param(prm), type(i_type) {}
		VAR_ATTR(const size_t& prm, const std::string& i_name, const std::string& std_name, const std::string& lg_name, const std::string& i_units, const double& hgt, const int& i_type)
								: name(i_name), standard_name(std_name), long_name(lg_name), units(i_units), height(hgt), param(prm), type(i_type) {}
		
		bool isUndef() const {return (type==-1);}
		std::string toString() const {std::ostringstream os; os << "["  << ((param<=ncpp::lastdimension)? getParameterName(param) : "Unknown") << " - " << name << " / " << standard_name << " / " << long_name << " , in " << units << " @ " << height << ", type=" << type << "]"; return os.str();}

		std::string name; ///< variable name (it is possible to retrieve a variable by name)
		std::string standard_name; ///< somehow human-friendly, standardized description of the name
		std::string long_name; ///< non-standard but often present, longer description of the variable
		std::string units; ///< unit string representation
		double height; ///< sensor height (currently unused)
		size_t param; ///< parameter index (from Dimensions or MeteoGrids::Parameters)
		int type; ///< contain NetCDF External Data Types, -1 for "none"
	} var_attr;

	/** This structure contains the metadata associated with a NetCDF variable that are file specific as well as contains the schema specific metadata */
	typedef struct NC_VARIABLE {
		NC_VARIABLE() : attributes(), dimids(), dimid_time(0.), dimid_X(0.), dimid_Y(0.), scale(1.), offset(0.), nodata(mio::IOUtils::nodata), varid(-1) {} //please do NOT use this constructor!
		NC_VARIABLE(const int& i_type) : attributes(i_type), dimids(), dimid_time(0.), dimid_X(0.), dimid_Y(0.), scale(1.), offset(0.), nodata(mio::IOUtils::nodata), varid(-1) {}
		NC_VARIABLE(const var_attr& attr, const double& i_nodata)
							: attributes(attr), dimids(), dimid_time(0.), dimid_X(0.), dimid_Y(0.), scale(1.), offset(0.), nodata(i_nodata), varid(-1) {}
		NC_VARIABLE(const var_attr& attr, const size_t& i_dimid_time, const size_t& i_dimid_X, const size_t& i_dimid_Y, const double& i_scale, const double& i_offset, const double& i_nodata, const int& i_varid)
							: attributes(attr), dimids(), dimid_time(i_dimid_time), dimid_X(i_dimid_X), dimid_Y(i_dimid_Y), scale(i_scale), offset(i_offset), nodata(i_nodata), varid(i_varid) {}
		
		bool isUndef() const {return (attributes.isUndef());}
		std::string toString() const {std::ostringstream os; os << "["; if(varid!=-1) os<<varid; else os<<"∅"; os << " - " << "\"" << attributes.name << "\" - packing( *" << scale << ", +" << offset << "), nodata=" << nodata << " - depends on ("; for(size_t ii=0; ii<dimids.size(); ii++) os << " " << dimids[ii]; os << ") ]"; return os.str();}
		
		var_attr attributes; ///< metadata about the variable
		std::vector<int> dimids;  ///< dimensions this variable depends on
		// Store the sequence in which dimensions vary:
		size_t dimid_time;	///< dimension sequence for time variable
		size_t dimid_X;	///< dimension sequence for longitude/easting
		size_t dimid_Y;	///< dimension sequence for latitude/northing
		double scale, offset, nodata; ///< scale and offset for data packing, nodata value
		int varid; ///< variable ID, set to -1 and then to a positive value after reading/writing to/from a file
	} nc_variable;
	
	/** This structure contains the metadata associated with a NetCDF dimension */
	typedef struct NC_DIMENSION {
		NC_DIMENSION() : name(), length(0), dimid(-1), param(mio::IOUtils::npos), isUnlimited(false) {}
		NC_DIMENSION(const size_t& i_param, const std::string& i_name)
					: name(i_name), length(0), dimid(-1), param(i_param), isUnlimited(false) {}
		NC_DIMENSION(const size_t& i_param, const std::string& i_name, const size_t& len, const int& i_dimid, const bool& unlimited)
					: name(i_name), length(len), dimid(i_dimid), param(i_param), isUnlimited(unlimited) {}
		std::string toString() const {std::ostringstream os; os << getParameterName(param) << " -> [ "; if(dimid!=-1) os<<dimid; else os<<"∅"; os << " - " << name << ", length " << length; if (isUnlimited) os << ", unlimited"; os << "]"; return os.str();}

		std::string name; ///< dimension name
		size_t length; ///< dimension length (irrelevant when the dimension is "unlimited")
		int dimid; ///< dimension ID, set to -1 and then to a positive value after reading/writing to/from a file
		size_t param; ///< parameter index (from Dimensions or MeteoGrids::Parameters)
		bool isUnlimited; ///< at most, one dimension can be "unlimited"
	} nc_dimension;
	
	void open_file(const std::string& filename, const int& omode, int& ncid);
	void create_file(const std::string& filename, const int& cmode, int& ncid);
	void file_redef(const std::string& filename, const int& ncid);
	void create_variable(const int& ncid, ncpp::nc_variable& var);
	void end_definitions(const std::string& filename, const int& ncid);
	void close_file(const std::string& filename, const int& ncid);
	
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const float& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const int& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value, const int& data_type);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
	void writeACDDAttributes(const int& ncid, const mio::ACDD& acdd);
	bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name);
	void getGlobalAttribute(const int& ncid, const std::string& attr_name, std::string& attr_value);
	void getGlobalAttribute(const int& ncid, const std::string& attr_name, int& attr_value);
	void getAttribute(const int& ncid, const nc_variable& var, const std::string& attr_name, std::string& attr_value);
	void getAttribute(const int& ncid, const nc_variable& var, const std::string& attr_name, double& attr_value);
	
	void read_data(const int& ncid, const nc_variable& var, const size_t& pos, const size_t& nrows, const size_t& ncols, const size_t& pos_i, const size_t& row_i, const size_t& col_i, double* data);
	void read_data_point(const int& ncid, const nc_variable& var, const size_t& row, const size_t& col, const size_t& row_i, const size_t& col_i, double* data);
	void read_data_point(const int& ncid, const nc_variable& var, const size_t& pos, const size_t& row, const size_t& col, const size_t& pos_i, const size_t& row_i, const size_t& col_i, double* data);
	void read_data(const int& ncid, const nc_variable& var, double* data);
	void read_data(const int& ncid, const nc_variable& var, int* data);
	void readVariableMetadata(const int& ncid, ncpp::nc_variable& var, const bool& readTimeTransform=false, const double& TZ=0.);
	void write_data(const int& ncid, const nc_variable& var, const size_t& pos, const size_t& nrows, const size_t& ncols, const double * const data);
	void write_1Ddata(const int& ncid, const nc_variable& var, const std::vector<double>& data, const bool& isUnlimited=false);
	void write_1Ddata(const int& ncid, const nc_variable& var, const std::vector<std::string>& data, const int& strMaxLen);

	void fill2DGrid(mio::Grid2DObject& grid, const double data[], const double& nodata, const bool& normal_Xorder=true, const bool& normal_Yorder=true);
	void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier);
	void createDimension(const int& ncid, nc_dimension& dimension, const size_t& length);
	std::string generateHistoryAttribute();
} // end namespace

/**
 * @class NC_SCHEMA
 * @brief This class contains and handles NetCDF schemas.
 * @details As many applications have their own naming schemes, data types, units, etc we define here schemas that contain all this information
 * so generic NetCDF methods can find here all the application specific information they need to successfully read and interpret the structure and
 * data contained in a NetCDF file.
*/
//TODO: redo the whole user_schema thing: we should fill vars / dimensions with the schema, then add/overwrite with the user schema
class NC_SCHEMA {
	public:
		NC_SCHEMA(const mio::Config& cfg, const std::string& schema);
		
		void initFromSchema(std::map<size_t, ncpp::nc_variable> &vars, std::map<size_t, ncpp::nc_dimension> &dimensions_map);
		const ncpp::var_attr getSchemaAttributes(const std::string& var) const;
		const ncpp::var_attr getSchemaAttributes(const size_t& param) const;
		const ncpp::nc_dimension getSchemaDimension(const std::string& dimname) const;
		
	private:
		static std::map< std::string, std::vector<ncpp::var_attr> > initSchemasVars();
		static std::map< std::string, std::vector<ncpp::nc_dimension> > initSchemasDims();
		static std::vector<ncpp::var_attr> initUserSchemas(const mio::Config& i_cfg);
		static std::vector<ncpp::nc_dimension> initUserDimensions(const mio::Config& i_cfg);
		void initSchemaCst(const std::string& schema);
		
		static std::map< std::string, std::vector<ncpp::var_attr> > schemas_vars; ///< all the variables' attributes for all schemas
		static std::map< std::string, std::vector<ncpp::nc_dimension> > schemas_dims; ///< all the dimensions' attributes for all schemas
		
		std::vector<ncpp::var_attr> user_schemas; ///< all the variables' attributes for the user defined schema
		std::vector<ncpp::nc_dimension> user_dimensions; ///< all the variables' attributes for the user defined schema
		
	public:
		std::string name; ///< name of the current schema
		double nodata; ///< nodata value as defined in the schema
		int dflt_type; ///< default data type as defined in the schema
		bool force_station_dimension; ///< force writing a station dimension even if only one station is present
};

#endif
