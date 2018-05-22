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
#ifndef LIBNCPP_H
#define LIBNCPP_H

#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>

#include <netcdf.h>
#include <string>
#include <vector>

//https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_set_components.html

namespace ncpp {
	enum Dimensions {firstdimension=mio::MeteoGrids::lastparam+10, NONE=firstdimension, TIME, LATITUDE, LONGITUDE, NORTHING, EASTING, STATION, STATSTRLEN, lastdimension=STATSTRLEN};
	enum types {nc_none=0, nc_byte, nc_char, nc_short, nc_int, nc_long=nc_int, nc_float, nc_double, nc_ubyte, nc_ushort, nc_uint};
	
	std::string getParameterName(const size_t& param);
	size_t getParameterIndex(const std::string& param);
	
	typedef struct VAR_ATTR {
		VAR_ATTR() : name(), standard_name(), long_name(), units(), height(mio::IOUtils::nodata), param(mio::IOUtils::npos), type(nc_none) {};
		VAR_ATTR(const std::string& i_name) : name(i_name), standard_name(), long_name(), units(), height(mio::IOUtils::nodata), param(mio::IOUtils::npos), type(nc_none) {};
		VAR_ATTR(const size_t& prm, const std::string& str1, const double& hgt, const types& i_type)
								: name(str1), standard_name(), long_name(), units(), height(hgt), param(prm), type(i_type) {};
		VAR_ATTR(const size_t& prm, const std::string& str1, const std::string& str2, const std::string& str3, const std::string& str4, const double& hgt, const types& i_type)
								: name(str1), standard_name(str2), long_name(str3), units(str4), height(hgt), param(prm), type(i_type) {};
		std::string toString() const {std::ostringstream os; os << "["  << getParameterName(param) << " - " << name << " / " << standard_name << " / " << long_name << " , in " << units << " @ " << height << ", type=" << type << "]"; return os.str();};

		std::string name;
		std::string standard_name;
		std::string long_name;
		std::string units;
		double height;
		size_t param; //mapping to our MeteoGrids::Parameters or Dimensions
		types type;
	} var_attr;

	typedef struct NC_VARIABLE {
		NC_VARIABLE() : attributes(), dimids(), scale(1.), offset(0.), nodata(mio::IOUtils::nodata), varid(-1) {};
		NC_VARIABLE(const var_attr& attr)
							: attributes(attr), dimids(), scale(1.), offset(0.), nodata(mio::IOUtils::nodata), varid(-1) {};
		NC_VARIABLE(const var_attr& attr, const double& i_scale, const double& i_offset, const double& i_nodata, const int& i_varid)
							: attributes(attr), dimids(), scale(i_scale), offset(i_offset), nodata(i_nodata), varid(i_varid) {};
		std::string toString() const {std::ostringstream os; os << "[" << varid << " - " << "\"" << attributes.name << "\" - packing( *" << scale << ", +" << offset << "), nodata=" << nodata << " - depends on ("; for(size_t ii=0; ii<dimids.size(); ii++) os << " " << dimids[ii]; os << ") ]"; return os.str();};
		
		var_attr attributes;
		std::vector<int> dimids;  //dimensions this variable depends on
		double scale, offset, nodata;
		int varid;
	} nc_variable;
	
	typedef struct NC_DIMENSION {
			NC_DIMENSION() : name(), length(0), dimid(-1), type(mio::IOUtils::npos), isUnlimited(false) {};
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
	
	void open_file(const std::string& filename, const int& omode, int& ncid);
	void create_file(const std::string& filename, const int& cmode, int& ncid);
	void file_redef(const std::string& filename, const int& ncid);
	void end_definitions(const std::string& filename, const int& ncid);
	void close_file(const std::string& filename, const int& ncid);
	
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const int& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
	bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name);
	void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, std::string& attr_value);
	void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value);
	
	void read_data(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data);
	void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
	void readVariableMetadata(const int& ncid, ncpp::nc_variable& var, const bool& readTimeTransform=false, const double& TZ=0.);
	void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols, const size_t& pos_start, const double * const data);
	void write_data(const int& ncid, const nc_variable& var, const std::vector<double>& data, const bool& isUnlimited);
	void write_data(const int& ncid, const nc_variable& var, const std::vector<std::string>& data, const int& strMaxLen);

	double calculate_cellsize(double& factor_x, double& factor_y, const std::vector<double>& vecX, const std::vector<double>& vecY);
	double calculate_XYcellsize(double& factor_x, double& factor_y, const std::vector<double>& vecX, const std::vector<double>& vecY);
	void fill2DGrid(mio::Grid2DObject& grid, const double data[], const double& nodata, const bool& normal_Xorder=true, const bool& normal_Yorder=true);
	void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier);
	void createDimension(const int& ncid, nc_dimension& dimension, const size_t& length);
	
	//std::vector<double> read_1Dvariable(const int& ncid, const size_t& param, std::map<size_t, ncpp::nc_variable> vars, const std::map<size_t, ncpp::nc_dimension>& dimensions_map, const std::string& file_and_path);
	//size_t read_1DvariableLength(const ncpp::nc_variable& var, const std::map<size_t, ncpp::nc_dimension>& dimensions_map, const std::string& file_and_path);
} // end namespace

#endif
