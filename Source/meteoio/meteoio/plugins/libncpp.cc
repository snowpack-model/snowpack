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
#include <meteoio/plugins/libncpp.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOExceptions.h>

#include <algorithm>
#include <netcdf.h>
#include <cstring>
#include <cerrno>

using namespace std;

namespace ncpp {

/**
* @brief Set the names of the dimensions
* @return vector of names that should be in the same order as the enum
*/
inline std::vector<std::string> initDimensionNames()
{
	//the order must be the same as in the enum Dimensions
	std::vector<std::string> tmp;
	tmp.push_back("NONE"); tmp.push_back("TIME"); 
	tmp.push_back("LATITUDE"); tmp.push_back("LONGITUDE");
	tmp.push_back("NORTHING"); tmp.push_back("EASTING"); 
	tmp.push_back("STATION"); tmp.push_back("STATSTRLEN");
	tmp.push_back("DATESTRLEN");
	tmp.push_back("ZREF"); tmp.push_back("UREF");
	
	return tmp;
}
const std::vector<std::string> dimnames( initDimensionNames() );

void open_file(const std::string& filename, const int& omode, int& ncid)
{	
	const int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void create_file(const std::string& filename, const int& cmode, int& ncid)
{
	mio::ofilestream::createDirectoriesOfFile(filename.c_str());
	const int status = nc_create(filename.c_str(), cmode, &ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not create netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

/**
* @brief Add an attribute to the file pointed to by ncid.
* @details The provided attribute value will be casted to the data type that is provided as argument.
* @param[in] ncid file ID
* @param[in] varid ID of the variable this attribute belongs to
* @param[in] attr_name name of the attribute
* @param[in] attr_value value of the attribute (represented as a double)
* @param[in] data_type data type to cast the value to (according to NetCDF's <a href="https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html">external data types</a>)
*/
void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value, const int& data_type)
{
	const int status = nc_put_att_double(ncid, varid, attr_name.c_str(), data_type, 1, &attr_value);
	if (status != NC_NOERR)
		throw mio::IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

/**
* @brief Add an attribute to the file pointed to by ncid.
* @details In the target NetCDF file, the attribute will have the same type as the provided attribute value argument provided in this call
* @param[in] ncid file ID
* @param[in] varid ID of the variable this attribute belongs to
* @param[in] attr_name name of the attribute
* @param[in] attr_value value of the attribute
*/
void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value)
{
	const int status = nc_put_att_double(ncid, varid, attr_name.c_str(), NC_DOUBLE, 1, &attr_value);
	if (status != NC_NOERR)
		throw mio::IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const float& attr_value)
{
	const int status = nc_put_att_float(ncid, varid, attr_name.c_str(), NC_FLOAT, 1, &attr_value);
	if (status != NC_NOERR)
		throw mio::IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const int& attr_value)
{
	const int status = nc_put_att_int(ncid, varid, attr_name.c_str(), NC_INT, 1, &attr_value);
	if (status != NC_NOERR)
		throw mio::IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	const int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw mio::IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

/**
* @brief Add all the ACDD attributes contained in the acdd object to the file pointed to by ncid.
* @param[in] ncid file ID
* @param[in] acdd ACDD object containing all ACDD attributes (only those that have a value will be written)
*/
void writeACDDAttributes(const int& ncid, const mio::ACDD& acdd)
{
	for (auto it = acdd.cbegin(); it!=acdd.cend(); ++it) {
		const std::string name( it->second.getName() );
		const std::string value( it->second.getValue() );
		if (name.empty() || value.empty()) continue;
		
		ncpp::add_attribute(ncid, NC_GLOBAL, name, value);
	}
}

/**
* @brief Check if a variable has a given attribute
* @param[in] ncid file ID
* @param[in] varid ID of the variable those attributes should be checked
* @param[in] attr_name name of the attribute to check
* @return true if the variable has the given attribute, false otherwise
*/
bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name)
{
	size_t attr_len;
	const int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);

	if (status != NC_NOERR) return false;

	return true;
}

/**
* @brief Add attributes to non-standard variables
* @param[in,out] var variable whose attributes should be set.
*/
inline void overload_attributes(ncpp::nc_variable& var)
{
	if (var.attributes.name == "MS_SNOWPACK_RUNOFF" ) {
		var.attributes.standard_name = "Snowpack_runoff";
		var.attributes.long_name = "Runoff at bottom of snowpack (positive = downward flux)";
		var.attributes.units = "kg/m2";
	} else if (var.attributes.name == "MS_SOIL_RUNOFF" ) {
		var.attributes.standard_name = "Soil_runoff";
		var.attributes.long_name = "Runoff at bottom of soil (positive = downward flux)";
		var.attributes.units = "kg/m2";
	} else if (var.attributes.name == "RSNO") {
		var.attributes.standard_name = "Mean_snow_density";
		var.attributes.long_name = "Mean snow density";
		var.attributes.units = "kg/m3";
	} else if (var.attributes.name == "SFC_SUBL") {
		var.attributes.standard_name = "Surface_evaporation_sublimation";
		var.attributes.long_name = "Evaporation and sublimated mass from surface (positive = downward flux)";
		var.attributes.units = "kg/m2";
	}
}

/**
* @brief Write a pre-defined set of attributes for the given variable.
* @details Please note that during this call, a variable will be created, therefore the nc_variable structure will get a positive varid.
* If the variable already exists, it will return without doing anything.
* @param[in] ncid file ID
* @param[in,out] var variable whose attributes should be set.
*/
void create_variable(const int& ncid, ncpp::nc_variable& var)
{
	if (var.varid != -1) return; //the variable already exists
	const int ndims = static_cast<int>( var.dimids.size() );
	if (var.attributes.type==-1) throw mio::InvalidArgumentException("Undefined data type for variable '"+var.attributes.standard_name+"'", AT);
	const int status = nc_def_var(ncid, var.attributes.name.c_str(), var.attributes.type, ndims, &var.dimids[0], &var.varid);
	if (status != NC_NOERR) throw mio::IOException("Could not define variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	
	if (var.attributes.standard_name.empty()) overload_attributes(var);		// If no attributes, check if variable is defined in function overload_attributes.
	if (!var.attributes.standard_name.empty()) ncpp::add_attribute(ncid, var.varid, "standard_name", var.attributes.standard_name);
	if (!var.attributes.long_name.empty()) ncpp::add_attribute(ncid, var.varid, "long_name", var.attributes.long_name);
	if (!var.attributes.units.empty()) ncpp::add_attribute(ncid, var.varid, "units", var.attributes.units);
	if (var.attributes.param!=ncpp::STATION && var.attributes.param!=ncpp::TIME) ncpp::add_attribute(ncid, var.varid, "_FillValue", var.nodata, var.attributes.type);

	if (var.attributes.param==ncpp::STATION) ncpp::add_attribute(ncid, var.varid, "cf_role", "timeseries_id");
	if (var.attributes.param==ncpp::EASTING || var.attributes.param==ncpp::LONGITUDE) ncpp::add_attribute(ncid, var.varid, "axis", "X");
	if (var.attributes.param==ncpp::NORTHING || var.attributes.param==ncpp::LATITUDE) ncpp::add_attribute(ncid, var.varid, "axis", "Y");
	if (var.attributes.param==ncpp::TIME) {
		ncpp::add_attribute(ncid, var.varid, "calendar", "gregorian");
		ncpp::add_attribute(ncid, var.varid, "axis", "T");
	}
	if (var.attributes.param==mio::MeteoGrids::DEM) {
		ncpp::add_attribute(ncid, var.varid, "positive", "up");
		ncpp::add_attribute(ncid, var.varid, "axis", "Z");
	}
}

/**
* @brief Re-open the file in "definition" mode
* @param[in] filename filename to use when reporting errors
* @param[in] ncid file ID
*/
void file_redef(const std::string& filename, const int& ncid)
{
	const int status = nc_redef(ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not open define mode for file '" + filename + "': " + nc_strerror(status), AT);
}

void end_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_enddef(ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not close define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void close_file(const std::string& filename, const int& ncid)
{
	const int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

/**
* @brief Read 2D gridded data at the provided time position for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to read
* @param[in] pos time index in the file
* @param[in] nrows number of rows
* @param[in] ncols number of columns
* @param[in] pos_i dimension index of time
* @param[in] row_i dimension index of the rows
* @param[in] col_i dimension index of the columns
* @param[out] data data extracted from the file
*/
void read_data(const int& ncid, const nc_variable& var,
               const size_t& pos, const size_t& nrows, const size_t& ncols, const size_t& pos_i, const size_t& row_i, const size_t& col_i, double* data)
{
	//map dimensions for variable (start)
	size_t start[NC_MAX_VAR_DIMS] = { 0 };
	start[pos_i] = pos;
	start[col_i] = 0;
	start[row_i] = 0;

	//map dimensions for variable (count)
	size_t count[NC_MAX_VAR_DIMS];
	std::fill(count, count + NC_MAX_VAR_DIMS, 1);
	count[pos_i] = 1;		// read only one pos
	count[col_i] = ncols;
	count[row_i] = nrows;

	const int status = nc_get_vara_double(ncid, var.varid, start, count, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
* @brief Read grid point in 2D gridded data for non time dependent data
* @param[in] ncid file ID
* @param[in] var variable to read
* @param[in] row row index to read
* @param[in] col column index to read
* @param[in] row_i dimension index of the rows
* @param[in] col_i dimension index of the columns
* @param[out] data data extracted from the file
*/
void read_data_point(const int& ncid, const nc_variable& var,
               const size_t& row, const size_t& col, const size_t& row_i, const size_t& col_i, double* data)
{
	//map dimensions for variable
	size_t c[NC_MAX_VAR_DIMS] = { 0 };
	c[col_i] = col;
	c[row_i] = row;

	const int status = nc_get_var1_double(ncid, var.varid, c, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
* @brief Read grid point in 2D gridded data at the provided time position for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to read
* @param[in] pos time index in the file
* @param[in] row row index to read
* @param[in] col column index to read
* @param[in] pos_i dimension index of time
* @param[in] row_i dimension index of the rows
* @param[in] col_i dimension index of the columns
* @param[out] data data extracted from the file
*/
void read_data_point(const int& ncid, const nc_variable& var,
               const size_t& pos, const size_t& row, const size_t& col, const size_t& pos_i, const size_t& row_i, const size_t& col_i, double* data)
{
	//map dimensions for variable
	size_t c[NC_MAX_VAR_DIMS] = { 0 };
	c[pos_i] = pos;
	c[col_i] = col;
	c[row_i] = row;

	const int status = nc_get_var1_double(ncid, var.varid, c, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
* @brief Read all the data for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to read
* @param[out] data data extracted from the file
*/
void read_data(const int& ncid, const nc_variable& var, double* data)
{
	const int status = nc_get_var_double(ncid, var.varid, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
 * @brief Read all the data for a specific variable
 * @param[in] ncid file ID
 * @param[in] var variable to read
 * @param[out] data data extracted from the file
 */
void read_data(const int& ncid, const nc_variable& var, int* data)
{
	const int status = nc_get_var_int(ncid, var.varid, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
* @brief Read a pre-defined set of attributes for the given variable, from the provided file.
* @details Please note that during this call, the nc_variable structure will get a positive varid.
* If any of the predefined set of attribute does not exist, it will be silently skipped. The attributes structure
* of the variable var will then be populated by what has been read.
* @param[in] ncid file ID
* @param[in,out] var variable whose attributes should be read.
* @param[in] readTimeTransform should the time-parsing arguments (offset and scale from a reference date and units) be read? default=false
* @param[in] TZ timezone to use when/if reading a date
*/
void readVariableMetadata(const int& ncid, ncpp::nc_variable& var, const bool& readTimeTransform, const double& TZ)
{
	int nrdims;
	int dimids[NC_MAX_VAR_DIMS];

	int status = nc_inq_varid (ncid, var.attributes.name.c_str(), &var.varid);
	if (status != NC_NOERR) throw mio::IOException(nc_strerror(status), AT);
	status = nc_inq_var(ncid, var.varid, nullptr, nullptr, &nrdims, dimids, nullptr);
	if (status != NC_NOERR) throw mio::IOException(nc_strerror(status), AT);
	var.dimids.assign(dimids, dimids+nrdims);
	
	ncpp::getAttribute(ncid, var, "_FillValue", var.nodata);
	ncpp::getAttribute(ncid, var, "missing_value", var.nodata);
	ncpp::getAttribute(ncid, var, "scale_factor", var.scale);
	ncpp::getAttribute(ncid, var, "add_offset", var.offset);
	ncpp::getAttribute(ncid, var, "units", var.attributes.units);
	status = nc_inq_vartype(ncid, var.varid, &var.attributes.type);
	if (status != NC_NOERR) throw mio::IOException(nc_strerror(status), AT);
	
	if (readTimeTransform) {
		//trying to be more robust as some might misplace the time units
		if (var.attributes.units.empty()) 
			ncpp::getAttribute(ncid, var, "description", var.attributes.units);
		ncpp::getTimeTransform(var.attributes.units, TZ, var.offset, var.scale);
	}
}

/**
* @brief Write 2D gridded data at the provided time position for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to write out
* @param[in] nrows number of rows
* @param[in] ncols number of columns
* @param[in] pos time index in the file (IOUtils::npos for a variable that is not time dependent)
* @param[in] data data to write to the file
*/
void write_data(const int& ncid, const nc_variable& var, const size_t& pos, const size_t& nrows, const size_t& ncols,
                const double * const data)
{
	if (pos==mio::IOUtils::npos) { //no time dependency
		const size_t start[] = {0, 0};
		const size_t count[] = {nrows, ncols};

		const int status = nc_put_vara_double(ncid, var.varid, start, count, data);
		if (status != NC_NOERR)
			throw mio::IOException("Could not write variable '" + var.attributes.name + "': " + string(nc_strerror(status)), AT);
	} else { //time dependent variable
		const size_t start[] = {pos, 0, 0};
		const size_t count[] = {1, nrows, ncols};

		const int status = nc_put_vara_double(ncid, var.varid, start, count, data);
		if (status != NC_NOERR)
			throw mio::IOException("Could not write variable '" + var.attributes.name + "': " + string(nc_strerror(status)), AT);
	}
}

/**
* @brief Write a vector of data for a given 1D variable
* @param[in] ncid file ID
* @param[in] var variable to write out
* @param[in] data vector that has to be written
* @param[in] isUnlimited Is the variable the associated variable of an unlimited dimension? (default: false)
*/
void write_1Ddata(const int& ncid, const nc_variable& var, const std::vector<double>& data, const bool& isUnlimited)
{
	if (!isUnlimited) {
		const int status = nc_put_var_double(ncid, var.varid, &data[0]); //this call seems to only work properly with 1D data
		if (status != NC_NOERR) throw mio::IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	} else {
		//because nc_put_var_double does not work for unlimited dimensions
		const size_t start[] = {0};
		const size_t count[] = {data.size()};
		const int status = nc_put_vara_double(ncid, var.varid, start, count, &data[0]); 
		if (status != NC_NOERR) throw mio::IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	}
}

/**
* @brief Write a vector of strings for a given 1D variable
* @param[in] ncid file ID
* @param[in] var variable properties
* @param[in] data vector that has to be written
* @param[in] strMaxLen maximum length of the strings in the vector (this MUST have been defined as a dimension before)
*/
void write_1Ddata(const int& ncid, const nc_variable& var, const std::vector<std::string>& data, const int& strMaxLen)
{
	//the handling of arrays of strings is half broken in netcdf<4, therefore this hacky code below...
	for (size_t ii=0; ii<data.size(); ii++) {
		const std::string text(data[ii], 0, strMaxLen);
		const size_t start[] = {ii, 0};
		const size_t count[] = {1, text.size() + 1}; //only one record, and that many chars to write
		const int status = nc_put_vara_text(ncid, var.varid, start, count, text.c_str());
		if (status != NC_NOERR) throw mio::IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	}
}

/**
* @brief Read a given attribute from a variable (if not found, an empty string is returned)
* @param[in] ncid file ID
* @param[in] var variable properties
* @param[in] attr_name attribute name
* @param[out] attr_value attribute value as read
*/
void getAttribute(const int& ncid, const nc_variable& var, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, var.varid, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		char* value = new char[attr_len + 1]; // +1 for trailing null
		status = nc_get_att_text(ncid, var.varid, attr_name.c_str(), value);
		if (status != NC_NOERR) throw mio::IOException("Could not read attribute '" + attr_name + "' for '" + var.attributes.name + "': " + nc_strerror(status), AT);

		value[attr_len] = '\0';
		attr_value = value;
		delete[] value;
	}
}

/**
* @brief Read a given attribute from a variable (if not found, attr_value is left unchanged)
* @param[in] ncid file ID
* @param[in] var variable properties
* @param[in] attr_name attribute name
* @param[out] attr_value attribute value as read
*/
void getAttribute(const int& ncid, const nc_variable& var, const std::string& attr_name, double& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, var.varid, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		status = nc_get_att_double(ncid, var.varid, attr_name.c_str(), &attr_value);
		if (status != NC_NOERR) throw mio::IOException("Could not read attribute '" + attr_name + "' for '" + var.attributes.name + "': " + nc_strerror(status), AT);
	}
}

/**
* @brief Read a given global attribute (if not found, an empty string is returned)
* @param[in] ncid file ID
* @param[in] attr_name attribute name
* @param[out] attr_value attribute value as read
*/
void getGlobalAttribute(const int& ncid, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, NC_GLOBAL, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		char* value = new char[attr_len + 1]; // +1 for trailing null
		status = nc_get_att_text(ncid, NC_GLOBAL, attr_name.c_str(), value);
		if (status != NC_NOERR) throw mio::IOException("Could not read global attribute '" + attr_name + "': " + nc_strerror(status), AT);

		value[attr_len] = '\0';
		attr_value = value;
		delete[] value;
	}
}

void getGlobalAttribute(const int& ncid, const std::string& attr_name, int& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, NC_GLOBAL, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		status = nc_get_att_int(ncid, NC_GLOBAL, attr_name.c_str(), &attr_value);
		if (status != NC_NOERR) throw mio::IOException("Could not read global attribute '" + attr_name + "': " + nc_strerror(status), AT);
	}
}

/**
* @brief Parse a time unit specification
* @details Time is often defined as a number of intervals (hours, seconds, etc) from a reference date. This call parses such as specification string
* and return the necessary offset and multiplier compared to julian date.
* @param[in] time_units time specification string
* @param[in] i_TZ timezone to use to interpret the reference date
* @param[out] o_time_offset offset to apply to convert the packed values to julian date
* @param[out] o_time_divisor multiplier to apply to convert the packed values to julian date
*/
void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_divisor)
{
	static const double equinox_year = 365.242198781; //definition used by the NetCDF Udunits package
	
	std::vector<std::string> vecString;
	const size_t nrWords = mio::IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>5) throw mio::InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	if (vecString[0]=="years") o_time_divisor = 1./equinox_year;
	else if (vecString[0]=="months") o_time_divisor = 12./equinox_year;
	else if (vecString[0]=="days") o_time_divisor = 1.;
	else if (vecString[0]=="hours") o_time_divisor = 24.;
	else if (vecString[0]=="minutes") o_time_divisor = (24.*60.);
	else if (vecString[0]=="seconds") o_time_divisor = (24.*3600);
	else throw mio::InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	std::string ref_date_str = (nrWords>3)? vecString[2]+"T"+vecString[3] : vecString[2];
	if (nrWords==5) {
		const char first = vecString[4][0];
		if (first>=48 && first<=57) //ie starts with a number
			ref_date_str += (std::string("+")+vecString[4]);
		else 
			ref_date_str += vecString[4];
	}
	mio::Date refDate;
	if (!mio::IOUtils::convertString(refDate, ref_date_str, i_TZ))
		throw mio::InvalidArgumentException("Invalid reference date \'"+ref_date_str+"\'", AT);
	
	o_time_offset = refDate.getJulian();
}

/**
* @brief Create a new dimension
* @details If the requested dimension already exists, nothing is done
* @param[in] ncid file ID
* @param[in] dimension dimension to create
* @param[in] length length to set for this dimension (for the unlimited dimension, this will be ignored)
*/
void createDimension(const int& ncid, ncpp::nc_dimension& dimension, const size_t& length)
{
	if (dimension.dimid == -1) {
		dimension.length = length;
		const nc_type len = (dimension.isUnlimited)? NC_UNLIMITED : static_cast<int>(dimension.length);
		const int status = nc_def_dim(ncid, dimension.name.c_str(), len, &dimension.dimid);
		if (status != NC_NOERR) throw mio::IOException("Could not define dimension '" + dimension.name + "': " + nc_strerror(status), AT);
	} else {
		if (!dimension.isUnlimited && dimension.length != length)
			throw mio::InvalidArgumentException("Attempting to write an inconsistent lenght for dimension '" + dimension.name+"'", AT);
	}
}

/**
* @brief Fill a Grid2DObject with 2D gridded data as read from a NetCDF file
* @details The provided Grid2DObject must have been properly initialized before (ie proper Nx, Ny). Grids whose llcorner/urcorner have been reversed
* are properly handled by providing the normal_Xorder and/or normal_Yorder booleans (as well as any combination).
* @param[out] grid grid to populate
* @param[in] data serialized data, as read from the NetCDF file
* @param[in] nodata value that indicates nodata
* @param[in] normal_Xorder set to false if the X coordinate is reversed
* @param[in] normal_Yorder set to false if the Y coordinate is reversed
*/
void fill2DGrid(mio::Grid2DObject& grid, const double data[], const double& nodata, const bool& normal_Xorder, const bool& normal_Yorder)
{
	const size_t ncols = grid.getNx();
	const size_t nrows = grid.getNy();

	if (normal_Yorder) {
		for (size_t kk=0; kk < nrows; kk++) {
			const size_t row = kk*ncols;
			if (normal_Xorder) {
				for (size_t ll=0; ll < ncols; ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < ncols; ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (ncols -1) - ll], nodata);
			}
		}
	} else {
		for (size_t kk=0; kk < nrows; kk++) {
			const size_t row = ((nrows-1) - kk)*ncols;
			if (normal_Xorder) {
				for (size_t ll=0; ll < ncols; ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < ncols; ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (ncols -1) - ll], nodata);
			}
		}
	}
}

/**
* @brief Given a parameter index, return its associated name
* @details Since the MeteoGrids::Parameters have been extended inncpp, this method had to be redefined.
* @param[in] param parameter index to get the name for
* @return parameter name
*/
std::string getParameterName(const size_t& param)
{
	if (param==mio::IOUtils::npos) return "";
	
	if (param>=NONE) {
		if (param>lastdimension) 
			throw mio::IndexOutOfBoundsException("Trying to get name for a dimension that does not exist", AT);
		return dimnames[ param - firstdimension ];
	}
	
	return mio::MeteoGrids::getParameterName( param );
}

std::string getParameterDescription(const size_t& param)
{
	if (param==mio::IOUtils::npos || param>=NONE) return "";
	return mio::MeteoGrids::getParameterDescription( param );
}

std::string getParameterUnits(const size_t& param)
{
	if (param==mio::IOUtils::npos || param>=NONE) return "";
	return mio::MeteoGrids::getParameterUnits( param );
}

/**
* @brief Given a parameter name, return its associated index
* @details Since the MeteoGrids::Parameters have been extended inncpp, this method had to be redefined.
* @param[in] param parameter name to get the index for
* @return parameter index
*/
size_t getParameterIndex(const std::string& param)
{
	for (size_t ii=firstdimension; ii<=lastdimension; ii++) {
		if (dimnames[ii - firstdimension]==param) return ii;
	}
	
	return mio::MeteoGrids::getParameterIndex( param );
}

/**
* @brief Build a CF-1 history string (date of creation, creator, software version)
* @return string describing the file creation metadata
*/
std::string generateHistoryAttribute()
{
	mio::Date now;
	now.setFromSys();
	return now.toString(mio::Date::ISO_Z) + ", " + mio::IOUtils::getLogName() + "@" + mio::IOUtils::getHostName() + ", MeteoIO-" + mio::getLibVersion(true);
}

} //end namespace


///////////////////////////////////////////////////// Now the NC_SCHEMA class starts //////////////////////////////////////////
std::map< std::string, std::vector<ncpp::nc_dimension> > NC_SCHEMA::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncpp::var_attr> > NC_SCHEMA::schemas_vars( initSchemasVars() );

NC_SCHEMA::NC_SCHEMA(const mio::Config& cfg, const std::string& schema) 
          : user_schemas( initUserSchemas(cfg) ), user_dimensions( initUserDimensions(cfg) ), name(schema), 
            nodata(mio::IOUtils::nodata), dflt_type(NC_DOUBLE), force_station_dimension(false)
{
	initSchemaCst(schema);
}

void NC_SCHEMA::initSchemaCst(const std::string& schema)
{
	if (schema=="CF-1.6") {
		dflt_type = NC_FLOAT;
	} else if (schema=="CROCUS") {
		dflt_type = NC_FLOAT;
		nodata =  -9999999.; //CNRM-GAME nodata value
		force_station_dimension = true;
	} else if (schema=="ERA-INTERIM" || schema=="ERA5") {
		dflt_type = NC_DOUBLE;
	} else if (schema=="WRF") {
		dflt_type = NC_DOUBLE;
	} else if (schema=="AMUNDSEN") {
		dflt_type = NC_FLOAT;
	} else if (schema=="METEOCH") {
		dflt_type = NC_FLOAT;
	} else if (schema=="ECMWF") {
		throw mio::InvalidArgumentException("The ECMWF schema has been replaced by the ERA-INTERIM and the ERA5 schemas, please update your configuration file", AT);
	} else 
		throw mio::InvalidArgumentException("Unsupported NetCDF schema "+schema, AT);
}

std::map< std::string, std::vector<ncpp::nc_dimension> > NC_SCHEMA::initSchemasDims()
{
	std::map< std::string, std::vector<ncpp::nc_dimension> > results;
	std::vector<ncpp::nc_dimension> tmp;
	
	//CF-1.6 schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "surface_altitude") );
	results["CF-1.6"] = tmp;
	
	//CROCUS schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "Number_of_points") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "ZS") );
	results["CROCUS"] = tmp;
	
	//AMUNDSEN schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "lat") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "lon") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "x") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "y") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "alt") );
	results["AMUNDSEN"] = tmp;
	
	//ERA-Interim and ERA5 schemas
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "geopotential_height") );
	results["ERA-INTERIM"] = tmp;
	results["ERA5"] = tmp;
	
	//WRF schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "Time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "south_north") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "west_east") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::DATESTRLEN, "DateStrLen") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "HGT") );
	results["WRF"] = tmp;
	
	//METEOCH schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "REFERENCE_TS") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "x") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "y") );
	tmp.push_back( ncpp::nc_dimension(mio::MeteoGrids::DEM, "alt") );
	results["METEOCH"] = tmp;
	
	return results;
}

std::map< std::string, std::vector<ncpp::var_attr> > NC_SCHEMA::initSchemasVars()
{ //HACK: vars/dims should be identified based on standard_name, not name (cf1)
	std::map< std::string, std::vector<ncpp::var_attr> > results;
	std::vector<ncpp::var_attr> tmp;

	//CF-1.6 schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "", "min", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "", "degree_north", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "", "degree_east", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "projection_x_coordinate", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "projection_y_coordinate", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DEM, "orog", "surface_altitude", "height above mean sea level", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SLOPE, "slope", "slope_angle", "slope angle", "degrees from horizontal", mio::IOUtils::nodata, NC_FLOAT) ); //HACK this is not official cf!
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::AZI, "azimuth", "slope_azimuth", "slope azimuth", "degrees from north", mio::IOUtils::nodata, NC_FLOAT) ); //HACK this is not official cf!
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TA, "ta", "air_temperature", "near surface air temperature", "K", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RH, "hur", "relative_humidity", "", "1", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DW, "dw", "wind_from_direction", "", "degree", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::VW, "ws", "wind_speed", "", "m/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P, "ps", "surface_air_pressure", "near surface air pressure", "Pa", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P_SEA, "psl", "air_pressure_at_mean_sea_level", "", "Pa", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR_DIR, "iswr_dir", "direct_downwelling_shortwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR_DIFF, "iswr_diff", "diffuse_downwelling_shortwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "rsds", "surface_downwelling_shortwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RSWR, "rsus", "surface_upwelling_shortwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ILWR, "rlds", "surface_downwelling_longwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::OLWR, "rlus", "surface_upwelling_longwave_flux_in_air", "", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::HS, "snd", "surface_snow_thickness", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RSNO, "snow_density", "snow_density", "", "kg/m3", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SWE, "swe", "lwe_thickness_of_surface_snow_amount", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM, "pr", "precipitation_flux", "", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_S, "solid_precipitation_flux", "solid_precipitation_flux", "", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TSS, "ts", "surface_temperature", "", "K", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::VW_MAX, "ws_max", "wind_speed_of_gust", "", "m/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ALB, "surface_albedo", "surface_albedo", "", "1", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ROT, "ro", "surface_runoff_amount", "", "kg/m2", mio::IOUtils::nodata, NC_FLOAT) ); //proper variable name?
	results["CF-1.6"] = tmp;

	//CROCUS schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "min", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "LAT", "latitude", "latitude", "degrees_north", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "LON", "longitude", "longitude", "degrees_east", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::ZREF, "ZREF", "", "Reference_Height", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::UREF, "UREF", "", "Reference_Height_for_Wind", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DEM, "ZS", "altitude", "altitude", "m", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TA, "Tair", "", "Near Surface Air Temperature", "K", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RH, "HUMREL", "", "Relative Humidity", "%", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::QI, "Qair", "", "Near Surface Specific Humidity", "kg/kg", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::VW, "Wind", "", "Wind Speed", "m/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DW, "Wind_DIR", "", "Wind Direction", "deg", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_L, "Rainf", "", "Rainfall Rate", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_S, "Snowf", "", "Snowfall Rate", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "rsds", "", "Surface Incident total Shortwave radiation", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR_DIR, "DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR_DIFF, "SCA_SWdown", "", "Surface Incident Diffuse Shortwave Radiation", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P, "PSurf", "", "Surface Pressure", "Pa", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ILWR, "LWdown", "", "Surface Incident Longwave Radiation", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	results["CROCUS"] = tmp;
	
	//AMUNDSEN schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "h", mio::IOUtils::nodata, NC_INT) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "lat", "latitude", "latitude", "degrees_north", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "lon", "longitude", "longitude", "degrees_east", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "x", "projection_x_coordinate", "x coordinate of projection", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "y", "projection_y_coordinate", "x coordinate of projection", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DEM, "alt", "surface_altitude", "", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TA, "tas", "", "air_temperature", "K", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RH, "hurs", "", "relative_humidity", "%", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::VW, "wss", "", "wind_speed", "m/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM, "pr", "", "precipitation_flux", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "rsds", "", "surface_downwelling_shortwave_flux_in_air", "W/m2", mio::IOUtils::nodata, NC_FLOAT) );
	results["AMUNDSEN"] = tmp;

	//ERA-INTERIM and ERA5 schemas
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "h", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "latitude", "degrees", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "longitude", "degrees", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DEM, "z", "geopotential_height", "geopotential_height", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TA, "t2m", "", "2 metre temperature", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TD, "d2m", "", "2 metre dewpoint temperature", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P, "sp", "surface_air_pressure", "Surface pressure", "Pa", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P_SEA, "msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM, "tp", "", "Total precipitation", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::U, "u10", "", "10 metre U wind component", "m/s", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::V, "v10", "", "10 metre V wind component", "m/s", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SWE, "sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TSS, "skt", "", "Skin temperature", "K", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TSG, "stl1", "surface_temperature", "Soil temperature level 1", "K", mio::IOUtils::nodata, NC_DOUBLE) ); //this is from 0 to -7cm
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ALB, "al", "surface_albedo", "Albedo", "(0 - 1)", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ALB, "fal", "", "Forecast albedo", "(0 - 1)", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::RSNO, "rsn", "", "Snow density", "kg/m3", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ROT, "ro", "", "Runoff", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	results["ERA-INTERIM"] = tmp;
	//ERA-INTERIM schema
	results["ERA-INTERIM"].push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	results["ERA-INTERIM"].push_back( ncpp::var_attr(mio::MeteoGrids::ISWR_DIR, "fdir", "", "DIRect solar radiation at the surface", "J/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	results["ERA-INTERIM"].push_back( ncpp::var_attr(mio::MeteoGrids::ILWR, "strd", "", "Surface thermal radiation downwards", "J/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	//ERA5 schema
	results["ERA5"] = tmp;
	results["ERA5"].push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "msdwswrf", "", "Mean surface downward short-wave radiation flux", "W/m2", mio::IOUtils::nodata, NC_DOUBLE));
	results["ERA5"].push_back( ncpp::var_attr(mio::MeteoGrids::ILWR, "msdwlwrf", "", "Mean surface downward long-wave radiation flux", "W/m2", mio::IOUtils::nodata, NC_DOUBLE));
	

	//WRF schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "Times", "Times", "Times", "h", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "XLAT", "latitude", "LATITUDE, SOUTH IS NEGATIVE", "degree_north", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "XLONG", "longitude", "LONGITUDE, WEST IS NEGATIVE", "degree_east", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::DEM, "HGT", "Terrain Height", "Terrain Height", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::P, "PSFC", "Surface pressure", "SFC PRESSURE", "Pa", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TA, "T2", "2-meter temperature", "TEMP at 2 M", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TSG, "TSK", "surface skin temperature", "SURFACE SKIN TEMPERATURE", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::QI, "Q2", "2-meter specific humidity", "QV at 2 M", "kg/kg", 2, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ISWR, "SWDOWN", "Downward Short Wave flux at ground surface", "Downward Short Wave Flux at ground surface", "W/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ILWR, "GLW", "Downward Long Wave flux at ground surface", "Downward Long Wave flux at ground surface", "W/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ALB, "ALBEDO", "surface albedo", "ALBEDO", "-", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::ROT, "SFROFF", "Surface runoff", "SURFACE RUNOFF", "kg*m2/s", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::SWE, "SNOW", "snow water equivalent", "SNOW WATER EQUIVALENT", "kg/m2", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_PH, "SR", "", "fraction of frozen precipitation", "-", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_L, "RAINNC", "", "ACCUMULATED TOTAL GRID SCALE PRECIPITATION", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM_S, "SNOWNC", "", "ACCUMULATED TOTAL GRID SCALE SNOW AND ICE", "kg/m2/s", mio::IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::HS, "SNOWH", "Snow depth", "PHYSICAL SNOW DEPTH", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::TSS, "TSK", "Surface skin temperature", "SURFACE SKIN TEMPERATURE", "K", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::U, "U10", "10-meter wind speed", "U at 10 M", "m/s", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::V, "V10", "10-meter wind speed", "V at 10 M", "m/s", 10., NC_DOUBLE) );
	results["WRF"] = tmp;
	
	//METEOCH schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "REFERENCE_TS", "REFERENCE_TS", "REFERENCE_TS", "d", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "LAT", "latitude", "latitude", "degrees", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "LONG", "longitude", "longitude", "degrees", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", mio::IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "x", "projection_x_coordinate", "x coordinate of projection", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "y", "projection_y_coordinate", "y coordinate of projection", "m", mio::IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(mio::MeteoGrids::PSUM, "RhiresD", "", "", "mm", mio::IOUtils::nodata, NC_FLOAT) );
	results["METEOCH"] = tmp;
	
	return results;
}

//The user can provide his own variables properties as NETCDF_VAR::{param} = {name}
std::vector<ncpp::var_attr> NC_SCHEMA::initUserSchemas(const mio::Config& i_cfg)
{
	std::vector<ncpp::var_attr> results;
	const std::string user_type_str = i_cfg.get("NC_TYPE", "Input", "");
	char user_type = -1;
	if (!user_type_str.empty()) {
		if (user_type_str=="DOUBLE") user_type = NC_DOUBLE;
		else if (user_type_str=="FLOAT") user_type = NC_FLOAT;
		else if (user_type_str=="INT") user_type = NC_INT;
		else
			throw mio::InvalidArgumentException("Unknown NC_TYPE value "+user_type_str, AT);
	}
	
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_VAR::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string meteo_grid( custom_attr[ii].substr(found+1) );
		const std::string netcdf_param = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = ncpp::getParameterIndex(meteo_grid);
		if (param_index==mio::IOUtils::npos)
			throw mio::InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( ncpp::var_attr(param_index, netcdf_param, mio::IOUtils::nodata, user_type) );
	}
	
	return results;
}

//The user can provide his own dimensions properties as NETCDF_DIM::{dimension_param} = {name_in_current_file}
std::vector<ncpp::nc_dimension> NC_SCHEMA::initUserDimensions(const mio::Config& i_cfg)
{
	std::vector<ncpp::nc_dimension> results;
	
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_DIM::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string dim_str( custom_attr[ii].substr(found+1) );
		const std::string netcdf_dim = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = ncpp::getParameterIndex(dim_str);
		if (param_index==mio::IOUtils::npos || param_index<ncpp::firstdimension || param_index>ncpp::lastdimension)
			throw mio::InvalidArgumentException("Dimension '"+dim_str+"' is not a valid dimension! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( ncpp::nc_dimension( static_cast<ncpp::Dimensions>(param_index), netcdf_dim) );
	}
	
	return results;
}

//populate the dimensions_map from the selected schema
void NC_SCHEMA::initFromSchema(std::map<size_t, ncpp::nc_variable> &vars, std::map<size_t, ncpp::nc_dimension> &dimensions_map)
{
	for (size_t ii=0; ii<schemas_dims[name].size(); ii++) {
		dimensions_map[ schemas_dims[name][ii].param ] = schemas_dims[name][ii];
	}
	if (dimensions_map.count(ncpp::TIME)==0) throw mio::IOException("No TIME dimension in schema '"+name+"'", AT);
	dimensions_map[ ncpp::TIME ].isUnlimited = true;
	
	for (size_t ii=0; ii<schemas_vars[name].size(); ii++) {
		vars[ schemas_vars[name][ii].param ] = ncpp::nc_variable( schemas_vars[name][ii], nodata );
	}
}

const ncpp::var_attr NC_SCHEMA::getSchemaAttributes(const std::string& var) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].name==var) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncpp::var_attr> >::const_iterator it = schemas_vars.find( name );
	if (it==schemas_vars.end())
		throw mio::InvalidArgumentException("Invalid schema selected for NetCDF: \""+name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==var) return it->second[ii];
	}
	
	return ncpp::var_attr(var, dflt_type);
}

const ncpp::var_attr NC_SCHEMA::getSchemaAttributes(const size_t& param) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].param==param) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncpp::var_attr> >::const_iterator it = schemas_vars.find( name );
	if (it==schemas_vars.end())
		throw mio::InvalidArgumentException("Invalid schema selected for NetCDF: \""+name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].param==param) return it->second[ii];
	}
	
	return ncpp::var_attr(dflt_type);
}

const ncpp::nc_dimension NC_SCHEMA::getSchemaDimension(const std::string& dimname) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_dimensions.size(); ii++) {
		if (user_dimensions[ii].name==dimname) return user_dimensions[ii];
	}
	
	std::map< std::string, std::vector<ncpp::nc_dimension> >::const_iterator it = schemas_dims.find( name );
	if (it==schemas_dims.end())
		throw mio::InvalidArgumentException("Invalid schema selected for NetCDF: \""+name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==dimname) return it->second[ii];
	}
	
	return ncpp::nc_dimension();
}

