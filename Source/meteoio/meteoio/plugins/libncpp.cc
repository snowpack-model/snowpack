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
#include <meteoio/IOExceptions.h>

#include <algorithm>

using namespace std;

namespace ncpp {

/**
* @brief Set the names of the dimensions
* @return vector of names that should be in the same order as the enum
*/
std::vector<std::string> initDimensionNames()
{
	//the order must be the same as in the enum Dimensions
	std::vector<std::string> tmp;
	tmp.push_back("NONE"); tmp.push_back("TIME"); 
	tmp.push_back("LATITUDE"); tmp.push_back("LONGITUDE");
	tmp.push_back("NORTHING"); tmp.push_back("EASTING"); 
	tmp.push_back("STATION"); tmp.push_back("STATSTRLEN");
	
	return tmp;
}
std::vector<std::string> dimnames( initDimensionNames() );

void open_file(const std::string& filename, const int& omode, int& ncid)
{
	const int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw mio::IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void create_file(const std::string& filename, const int& cmode, int& ncid)
{
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
	
	if (!var.attributes.standard_name.empty()) ncpp::add_attribute(ncid, var.varid, "standard_name", var.attributes.standard_name);
	if (!var.attributes.long_name.empty()) ncpp::add_attribute(ncid, var.varid, "long_name", var.attributes.long_name);
	if (!var.attributes.units.empty()) ncpp::add_attribute(ncid, var.varid, "units", var.attributes.units);
	if (var.attributes.param!=ncpp::STATION) ncpp::add_attribute(ncid, var.varid, "_FillValue", var.nodata, var.attributes.type);
	
	if (var.attributes.param==ncpp::TIME) ncpp::add_attribute(ncid, var.varid, "calendar", "gregorian");
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
* @param[in] ncols number of longitudes
* @param[out] data data extracted from the file
*/
void read_data(const int& ncid, const nc_variable& var,
               const size_t& pos, const size_t& nrows, const size_t& ncols, double*& data)
{
	const size_t start[] = {pos, 0, 0};
	const size_t count[] = {1, nrows, ncols};

	const int status = nc_get_vara_double(ncid, var.varid, start, count, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not retrieve data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
}

/**
* @brief Read all the data for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to read
* @param[out] data data extracted from the file
*/
void read_data(const int& ncid, const nc_variable& var, double*& data)
{
	const int status = nc_get_var_double(ncid, var.varid, data);
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
	status = nc_inq_var(ncid, var.varid, NULL, NULL, &nrdims, dimids, NULL);
	if (status != NC_NOERR) throw mio::IOException(nc_strerror(status), AT);
	var.dimids.assign(dimids, dimids+nrdims);
	
	ncpp::getAttribute(ncid, var, "_FillValue", var.nodata);
	ncpp::getAttribute(ncid, var, "missing_value", var.nodata);
	ncpp::getAttribute(ncid, var, "scale_factor", var.scale);
	ncpp::getAttribute(ncid, var, "add_offset", var.offset);
	ncpp::getAttribute(ncid, var, "units", var.attributes.units);
	status = nc_inq_vartype(ncid, var.varid, &var.attributes.type);
	if (status != NC_NOERR) throw mio::IOException(nc_strerror(status), AT);
	
	if (readTimeTransform)
		ncpp::getTimeTransform(var.attributes.units, TZ, var.offset, var.scale);
}

/**
* @brief Write 2D gridded data at the provided time position for a specific variable
* @param[in] ncid file ID
* @param[in] var variable to write out
* @param[in] nrows number of rows
* @param[in] ncols number of columns
* @param[in] pos time index in the file
* @param[in] data data to write to the file
*/
void write_data(const int& ncid, const nc_variable& var, const size_t& pos, const size_t& nrows, const size_t& ncols,
                const double * const data)
{
	const size_t start[] = {pos, 0, 0};
	const size_t count[] = {1, nrows, ncols};

	const int status = nc_put_vara_double(ncid, var.varid, start, count, data);
	if (status != NC_NOERR)
		throw mio::IOException("Could not write variable '" + var.attributes.name + "': " + string(nc_strerror(status)), AT);
}

/**
* @brief Write a vector of data for a given 1D variable
* @param[in] ncid file ID
* @param[in] var variable to write out
* @param[in] data vector that has to be written
* @param[in] isUnlimited Is the variable the associated variable of an unlimited dimension?
*/
void write_data(const int& ncid, const nc_variable& var, const std::vector<double>& data, const bool& isUnlimited)
{
	if (!isUnlimited) {
		const int status = nc_put_var_double(ncid, var.varid, &data[0]);
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
void write_data(const int& ncid, const nc_variable& var, const std::vector<std::string>& data, const int& strMaxLen)
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
* @brief Parse a time unit specification
* @details Time is often defined as a number of intervals (hours, seconds, etc) from a reference date. This call parses such as specification string
* and return the necessary offset and multiplier compared to julian date.
* @param[in] time_units time specification string
* @param[in] i_TZ timezone to use to interpret the reference date
* @param[out] o_time_offset offset to apply to convert the read values to julian date
* @param[out] o_time_multiplier multiplier to apply to convert the read values to julian date
*/
void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier)
{
	static const double equinox_year = 365.242198781; //definition used by the NetCDF Udunits package
	
	std::vector<std::string> vecString;
	const size_t nrWords = mio::IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>4) throw mio::InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	if (vecString[0]=="years") o_time_multiplier = equinox_year;
	else if (vecString[0]=="months") o_time_multiplier = equinox_year/12.;
	else if (vecString[0]=="days") o_time_multiplier = 1.;
	else if (vecString[0]=="hours") o_time_multiplier = 1./24.;
	else if (vecString[0]=="minutes") o_time_multiplier = 1./(24.*60.);
	else if (vecString[0]=="seconds") o_time_multiplier = 1./(24.*3600);
	else throw mio::InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	const std::string ref_date_str = (nrWords==3)? vecString[2] : vecString[2]+"T"+vecString[3];
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
		if (dimension.length != length)
			throw mio::InvalidArgumentException("Attempting to write an inconsistent lenght for dimension '" + dimension.name+"'", AT);
	}
}

double calculate_cellsize(double& factor_x, double& factor_y, const std::vector<double>& vecX, const std::vector<double>& vecY)
{
	//in order to handle swapped llcorner/urcorner, we use "fabs" everywhere
	double alpha;
	const double cntr_lat = .5*fabs(vecY.front()+vecY.back());
	const double cntr_lon = .5*fabs(vecX.front()+vecX.back());
	const double distanceX = mio::CoordsAlgorithms::VincentyDistance(cntr_lat, vecX.front(), cntr_lat, vecX.back(), alpha);
	const double distanceY = mio::CoordsAlgorithms::VincentyDistance(vecY.front(), cntr_lon, vecY.back(), cntr_lon, alpha);

	//round to 1cm precision for numerical stability
	const double cellsize_x = static_cast<double>(mio::Optim::round( distanceX / static_cast<double>(vecX.size())*100. )) / 100.;
	const double cellsize_y = static_cast<double>(mio::Optim::round( distanceY / static_cast<double>(vecY.size())*100. )) / 100.;
	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		const double cellsize = std::min(cellsize_x, cellsize_y);
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;
		return cellsize;
	}
}

double calculate_XYcellsize(double& factor_x, double& factor_y, const std::vector<double>& vecX, const std::vector<double>& vecY)
{
	const double distanceX = fabs(vecX.front() - vecX.back());
	const double distanceY = fabs(vecY.front() - vecY.back());

	//round to 1cm precision for numerical stability
	const double cellsize_x = static_cast<double>(mio::Optim::round( distanceX / static_cast<double>(vecX.size())*100. )) / 100.;
	const double cellsize_y = static_cast<double>(mio::Optim::round( distanceY / static_cast<double>(vecY.size())*100. )) / 100.;

	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		const double cellsize = std::min(cellsize_x, cellsize_y);
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;
		return cellsize;
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

/**
* @brief Given a parameter name, return its associated index
* @details Since the MeteoGrids::Parameters have been extended inncpp, this method had to be redefined.
* @param[in] param parameter name to get the index for
* @return parameter index
*/
size_t getParameterIndex(const std::string& param)
{
	for (size_t ii=firstdimension; ii<=lastdimension; ii++) {
		if (dimnames[ii]==param) return ii;
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

