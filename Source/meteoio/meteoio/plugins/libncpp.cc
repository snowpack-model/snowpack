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
#include <meteoio/plugins/libncpp.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <algorithm>

using namespace std;
using namespace mio;  // for the IOExceptions and IOUtils

namespace ncpp {

std::vector<std::string> initDimensionNames()
{
	//the order must be the same as in the enum Dimensions
	std::vector<std::string> tmp;
	tmp.push_back("NONE"); tmp.push_back("TIME"); 
	tmp.push_back("LATITUDE"); tmp.push_back("LONGITUDE");
	tmp.push_back("NORTHING"); tmp.push_back("EASTING"); tmp.push_back("STATION"); tmp.push_back("STATSTRLEN");
	
	return tmp;
}
std::vector<std::string> dimnames( initDimensionNames() );

void open_file(const std::string& filename, const int& omode, int& ncid)
{
	const int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void create_file(const std::string& filename, const int& cmode, int& ncid)
{
	const int status = nc_create(filename.c_str(), cmode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not create netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value)
{
	const int status = nc_put_att_double(ncid, varid, attr_name.c_str(), NC_DOUBLE, 1, &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const int& attr_value)
{
	const int status = nc_put_att_int(ncid, varid, attr_name.c_str(), NC_INT, 1, &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}


void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	const int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name)
{
	size_t attr_len;
	const int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);

	if (status != NC_NOERR) return false;

	return true;
}

void file_redef(const std::string& filename, const int& ncid)
{
	const int status = nc_redef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open define mode for file '" + filename + "': " + nc_strerror(status), AT);
}

void end_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_enddef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void close_file(const std::string& filename, const int& ncid)
{
	const int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

void read_data(const int& ncid, const std::string& varname, const int& varid,
               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data)
{
	const size_t start[] = {pos, 0, 0};
	const size_t count[] = {1, latlen, lonlen};

	const int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	const int status = nc_get_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

//attributes.name is used as a handle to get all the metadata from the file
//all other attributes are ignored, attributes.units is overwritten
void readVariableMetadata(const int& ncid, ncpp::nc_variable& var, const bool& readTimeTransform, const double& TZ)
{
	int nrdims;
	int dimids[NC_MAX_VAR_DIMS];

	int status = nc_inq_varid (ncid, var.attributes.name.c_str(), &var.varid);
	if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
	status = nc_inq_var(ncid, var.varid, NULL, NULL, &nrdims, dimids, NULL);
	if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
	var.dimids.assign(dimids, dimids+nrdims);
	
	ncpp::getAttribute(ncid, var.varid, var.attributes.name, "_FillValue", var.nodata);
	ncpp::getAttribute(ncid, var.varid, var.attributes.name, "missing_value", var.nodata);
	ncpp::getAttribute(ncid, var.varid, var.attributes.name, "scale_factor", var.scale);
	ncpp::getAttribute(ncid, var.varid, var.attributes.name, "add_offset", var.offset);
	ncpp::getAttribute(ncid, var.varid, var.attributes.name, "units", var.attributes.units);
	nc_type type;
	status = nc_inq_vartype(ncid, var.varid, &type);
	if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
	var.attributes.type = static_cast<ncpp::types>(type);
	
	if (readTimeTransform)
		ncpp::getTimeTransform(var.attributes.units, TZ, var.offset, var.scale);
}

void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols,
                const size_t& pos_start, const double * const data)
{
	const size_t start[] = {pos_start, 0, 0};
	const size_t count[] = {1, nrows, ncols};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR) {
		throw IOException("Could not write variable '" + varname + "': " + string(nc_strerror(status)), AT);
	}
}

void write_data(const int& ncid, const nc_variable& var, const std::vector<double>& data, const bool& isUnlimited)
{
	if (!isUnlimited) {
		const int status = nc_put_var_double(ncid, var.varid, &data[0]);
		if (status != NC_NOERR) throw IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	} else {
		//because nc_put_var_double does not work for unlimited dimensions
		const size_t start[] = {0};
		const size_t count[] = {data.size()};
		const int status = nc_put_vara_double(ncid, var.varid, start, count, &data[0]); 
		if (status != NC_NOERR) throw IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	}
}

void write_data(const int& ncid, const nc_variable& var, const std::vector<std::string>& data, const int& strMaxLen)
{
	//the handling of arrays of strings is half broken in netcdf<4, therefore this hacky code below...
	for (size_t ii=0; ii<data.size(); ii++) {
		const std::string text(data[ii], 0, strMaxLen);
		const size_t start[] = {ii, 0};
		const size_t count[] = {1, text.size() + 1}; //only one record, and that many chars to write
		const int status = nc_put_vara_text(ncid, var.varid, start, count, text.c_str());
		if (status != NC_NOERR) throw IOException("Could not write data for variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	}
}

//if the attribute is not found, an empty string is returned
void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		char* value = new char[attr_len + 1]; // +1 for trailing null
		status = nc_get_att_text(ncid, value_id, attr_name.c_str(), value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);

		value[attr_len] = '\0';
		attr_value = value;
		delete[] value;
	}
}

//if the attribute is not found, the attr_value is not changed
void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		status = nc_get_att_double(ncid, value_id, attr_name.c_str(), &attr_value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);
	}
}

void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier)
{
	static const double equinox_year = 365.242198781; //definition used by the NetCDF Udunits package
	
	std::vector<std::string> vecString;
	const size_t nrWords = IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>4) throw InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	if (vecString[0]=="years") o_time_multiplier = equinox_year;
	else if (vecString[0]=="months") o_time_multiplier = equinox_year/12.;
	else if (vecString[0]=="days") o_time_multiplier = 1.;
	else if (vecString[0]=="hours") o_time_multiplier = 1./24.;
	else if (vecString[0]=="minutes") o_time_multiplier = 1./(24.*60.);
	else if (vecString[0]=="seconds") o_time_multiplier = 1./(24.*3600);
	else throw InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	const std::string ref_date_str = (nrWords==3)? vecString[2] : vecString[2]+"T"+vecString[3];
	Date refDate;
	if (!IOUtils::convertString(refDate, ref_date_str, i_TZ))
		throw InvalidArgumentException("Invalid reference date \'"+ref_date_str+"\'", AT);
	
	o_time_offset = refDate.getJulian();
}

//check/add dimension as necessary
void createDimension(const int& ncid, ncpp::nc_dimension& dimension, const size_t& length)
{
	if (dimension.dimid == -1) {
		dimension.length = length;
		const nc_type len = (dimension.isUnlimited)? NC_UNLIMITED : static_cast<int>(dimension.length);
		const int status = nc_def_dim(ncid, dimension.name.c_str(), len, &dimension.dimid);
		if (status != NC_NOERR) throw IOException("Could not define dimension '" + dimension.name + "': " + nc_strerror(status), AT);
	} else {
		if (dimension.length != length)
			throw InvalidArgumentException("Attempting to write an inconsistent lenght for dimension '" + dimension.name+"'", AT);
	}
}

double calculate_cellsize(double& factor_x, double& factor_y, const std::vector<double>& vecX, const std::vector<double>& vecY)
{
	//in order to handle swapped llcorner/urcorner, we use "fabs" everywhere
	double alpha;
	const double cntr_lat = .5*fabs(vecY.front()+vecY.back());
	const double cntr_lon = .5*fabs(vecX.front()+vecX.back());
	const double distanceX = CoordsAlgorithms::VincentyDistance(cntr_lat, vecX.front(), cntr_lat, vecX.back(), alpha);
	const double distanceY = CoordsAlgorithms::VincentyDistance(vecY.front(), cntr_lon, vecY.back(), cntr_lon, alpha);

	//round to 1cm precision for numerical stability
	const double cellsize_x = static_cast<double>(Optim::round( distanceX / static_cast<double>(vecX.size())*100. )) / 100.;
	const double cellsize_y = static_cast<double>(Optim::round( distanceY / static_cast<double>(vecY.size())*100. )) / 100.;
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
	const double cellsize_x = static_cast<double>(Optim::round( distanceX / static_cast<double>(vecX.size())*100. )) / 100.;
	const double cellsize_y = static_cast<double>(Optim::round( distanceY / static_cast<double>(vecY.size())*100. )) / 100.;

	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		const double cellsize = std::min(cellsize_x, cellsize_y);
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;
		return cellsize;
	}
}

//populate the results grid and handle the case of llcorner/urcorner swapped
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

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
std::string getParameterName(const size_t& param)
{
	if (param==IOUtils::npos) return "";
	
	if (param>=NONE) {
		if (param>lastdimension) 
			throw IndexOutOfBoundsException("Trying to get name for a dimension that does not exist", AT);
		return dimnames[ param - firstdimension ];
	}
	
	return mio::MeteoGrids::getParameterName( param );
}

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
size_t getParameterIndex(const std::string& param)
{
	for (size_t ii=firstdimension; ii<=lastdimension; ii++) {
		if (dimnames[ii]==param) return ii;
	}
	
	return mio::MeteoGrids::getParameterIndex( param );
}

/*std::vector<double> read_1Dvariable(const int& ncid, const size_t& param, std::map<size_t, ncpp::nc_variable> vars, const std::map<size_t, ncpp::nc_dimension>& dimensions_map, const std::string& file_and_path)
{
	const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end() || it->second.varid==-1) throw InvalidArgumentException("Could not find parameter \""+ncpp::getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	const size_t length = read_1DvariableLength(it->second, dimensions_map, file_and_path);
	
	std::vector<double> results( length );
	double *data = new double[ length ];
	ncpp::read_data(ncid, it->second.attributes.name, it->second.varid, data);
	std::copy(data, data+length, results.begin());
	delete[] data;
	return results;
}

size_t read_1DvariableLength(const ncpp::nc_variable& var, const std::map<size_t, ncpp::nc_dimension>& dimensions_map, const std::string& file_and_path)
{
	if (var.dimids.size()!=1) throw InvalidArgumentException("Parameter \""+ncpp::getParameterName(var.attributes.param)+"\" in file \""+file_and_path+"\" is not a 1D variable", AT);
	
	const int dimid = var.dimids[0];
	std::map<size_t, ncpp::nc_dimension>::const_iterator it = dimensions_map.begin();
	for (; it!=dimensions_map.end(); ++it) {
		if (it->second.dimid==dimid) break;
	}
	if (it==dimensions_map.end()) throw InvalidArgumentException("Could not find a dimension in file \""+file_and_path+"\"", AT);
	
	return it->second.length;
}*/

} //end namespace

