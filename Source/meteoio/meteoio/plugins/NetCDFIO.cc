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
#include <meteoio/plugins/NetCDFIO.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/FileUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/plugins/libncpp.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page netcdf NetCDF
 * @section netcdf_format Format
 * In order to promote creation, access and sharing of scientific data, the NetCDF format has been
 * created as a machine-independent format. NetCDF (network Common Data Form) is therefore an interface
 * for array-oriented data access and a library that provides an implementation of the interface. The
 * <A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF software</A> was developed
 * at the <A HREF="http://www.unidata.ucar.edu/">Unidata Program Center</A> in Boulder, Colorado.
 * In order to graphicaly explore the content and structure of NetCDF files, you can use the
 * <A HREF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse</A> java software or
 * <A HREF="http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</A>. It is also possible to run *ncdump* on a given
 * file in order to have a look at its structure (such as *ncdump {my_netcdf_file} | more*) and specially the parameters names
 * (this is useful if remapping is needed, see below for in the \ref netcdf_keywords "keywords" section).
 *
 * The NetCDF format does not impose a specific set of metadata and therefore in order to easily exchange data
 * within a given field, it is a good idea to standardize the metadata. Several such metadata schema can be used
 * by this plugin:
 * - CF1 - the <A HREF="http://cfconventions.org">conventions</A> for climate and forecast (CF) metadata;
 * - ECMWF - from the <A HREF="http://www.ecmwf.int/">European Centre for Medium-Range Weather Forecasts</A>, see the <A HREF="https://software.ecmwf.int/wiki/display/TIGGE/Soil+temperature">ECMWF Wiki</A> for a description of the available fields;
 * - CNRM - from the <A HREF="http://www.cnrm.meteo.fr/">National Centre for Meteorological Research</A>;
 * - WRF - the <A HREF="http://www.wrf-model.org/index.php">Weather Research & Forecasting</A> model.
 *
 * If you want to better understand the structure of the NetCDF file format, you are highly encouraged to read about
 * its <A HREF="https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_set_components.html">components</A>.
 *
 * @section netcdf_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 *
 * @section netcdf_keywords Keywords
 * This plugin uses the following keywords:
 * - TIME_ZONE: the time zone to use when interpreting date/time information; [Input] and [Output] section
 * - DEMFILE: The filename of the file containing the DEM; [Input] section
 * - DEMVAR: The variable name of the DEM within the DEMFILE; [Input] section
 * - GRID2DPATH: if this directory contains files, they will be used for reading the input from; [Input] and [Output] section
 * - NC_EXT: only the files containing this pattern in their filename will be used; [Input] section (default: .nc)
 * - GRID2DFILE: if GRID2DPATH has not been defined or if it does not contain files matching the NC_EXT extension, provides
 * the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 * - NETCDF_SCHEMA: the schema to use (either CF1 or CNRM or ECMWF or WRF); [Input] and [Output] section (default: ECMWF)
 * - NETCDF_VAR::{MeteoGrids::Parameters} = {netcdf_param_name} : this allows to remap the names as found in the NetCDF file to the MeteoIO grid parameters; [Input] section;
 * - NETCDF_DIM::{MeteoGrids::Parameters} = {netcdf_dimension_name} : this allows to remap the names as found in the NetCDF file to the ncParameters Dimensions; [Input] section;
 * - DEM_FROM_PRESSURE: if no dem is found but local and sea level pressure grids are found, use them to rebuild a DEM; [Input] section
 *
 * When providing multiple files in one directory, in case of overlapping files (because each file can provide multiple timestamps), the file containing the newest data has priority. This is
 * convenient when using forecats data to automatically use the most short-term forecast.
 *
 * @section netcdf_example Example use
 * Using this plugin to build downscaled time series at virtual stations, with the ECMWF Era Interim data set (see section below):
 * @code
 * [Input]
 * GRID2D    = NETCDF
 * GRID2DPATH =  /data/meteo_reanalysis
 * NETCDF_SCHEMA = ECMWF
 *
 * DEM = NETCDF
 * DEMFILE = /data/meteo_reanalysis/ECMWF_Europe_20150101-20150701.nc
 * DEM_FROM_PRESSURE = true
 *
 * #The lines below have nothing to do with this plugin
 * Downscaling = true
 * VSTATION1 = 46.793029 9.821343 ;this is Davos
 * Virtual_parameters = TA RH PSUM ISWR ILWR P VW DW TSS HS RSWR TSG ;this has to fit the parameter set in the data files
 * @endcode
 *
 * Another example, to extract precipitation from the MeteoSwiss daily precipitation reanalysis, RhiresD
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/ch02_lonlat.nc
 *
 * GRID2D    = NETCDF
 * GRID2DPATH =  /data/meteo_reanalysis
 * NC_EXT = .nc
 * NETCDF_VAR::PSUM = RhiresD               ;overwrite the PSUM parameter with "RhiresD", for example for MeteoCH reanalysis
 *
 * #The lines below have nothing to do with this plugin
 * Downscaling = true
 * VSTATION1 = 46.793029 9.821343 ;this is Davos
 * Virtual_parameters = PSUM ;this has to fit the parameter set in the data files
 * @endcode
 *
 * @section netcdf_meteoch MeteoCH RhiresD & similar products
 * <A HREF="http://www.meteoswiss.admin.ch/home.html?tab=overview">MeteoSwiss</A> provides <A HREF="http://www.ifu.ethz.ch/hydrologie/research/research_data/proddoc.pdf">reanalysis</A> of precipitation and other meteo fields from 1961 to present over Switzerland for different time scales: daily, monthly, yearly, as well as hourly (CPC dataset). The DEM are also provided, either in lat/lon,
 * Swiss coordinates, rotated lat/lon, ... These data sets must be requested from MeteoSwiss and are available with a specific license for research.
 *
 * @section netcdf_wrf WRF output files
 * While <A HREF="http://www.wrf-model.org/index.php">WRF</A> can write its <A HREF="http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap5.htm#fields">outputs in NetCDF</A>, unfortunately
 * it does not follow the CF1 convention and relies on lots of idiosyncracies (see http://www.ncl.ucar.edu/Applications/wrfnetcdf.shtml) that break lots of
 * applications dealing with NetCDF. If some fields are not read by MeteoIO,
 * please follow the tips given \ref netcdf_tricks "below". Moreover, WRF assumes that latitudes / longitudes are given on an ideal sphere while standard
 * coordinates systems assume an ellipsoid. This may lead to trouble when converting model coordinates to real world coordinates (see
 * http://www.pkrc.net/wrf-lambert.html).
 *
 * @section netcdf_ecmwf ECMWF Era Interim
 * The Era Interim data can be downloaded on the <A HREF="http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/">ECMWF dataserver</A>
 * after creating an account and login in.
 *
 * It is recommended to extract data at 00:00, and 12:00 for all steps 3, 6, 9, 12. The select the following fields:
 * 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, Forecast albedo, Mean sea level pressure, Skin temperature, Snow density, Snow depth, Soil temperature level 1, Surface pressure, Surface solar radiation downwards, Surface thermal radiation downwards, Total precipitation
 *
 * Here we have included the *forecast albedo* so the RSWR can be computed from ISWR and the *mean sea level pressure* and *surface pressure*
 * as proxies to compute the elevation. If you have the altitude in a separate file, it can be declared as DEM and there would be no need for the sea
 *level pressure (this would also be much more precise).
 *
 * You should therefore have the following request:
 * @code
 * Parameter: 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, Forecast albedo,
 *            Mean sea level pressure, Skin temperature, Snow density, Snow depth, Soil temperature level 1, Surface pressure,
 *            Surface solar radiation downwards, Surface thermal radiation downwards, Total precipitation
 *      Step: 3 to 12 by 3
 *      Type: Forecast
 *      Time: 00:00:00, 12:00:00
 * @endcode
 *
 * With the <A HREF="https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets">ECMWF Python Library</A>, the request
 * would be for example (the area is defined as North/West/South/East, see the
 * <A HREF="https://software.ecmwf.int/wiki/display/UDOC/Post-processing+keywords#Post-processingkeywords-area">WEBAPI</a> documentation):
 * @code
 * #!/usr/bin/env python
 * from ecmwfapi import ECMWFDataServer
 * server = ECMWFDataServer()
 * server.retrieve({
 * "class": "ei",
 * "dataset": "interim",
 * "date": "2015-01-01/to/2015-01-31",
 * "expver": "1",
 * "grid": "0.75/0.75",
 * "levtype": "sfc",
 * "param": "33.128/134.128/139.128/141.128/151.128/165.128/166.128/167.128/168.128/169.128/175.128/205.128/228.128/235.128/243.128",
 * "step": "3/6/9/12",
 * "area":"42.2/-1.5/51.7/15.7",
 * "stream": "oper",
 * "format":"netcdf",
 * "target": "my-era-interim.nc",
 * "time": "00/12",
 * "type": "fc",
 * })
 * @endcode
 *
 * @section netcdf_tricks Saving the day when a file is not standard compliant
 * Unfortunatelly, the naming of the parameters and dimensions within the files is not always standard nor consistent. In order to handle the parameters names,
 * simply run *ncdump {my_netcdf_file} | more* and use the name mapping facility of this plugin to map the non-standard parameters to our internal names
 * (see the \ref netcdf_keywords "plugin keywords"). When the dimensions are not standard (for example the time axis being called "TIME_T"),
 * use first the <A HREF="http://linux.die.net/man/1/ncrename">ncrename</A> tool that is part of the
 * <A HREF="http://nco.sourceforge.net/">NCO utilities</A> to rename both the dimension (-d) and the variable (-v):
 * @code
 * ncrename -d TIME_T,time -v TIME_T,time {my_netcdf_file}
 * @endcode
 */

//helper function to sort the cache of grid files
inline bool sort_cache_grids(const std::pair<std::pair<Date,Date>,ncParameters> &left, const std::pair<std::pair<Date,Date>,ncParameters> &right) {
	if (left.first.first < right.first.first) return true;
	if (left.first.first > right.first.first) return false;
	return left.first.second < right.first.second; //date_start equallity case
}

NetCDFIO::NetCDFIO(const std::string& configfile) 
         : cfg(configfile), cache_grid_files(), available_params(), in_schema("ECMWF"), out_schema("ECMWF"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), grid2d_out_file(), 
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) 
         : cfg(cfgreader), cache_grid_files(), available_params(), in_schema("ECMWF"), out_schema("ECMWF"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), grid2d_out_file(), 
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

void NetCDFIO::parseInputOutputSection()
{
	std::string in_grid2d, out_grid2d;
	cfg.getValue("GRID2D", "Input", in_grid2d, IOUtils::nothrow);
	cfg.getValue("GRID2D", "Output", out_grid2d, IOUtils::nothrow);
	if (in_grid2d=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Input", in_schema, IOUtils::nothrow); IOUtils::toUpper(in_schema);
		cfg.getValue("GRID2DPATH", "Input", in_grid2d_path);
		cfg.getValue("NC_EXT", "INPUT", in_nc_ext, IOUtils::nothrow);
		cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow);
		cfg.getValue("NC_DEBUG", "INPUT", debug, IOUtils::nothrow);
	}
	
	if (out_grid2d=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Output", out_schema, IOUtils::nothrow); IOUtils::toUpper(out_schema);
		cfg.getValue("GRID2DPATH", "Output", out_grid2d_path);
		cfg.getValue("GRID2DFILE", "Output", grid2d_out_file);
	}
	
	std::string in_meteo, out_meteo;
	cfg.getValue("METEO", "Input", in_meteo, IOUtils::nothrow);
	cfg.getValue("METEO", "Output", out_meteo, IOUtils::nothrow);
	if (in_meteo=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Input", out_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Input", in_schema, IOUtils::nothrow); IOUtils::toUpper(in_schema);
		cfg.getValue("METEOPATH", "Input", in_meteo_path);
	}
	if (out_meteo=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Output", out_schema, IOUtils::nothrow); IOUtils::toUpper(out_schema);
		cfg.getValue("METEOPATH", "Output", out_meteo_path);
		cfg.getValue("METEOFILE", "Output", out_meteo_file);
	}
}

void NetCDFIO::scanMeteoPath(const std::string& in_path, const std::string& nc_ext, std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files)
{
	meteo_files.clear();
	std::list<std::string> dirlist( FileUtils::readDirectory(in_path, nc_ext) );
	if (dirlist.empty()) return; //nothing to do if the directory is empty, we will transparently swap to using GRID2DFILE
	dirlist.sort();

	//Check date range in every filename and cache it
	std::list<std::string>::const_iterator it = dirlist.begin();
	while ((it != dirlist.end())) {
		const std::string filename( in_path + "/" + *it );
		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		const ncParameters ncFile = ncParameters(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		meteo_files.push_back( make_pair(ncFile.getDateRange(), ncFile) );
		it++;
	}
	std::sort(meteo_files.begin(), meteo_files.end(), &sort_cache_grids);

	//now handle overlaping files: truncate the end date of the file starting earlier
	for (size_t ii=0; ii<(meteo_files.size()-1); ii++) {
		if (meteo_files[ii].first.second > meteo_files[ii+1].first.first)
			meteo_files[ii].first.second = meteo_files[ii+1].first.first;
	}
}

bool NetCDFIO::list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list)
{
	if (cache_grid_files.empty()) scanMeteoPath(in_grid2d_path, in_nc_ext, cache_grid_files);
	if (cache_grid_files.empty()) return true; //there are no grids to read
	
	//HACK handle the case of file_start & file_end are undef() (example: DEM)
	for (size_t ii=0; ii<cache_grid_files.size(); ii++) {
		const Date file_start( cache_grid_files[ii].first.first );
		const Date file_end( cache_grid_files[ii].first.second );
		
		if (file_start > end) return true; //no more files to process (since the files are sorted in cache_grid_files)
		if (file_end < start) continue;
		
		//we consider that the exact same parameters are available at all time steps in the current file
		const std::set<size_t> params_set( cache_grid_files[ii].second.getParams() );
		const std::vector<Date> ts( cache_grid_files[ii].second.getTimestamps() );
		for (size_t jj=0; jj<ts.size(); jj++) {
			if (ts[jj]>end) break; //no more timestamps in the right range
			if (ts[jj]<start) continue;
			list[ ts[jj] ] = params_set;
		}
	}

	return true;
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	std::vector<std::string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		const ncParameters ncFile(vec_argument[0], ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::read2DGrid is filename:varname", AT);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (cache_grid_files.empty()) scanMeteoPath(in_grid2d_path, in_nc_ext, cache_grid_files);
	
	if (!cache_grid_files.empty()) {
		for (size_t ii=0; ii<cache_grid_files.size(); ii++) {
			const Date file_start( cache_grid_files[ii].first.first );
			const Date file_end( cache_grid_files[ii].first.second );
			if (file_start > date) return;
			if (file_end < date) continue;
			
			if (date>=file_start && date<=file_end) {
				grid_out = cache_grid_files[ii].second.read2DGrid(parameter, date);
				return;
			}
		}
		//the date was not found
		throw InvalidArgumentException("No Gridded data found for "+date.toString(Date::ISO)+" in '"+in_grid2d_path+"'", AT);
	} else {
		const std::string filename = cfg.get("GRID2DFILE", "Input");
		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		const ncParameters ncFile(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(parameter, date);
	}
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	const std::string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	
	const ncParameters ncFile(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
	const Grid2DObject grid = (varname.empty())? ncFile.read2DGrid(MeteoGrids::DEM, Date()) : ncFile.read2DGrid(varname);
	dem_out = DEMObject( grid ); //we can not directly assign a Grid2DObject to a DEMObject
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const std::string& arguments)
{
	// arguments is a string of the format filename:varname
	std::vector<std::string> vec_argument;
	if (IOUtils::readLineToVec(arguments, vec_argument, ':')  != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	ncParameters ncFile(out_grid2d_path+"/"+vec_argument[0], ncParameters::WRITE, cfg, out_schema, out_dflt_TZ, debug);
	const ncParameters::var_attr attr(-1, vec_argument[1], IOUtils::nodata);
	ncParameters::nc_variable tmp_var(attr);
	ncFile.write2DGrid(grid_in, tmp_var, Date());
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const std::string file_and_path( out_grid2d_path + "/" + grid2d_out_file );
	ncParameters ncFile(file_and_path, ncParameters::WRITE, cfg, out_schema, out_dflt_TZ, debug);
	if (parameter==MeteoGrids::DEM || parameter==MeteoGrids::SHADE || parameter==MeteoGrids::SLOPE || parameter==MeteoGrids::AZI)
		ncFile.write2DGrid(grid_in, parameter, Date()); //do not assign a date to a DEM?
	else
		ncFile.write2DGrid(grid_in, parameter, date);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	throw IOException("Not implemented yet! Be patient!", AT);
	const std::string file_and_path( out_meteo_path + "/" + out_meteo_file );
	ncParameters ncFile(file_and_path, ncParameters::WRITE, cfg, out_schema, out_dflt_TZ, debug);
	ncFile.writeMeteo(vecMeteo);
}


///////////////////////////////////////////////////// Now the ncParameters class starts //////////////////////////////////////////
std::vector<std::string> ncParameters::dimnames( initDimensionNames() );
std::map< std::string, std::vector<ncParameters::nc_dimension> > ncParameters::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::schemas_vars( initSchemasVars() );

std::vector<std::string> ncParameters::initDimensionNames()
{
	//the order must be the same as in the enum Dimensions
	std::vector<std::string> tmp;
	tmp.push_back("NONE"); tmp.push_back("TIME"); 
	tmp.push_back("LATITUDE"); tmp.push_back("LONGITUDE");
	tmp.push_back("NORTHING"); tmp.push_back("EASTING"); tmp.push_back("STATION");
	
	return tmp;
}

std::map< std::string, std::vector<ncParameters::nc_dimension> > ncParameters::initSchemasDims()
{
	std::map< std::string, std::vector<ncParameters::nc_dimension> > results;
	std::vector<ncParameters::nc_dimension> tmp;
	
	//CF1 schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "latitude") );
	tmp.push_back( nc_dimension(LONGITUDE, "longitude") );
	tmp.push_back( nc_dimension(STATION, "station") );
	tmp.push_back( nc_dimension(EASTING, "easting") );
	tmp.push_back( nc_dimension(NORTHING, "northing") );
	tmp.push_back( nc_dimension(MeteoGrids::DEM, "surface_altitude") );
	results["CF1"] = tmp;
	
	//CNRM schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "latitude") );
	tmp.push_back( nc_dimension(LONGITUDE, "longitude") );
	tmp.push_back( nc_dimension(STATION, "station") );
	tmp.push_back( nc_dimension(EASTING, "easting") );
	tmp.push_back( nc_dimension(NORTHING, "northing") );
	tmp.push_back( nc_dimension(MeteoGrids::DEM, "ZS") );
	results["CNRM"] = tmp;
	
	//ECMWF schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "latitude") );
	tmp.push_back( nc_dimension(LONGITUDE, "longitude") );
	tmp.push_back( nc_dimension(STATION, "station") );
	tmp.push_back( nc_dimension(EASTING, "easting") );
	tmp.push_back( nc_dimension(NORTHING, "northing") );
	tmp.push_back( nc_dimension(MeteoGrids::DEM, "geopotential_height") );
	results["ECMWF"] = tmp;
	
	//WRF schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "Time") );
	tmp.push_back( nc_dimension(LATITUDE, "south_north") );
	tmp.push_back( nc_dimension(LONGITUDE, "west_east") );
	tmp.push_back( nc_dimension(STATION, "station") );
	tmp.push_back( nc_dimension(EASTING, "easting") );
	tmp.push_back( nc_dimension(NORTHING, "northing") );
	tmp.push_back( nc_dimension(MeteoGrids::DEM, "HGT") );
	results["WRF"] = tmp;
	
	return results;
}

std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::initSchemasVars()
{ //HACK: vars/dims should be identified based on standard_name, not name (cf1)
	std::map< std::string, std::vector<ncParameters::var_attr> > results;
	std::vector<ncParameters::var_attr> tmp;

	//CF1 schema -> to be checked and improved from CF1 documentation
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "", "s", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "latitude", "latitude", "", "degree_north", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "longitude", "longitude", "", "degree_east", IOUtils::nodata) );
	tmp.push_back( var_attr(STATION, "station", "timeseries_id", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(EASTING, "easting", "projection_x_coordinate", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(NORTHING, "northing", "projection_y_coordinate", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "surface_altitude", "surface_altitude", "height above mean sea level", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "air_temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RH, "relative_humidity", "relative_humidity", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DW, "wind_from_direction", "wind_from_direction", "", "degree", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::VW, "wind_speed", "wind_speed", "", "m/s", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "air_pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR, "downwelling_shortwave_flux_in_air", "downwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RSWR, "upwelling_shortwave_flux_in_air", "upwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "downwelling_longwave_flux_in_air", "downwelling_longwave_flux_in_air", "", "W/m2", IOUtils::nodata) );
	results["CF1"] = tmp;

	//CNRM schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "time", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(STATION, "station", "timeseries_id", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(EASTING, "easting", "easting", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(NORTHING, "northing", "northing", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "ZS", "", "altitude", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RH, "HUMREL", "", "Relative Humidity", "%", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::QI, "Qair", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::VW, "Wind", "", "Wind Speed", "m/s", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DW, "Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM_L, "Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM_S, "Snowf", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR_DIR, "DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR_DIFF, "SCA_SWdown", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata) );
	results["CNRM"] = tmp;

	//ECMWF schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "time", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(STATION, "station", "timeseries_id", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(EASTING, "easting", "easting", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(NORTHING, "northing", "northing", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "z", "geopotential_height", "geopotential_height", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "t2m", "", "2 metre temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::TD, "d2m", "", "2 metre dewpoint temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::P, "sp", "surface_air_pressure", "Surface pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P_SEA, "msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR, "ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "strd", "", "Surface thermal radiation downwards", "J m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM, "tp", "", "Total precipitation", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::U, "u10", "", "10 metre U wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::V, "v10", "", "10 metre V wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::SWE, "sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSS, "skt", "", "Skin temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSG, "stl1", "surface_temperature", "Soil temperature level 1", "K", IOUtils::nodata) ); //this is from 0 to -7cm
	tmp.push_back( var_attr(MeteoGrids::ALB, "al", "surface_albedo", "Albedo", "(0 - 1)", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ALB, "fal", "", "Forecast albedo", "(0 - 1)", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RSNO, "rsn", "", "Snow density", "kg m**-3", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ROT, "ro", "", "Runoff", "m", IOUtils::nodata) );
	results["ECMWF"] = tmp;

	//WRF schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "Times", "Times", "Times", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "XLAT", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "XLONG", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(STATION, "station", "timeseries_id", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(EASTING, "easting", "easting", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(NORTHING, "northing", "northing", "", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "HGT", "Terrain Height", "Terrain Height", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "PSFC", "Surface pressure", "Surface pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "T2", "2-meter temperature", "2-meter temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::QI, "Q2", "2-meter specific humidity", "2-meter specific humidity", "kg kg-1", 2) );
	tmp.push_back( var_attr(MeteoGrids::ISWR, "ACSWDNB", "Downward SW surface radiation", "Downward SW surface radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RSWR, "ACSWUPB", "Upwelling Surface Shortwave Radiation", "Upwelling Surface Shortwave Radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "ACLWDNB", "Downward LW surface radiation", "Downward LW surface radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ROT, "SFROFF", "Surface runoff ", "Surface runoff ", "kg*m2*s-1", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::HS, "SNOWH", "Snow depth", "Snow depth", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSS, "TSK", "Surface skin temperature", "Surface skin temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::U, "U10", "10-meter wind speed", "10 metre U wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::V, "V10", "10-meter wind speed", "10 metre V wind component", "m s**-1", 10.) );
	results["WRF"] = tmp;
	
	return results;
}

//The user can provide his own variables properties as NETCDF_VAR::{param} = {name}
std::vector<ncParameters::var_attr> ncParameters::initUserSchemas(const Config& i_cfg)
{
	std::vector<ncParameters::var_attr> results;
	
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_VAR::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string meteo_grid( custom_attr[ii].substr(found+1) );
		const std::string netcdf_param = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = getParameterIndex(meteo_grid);
		if (param_index==IOUtils::npos)
			throw InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( var_attr(param_index, netcdf_param, IOUtils::nodata) );
	}
	
	return results;
}

//The user can provide his own dimensions properties as NETCDF_DIM::{dimension_param} = {name_in_current_file}
std::vector<ncParameters::nc_dimension> ncParameters::initUserDimensions(const Config& i_cfg)
{
	std::vector<ncParameters::nc_dimension> results;
	
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_DIM::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string dim_str( custom_attr[ii].substr(found+1) );
		const std::string netcdf_dim = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = getParameterIndex(dim_str);
		if (param_index==IOUtils::npos || param_index<firstdimension || param_index>lastdimension)
			throw InvalidArgumentException("Dimension '"+dim_str+"' is not a valid dimension! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( nc_dimension( static_cast<Dimensions>(param_index), netcdf_dim) );
	}
	
	return results;
}

ncParameters::ncParameters(const std::string& filename, const Mode& mode, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug)
             : user_schemas( initUserSchemas(cfg) ), user_dimensions( initUserDimensions(cfg) ), vars(), unknown_vars(), vecTime(), vecX(), vecY(), dimensions_map(), file_and_path(filename), coord_sys(), coord_param(), TZ(tz_in), wrf_hacks(schema=="WRF"), debug(i_debug), isLatLon(false)
{
	IOUtils::getProjectionParameters(cfg, coord_sys, coord_param);
	
	if (mode==WRITE) {
		initFromSchema(schema);
		if (FileUtils::fileExists(filename)) initFromFile(filename, schema);
	} else if (mode==READ)
		initFromFile(filename, schema);
	
	if (debug) {
		std::cout << filename << ":\n";
		std::cout << "\tDimensions:\n";
		for (std::map<size_t, nc_dimension>::const_iterator it = dimensions_map.begin(); it!=dimensions_map.end(); ++it)
			std::cout << "\t\t" << it->second.toString() << "\n";
		if (!vecTime.empty()) std::cout << "\ttime range: [" << vecTime.front().toString(Date::ISO) << " - " << vecTime.back().toString(Date::ISO) << "]\n";
		std::cout << "\tVariables:\n";
		for (std::map<size_t, nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
			std::cout << "\t\t" << getParameterName( it->first ) << " -> " << it->second.toString() << "\n";
		std::cout << "\tUnrecognized variables:\n";
		for (std::map<std::string, nc_variable>::const_iterator it=unknown_vars.begin(); it!=unknown_vars.end(); ++it)
			std::cout << "\t\t" << it->first << " -> " << it->second.toString() << "\n";
	}
}

//populate the dimensions_map from the selected schema
void ncParameters::initFromSchema(const std::string& schema)
{
	for (size_t ii=0; ii<schemas_dims[schema].size(); ii++) {
		dimensions_map[ schemas_dims[schema][ii].type ] = schemas_dims[schema][ii];
	}
	if (dimensions_map.count(TIME)==0) throw IOException("No TIME dimension in schema '"+schema+"'", AT);
	dimensions_map[ TIME ].isUnlimited = true;
	
	for (size_t ii=0; ii<schemas_vars[schema].size(); ii++) {
		vars[ schemas_vars[schema][ii].param ] = nc_variable( schemas_vars[schema][ii] );
	}
}

//populate the dimensions_map and vars and unknown_vars from the file
void ncParameters::initFromFile(const std::string& filename, const std::string& schema)
{
	if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
	
	int ncid;
	ncpp::open_file(filename, NC_NOWRITE, ncid);
	
	//read the dimensions and variables
	initDimensionsFromFile(ncid, schema);
	initVariablesFromFile(ncid, schema);

	isLatLon = hasDimension(LATITUDE) && hasDimension(LONGITUDE);
	const bool isXY = hasDimension(EASTING) && hasDimension(NORTHING);
	if (isLatLon) {
		vecY = read_1Dvariable(ncid, LATITUDE);
		vecX = read_1Dvariable(ncid, LONGITUDE);
	} else if (isXY) {
		vecX = read_1Dvariable(ncid, EASTING);
		vecY = read_1Dvariable(ncid, NORTHING);
	}
	if (hasDimension(TIME)) vecTime = read_1Dvariable(ncid);
	
	if (ncpp::check_attribute(ncid, NC_GLOBAL, "epsg")) {
		int epsg = IOUtils::inodata;
		const int status = nc_get_att_int(ncid, NC_GLOBAL, "epsg", &epsg);
		if (status != NC_NOERR) throw IOException("Could not read attribute 'epsg': " + std::string(nc_strerror(status)), AT);
		if (epsg!=IOUtils::inodata) CoordsAlgorithms::EPSG_to_str(epsg, coord_sys, coord_param);
	}
	
	ncpp::close_file(filename, ncid);
}

std::pair<Date, Date> ncParameters::getDateRange() const
{
	if (vecTime.empty()) return make_pair( Date(), Date() );
	return make_pair( vecTime.front(), vecTime.back() );
}

std::set<size_t> ncParameters::getParams() const 
{
	std::set<size_t> available_params;
	for (std::map<size_t, nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
		available_params.insert( it->first );
	
	return available_params;
}

Grid2DObject ncParameters::read2DGrid(const std::string& varname) const
{
	for (std::map<size_t, nc_variable>::const_iterator it = vars.begin(); it!=vars.end(); ++it) {
		if (it->second.attributes.name==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	for (std::map<std::string, nc_variable>::const_iterator it = unknown_vars.begin(); it!=unknown_vars.end(); ++it) {
		if (it->first==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	throw NotFoundException("The variable '"+varname+"' could not be found in file '"+file_and_path+"'", AT);
}

Grid2DObject ncParameters::read2DGrid(const size_t& param, const Date& date) const
{
	const std::map <size_t, nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end()) NoDataException("No "+MeteoGrids::getParameterName( param )+" grid in file "+file_and_path, AT);
	
	size_t time_pos = IOUtils::npos;
	if (!date.isUndef()) {
		const std::vector<Date>::const_iterator low = std::lower_bound(vecTime.begin(), vecTime.end(), date);
		if (*low!=date) throw NoDataException("No data at "+date.toString(Date::ISO)+" in file "+file_and_path, AT);
		time_pos = static_cast<size_t>( std::distance(vecTime.begin(), low) );
	} else { 
		//no date has been provided, check if this parameter depends on time
		const std::map<size_t, nc_dimension>::const_iterator it2 = dimensions_map.find(TIME);
		const int time_id = (it2!=dimensions_map.end())? it2->second.dimid : -1;
		const bool depend_on_time = (std::find(it->second.dimids.begin(), it->second.dimids.end(), time_id) != it->second.dimids.end());
		if (depend_on_time && !vecTime.empty() && (param!=MeteoGrids::DEM && param!=MeteoGrids::SLOPE && param!=MeteoGrids::AZI))
			throw InvalidFormatException("No time requirement has been provided for a file that contains multiple timestamps", AT);
	}
	
	const bool isPrecip = (param==MeteoGrids::PSUM || param==MeteoGrids::PSUM_L || param==MeteoGrids::PSUM_S);
	const bool isRad = (param==MeteoGrids::ISWR || param==MeteoGrids::RSWR || param==MeteoGrids::ISWR_DIFF || param==MeteoGrids::ISWR_DIR);
	return read2DGrid(it->second, time_pos, isPrecip, (isPrecip || isRad));
}

Grid2DObject ncParameters::read2DGrid(const nc_variable& var, const size_t& time_pos, const bool& m2mm, const bool& reZero) const
{
	if (isLatLon && (!hasDimension(LATITUDE) || !hasDimension(LONGITUDE))) throw IOException("No latitude / longitude could be identified in file "+file_and_path, AT);
	if (!isLatLon && (!hasDimension(EASTING) || !hasDimension(NORTHING))) throw IOException("No easting / northing could be identified in file "+file_and_path, AT);
		
	//define the results grid
	mio::Coords llcorner(coord_sys, coord_param); //if an EPSG was provided, this has been converted to coord_sys/coord_param
	if (isLatLon)
		llcorner.setLatLon( std::min(vecY.front(), vecY.back()), std::min(vecX.front(), vecX.back()), IOUtils::nodata);
	else
		llcorner.setXY( std::min(vecY.front(), vecY.back()), std::min(vecX.front(), vecX.back()), IOUtils::nodata);
	//HACK expand the definition of Grid2DObject to support lat/lon grids and reproject in GridsManager
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = (isLatLon)? calculate_cellsize(resampling_factor_x, resampling_factor_y) : calculate_XYcellsize(resampling_factor_x, resampling_factor_y);
	Grid2DObject grid(vecX.size(), vecY.size(), cellsize, llcorner);
	
	//read the raw data, copy it into the Grid2DObject
	int ncid;
	ncpp::open_file(file_and_path, NC_NOWRITE, ncid);
	
	double *data = new double[vecY.size()*vecX.size()];
	if (time_pos!=IOUtils::npos)
		ncpp::read_data(ncid, var.attributes.name, var.varid, time_pos, vecY.size(), vecX.size(), data);
	else
		ncpp::read_data(ncid, var.attributes.name, var.varid, data);
	fill2DGrid(grid, data, var.nodata);
	delete[] data;
	ncpp::close_file(file_and_path, ncid);
	
	//handle data packing and units, if necessary
	if (var.scale!=1.) grid *= var.scale;
	if (var.offset!=0.) grid += var.offset;
	const std::string units( var.attributes.units );
	if (!units.empty()) {
		if (units=="m**2 s**-2") grid /= Cst::gravity;
		else if (units=="%") grid /= 100.;
		else if (units=="J m**-2") {
			if (vecTime.size()>1 && time_pos!=IOUtils::npos) {
				const Date integration_period = (time_pos>0)? (vecTime[time_pos] - vecTime[time_pos-1]) : (vecTime[time_pos+1] - vecTime[time_pos]);
				grid /= (integration_period.getJulian()*24.*3600.); //converting back to W/m2
			}
		}
		else if (m2mm && units=="m") grid *= 1000.;
		
		if (reZero) {//reset very low values to zero
			for (size_t ii=0; ii<grid.size(); ii++)
				if (grid(ii)<1e-6 && grid(ii)!=mio::IOUtils::nodata) grid(ii)=0.;
		}
	}
	
	return grid;
}

void ncParameters::write2DGrid(Grid2DObject grid_in, const size_t& param, const Date& date)
{
	if (vars.count(param)>0) {
		write2DGrid(grid_in, vars[param], date);
	} else {
		const std::string param_name( MeteoGrids::getParameterName(param) );
		const var_attr tmp_attr(param, param_name, IOUtils::nodata);
		nc_variable tmp_var(tmp_attr);
		if (unknown_vars.count( param_name )>0) tmp_var = unknown_vars[ param_name ];
		write2DGrid(grid_in, tmp_var, date);
		unknown_vars[ param_name ] = tmp_var; //save varid, etc for next calls
	}
}

void ncParameters::write2DGrid(const Grid2DObject& grid_in, nc_variable& var, const Date& date)
{
	//check that the necessary dimensions are available in the maps
	const bool hasLatLon = (dimensions_map.count(LATITUDE)!=0 && dimensions_map.count(LONGITUDE)!=0) && (vars.count(LATITUDE)!=0 && vars.count(LONGITUDE)!=0);
	const bool hasEastNorth = (dimensions_map.count(EASTING)!=0 && dimensions_map.count(NORTHING)!=0) && (vars.count(EASTING)!=0 && vars.count(NORTHING)!=0);
	const bool hasTime = dimensions_map.count(TIME)!=0 && vars.count(TIME)!=0;
	if ((isLatLon && !hasLatLon) || (!isLatLon && !hasEastNorth) || !hasTime)
		throw IOException("Error in the schema definition, some basic quantities are not defined!", AT);
	
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		const int status = nc_redef(ncid);
		if (status != NC_NOERR) throw IOException("Could not open define mode for file '" + file_and_path + "': " + nc_strerror(status), AT);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
		ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO "+getLibVersion());
		if (!isLatLon) ncpp::add_attribute(ncid, NC_GLOBAL, "epsg", grid_in.llcorner.getEPSG());
	}
	
	//create any potentially missing definition, otherwise check that everything is consistent
	const bool is_record = (!date.isUndef());
	if (is_record) create_TimeDimension(ncid, date, 0);
	bool fill_spatial_vars = true;
	if (isLatLon) {
		fill_spatial_vars = create_Dimension(ncid, LATITUDE, grid_in.getNy());
		fill_spatial_vars = create_Dimension(ncid, LONGITUDE, grid_in.getNx());
	} else {
		fill_spatial_vars = create_Dimension(ncid, NORTHING, grid_in.getNy());
		fill_spatial_vars = create_Dimension(ncid, EASTING, grid_in.getNx());
	}
	if (var.varid == -1) { //the variable must be created
		if (is_record) var.dimids.push_back(dimensions_map[TIME].dimid);
		if (isLatLon) {
			var.dimids.push_back(dimensions_map[LATITUDE].dimid);
			var.dimids.push_back(dimensions_map[LONGITUDE].dimid);
		} else {
			var.dimids.push_back(dimensions_map[NORTHING].dimid);
			var.dimids.push_back(dimensions_map[EASTING].dimid);
		}
		ncParameters::create_variable(ncid, var);
	}
	ncpp::end_definitions(file_and_path, ncid);
	
	//now write the data
	if (fill_spatial_vars) fill_SpatialDimensions(ncid, grid_in);
	
	double *data = new double[grid_in.size()];
	ncpp::fill_grid_data(grid_in, data);
	if (is_record) {
		const size_t time_pos = addTimestamp(ncid, date);
		ncpp::write_data(ncid, var.attributes.name, var.varid, grid_in.getNy(), grid_in.getNx(), time_pos, data);
	} else {
		ncpp::write_data(ncid, var.attributes.name, var.varid, data);
	}
	delete[] data;
	
	ncpp::close_file(file_and_path, ncid);
}

void ncParameters::fill_SpatialDimensions(const int& ncid, const Grid2DObject& grid_in)
{
	//write the dimensions' data
	double *Y_array = new double[grid_in.getNy()];
	double *X_array = new double[grid_in.getNx()];
	if (isLatLon) {
		ncpp::calculate_dimensions(grid_in, Y_array, X_array);
		ncpp::write_data(ncid, vars[LATITUDE].attributes.name, vars[LATITUDE].varid, Y_array);
		ncpp::write_data(ncid, vars[LONGITUDE].attributes.name, vars[LONGITUDE].varid, X_array);
	} else {
		const double cellsize = grid_in.cellsize;
		X_array[0] = grid_in.llcorner.getEasting();
		for (size_t ii=1; ii<grid_in.getNx(); ii++) X_array[ii] = X_array[ii-1] + cellsize;
		Y_array[0] = grid_in.llcorner.getNorthing();
		for (size_t ii=1; ii<grid_in.getNy(); ii++) Y_array[ii] = Y_array[ii-1] + cellsize;
		ncpp::write_data(ncid, vars[NORTHING].attributes.name, vars[NORTHING].varid, Y_array);
		ncpp::write_data(ncid, vars[EASTING].attributes.name, vars[EASTING].varid, X_array);
	}
	delete[] Y_array; delete[] X_array;
}

void ncParameters::fill_SpatialDimensions(const int& ncid, const std::vector< std::vector<MeteoData> >& vecMeteo)
{
	std::vector<double> vecLat, vecLon, vecAlt;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) { //we consider that all the metadata is available at vecMete[0]
		vecLat.push_back( vecMeteo[ii].front().meta.position.getLat() );
		vecLon.push_back( vecMeteo[ii].front().meta.position.getLon() );
		vecAlt.push_back( vecMeteo[ii].front().meta.position.getAltitude() );
	}

	ncpp::write_data(ncid, vars[LATITUDE].attributes.name, vars[LATITUDE].varid, &vecLat[0]);
	ncpp::write_data(ncid, vars[LONGITUDE].attributes.name, vars[LONGITUDE].varid, &vecLon[0]);
	ncpp::write_data(ncid, vars[MeteoGrids::DEM].attributes.name, vars[MeteoGrids::DEM].varid, &vecAlt[0]);
}

void ncParameters::writeMeteo(const std::vector< std::vector<MeteoData> >& vecMeteo)
{//HACK instead, have a vector of which dims and their associated vars should be writen. The call "writeDims" and "writeVars" on these.
	if (vecMeteo.empty()) return;
	isLatLon = true;
	//check that the necessary dimensions are available in the maps
	const bool hasLatLon = (dimensions_map.count(LATITUDE)!=0 && dimensions_map.count(LONGITUDE)!=0) && (vars.count(LATITUDE)!=0 && vars.count(LONGITUDE)!=0);
	const bool hasEastNorth = (dimensions_map.count(EASTING)!=0 && dimensions_map.count(NORTHING)!=0) && (vars.count(EASTING)!=0 && vars.count(NORTHING)!=0);
	const bool hasTime = dimensions_map.count(TIME)!=0 && vars.count(TIME)!=0;
	if ((isLatLon && !hasLatLon) || (!isLatLon && !hasEastNorth) || !hasTime)
		throw IOException("Error in the schema definition, some basic quantities are not defined!", AT);
	
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		const int status = nc_redef(ncid);
		if (status != NC_NOERR) throw IOException("Could not open define mode for file '" + file_and_path + "': " + nc_strerror(status), AT);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
		ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO "+getLibVersion());
		if (!isLatLon) ncpp::add_attribute(ncid, NC_GLOBAL, "epsg", vecMeteo.front().front().meta.position.getEPSG()); //HACK
	}
	
	//create any potentially missing definition, otherwise check that everything is consistent
	create_TimeDimension(ncid, vecMeteo.front().front().date, 0); //HACK
	create_Dimension(ncid, STATION, vecMeteo.size());
	bool fill_spatial_vars = create_Dimension(ncid, MeteoGrids::DEM, vecMeteo.size());
	if (isLatLon) {
		fill_spatial_vars = create_Dimension(ncid, LATITUDE, vecMeteo.size());
		fill_spatial_vars = create_Dimension(ncid, LONGITUDE, vecMeteo.size());
	} else {
		fill_spatial_vars = create_Dimension(ncid, NORTHING, vecMeteo.size());
		fill_spatial_vars = create_Dimension(ncid, EASTING, vecMeteo.size());
	}
	ncpp::end_definitions(file_and_path, ncid);
	
	//now write the data
	if (fill_spatial_vars) fill_SpatialDimensions(ncid, vecMeteo);
	
	//write data
	
	ncpp::close_file(file_and_path, ncid);
}

//this returns the index where to insert the new grid
size_t ncParameters::addTimestamp(const int& ncid, const Date& date)
{
	size_t time_pos = vecTime.size();
	bool create_timestamp = true;

	//search for the proper insertion position
	if (time_pos>0) {
		//we keep vecTime in sync with the file, so we can do the searches into vecTime
		if (date==vecTime.back()) {
			create_timestamp = false;
			time_pos--;
		} else if (date<vecTime.back()) {
			const std::vector<Date>::const_iterator low = std::lower_bound(vecTime.begin(), vecTime.end(), date);
			if (*low==date) create_timestamp = false;
			time_pos = static_cast<size_t>( std::distance<std::vector<Date>::const_iterator>(vecTime.begin(), low) ); //weird syntax to take the proper template
		}
	}
	
	if (create_timestamp) {
		const double dt = (date.getJulian() - vars[TIME].offset) / vars[TIME].scale;
		const size_t start[] = {time_pos};
		const size_t count[] = {1};
		const int status = nc_put_vara_double(ncid, vars[TIME].varid, start, count, &dt);
		if (status != NC_NOERR) throw IOException("Could not write data for record variable '" + vars[TIME].attributes.name + "': " + nc_strerror(status), AT);
		if (time_pos==vecTime.size())
			vecTime.push_back( date );
		else
			vecTime.insert(vecTime.begin()+time_pos, date);
	}
	return time_pos;
}

//this creates what must be created, otherwise checks that what is present is consistent. Returns false if nothing had to be created
bool ncParameters::create_Dimension(const int& ncid, const size_t& param, const size_t& length)
{
	if (dimensions_map[ param ].dimid == -1) {
		dimensions_map[ param ].length = (dimensions_map[ param ].isUnlimited)? 0 : length;
		const nc_type len = (dimensions_map[ param ].isUnlimited)? NC_UNLIMITED : static_cast<int>(length);
		const int status = nc_def_dim(ncid, dimensions_map[ param ].name.c_str(), len, &dimensions_map[ param ].dimid);
		if (status != NC_NOERR) throw IOException("Could not define dimension '" + dimensions_map[ param ].name + "': " + nc_strerror(status), AT);
	} else {
		if (dimensions_map[ param ].length != length)
			throw InvalidArgumentException("Attempting to write an inconsistent lenght into dimension '" + getParameterName(param)+"' into file '"+file_and_path+"'", AT);
	}
	
	if (vars[ param ].varid == -1) {
		vars[ param ].dimids.push_back( dimensions_map[ param ].dimid );
		ncParameters::create_variable(ncid, vars[ param ]);
		return true;
	}
	
	return false;
}

bool ncParameters::create_TimeDimension(const int& ncid, const Date& date, const size_t& length)
{
	if (dimensions_map[ TIME ].dimid == -1) {
		dimensions_map[ TIME ].length = (dimensions_map[ TIME ].isUnlimited)? 0 : length;
		const nc_type len = (dimensions_map[ TIME ].isUnlimited)? NC_UNLIMITED : static_cast<int>(length);
		const int status = nc_def_dim(ncid, dimensions_map[ TIME ].name.c_str(), len, &dimensions_map[ TIME ].dimid);
		if (status != NC_NOERR) throw IOException("Could not define dimension '" + dimensions_map[ TIME ].name + "': " + nc_strerror(status), AT);
	}
	
	if (vars[ TIME ].varid == -1) {
		vars[ TIME ].dimids.push_back( dimensions_map[ TIME ].dimid );
		const Date ref_date(date.getYear(), 1, 1, 0, 0, TZ); //HACK move all to GMT or keep local?
		std::string date_str( ref_date.toString(Date::ISO) );
		date_str[ 10 ] = ' '; //replace "T" by " "
		vars[TIME].attributes.units = "hours since " + date_str;
		vars[TIME].offset = ref_date.getJulian();
		vars[TIME].scale = 1./24.;
		
		ncParameters::create_variable(ncid, vars[TIME]);
		return true;
	}
	
	return false;
}

//write the variable's attributes into the file
void ncParameters::create_variable(const int& ncid, nc_variable& var)
{
	const int ndims = static_cast<int>( var.dimids.size() );
	const int status = nc_def_var(ncid, var.attributes.name.c_str(), NC_DOUBLE, ndims, &var.dimids[0], &var.varid);
	if (status != NC_NOERR) throw IOException("Could not define variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	
	if (!var.attributes.standard_name.empty()) ncpp::add_attribute(ncid, var.varid, "standard_name", var.attributes.standard_name);
	if (!var.attributes.long_name.empty()) ncpp::add_attribute(ncid, var.varid, "long_name", var.attributes.long_name);
	if (!var.attributes.units.empty()) ncpp::add_attribute(ncid, var.varid, "units", var.attributes.units);
	ncpp::add_attribute(ncid, var.varid, "_FillValue", var.nodata);
	
	if (var.attributes.param==TIME) ncpp::add_attribute(ncid, var.varid, "calendar", "gregorian");
	if (var.attributes.param==MeteoGrids::DEM) {
		ncpp::add_attribute(ncid, var.varid, "positive", "up");
		ncpp::add_attribute(ncid, var.varid, "axis", "Z");
	}
}

void ncParameters::initDimensionsFromFile(const int& ncid, const std::string& schema_name)
{
	int status;
	int ndims;
	status = nc_inq_ndims(ncid, &ndims);
	if (status != NC_NOERR) throw IOException("Could not retrieve number of dimensions: " + std::string(nc_strerror(status)), AT);
	
	int unlim_id;
	status = nc_inq_unlimdim(ncid, &unlim_id);
	if (status != NC_NOERR) throw IOException("Could not retrieve unlimited dimension: " + std::string(nc_strerror(status)), AT);
	
	int *dimids = (int*)malloc(ndims * sizeof(int));
	status = nc_inq_dimids(ncid, &ndims, dimids, 0);
	if (status != NC_NOERR) throw IOException("Could not retrieve dimensions IDs: " + std::string(nc_strerror(status)), AT);
	
	for (int idx=0; idx<ndims; idx++) {
		char name[NC_MAX_NAME+1];
		status = nc_inq_dimname(ncid, dimids[idx], name);
		const std::string dimname( name );
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension name: " + std::string(nc_strerror(status)), AT);
		
		nc_dimension tmp_dim( getSchemaDimension(dimname, schema_name) ); //set name and type
		if (tmp_dim.type==NONE) { //unrecognized dimension -> try harder with some typical names that are not in the schema
			if (dimname=="lat")
				tmp_dim.type = LATITUDE;
			else if (dimname=="lon")
				tmp_dim.type = LONGITUDE;
			else continue;
			if (dimensions_map.count( tmp_dim.type )>0) continue; //if this parameter has already been read, skip it (so the schema naming has priority)
		}

		tmp_dim.dimid = idx;
		tmp_dim.isUnlimited = (idx==unlim_id);
		status = nc_inq_dimlen(ncid, dimids[idx], &tmp_dim.length);
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension lenght: " + std::string(nc_strerror(status)), AT);
		
		dimensions_map[ tmp_dim.type ] = tmp_dim;
	}
	
	free( dimids );
}

//attributes.name is used as a handle to get all the metadata from the file
//all other attributes are ignored, attributes.units is overwritten
void ncParameters::initVariableFromFile(const int& ncid, nc_variable& var) const
{
	int nrdims;
	int dimids[NC_MAX_VAR_DIMS];

	int status = nc_inq_varid (ncid, var.attributes.name.c_str(), &var.varid);
	if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
	status = nc_inq_var(ncid, var.varid, NULL, NULL, &nrdims, dimids, NULL);
	if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
	
	getAttribute(ncid, var.varid, var.attributes.name, "_FillValue", var.nodata);
	getAttribute(ncid, var.varid, var.attributes.name, "missing_value", var.nodata);
	getAttribute(ncid, var.varid, var.attributes.name, "scale_factor", var.scale);
	getAttribute(ncid, var.varid, var.attributes.name, "add_offset", var.offset);
	getAttribute(ncid, var.varid, var.attributes.name, "units", var.attributes.units);
	var.dimids.assign(dimids, dimids+nrdims);
	
	if (var.attributes.param==TIME && !wrf_hacks) 
		getTimeTransform(var.attributes.units, TZ, var.offset, var.scale);
}

void ncParameters::initVariablesFromFile(const int& ncid, const std::string& schema_name)
{
	int nr_of_variables = -1;
	int status = nc_inq_nvars(ncid, &nr_of_variables);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve variables for dataset: " + string(nc_strerror(status)), AT);
	
	// Variable IDs in a NetCDF file are consecutive integers starting with 0
	for (int ii=0; ii<nr_of_variables; ++ii) {
		char name[NC_MAX_NAME+1];
		status = nc_inq_varname(ncid, ii, name);
		const std::string varname( name );
		var_attr tmp_attr( getSchemaAttributes(varname, schema_name) );

		//try harder: we try some typical names that are NOT part of the schema (like for DEM)
		if (tmp_attr.param==IOUtils::npos) {
			if (varname=="Band1" || varname=="z" || varname=="height" || varname=="HGT")
				tmp_attr.param = MeteoGrids::DEM;
			else if (varname=="lat")
				tmp_attr.param = LATITUDE;
			else if (varname=="lon")
				tmp_attr.param = LONGITUDE;
			if (vars.count( tmp_attr.param )>0) continue; //if this parameter has already been read, skip it (so the schema naming has priority)
		}

		nc_variable tmp_var( tmp_attr );
		initVariableFromFile(ncid, tmp_var);
		
		if (tmp_var.attributes.param!=IOUtils::npos) {
			vars[ tmp_var.attributes.param ] = tmp_var;
		} else {
			unknown_vars[ tmp_var.attributes.name ] = tmp_var;
		}
	}
}

std::vector<Date> ncParameters::read_1Dvariable(const int& ncid) const
{
	if (!wrf_hacks) {
		const std::vector<double> tmp_results( read_1Dvariable(ncid, TIME) );
		const std::map<size_t, nc_variable>::const_iterator it = vars.find( TIME ); //it exists since it has been read above
		std::vector<Date> results(tmp_results.size());
		for (size_t ii=0; ii<tmp_results.size(); ii++)
			results[ii].setDate(tmp_results[ii]*it->second.scale + it->second.offset, TZ);
		return results;
	} else {
		static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions
		const std::map<size_t, nc_variable>::const_iterator it = vars.find( TIME );
		if (it==vars.end()) throw InvalidArgumentException("Could not find parameter \"TIME\" in file \""+file_and_path+"\"", AT);
		const size_t length = read_1DvariableLength(it->second);
		
		char *data = (char*)calloc(length, sizeof(char)*DateStrLen);
		const int status = nc_get_var_text(ncid, it->second.varid, data);
		if (status != NC_NOERR) throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);

		std::vector<Date> results(length);
		for(size_t ii=0; ii<length; ii++) {
			std::string tmp(DateStrLen, '\0');
			for(size_t jj=0; jj<DateStrLen; jj++)
				tmp[jj] = data[ii*DateStrLen+jj];
			tmp[ 10 ] = 'T';
			IOUtils::convertString(results[ii], tmp, TZ);
		}
		free( data );
		return results;
	}
}

std::vector<double> ncParameters::read_1Dvariable(const int& ncid, const size_t& param) const
{
	const std::map<size_t, nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end()) throw InvalidArgumentException("Could not find parameter \""+getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	const size_t length = read_1DvariableLength(it->second);
	
	std::vector<double> results( length );
	double *data = new double[ length ];
	ncpp::read_data(ncid, it->second.attributes.name, it->second.varid, data);
	std::copy(data, data+length, results.begin());
	delete[] data;
	return results;
}

size_t ncParameters::read_1DvariableLength(const nc_variable& var) const
{
	if (var.dimids.size()!=1) throw InvalidArgumentException("Parameter \""+getParameterName(var.attributes.param)+"\" in file \""+file_and_path+"\" is not a 1D variable", AT);
	
	const int dimid = var.dimids[0];
	std::map<size_t, nc_dimension>::const_iterator it = dimensions_map.begin();
	for (; it!=dimensions_map.end(); ++it) {
		if (it->second.dimid==dimid) break;
	}
	if (it==dimensions_map.end()) throw InvalidArgumentException("Could not find a dimension in file \""+file_and_path+"\"", AT);
	
	return it->second.length;
}

//check that a given dimension exists, has a dimid and an associated variable with a varid
bool ncParameters::hasDimension(const size_t& dim) const
{
	const std::map<size_t, nc_dimension>::const_iterator it_dim = dimensions_map.find( dim );
	if (it_dim==dimensions_map.end()) return false;
	if (it_dim->second.dimid==-1) return false;
	
	const std::map<size_t, nc_variable>::const_iterator it_var = vars.find( dim );
	if (it_var==vars.end()) return false;
	if (it_var->second.varid==-1) return false;
	
	return true;
}

double ncParameters::calculate_cellsize(double& factor_x, double& factor_y) const
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

double ncParameters::calculate_XYcellsize(double& factor_x, double& factor_y) const
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
void ncParameters::fill2DGrid(Grid2DObject& grid, const double data[], const double& nodata) const
{
	if (vecY.front()<=vecY.back()) {
		for (size_t kk=0; kk < vecY.size(); kk++) {
			const size_t row = kk*vecX.size();
			if (vecX.front()<=vecX.back()) {
				for (size_t ll=0; ll < vecX.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < vecX.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (vecX.size() -1) - ll], nodata);
			}
		}
	} else {
		for (size_t kk=0; kk < vecY.size(); kk++) {
			const size_t row = ((vecY.size()-1) - kk)*vecX.size();
			if (vecX.front()<=vecX.back()) {
				for (size_t ll=0; ll < vecX.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < vecX.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (vecX.size() -1) - ll], nodata);
			}
		}
	}
}

void ncParameters::getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier)
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

const ncParameters::var_attr ncParameters::getSchemaAttributes(const std::string& var, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].name==var) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncParameters::var_attr> >::const_iterator it = schemas_vars.find( schema_name );
	if (it==schemas_vars.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==var) return it->second[ii];
	}
	
	return var_attr(var);
}

const ncParameters::nc_dimension ncParameters::getSchemaDimension(const std::string& dimname, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_dimensions.size(); ii++) {
		if (user_dimensions[ii].name==dimname) return user_dimensions[ii];
	}
	
	std::map< std::string, std::vector<ncParameters::nc_dimension> >::const_iterator it = schemas_dims.find( schema_name );
	if (it==schemas_dims.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==dimname) return it->second[ii];
	}
	
	return nc_dimension();
}

//if the attribute is not found, an empty string is returned
void ncParameters::getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, std::string& attr_value)
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
void ncParameters::getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		status = nc_get_att_double(ncid, value_id, attr_name.c_str(), &attr_value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);
	}
}

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
std::string ncParameters::getParameterName(const size_t& param)
{
	if (param==IOUtils::npos) return "";
	
	if (param>=NONE) {
		if (param>lastdimension) 
			throw IndexOutOfBoundsException("Trying to get name for a dimension that does not exist", AT);
		return dimnames[ param - firstdimension ];
	}
	
	return MeteoGrids::getParameterName( param );
}

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
size_t ncParameters::getParameterIndex(const std::string& param)
{
	for (size_t ii=firstdimension; ii<=lastdimension; ii++) {
		if (dimnames[ii]==param) return ii;
	}
	
	return MeteoGrids::getParameterIndex( param );
}

} //namespace
