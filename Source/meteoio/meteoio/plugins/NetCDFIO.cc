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

#define DFLT_STAT_STR_LEN 16

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
 * - NETCDF_SCHEMA: the schema to use (either CF-1 or CNRM or ECMWF or WRF); [Input] and [Output] section (default: ECMWF)
 * - NETCDF_VAR::{MeteoGrids::Parameters} = {netcdf_param_name} : this allows to remap the names as found in the NetCDF file to the MeteoIO grid parameters; [Input] section;
 * - NETCDF_DIM::{MeteoGrids::Parameters} = {netcdf_dimension_name} : this allows to remap the names as found in the NetCDF file to the ncParameters Dimensions; [Input] section;
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
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) 
         : cfg(cfgreader), cache_grid_files(), available_params(), in_schema("ECMWF"), out_schema("ECMWF"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), grid2d_out_file(), 
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), debug(false)
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
		if (!FileUtils::fileExists(filename)) throw NotFoundException(filename, AT); //prevent invalid filenames
		const ncParameters ncFile(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(parameter, date);
	}
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	const std::string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	
	if (!FileUtils::fileExists(filename)) throw NotFoundException(filename, AT);
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

	const ncParameters::Mode file_mode = (FileUtils::fileExists(vec_argument[0]))? ncParameters::READ : ncParameters::WRITE;
	ncParameters ncFile(out_grid2d_path+"/"+vec_argument[0], file_mode, cfg, out_schema, out_dflt_TZ, debug);
	ncFile.write2DGrid(grid_in, IOUtils::npos, vec_argument[1], Date());
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const std::string file_and_path( out_grid2d_path + "/" + grid2d_out_file );
	if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException("Invalid output file name '"+file_and_path+"'", AT);
	
	const ncParameters::Mode file_mode = (FileUtils::fileExists(file_and_path))? ncParameters::READ : ncParameters::WRITE;
	ncParameters ncFile(file_and_path, file_mode, cfg, out_schema, out_dflt_TZ, debug);
	if (parameter==MeteoGrids::DEM || parameter==MeteoGrids::SHADE || parameter==MeteoGrids::SLOPE || parameter==MeteoGrids::AZI)
		ncFile.write2DGrid(grid_in, parameter, std::string(), Date()); //do not assign a date to a DEM?
	else
		ncFile.write2DGrid(grid_in, parameter, std::string(), date);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	const std::string file_and_path( out_meteo_path + "/" + out_meteo_file );
	if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException("Invalid output file name '"+file_and_path+"'", AT);
	
	const ncParameters::Mode file_mode = (FileUtils::fileExists(file_and_path))? ncParameters::READ : ncParameters::WRITE;
	ncParameters ncFile(file_and_path, file_mode, cfg, out_schema, out_dflt_TZ, debug);
	ncFile.writeMeteo(vecMeteo);
}


///////////////////////////////////////////////////// Now the ncParameters class starts //////////////////////////////////////////
std::map< std::string, std::vector<ncpp::nc_dimension> > ncParameters::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncpp::var_attr> > ncParameters::schemas_vars( initSchemasVars() );

std::map< std::string, std::vector<ncpp::nc_dimension> > ncParameters::initSchemasDims()
{
	std::map< std::string, std::vector<ncpp::nc_dimension> > results;
	std::vector<ncpp::nc_dimension> tmp;
	
	//CF-1 schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "surface_altitude") );
	results["CF-1"] = tmp;
	
	//CNRM schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "ZS") );
	results["CNRM"] = tmp;
	
	//ECMWF schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "latitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "longitude") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "geopotential_height") );
	results["ECMWF"] = tmp;
	
	//WRF schema
	tmp.clear();
	tmp.push_back( ncpp::nc_dimension(ncpp::TIME, "Time") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LATITUDE, "south_north") );
	tmp.push_back( ncpp::nc_dimension(ncpp::LONGITUDE, "west_east") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATION, "station") );
	tmp.push_back( ncpp::nc_dimension(ncpp::STATSTRLEN, "station_str_len") );
	tmp.push_back( ncpp::nc_dimension(ncpp::EASTING, "easting") );
	tmp.push_back( ncpp::nc_dimension(ncpp::NORTHING, "northing") );
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "HGT") );
	results["WRF"] = tmp;
	
	return results;
}

std::map< std::string, std::vector<ncpp::var_attr> > ncParameters::initSchemasVars()
{ //HACK: vars/dims should be identified based on standard_name, not name (cf1)
	std::map< std::string, std::vector<ncpp::var_attr> > results;
	std::vector<ncpp::var_attr> tmp;

	//CF-1 schema -> to be checked and improved from CF-1 documentation
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "", "s", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "", "degree_north", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "", "degree_east", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, ncpp::nc_char) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "projection_x_coordinate", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "projection_y_coordinate", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "surface_altitude", "surface_altitude", "height above mean sea level", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "air_temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RH, "relative_humidity", "relative_humidity", "", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DW, "wind_from_direction", "wind_from_direction", "", "degree", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW, "wind_speed", "wind_speed", "", "m/s", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "air_pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "downwelling_shortwave_flux_in_air", "downwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSWR, "upwelling_shortwave_flux_in_air", "upwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "downwelling_longwave_flux_in_air", "downwelling_longwave_flux_in_air", "", "W/m2", IOUtils::nodata, ncpp::nc_double) );
	results["CF-1"] = tmp;

	//CNRM schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, ncpp::nc_char) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "ZS", "", "altitude", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RH, "HUMREL", "", "Relative Humidity", "%", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::QI, "Qair", "", "", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW, "Wind", "", "Wind Speed", "m/s", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DW, "Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM_L, "Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM_S, "Snowf", "", "", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIR, "DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIFF, "SCA_SWdown", "", "", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata, ncpp::nc_double) );
	results["CNRM"] = tmp;

	//ECMWF schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, ncpp::nc_char) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "z", "geopotential_height", "geopotential_height", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "t2m", "", "2 metre temperature", "K", 2., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TD, "d2m", "", "2 metre dewpoint temperature", "K", 2., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "sp", "surface_air_pressure", "Surface pressure", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P_SEA, "msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J m**-2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "strd", "", "Surface thermal radiation downwards", "J m**-2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM, "tp", "", "Total precipitation", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::U, "u10", "", "10 metre U wind component", "m s**-1", 10., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::V, "v10", "", "10 metre V wind component", "m s**-1", 10., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SWE, "sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSS, "skt", "", "Skin temperature", "K", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSG, "stl1", "surface_temperature", "Soil temperature level 1", "K", IOUtils::nodata, ncpp::nc_double) ); //this is from 0 to -7cm
	tmp.push_back( ncpp::var_attr(MeteoGrids::ALB, "al", "surface_albedo", "Albedo", "(0 - 1)", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ALB, "fal", "", "Forecast albedo", "(0 - 1)", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSNO, "rsn", "", "Snow density", "kg m**-3", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ROT, "ro", "", "Runoff", "m", IOUtils::nodata, ncpp::nc_double) );
	results["ECMWF"] = tmp;

	//WRF schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "Times", "Times", "Times", "", IOUtils::nodata, ncpp::nc_char) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "XLAT", "latitude", "latitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "XLONG", "longitude", "longitude", "degrees", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, ncpp::nc_char) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "HGT", "Terrain Height", "Terrain Height", "m", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "PSFC", "Surface pressure", "Surface pressure", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "T2", "2-meter temperature", "2-meter temperature", "K", 2., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::QI, "Q2", "2-meter specific humidity", "2-meter specific humidity", "kg kg-1", 2, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "ACSWDNB", "Downward SW surface radiation", "Downward SW surface radiation", "W m**-2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSWR, "ACSWUPB", "Upwelling Surface Shortwave Radiation", "Upwelling Surface Shortwave Radiation", "W m**-2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "ACLWDNB", "Downward LW surface radiation", "Downward LW surface radiation", "W m**-2", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ROT, "SFROFF", "Surface runoff ", "Surface runoff ", "kg*m2*s-1", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::HS, "SNOWH", "Snow depth", "Snow depth", "Pa", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSS, "TSK", "Surface skin temperature", "Surface skin temperature", "K", IOUtils::nodata, ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::U, "U10", "10-meter wind speed", "10 metre U wind component", "m s**-1", 10., ncpp::nc_double) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::V, "V10", "10-meter wind speed", "10 metre V wind component", "m s**-1", 10., ncpp::nc_double) );
	results["WRF"] = tmp;
	
	return results;
}

//The user can provide his own variables properties as NETCDF_VAR::{param} = {name}
std::vector<ncpp::var_attr> ncParameters::initUserSchemas(const Config& i_cfg)
{
	std::vector<ncpp::var_attr> results;
	const std::string user_type_str = i_cfg.get("NC_TYPE", "Input", IOUtils::nothrow);
	ncpp::types user_type = ncpp::nc_none;
	if (!user_type_str.empty()) {
		if (user_type_str=="DOUBLE") user_type = ncpp::nc_double;
		else if (user_type_str=="FLOAT") user_type = ncpp::nc_float;
		else if (user_type_str=="INT") user_type = ncpp::nc_int;
		else
			throw InvalidArgumentException("Unknown NC_TYPE value "+user_type_str, AT);
	}
	
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_VAR::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string meteo_grid( custom_attr[ii].substr(found+1) );
		const std::string netcdf_param = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = ncpp::getParameterIndex(meteo_grid);
		if (param_index==IOUtils::npos)
			throw InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( ncpp::var_attr(param_index, netcdf_param, IOUtils::nodata, user_type) );
	}
	
	return results;
}

//The user can provide his own dimensions properties as NETCDF_DIM::{dimension_param} = {name_in_current_file}
std::vector<ncpp::nc_dimension> ncParameters::initUserDimensions(const Config& i_cfg)
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
		if (param_index==IOUtils::npos || param_index<ncpp::firstdimension || param_index>ncpp::lastdimension)
			throw InvalidArgumentException("Dimension '"+dim_str+"' is not a valid dimension! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( ncpp::nc_dimension( static_cast<ncpp::Dimensions>(param_index), netcdf_dim) );
	}
	
	return results;
}

//TODO: redo the whole user_schema thing: we should fill vars / dimensions with the schema, then add/overwrite with the user schema
ncParameters::ncParameters(const std::string& filename, const Mode& mode, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug)
             : user_schemas( initUserSchemas(cfg) ), user_dimensions( initUserDimensions(cfg) ), vars(), unknown_vars(), vecTime(), vecX(), vecY(), dimensions_map(), file_and_path(filename), current_schema(schema), coord_sys(), coord_param(), TZ(tz_in), wrf_hacks(schema=="WRF"), debug(i_debug), isLatLon(false)
{
	IOUtils::getProjectionParameters(cfg, coord_sys, coord_param);
	initFromSchema(schema);
	
	if (mode==WRITE) {
		if (FileUtils::fileExists(filename)) initFromFile(filename, schema);
	} else if (mode==READ) {
		initFromFile(filename, schema);
	}
	
	if (debug) {
		std::cout << filename << ":\n";
		std::cout << "\tDimensions:\n";
		for (std::map<size_t, ncpp::nc_dimension>::const_iterator it = dimensions_map.begin(); it!=dimensions_map.end(); ++it)
			std::cout << "\t\t" << it->second.toString() << "\n";
		if (!vecTime.empty()) std::cout << "\ttime range: [" << vecTime.front().toString(Date::ISO) << " - " << vecTime.back().toString(Date::ISO) << "]\n";
		std::cout << "\tVariables:\n";
		for (std::map<size_t, ncpp::nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
			std::cout << "\t\t" << ncpp::getParameterName( it->first ) << " -> " << it->second.toString() << "\n";
		std::cout << "\tUnrecognized variables:\n";
		for (std::map<std::string, ncpp::nc_variable>::const_iterator it=unknown_vars.begin(); it!=unknown_vars.end(); ++it)
			std::cout << "\t\t" << it->first << " -> " << it->second.toString() << "\n";
	}
	
	//check that the used schema has declared the minimum required dimensions and potentially, associated variables (if based on schema)
	const bool hasLatLon = ((dimensions_map.count(ncpp::LATITUDE)!=0 && dimensions_map.count(ncpp::LONGITUDE)!=0) && (mode!=WRITE || (vars.count(ncpp::LATITUDE)!=0 && vars.count(ncpp::LONGITUDE)!=0)));
	const bool hasEastNorth = ((dimensions_map.count(ncpp::EASTING)!=0 && dimensions_map.count(ncpp::NORTHING)!=0) && (mode!=WRITE || (vars.count(ncpp::EASTING)!=0 && vars.count(ncpp::NORTHING)!=0)));
	const bool hasTime = dimensions_map.count(ncpp::TIME)!=0 && (mode!=WRITE || vars.count(ncpp::TIME)!=0);
	if (!hasLatLon || !hasEastNorth || !hasTime) throw IOException("Error in the schema definition, some basic quantities are not defined!", AT);
}

//populate the dimensions_map from the selected schema
void ncParameters::initFromSchema(const std::string& schema)
{
	for (size_t ii=0; ii<schemas_dims[schema].size(); ii++) {
		dimensions_map[ schemas_dims[schema][ii].type ] = schemas_dims[schema][ii];
	}
	if (dimensions_map.count(ncpp::TIME)==0) throw IOException("No TIME dimension in schema '"+schema+"'", AT);
	dimensions_map[ ncpp::TIME ].isUnlimited = true;
	
	for (size_t ii=0; ii<schemas_vars[schema].size(); ii++) {
		vars[ schemas_vars[schema][ii].param ] = ncpp::nc_variable( schemas_vars[schema][ii] );
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

	isLatLon = hasDimension(ncpp::LATITUDE) && hasDimension(ncpp::LONGITUDE);
	const bool isXY = hasDimension(ncpp::EASTING) && hasDimension(ncpp::NORTHING);
	if (isLatLon) {
		vecY = read_1Dvariable(ncid, ncpp::LATITUDE);
		vecX = read_1Dvariable(ncid, ncpp::LONGITUDE);
	} else if (isXY) {
		vecX = read_1Dvariable(ncid, ncpp::EASTING);
		vecY = read_1Dvariable(ncid, ncpp::NORTHING);
	}
	if (hasDimension(ncpp::TIME)) vecTime = read_1Dvariable(ncid);
	
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
	for (std::map<size_t, ncpp::nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
		available_params.insert( it->first );
	
	return available_params;
}

Grid2DObject ncParameters::read2DGrid(const std::string& varname) const
{
	for (std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.begin(); it!=vars.end(); ++it) {
		if (it->second.attributes.name==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	for (std::map<std::string, ncpp::nc_variable>::const_iterator it = unknown_vars.begin(); it!=unknown_vars.end(); ++it) {
		if (it->first==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	throw NotFoundException("The variable '"+varname+"' could not be found in file '"+file_and_path+"'", AT);
}

Grid2DObject ncParameters::read2DGrid(const size_t& param, const Date& date) const
{
	const std::map <size_t, ncpp::nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end() || it->second.varid==-1) 
		NoDataException("No "+MeteoGrids::getParameterName( param )+" grid in file "+file_and_path, AT);
	
	size_t time_pos = IOUtils::npos;
	if (!date.isUndef()) {
		const std::vector<Date>::const_iterator low = std::lower_bound(vecTime.begin(), vecTime.end(), date);
		if (*low!=date) throw NoDataException("No data at "+date.toString(Date::ISO)+" in file "+file_and_path, AT);
		time_pos = static_cast<size_t>( std::distance(vecTime.begin(), low) );
	} else { 
		//no date has been provided, check if this parameter depends on time
		const std::map<size_t, ncpp::nc_dimension>::const_iterator it2 = dimensions_map.find(ncpp::TIME);
		const int time_id = (it2!=dimensions_map.end())? it2->second.dimid : -1;
		const bool depend_on_time = (std::find(it->second.dimids.begin(), it->second.dimids.end(), time_id) != it->second.dimids.end());
		if (depend_on_time && vecTime.size()>1) //if only one timestep is present, we take it
			throw InvalidFormatException("No time requirement has been provided for a file that contains multiple timestamps", AT);
	}
	
	const bool isPrecip = (param==MeteoGrids::PSUM || param==MeteoGrids::PSUM_L || param==MeteoGrids::PSUM_S);
	const bool isRad = (param==MeteoGrids::ISWR || param==MeteoGrids::RSWR || param==MeteoGrids::ISWR_DIFF || param==MeteoGrids::ISWR_DIR);
	return read2DGrid(it->second, time_pos, isPrecip, (isPrecip || isRad));
}

Grid2DObject ncParameters::read2DGrid(const ncpp::nc_variable& var, const size_t& time_pos, const bool& m2mm, const bool& reZero) const
{
	if (isLatLon && (!hasDimension(ncpp::LATITUDE) || !hasDimension(ncpp::LONGITUDE))) throw IOException("No latitude / longitude could be identified in file "+file_and_path, AT);
	if (!isLatLon && (!hasDimension(ncpp::EASTING) || !hasDimension(ncpp::NORTHING))) throw IOException("No easting / northing could be identified in file "+file_and_path, AT);
	
	//define the results grid
	mio::Coords llcorner(coord_sys, coord_param); //if an EPSG was provided, this has been converted to coord_sys/coord_param
	if (isLatLon)
		llcorner.setLatLon( std::min(vecY.front(), vecY.back()), std::min(vecX.front(), vecX.back()), IOUtils::nodata);
	else
		llcorner.setXY( std::min(vecY.front(), vecY.back()), std::min(vecX.front(), vecX.back()), IOUtils::nodata);
	//HACK expand the definition of Grid2DObject to support lat/lon grids and reproject in GridsManager
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = (isLatLon)? ncpp::calculate_cellsize(resampling_factor_x, resampling_factor_y, vecX, vecY) : ncpp::calculate_XYcellsize(resampling_factor_x, resampling_factor_y, vecX, vecY);
	Grid2DObject grid(vecX.size(), vecY.size(), cellsize, llcorner);
	
	//read the raw data, copy it into the Grid2DObject
	int ncid;
	ncpp::open_file(file_and_path, NC_NOWRITE, ncid);
	
	double *data = new double[vecY.size()*vecX.size()];
	if (time_pos!=IOUtils::npos)
		ncpp::read_data(ncid, var.attributes.name, var.varid, time_pos, vecY.size(), vecX.size(), data);
	else
		ncpp::read_data(ncid, var.attributes.name, var.varid, data);
	ncpp::fill2DGrid(grid, data, var.nodata, (vecX.front()<=vecX.back()), (vecY.front()<=vecY.back()) );
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

//this should be most often used as wrapper, to select the proper parameters for a given param or param_name
//If both are provided, param has the priority
void ncParameters::write2DGrid(const Grid2DObject& grid_in, size_t param, std::string param_name, const Date& date)
{
	if (!param_name.empty() && param==IOUtils::npos) param = MeteoGrids::getParameterIndex( param_name );
	
	if (param!=IOUtils::npos && vars.count( param )>0) {
		write2DGrid(grid_in, vars[param], date);
	} else {
		if (param!=IOUtils::npos && param_name.empty()) param_name = MeteoGrids::getParameterName(param);
		ncpp::var_attr tmp_attr( getSchemaAttributes(param, current_schema) );
		if (tmp_attr.name.empty()) tmp_attr.name = param_name; //ie it was not found
		ncpp::nc_variable tmp_var(tmp_attr);
		if (unknown_vars.count( param_name )>0) tmp_var = unknown_vars[ param_name ];
		
		write2DGrid(grid_in, tmp_var, date);
		unknown_vars[ param_name ] = tmp_var; //save varid, etc for next calls
	}
}

void ncParameters::write2DGrid(const Grid2DObject& grid_in, ncpp::nc_variable& var, const Date& date)
{
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		ncpp::file_redef(file_and_path, ncid);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", current_schema);
		ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO "+getLibVersion());
		if (!isLatLon) ncpp::add_attribute(ncid, NC_GLOBAL, "epsg", grid_in.llcorner.getEPSG());
	}
	
	//create any potentially missing definition, otherwise check that everything is consistent
	std::vector<size_t> nc_variables, dimensions;
	if (!date.isUndef()) dimensions.push_back( ncpp::TIME );
	if (isLatLon) {
		dimensions.push_back( ncpp::LATITUDE );
		dimensions.push_back( ncpp::LONGITUDE );
	} else {
		dimensions.push_back( ncpp::NORTHING );
		dimensions.push_back( ncpp::EASTING );
	}
	
	for (size_t ii=0; ii<dimensions.size(); ii++) {
		const size_t param = dimensions[ii];
		size_t length = 0;
		if (param==ncpp::LATITUDE || param==ncpp::NORTHING) length = grid_in.getNy();
		if (param==ncpp::LONGITUDE || param==ncpp::EASTING) length = grid_in.getNx();
		ncpp::createDimension(ncid, dimensions_map[ param ], length);
		if (setAssociatedVariable(ncid, param, date)) nc_variables.push_back( dimensions[ii] ); //associated variable will have to be filled
		if (var.varid == -1) var.dimids.push_back( dimensions_map[param].dimid );
	}
	if (var.varid == -1) ncParameters::create_variable(ncid, var); //create the "main" variable if necessary
	nc_variables.push_back( var.attributes.param );
	vars[ var.attributes.param ] = var;
	
	ncpp::end_definitions(file_and_path, ncid);
	
	//now write the data
	const size_t time_pos = (!date.isUndef())? addTimestamp(ncid, date) : IOUtils::npos;
	for (size_t ii=0; ii<nc_variables.size(); ii++) {
		const size_t param = nc_variables[ii];
		if (param==ncpp::TIME) continue; //this was done above
		
		const std::vector<double> data( fillBufferForVar(grid_in, vars[ param ]) );
		if (data.empty()) continue;
		
		if (vars[ param ].dimids.size()>0 && vars[ param ].dimids.front()==ncpp::TIME) { //as unlimited dimension, TIME is always first
			ncpp::write_data(ncid, vars[ param ].attributes.name, vars[ param ].varid, grid_in.getNy(), grid_in.getNx(), time_pos, &data[0]);
		} else {
			ncpp::write_data(ncid, vars[ param ], data, false);
		}
	}
	
	ncpp::close_file(file_and_path, ncid);
}

//HACK this assumes that all stations have the same parameters, the same timestamps and the same coordinate system
void ncParameters::writeMeteo(const std::vector< std::vector<MeteoData> >& vecMeteo/*, const size_t& station_idx*/) //HACK
{
	if (vecMeteo.empty()) return;
	if (vecMeteo.front().empty()) return;
	isLatLon = true; //HACK
	
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		ncpp::file_redef(file_and_path, ncid);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", current_schema);
		ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO "+getLibVersion());
		if (!isLatLon) ncpp::add_attribute(ncid, NC_GLOBAL, "epsg", vecMeteo.front().front().meta.position.getEPSG());
	}
	
	std::vector<size_t> nc_variables, dimensions;
	dimensions.push_back( ncpp::TIME );
	dimensions.push_back( ncpp::STATSTRLEN ); //this MUST be before STATION
	dimensions.push_back( ncpp::STATION );
	const Date ref_date( vecMeteo.front().front().date ); //TODO: getRefDate(vecMeteo)
	for (size_t ii=0; ii<dimensions.size(); ii++) {
		const size_t param = dimensions[ii];
		const size_t length = (param==ncpp::TIME)? 0 : ((param==ncpp::STATSTRLEN)? DFLT_STAT_STR_LEN : vecMeteo.size()) ;
		ncpp::createDimension(ncid, dimensions_map[ param ], length);
		
		if (param==ncpp::STATSTRLEN) continue; //no associated variable for STATSTRLEN
		if (setAssociatedVariable(ncid, param, ref_date)) nc_variables.push_back( dimensions[ii] ); //associated variable will have to be filled
	}
	
	appendVariablesList(nc_variables, vecMeteo);
	
	//associate dimids to vars and create the variables
	for (size_t ii=0; ii<nc_variables.size(); ii++) {
		const size_t param = nc_variables[ii];
		if (vars[ param ].varid == -1) { //skip existing nc_variables
			if (param<ncpp::firstdimension && param!=MeteoGrids::DEM) vars[ param ].dimids.push_back( dimensions_map[ncpp::TIME].dimid );
			vars[ param ].dimids.push_back( dimensions_map[ncpp::STATION].dimid );
			if (param==ncpp::STATION) vars[ param ].dimids.push_back( dimensions_map[ncpp::STATSTRLEN].dimid );
			
			ncParameters::create_variable(ncid, vars[ param ]);
		}
	}
	
	ncpp::end_definitions(file_and_path, ncid);
	
	//write data: fill associated nc_variables and normal nc_variables
	for (size_t ii=0; ii<nc_variables.size(); ii++) {
		const size_t param = nc_variables[ii];
		
		if (param==ncpp::STATION) {
			std::vector<std::string> txtdata( vecMeteo.size() );
			for (size_t jj=0; jj<vecMeteo.size(); jj++) txtdata[jj] = vecMeteo[jj].front().meta.stationID;
			ncpp::write_data(ncid, vars[param], txtdata, DFLT_STAT_STR_LEN);
		} else {
			const std::vector<double> data( fillBufferForVar(vecMeteo, vars[ param ]) );
			if (data.empty()) continue;
			ncpp::write_data(ncid, vars[ param ], data, (param==ncpp::TIME));
		}
	}
	
	ncpp::close_file(file_and_path, ncid);
}

void ncParameters::appendVariablesList(std::vector<size_t> &nc_variables, const std::vector< std::vector<MeteoData> >& vecMeteo)
{
	//add metadata variables
	nc_variables.push_back( MeteoGrids::DEM );
	if (isLatLon) {
		nc_variables.push_back( ncpp::LATITUDE );
		nc_variables.push_back( ncpp::LONGITUDE );
	} else {
		nc_variables.push_back( ncpp::EASTING );
		nc_variables.push_back( ncpp::NORTHING );
	}
	
	//add all vars found in vecMeteo, assuming that the first timestep for any given station has all the parameters for that station
	for (size_t st=0; st<vecMeteo.size(); st++) {
		const MeteoData md( vecMeteo[st].front() );
		
		for (size_t ii=0; ii<md.getNrOfParameters(); ii++) { //TODO check which parameters do not have any values!
			const size_t param = MeteoGrids::getParameterIndex( md.getNameForParameter( ii ) );
			if (param>=MeteoGrids::nrOfParameters) continue;
			if (vars.count(param)==0) { //ie unrecognized in loaded schema, adding it
				const ncpp::var_attr tmp_attr(param, md.getNameForParameter( ii ), IOUtils::nodata, ncpp::nc_double); //using NC_DOUBLE as default
				vars[param] = ncpp::nc_variable(tmp_attr);
			}
			if (std::find(nc_variables.begin(), nc_variables.end(), param) == nc_variables.end())
				nc_variables.push_back( param );
		}
	}
}

bool ncParameters::setAssociatedVariable(const int& ncid, const size_t& param, const Date& ref_date)
{
	if (vars[ param ].varid == -1) {
		vars[ param ].dimids.push_back( dimensions_map[ param ].dimid );
		if (param==ncpp::STATION) vars[ param ].dimids.push_back( dimensions_map[ ncpp::STATSTRLEN ].dimid );
		if (param==ncpp::TIME) {
			const Date ref_date_simplified(ref_date.getYear(), 1, 1, 0, 0, 0.); //we force data into GMT
			std::string date_str( ref_date_simplified.toString(Date::ISO) );
			date_str[ 10 ] = ' '; //replace "T" by " "
			vars[param].attributes.units = "hours since " + date_str;
			vars[param].offset = ref_date_simplified.getJulian();
			vars[param].scale = 1./24.;
		}
		
		ncParameters::create_variable(ncid, vars[ param ]);
		return true;
	}
	return false;
}

const std::vector<double> ncParameters::fillBufferForVar(const std::vector< std::vector<MeteoData> >& vecMeteo, ncpp::nc_variable& var)
{
	const size_t param = var.attributes.param;
	
	if (param>=ncpp::firstdimension || param==MeteoGrids::DEM) { //associated nc_variables
		if (param==ncpp::TIME) {
			std::vector<double> data(vecMeteo.front().size());
			for (size_t ll=0; ll<vecMeteo.front().size(); ll++)
				data[ll] = (vecMeteo.front()[ll].date.getJulian() - var.offset) / var.scale;
			return data;
		} else {
			std::vector<double> data(vecMeteo.size());
			for (size_t jj=0; jj<vecMeteo.size(); jj++) {
				if (param==MeteoGrids::DEM) {
					data[jj] = vecMeteo[jj].front().meta.position.getAltitude();
				} else if (param==ncpp::EASTING) {
					data[jj] = vecMeteo[jj].front().meta.position.getEasting();
				} else if (param==ncpp::NORTHING) {
					data[jj] = vecMeteo[jj].front().meta.position.getNorthing();
				} else if (param==ncpp::LATITUDE) {
					data[jj] = vecMeteo[jj].front().meta.position.getLat();
				} else if (param==ncpp::LONGITUDE) {
					data[jj] = vecMeteo[jj].front().meta.position.getLon();
				} else
					throw UnknownValueException("Unknown dimension found when trying to write out netcdf file", AT);
			}
			return data;
		}
	} else { //normal nc_variables
		const MeteoData md( vecMeteo.front().front() );
		const size_t meteodata_param = md.getParameterIndex( MeteoGrids::getParameterName( param ) ); //retrieve the equivalent parameter in vecMeteo
		if (meteodata_param==IOUtils::npos) return std::vector<double>(); //this should not have happened...
		
		const size_t nrTimeSteps = vecMeteo.front().size();
		const size_t nrStations = vecMeteo.size();
		std::vector<bool> stationHasParameter(nrStations, false);
		for (size_t jj=0; jj<nrStations; jj++) {
			if(!vecMeteo[jj].empty()) stationHasParameter[jj] = (vecMeteo[jj].front().getNrOfParameters()>meteodata_param);
		}
		
		std::vector<double> data(nrTimeSteps*nrStations, IOUtils::nodata);
		for (size_t ll=0; ll<nrTimeSteps; ll++) {
			for (size_t jj=0; jj<nrStations; jj++) {
				if (stationHasParameter[jj]) data[ll*nrStations+jj] = vecMeteo[jj][ll]( meteodata_param );
			}
		}
		return data;
	}
}

const std::vector<double> ncParameters::fillBufferForVar(const Grid2DObject& grid, ncpp::nc_variable& var)
{
	const size_t param = var.attributes.param;
	
	if (param>=ncpp::firstdimension && param!=IOUtils::npos) { //associated nc_variables
		const double cellsize = grid.cellsize;
		const size_t nrPts = (param==ncpp::NORTHING)? grid.getNy() : grid.getNx();
		if (nrPts==0) return std::vector<double>(); //this should not have happened...
		std::vector<double> data( nrPts );
		
		if (param==ncpp::NORTHING || param==ncpp::EASTING) {
			data[0] = (param==ncpp::NORTHING)? grid.llcorner.getNorthing() : grid.llcorner.getEasting();
			for (size_t ii=1; ii<nrPts; ii++)
				data[ii] = data[0] + static_cast<double>(ii)*cellsize;
		} else if (param==ncpp::LATITUDE || param==ncpp::LONGITUDE) {
			//this is (very cheap) approximation of some kind of projection from x/y to lat/lon
			//There is a trick here: walking along a line of constant northing does NOT lead to a constant latitude. Both grids
			//are shifted (even if a little), which means that the center of lat/lon is != center of east./north..
			//So, in order to find the center of the domain, we do a few iteration to converge toward a reasonnable approximation
			double alpha;
			double lat_length, lon_length, cntr_lat=grid.llcorner.getLat(), cntr_lon=grid.llcorner.getLon();
			for(size_t ii=0; ii<5; ii++) {
				lat_length = CoordsAlgorithms::VincentyDistance(cntr_lat-.5, cntr_lon, cntr_lat+.5, cntr_lon, alpha);
				lon_length = CoordsAlgorithms::VincentyDistance(cntr_lat, cntr_lon-.5, cntr_lat, cntr_lon+.5, alpha);
				cntr_lat = (.5*static_cast<double>(grid.getNy())*grid.cellsize) / lat_length + grid.llcorner.getLat();
				cntr_lon = (.5*static_cast<double>(grid.getNx())*grid.cellsize) / lon_length + grid.llcorner.getLon();
			}
			
			const double center = (param==ncpp::NORTHING)? cntr_lat : cntr_lon;
			const double length = (param==ncpp::NORTHING)? lat_length : lon_length;
			const double min = center - (.5*static_cast<double>(nrPts)*cellsize) / length;
			const double max = center + (.5*static_cast<double>(nrPts)*cellsize) / length;
			const double interval =  abs(max - min);
			
			for (size_t ii=0; ii<nrPts; ii++)
				data[ii] = min + (interval * static_cast<double>(ii)) / (static_cast<double>(nrPts)-1.);
		} else
			throw UnknownValueException("Unsupported dimension '"+ncpp::getParameterName(param)+"' found when trying to write out netcdf file", AT);
		return data;
	} else { //normal grid variable
		std::vector<double> data( grid.size() );
		const size_t nrows = grid.getNy(), ncols = grid.getNx();
		for (size_t kk=0; kk<nrows; ++kk) {
			for (size_t ll=0; ll<ncols; ++ll) {
				data[kk*ncols + ll] = grid.grid2D(ll,kk);
			}
		}
		return data;
	}
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
		const double dt = (date.getJulian() - vars[ncpp::TIME].offset) / vars[ncpp::TIME].scale;
		const size_t start[] = {time_pos};
		const size_t count[] = {1};
		const int status = nc_put_vara_double(ncid, vars[ncpp::TIME].varid, start, count, &dt);
		if (status != NC_NOERR) throw IOException("Could not write data for record variable '" + vars[ncpp::TIME].attributes.name + "': " + nc_strerror(status), AT);
		if (time_pos==vecTime.size())
			vecTime.push_back( date );
		else
			vecTime.insert(vecTime.begin()+time_pos, date);
	}
	return time_pos;
}

//write the variable's attributes into the file (does nothing if the variable already exists)
void ncParameters::create_variable(const int& ncid, ncpp::nc_variable& var)
{
	if (var.varid != -1) return; //the variable already exists
	const int ndims = static_cast<int>( var.dimids.size() );
	const int status = nc_def_var(ncid, var.attributes.name.c_str(), var.attributes.type, ndims, &var.dimids[0], &var.varid);
	if (status != NC_NOERR) throw IOException("Could not define variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	
	if (!var.attributes.standard_name.empty()) ncpp::add_attribute(ncid, var.varid, "standard_name", var.attributes.standard_name);
	if (!var.attributes.long_name.empty()) ncpp::add_attribute(ncid, var.varid, "long_name", var.attributes.long_name);
	if (!var.attributes.units.empty()) ncpp::add_attribute(ncid, var.varid, "units", var.attributes.units);
	if (var.attributes.param!=ncpp::STATION) ncpp::add_attribute(ncid, var.varid, "_FillValue", var.nodata);
	
	if (var.attributes.param==ncpp::TIME) ncpp::add_attribute(ncid, var.varid, "calendar", "gregorian");
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
		
		ncpp::nc_dimension tmp_dim( getSchemaDimension(dimname, schema_name) ); //set name and type
		if (tmp_dim.type==IOUtils::npos) { //unrecognized dimension -> try harder with some typical names that are not in the schema
			if (dimname=="lat")
				tmp_dim.type = ncpp::LATITUDE;
			else if (dimname=="lon")
				tmp_dim.type = ncpp::LONGITUDE;
			else continue;
			if (dimensions_map.count( tmp_dim.type )>0) {
				//if this parameter has already been read, skip it (so the schema naming has priority)
				if (dimensions_map[ tmp_dim.type ].dimid != -1) continue;
			}
		}

		tmp_dim.dimid = idx;
		tmp_dim.isUnlimited = (idx==unlim_id);
		status = nc_inq_dimlen(ncid, dimids[idx], &tmp_dim.length);
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension lenght: " + std::string(nc_strerror(status)), AT);
		
		dimensions_map[ tmp_dim.type ] = tmp_dim;
	}
	
	free( dimids );
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
		ncpp::var_attr tmp_attr( getSchemaAttributes(varname, schema_name) );

		//try harder: we try some typical names that are NOT part of the schema (like for DEM)
		if (tmp_attr.param==IOUtils::npos) {
			if (varname=="Band1" || varname=="z" || varname=="height" || varname=="HGT")
				tmp_attr.param = MeteoGrids::DEM;
			else if (varname=="lat")
				tmp_attr.param = ncpp::LATITUDE;
			else if (varname=="lon")
				tmp_attr.param = ncpp::LONGITUDE;
			if (vars.count( tmp_attr.param )>0) {
				if (vars[ tmp_attr.param ].varid != -1)
					continue; //if this parameter has already been read, skip it (so the schema naming has priority)
			}
		}

		ncpp::nc_variable tmp_var( tmp_attr );
		const bool readTimeTransform = (tmp_var.attributes.param==ncpp::TIME && !wrf_hacks);
		ncpp::readVariableMetadata(ncid, tmp_var, readTimeTransform, TZ);
		
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
		const std::vector<double> tmp_results( read_1Dvariable(ncid, ncpp::TIME) );
		const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( ncpp::TIME ); //it exists since it has been read above
		std::vector<Date> results(tmp_results.size());
		for (size_t ii=0; ii<tmp_results.size(); ii++)
			results[ii].setDate(tmp_results[ii]*it->second.scale + it->second.offset, TZ);
		return results;
	} else {
		static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions
		const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( ncpp::TIME );
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
	const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end() || it->second.varid==-1) throw InvalidArgumentException("Could not find parameter \""+ncpp::getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	const size_t length = read_1DvariableLength(it->second);
	
	std::vector<double> results( length );
	double *data = new double[ length ];
	ncpp::read_data(ncid, it->second.attributes.name, it->second.varid, data);
	std::copy(data, data+length, results.begin());
	delete[] data;
	return results;
}

size_t ncParameters::read_1DvariableLength(const ncpp::nc_variable& var) const
{
	if (var.dimids.size()!=1) throw InvalidArgumentException("Parameter \""+ncpp::getParameterName(var.attributes.param)+"\" in file \""+file_and_path+"\" is not a 1D variable", AT);
	
	const int dimid = var.dimids[0];
	std::map<size_t, ncpp::nc_dimension>::const_iterator it = dimensions_map.begin();
	for (; it!=dimensions_map.end(); ++it) {
		if (it->second.dimid==dimid) break;
	}
	if (it==dimensions_map.end()) throw InvalidArgumentException("Could not find a dimension in file \""+file_and_path+"\"", AT);
	
	return it->second.length;
}

//check that a given dimension exists, has a dimid and an associated variable with a varid
bool ncParameters::hasDimension(const size_t& dim) const
{
	const std::map<size_t, ncpp::nc_dimension>::const_iterator it_dim = dimensions_map.find( dim );
	if (it_dim==dimensions_map.end()) return false;
	if (it_dim->second.dimid==-1) return false;
	
	const std::map<size_t, ncpp::nc_variable>::const_iterator it_var = vars.find( dim );
	if (it_var==vars.end()) return false;
	if (it_var->second.varid==-1) return false;
	
	return true;
}

const ncpp::var_attr ncParameters::getSchemaAttributes(const std::string& var, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].name==var) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncpp::var_attr> >::const_iterator it = schemas_vars.find( schema_name );
	if (it==schemas_vars.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==var) return it->second[ii];
	}
	
	return ncpp::var_attr(var);
}

const ncpp::var_attr ncParameters::getSchemaAttributes(const size_t& param, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].param==param) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncpp::var_attr> >::const_iterator it = schemas_vars.find( schema_name );
	if (it==schemas_vars.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].param==param) return it->second[ii];
	}
	
	return ncpp::var_attr();
}

const ncpp::nc_dimension ncParameters::getSchemaDimension(const std::string& dimname, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_dimensions.size(); ii++) {
		if (user_dimensions[ii].name==dimname) return user_dimensions[ii];
	}
	
	std::map< std::string, std::vector<ncpp::nc_dimension> >::const_iterator it = schemas_dims.find( schema_name );
	if (it==schemas_dims.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==dimname) return it->second[ii];
	}
	
	return ncpp::nc_dimension();
}

} //namespace
