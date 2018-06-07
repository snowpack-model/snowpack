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
#include <sstream>
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
 * - CROCUS - from the <A HREF="http://www.cnrm.meteo.fr/">National Centre for Meteorological Research</A>;
 * - AMUNDSEN - from the <A HREF="https://geographie.uibk.ac.at/blog/ahc/models/">Alpine, Hydro, climatology</A> group in Innsbruck;
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
 * - NETCDF_SCHEMA: the schema to use (either CF-1.6, CROCUS, AMUNDSEN,  ECMWF or WRF); [Input] and [Output] section (default: CF-1.6)
 * - NETCDF_VAR::{MeteoGrids::Parameters} = {netcdf_param_name} : this allows to remap the names as found in the NetCDF file to the MeteoIO grid parameters; [Input] section;
 * - NETCDF_DIM::{MeteoGrids::Parameters} = {netcdf_dimension_name} : this allows to remap the names as found in the NetCDF file to the ncParameters Dimensions; [Input] section;
 * - NC_SINGLE_FILE: when writing timeseries of station data, force all stations to be contained in a single file (default: false)
 *
 * @note When providing multiple files in one directory, in case of overlapping files (because each file can provide multiple timestamps), the file containing the newest data has priority. This is
 * convenient when using forecats data to automatically use the most short-term forecast.
 * @note When using the CROCUS schema, please note that the humidity should be provided as specific humidity, so please use a data 
 * creator / generator if needed to get a QI parameter (see HumidityGenerator).
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
 *
 * If you need to edit netCDF attributes (for example, to add your metadata), you can use <A href="http://nco.sourceforge.net/">ncatted</A> to do so. For example, to
 * replace the previous global attribute <i>summary</i> by yours or to append a line to the <i>history</i> global attribute:
 * @code
 * ncatted -a summary,global,o,c,"My summary" myFile.nc	#Overwrite the summary field
 * ncatted -a history,global,a,c,"Edited by me on 2018-06-06\n" myFile.nc	#Append a line to the history field
 * @endcode
 */

//helper function to sort the cache of grid files
inline bool sort_cache_grids(const std::pair<std::pair<Date,Date>,ncParameters> &left, const std::pair<std::pair<Date,Date>,ncParameters> &right) {
	if (left.first.first < right.first.first) return true;
	if (left.first.first > right.first.first) return false;
	return left.first.second < right.first.second; //date_start equallity case
}

NetCDFIO::NetCDFIO(const std::string& configfile) 
         : cfg(configfile), cache_grid_files(), available_params(), in_schema("CF-1.6"), out_schema("CF-1.6"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), grid2d_out_file(), 
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), debug(false), out_single_file(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) 
         : cfg(cfgreader), cache_grid_files(), available_params(), in_schema("CF-1.6"), out_schema("CF-1.6"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), grid2d_out_file(), 
         out_meteo_path(), in_meteo_path(), out_meteo_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), debug(false), out_single_file(false)
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
		cfg.getValue("NC_SINGLE_FILE", "Output", out_single_file, IOUtils::nothrow);
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
	
	//TODO handle the case of file_start & file_end are undef() (example: DEM)
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

	const std::string file_and_path( out_grid2d_path + "/" + vec_argument[0] );
	if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException("Invalid output file name '"+file_and_path+"'", AT);
	
	const ncParameters::Mode file_mode = (FileUtils::fileExists(file_and_path))? ncParameters::READ : ncParameters::WRITE;
	ncParameters ncFile(file_and_path, file_mode, cfg, out_schema, out_dflt_TZ, debug);
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
	if (out_single_file) {
		const std::string file_and_path( out_meteo_path + "/" + out_meteo_file );
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException("Invalid output file name '"+file_and_path+"'", AT);
		
		const ncParameters::Mode file_mode = (FileUtils::fileExists(file_and_path))? ncParameters::READ : ncParameters::WRITE;
		ncParameters ncFile(file_and_path, file_mode, cfg, out_schema, out_dflt_TZ, debug);
		ncFile.writeMeteo(vecMeteo);
	} else {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			const std::string file_and_path( out_meteo_path + "/" + vecMeteo[ii].front().meta.stationID + ".nc" );
			if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException("Invalid output file name '"+file_and_path+"'", AT);
			if (FileUtils::fileExists(file_and_path)) throw IOException("Appending data to timeseries is currently non-functional for NetCDF, please delete file "+file_and_path, AT);
			
			const ncParameters::Mode file_mode = (FileUtils::fileExists(file_and_path))? ncParameters::READ : ncParameters::WRITE;
			ncParameters ncFile(file_and_path, file_mode, cfg, out_schema, out_dflt_TZ, debug);
			ncFile.writeMeteo(vecMeteo, ii);
		}
	}
}


///////////////////////////////////////////////////// Now the ncParameters class starts //////////////////////////////////////////
std::map< std::string, std::vector<ncpp::nc_dimension> > ncParameters::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncpp::var_attr> > ncParameters::schemas_vars( initSchemasVars() );

void ncParameters::initSchemaCst(const std::string& schema)
{
	if (schema=="CF-1.6") {
		schema_dflt_type = NC_FLOAT;
	} else if (schema=="CROCUS") {
		schema_dflt_type = NC_DOUBLE;
		schema_nodata =  -9999999.; //CNRM-GAME nodata value
		//TODO uref, zref must be provided
		force_station_dimension = true;
	} else if (schema=="ECMWF") {
		schema_dflt_type = NC_DOUBLE;
	} else if (schema=="WRF") {
		schema_dflt_type = NC_DOUBLE;
	} else if (schema=="AMUNDSEN") {
		schema_dflt_type = NC_FLOAT;
	} else
		throw InvalidArgumentException("Unsupported NetCDF schema "+schema, AT);
}

std::map< std::string, std::vector<ncpp::nc_dimension> > ncParameters::initSchemasDims()
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
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "surface_altitude") );
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
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "ZS") );
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
	tmp.push_back( ncpp::nc_dimension(MeteoGrids::DEM, "alt") );
	results["AMUNDSEN"] = tmp;
	
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

	//CF-1.6 schema -> to be checked and improved from CF-1.6 documentation
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "", "min", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "", "degree_north", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "", "degree_east", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "projection_x_coordinate", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "projection_y_coordinate", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "orog", "surface_altitude", "height above mean sea level", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "ta", "air_temperature", "near surface air temperature", "K", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RH, "hur", "relative_humidity", "", "1", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DW, "dw", "wind_from_direction", "", "degree", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW, "ws", "wind_speed", "", "m/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "ps", "surface_air_pressure", "near surface air pressure", "Pa", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P_SEA, "psl", "air_pressure_at_mean_sea_level", "", "Pa", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIR, "iswr_dir", "direct_downwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIFF, "iswr_diff", "diffuse_downwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "rsds", "surface_downwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSWR, "rsus", "surface_upwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSWR, "rsus", "surface_upwelling_shortwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "rlds", "surface_downwelling_longwave_flux_in_air", "", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::HS, "snd", "surface_snow_thickness", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSNO, "snow_density", "snow_density", "", "kg/m3", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SWE, "swe", "lwe_thickness_of_surface_snow_amount", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM, "pr", "precipitation_flux", "", "kg/m2/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM_S, "solid_precipitation_flux", "solid_precipitation_flux", "", "kg/m2/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSS, "ts", "surface_temperature", "", "K", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW_MAX, "ws_max", "wind_speed_of_gust", "", "m/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ALB, "surface_albedo", "surface_albedo", "", "1", IOUtils::nodata, NC_FLOAT) );
	//tmp.push_back( ncpp::var_attr(MeteoGrids::TSG, "tsg", "soil_surface_temperature", "", "K", IOUtils::nodata, NC_FLOAT) ); //HACK this is non-standard!
	results["CF-1.6"] = tmp;

	//CROCUS schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "s", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "LAT", "latitude", "latitude", "degrees_north", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "LON", "longitude", "longitude", "degrees_east", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::ZREF, "ZREF", "", "Reference_Height", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(ncpp::UREF, "UREF", "", "Reference_Height_for_Wind", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "ZS", "altitude", "altitude", "m", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RH, "HUMREL", "", "Relative Humidity", "%", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::QI, "Qair", "", "Near Surface Specific Humidity", "kg/kg", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW, "Wind", "", "Wind Speed", "m/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DW, "Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM_L, "Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM_S, "Snowf", "", "Snowfall Rate", "kg/m2/s", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "rsds", "", "Surface Incident total Shortwave radiation", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIR, "DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR_DIFF, "SCA_SWdown", "", "Surface Incident Diffuse Shortwave Radiation", "W/m2", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata, NC_FLOAT) );
	results["CROCUS"] = tmp;
	
	//AMUNDSEN schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "h", IOUtils::nodata, NC_INT) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "lat", "latitude", "latitude", "degrees_north", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "lon", "longitude", "longitude", "degrees_east", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "x", "projection_x_coordinate", "x coordinate of projection", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "y", "projection_y_coordinate", "x coordinate of projection", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "alt", "surface_altitude", "", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "tas", "", "air_temperature", "K", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RH, "hurs", "", "relative_humidity", "%", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::VW, "wss", "", "wind_speed", "m s-1", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM, "pr", "", "precipitation_flux", "kg m-2 s-1", IOUtils::nodata, NC_FLOAT) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "rsds", "", "surface_downwelling_shortwave_flux_in_air", "W m-2", IOUtils::nodata, NC_FLOAT) );
	results["AMUNDSEN"] = tmp;

	//ECMWF schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "time", "time", "time", "h", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "z", "geopotential_height", "geopotential_height", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "t2m", "", "2 metre temperature", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TD, "d2m", "", "2 metre dewpoint temperature", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "sp", "surface_air_pressure", "Surface pressure", "Pa", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P_SEA, "msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J m**-2", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "strd", "", "Surface thermal radiation downwards", "J m**-2", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::PSUM, "tp", "", "Total precipitation", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::U, "u10", "", "10 metre U wind component", "m s**-1", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::V, "v10", "", "10 metre V wind component", "m s**-1", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::SWE, "sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSS, "skt", "", "Skin temperature", "K", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSG, "stl1", "surface_temperature", "Soil temperature level 1", "K", IOUtils::nodata, NC_DOUBLE) ); //this is from 0 to -7cm
	tmp.push_back( ncpp::var_attr(MeteoGrids::ALB, "al", "surface_albedo", "Albedo", "(0 - 1)", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ALB, "fal", "", "Forecast albedo", "(0 - 1)", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSNO, "rsn", "", "Snow density", "kg m**-3", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ROT, "ro", "", "Runoff", "m", IOUtils::nodata, NC_DOUBLE) );
	results["ECMWF"] = tmp;

	//WRF schema
	tmp.clear();
	tmp.push_back( ncpp::var_attr(ncpp::TIME, "Times", "Times", "Times", "h", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::LATITUDE, "XLAT", "latitude", "latitude", "degrees", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::LONGITUDE, "XLONG", "longitude", "longitude", "degrees", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::STATION, "station", "timeseries_id", "", "", IOUtils::nodata, NC_CHAR) );
	tmp.push_back( ncpp::var_attr(ncpp::EASTING, "easting", "easting", "", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(ncpp::NORTHING, "northing", "northing", "", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::DEM, "HGT", "Terrain Height", "Terrain Height", "m", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::P, "PSFC", "Surface pressure", "Surface pressure", "Pa", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TA, "T2", "2-meter temperature", "2-meter temperature", "K", 2., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::QI, "Q2", "2-meter specific humidity", "2-meter specific humidity", "kg kg-1", 2, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ISWR, "ACSWDNB", "Downward SW surface radiation", "Downward SW surface radiation", "W m**-2", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::RSWR, "ACSWUPB", "Upwelling Surface Shortwave Radiation", "Upwelling Surface Shortwave Radiation", "W m**-2", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ILWR, "ACLWDNB", "Downward LW surface radiation", "Downward LW surface radiation", "W m**-2", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::ROT, "SFROFF", "Surface runoff ", "Surface runoff ", "kg*m2*s-1", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::HS, "SNOWH", "Snow depth", "Snow depth", "Pa", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::TSS, "TSK", "Surface skin temperature", "Surface skin temperature", "K", IOUtils::nodata, NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::U, "U10", "10-meter wind speed", "10 metre U wind component", "m s**-1", 10., NC_DOUBLE) );
	tmp.push_back( ncpp::var_attr(MeteoGrids::V, "V10", "10-meter wind speed", "10 metre V wind component", "m s**-1", 10., NC_DOUBLE) );
	results["WRF"] = tmp;
	
	return results;
}

//The user can provide his own variables properties as NETCDF_VAR::{param} = {name}
std::vector<ncpp::var_attr> ncParameters::initUserSchemas(const Config& i_cfg)
{
	std::vector<ncpp::var_attr> results;
	const std::string user_type_str = i_cfg.get("NC_TYPE", "Input", IOUtils::nothrow);
	char user_type = -1;
	if (!user_type_str.empty()) {
		if (user_type_str=="DOUBLE") user_type = NC_DOUBLE;
		else if (user_type_str=="FLOAT") user_type = NC_FLOAT;
		else if (user_type_str=="INT") user_type = NC_INT;
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
             : user_schemas( initUserSchemas(cfg) ), user_dimensions( initUserDimensions(cfg) ), vars(), unknown_vars(), vecTime(), vecX(), vecY(), dimensions_map(), file_and_path(filename), current_schema(schema), coord_sys(), coord_param(), TZ(tz_in), dflt_zref(IOUtils::nodata), dflt_uref(IOUtils::nodata), dflt_slope(IOUtils::nodata), dflt_azi(IOUtils::nodata),
             schema_nodata(IOUtils::nodata), schema_dflt_type(NC_DOUBLE), debug(i_debug), isLatLon(false), force_station_dimension(false)
{
	IOUtils::getProjectionParameters(cfg, coord_sys, coord_param);
	cfg.getValue("ZREF", "Input", dflt_zref, IOUtils::nothrow);
	cfg.getValue("UREF", "Input", dflt_uref, IOUtils::nothrow);
	cfg.getValue("DEFAULT_SLOPE", "Input", dflt_slope, IOUtils::nothrow);
	cfg.getValue("DEFAULT_AZI", "Input", dflt_azi, IOUtils::nothrow);
	initSchemaCst(schema);
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
		vars[ schemas_vars[schema][ii].param ] = ncpp::nc_variable( schemas_vars[schema][ii], schema_nodata );
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
		if (it->second.attributes.name==varname) return read2DGrid(it->second, IOUtils::npos);
	}
	
	for (std::map<std::string, ncpp::nc_variable>::const_iterator it = unknown_vars.begin(); it!=unknown_vars.end(); ++it) {
		if (it->first==varname) return read2DGrid(it->second, IOUtils::npos);
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
		ncpp::read_data(ncid, var, time_pos, vecY.size(), vecX.size(), data);
	else
		ncpp::read_data(ncid, var, data);
	ncpp::fill2DGrid(grid, data, var.nodata, (vecX.front()<=vecX.back()), (vecY.front()<=vecY.back()) );
	delete[] data;
	ncpp::close_file(file_and_path, ncid);
	
	//handle data packing and units, if necessary
	if (var.scale!=1.) grid *= var.scale;
	if (var.offset!=0.) grid += var.offset;
	applyUnits(grid, var.attributes.units, time_pos, m2mm);
	if (reZero) {//reset very low values to zero
		for (size_t ii=0; ii<grid.size(); ii++)
			if (grid(ii)<1e-6 && grid(ii)!=mio::IOUtils::nodata) grid(ii)=0.;
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
		ncpp::nc_variable tmp_var(tmp_attr, schema_nodata);
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
		writeGridMetadataHeader(ncid, grid_in);
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
	if (var.varid == -1) ncpp::create_variable(ncid, var); //create the "main" variable if necessary
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
			ncpp::write_data(ncid, vars[ param ], time_pos, grid_in.getNy(), grid_in.getNx(), &data[0]);
		} else {
			ncpp::write_data(ncid, vars[ param ], data, false);
		}
	}
	
	ncpp::close_file(file_and_path, ncid);
}

//When writing multiple stations in one file, this assumes that they all have the same parameters, the same timestamps and the same coordinate system
void ncParameters::writeMeteo(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx)
{
	const bool station_dimension = (station_idx==IOUtils::npos) || force_station_dimension;
	const size_t ref_station_idx = (station_dimension)? 0 : station_idx;
	if (vecMeteo.empty()) return;
	if (vecMeteo[ref_station_idx].empty()) return;
	isLatLon = true; //for now, we force lat/lon coordinates for time series
	
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		ncpp::file_redef(file_and_path, ncid);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		writeMeteoMetadataHeader(ncid, vecMeteo, station_idx);
	}
	
	std::vector<size_t> nc_variables, dimensions;
	dimensions.push_back( ncpp::TIME );
	if (station_dimension) {
		dimensions.push_back( ncpp::STATSTRLEN ); //this MUST be before STATION
		dimensions.push_back( ncpp::STATION );
	}
	const Date ref_date( getRefDate(vecMeteo, station_idx) );
	const size_t nrStations = (station_dimension)? vecMeteo.size() : 1;
	for (size_t ii=0; ii<dimensions.size(); ii++) {
		const size_t param = dimensions[ii];
		const size_t length = (param==ncpp::TIME)? 0 : ((param==ncpp::STATSTRLEN)? DFLT_STAT_STR_LEN : nrStations) ;
		ncpp::createDimension(ncid, dimensions_map[ param ], length);
		
		if (param==ncpp::STATSTRLEN) continue; //no associated variable for STATSTRLEN
		if (!station_dimension && param==ncpp::STATION) continue;
		if (setAssociatedVariable(ncid, param, ref_date)) nc_variables.push_back( dimensions[ii] ); //associated variable will have to be filled
	}
	
	appendVariablesList(nc_variables, vecMeteo, station_idx);
	
	//associate dimids to vars and create the variables
	for (size_t ii=0; ii<nc_variables.size(); ii++) {
		const size_t param = nc_variables[ii];
		if (vars[ param ].varid == -1) { //skip existing nc_variables
			const bool varIsLocation = (param==MeteoGrids::DEM || param==MeteoGrids::SLOPE || param==MeteoGrids::AZI || param==ncpp::ZREF || param==ncpp::UREF);
			if (param<ncpp::firstdimension && !varIsLocation) vars[ param ].dimids.push_back( dimensions_map[ncpp::TIME].dimid );
			if (station_dimension) vars[ param ].dimids.push_back( dimensions_map[ncpp::STATION].dimid );
			if (param==ncpp::STATION) vars[ param ].dimids.push_back( dimensions_map[ncpp::STATSTRLEN].dimid );
			
			ncpp::create_variable(ncid, vars[ param ]);
		}
	}
	
	ncpp::end_definitions(file_and_path, ncid);
	
	//write data: fill associated nc_variables and normal nc_variables
	for (size_t ii=0; ii<nc_variables.size(); ii++) {
		const size_t param = nc_variables[ii];
		
		if (param==ncpp::STATION) { //this is not present if !station_dimension
			std::vector<std::string> txtdata( vecMeteo.size() );
			for (size_t jj=0; jj<vecMeteo.size(); jj++) txtdata[jj] = vecMeteo[jj].front().meta.stationID;
			ncpp::write_data(ncid, vars[param], txtdata, DFLT_STAT_STR_LEN);
		} else {
			const std::vector<double> data( fillBufferForVar(vecMeteo, station_idx, vars[ param ]) );
			if (data.empty()) continue;
			ncpp::write_data(ncid, vars[ param ], data, (param==ncpp::TIME));
		}
	}
	
	ncpp::close_file(file_and_path, ncid);
}

void ncParameters::writeGridMetadataHeader(const int& ncid, const Grid2DObject& grid_in) const
{
	ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", current_schema+",ACDD-1.3");
	if (current_schema=="CF-1.6") ncpp::add_attribute(ncid, NC_GLOBAL, "standard_name_vocabulary", "CF-1.6");
	ncpp::add_attribute(ncid, NC_GLOBAL, "cdm_data_type", "Grid");
	Date now; now.setFromSys();
	ncpp::add_attribute(ncid, NC_GLOBAL, "date_created", now.toString(Date::ISO_DATE));
	ncpp::add_attribute(ncid, NC_GLOBAL, "creator_name", IOUtils::getLogName());
	ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO-" + getLibVersion(true));
	ncpp::add_attribute(ncid, NC_GLOBAL, "history", ncpp::generateHistoryAttribute());
	ncpp::add_attribute(ncid, NC_GLOBAL, "keywords_vocabulary", "AGU Index Terms");
	ncpp::add_attribute(ncid, NC_GLOBAL, "keywords", "Cryosphere, Mass Balance, Energy Balance, Atmosphere, Land/atmosphere interactions, Climatology");
	ncpp::add_attribute(ncid, NC_GLOBAL, "title", "Gridded data for various parameters and timesteps");
	//The following are placeholders to help users know what has to be manually provided for ACDD compliance
	ncpp::add_attribute(ncid, NC_GLOBAL, "summary", "Please fill this field to be ACDD compliant");
	ncpp::add_attribute(ncid, NC_GLOBAL, "acknowledgement", "Please fill this field to be ACDD compliant");
	ncpp::add_attribute(ncid, NC_GLOBAL, "metadata_link", "Please fill this field (with DOI or URL) to be ACDD compliant");

	
	Coords urcorner(grid_in.llcorner);
	urcorner.moveByXY(static_cast<double>(grid_in.getNx())*grid_in.cellsize, static_cast<double>(grid_in.getNy())*grid_in.cellsize);
	
	std::string epsg_str = "4326";
	std::string geometry;
	if (isLatLon) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid_in.llcorner.getLon() << " " << grid_in.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << grid_in.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << urcorner.getLat() << ", ";
		ss << grid_in.llcorner.getLon() << " " << urcorner.getLat();
		geometry = ss.str();
	}else {
		std::ostringstream os;
		os << grid_in.llcorner.getEPSG();
		epsg_str = os.str();
		
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid_in.llcorner.getEasting() << " " << grid_in.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << grid_in.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << urcorner.getNorthing() << ", ";
		ss << grid_in.llcorner.getEasting() << " " << urcorner.getNorthing();
		geometry = ss.str();
	}
	ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_bounds_crs", "EPSG:"+epsg_str);
	ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_bounds", "Polygon (("+geometry+"))");
}

void ncParameters::writeMeteoMetadataHeader(const int& ncid, const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx) const
{
	ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", current_schema+",ACDD-1.3");
	if (current_schema=="CF-1.6") ncpp::add_attribute(ncid, NC_GLOBAL, "standard_name_vocabulary", "CF-1.6");
	ncpp::add_attribute(ncid, NC_GLOBAL, "cdm_data_type", "Station");
	Date now; now.setFromSys();
	ncpp::add_attribute(ncid, NC_GLOBAL, "date_created", now.toString(Date::ISO_DATE));
	ncpp::add_attribute(ncid, NC_GLOBAL, "creator_name", IOUtils::getLogName());
	ncpp::add_attribute(ncid, NC_GLOBAL, "source", "MeteoIO-" + getLibVersion(true));
	ncpp::add_attribute(ncid, NC_GLOBAL, "history", ncpp::generateHistoryAttribute());
	ncpp::add_attribute(ncid, NC_GLOBAL, "keywords_vocabulary", "AGU Index Terms");
	ncpp::add_attribute(ncid, NC_GLOBAL, "keywords", "Cryosphere, Mass Balance, Energy Balance, Atmosphere, Land/atmosphere interactions, Climatology, Time series analysis");
	ncpp::add_attribute(ncid, NC_GLOBAL, "institution", IOUtils::getDomainName());
	//The following are placeholders to help users know what has to be manually provided for ACDD compliance
	ncpp::add_attribute(ncid, NC_GLOBAL, "summary", "Please fill this field to be ACDD compliant");
	ncpp::add_attribute(ncid, NC_GLOBAL, "acknowledgement", "Please fill this field to be ACDD compliant");
	ncpp::add_attribute(ncid, NC_GLOBAL, "metadata_link", "Please fill this field (with DOI or URL) to be ACDD compliant");
	
	Date set_start, set_end;
	int sampling_period = -1;
	
	if (station_idx==IOUtils::npos) {
		ncpp::add_attribute(ncid, NC_GLOBAL, "title", "Meteorological data timeseries for multiple stations");
		
		double lat_min, lat_max, lon_min, lon_max;
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			const Date curr_start = vecMeteo[ii].front().date;
			const Date curr_end = vecMeteo[ii].back().date;
			const double curr_lat = vecMeteo[ii].front().meta.position.getLat();
			const double curr_lon = vecMeteo[ii].front().meta.position.getLon();
			
			if (set_start.isUndef()) {
				set_start = curr_start;
				set_end = curr_end;
				lat_min = lat_max = curr_lat;
				lon_min = lon_max = curr_lon;
			}
			if (set_start>curr_start) set_start = curr_start;
			if (set_end<curr_end) set_end = curr_end;
			if (lat_min>curr_lat) lat_min = curr_lat;
			if (lat_max<curr_lat) lat_max = curr_lat;
			if (lon_min>curr_lon) lon_min = curr_lon;
			if (lon_max<curr_lon) lon_max = curr_lon;
			
			if (vecMeteo[ii].size()==1) continue;
			const int curr_sampling = static_cast<int>( (curr_end.getJulian() - curr_start.getJulian()) / static_cast<double>(vecMeteo[ii].size() - 1) * 24.*3600. + .5);
			if (sampling_period<=0 || sampling_period>curr_sampling) sampling_period = curr_sampling;
		}
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_lat_min", lat_min);
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_lat_max", lat_max);
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_lon_min", lon_min);
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_lon_max", lon_max);
	} else {
		const std::string stationName = vecMeteo[station_idx].front().meta.stationName;
		const std::string name = (!stationName.empty())? stationName : vecMeteo[station_idx].front().meta.stationID;
		ncpp::add_attribute(ncid, NC_GLOBAL, "title", "Meteorological data timeseries for the "+name+" station");
		ncpp::add_attribute(ncid, NC_GLOBAL, "station_name", name);
		std::string epsg_str = "4326";
		std::string geometry;
		const Coords location = vecMeteo[station_idx].front().meta.position;
		if (isLatLon) {
			std::ostringstream ss;
			ss << std::fixed << std::setprecision(10) << location.getLon() << " " << location.getLat();
			geometry = ss.str();
		}else {
			std::ostringstream os;
			os << location.getEPSG();
			epsg_str = os.str();
			std::ostringstream ss;
			ss << std::fixed << std::setprecision(10) << location.getEasting() << " " << location.getNorthing();
		}
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_bounds_crs", "EPSG:"+epsg_str);
		ncpp::add_attribute(ncid, NC_GLOBAL, "geospatial_bounds", "Point ("+geometry+")");
		
		set_start = vecMeteo[station_idx].front().date;
		set_end = vecMeteo[station_idx].back().date;
		const size_t npts = vecMeteo[station_idx].size();
		if (npts>1) sampling_period = static_cast<int>( (set_end.getJulian() - set_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
	}
	ncpp::add_attribute(ncid, NC_GLOBAL, "time_coverage_start", set_start.toString(Date::ISO_TZ));
	ncpp::add_attribute(ncid, NC_GLOBAL, "time_coverage_end", set_end.toString(Date::ISO_TZ));
	
	if (sampling_period>0) {
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		ncpp::add_attribute(ncid, NC_GLOBAL, "time_coverage_resolution", os.str());
	}
}

Date ncParameters::getRefDate(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx)
{
	if (station_idx==IOUtils::npos) {
		Date refDate;
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (refDate.isUndef()) refDate = vecMeteo[ii].front().date;
			if (vecMeteo[ii].front().date<refDate) refDate = vecMeteo[ii].front().date;
		}
		return refDate;
	} else {
		if (vecMeteo[station_idx].empty()) return Date();
		return vecMeteo[station_idx].front().date;
	}
}

void ncParameters::pushVar(std::vector<size_t> &nc_variables, const size_t& param)
{
	if (std::find(nc_variables.begin(), nc_variables.end(), param) == nc_variables.end())
		nc_variables.push_back( param );
}

//add an out-of-schema parameter to the vars map (if necessary)
void ncParameters::addToVars(const size_t& param)
{
	if (vars.count(param)==0) { //ie unrecognized in loaded schema, adding it
		const std::string varname( ncpp::getParameterName(param) );
		const ncpp::var_attr tmp_attr(param, varname, IOUtils::nodata, schema_dflt_type);
		vars[param] = ncpp::nc_variable(tmp_attr, schema_nodata);
	}
}

void ncParameters::appendVariablesList(std::vector<size_t> &nc_variables, const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx)
{
	//add metadata variables
	pushVar(nc_variables, MeteoGrids::DEM);
	if (isLatLon) {
		pushVar(nc_variables, ncpp::LATITUDE );
		pushVar(nc_variables, ncpp::LONGITUDE );
	} else {
		pushVar(nc_variables, ncpp::EASTING );
		pushVar(nc_variables, ncpp::NORTHING );
	}
	addToVars( MeteoGrids::SLOPE); pushVar(nc_variables, MeteoGrids::SLOPE );
	addToVars( MeteoGrids::AZI ); pushVar(nc_variables, MeteoGrids::AZI );
	if (dflt_zref!=IOUtils::nodata || dflt_uref!=IOUtils::nodata) { //This is required by Crocus and we don't have any better solution for now...
		addToVars( ncpp::ZREF ); pushVar(nc_variables, ncpp::ZREF );
		addToVars( ncpp::UREF ); pushVar(nc_variables, ncpp::UREF );
	}
	
	//add all vars found in vecMeteo
	const size_t st_start = (station_idx==IOUtils::npos)? 0 : station_idx;
	const size_t st_end = (station_idx==IOUtils::npos)? vecMeteo.size() : station_idx+1;
	for (size_t st=st_start; st<st_end; st++) {
		const std::set<std::string> parameters( MeteoData::listAvailableParameters( vecMeteo[st] ));
		
		for (std::set<std::string>::const_iterator it=parameters.begin(); it!=parameters.end(); ++it) {
			const size_t param = MeteoGrids::getParameterIndex( *it );
			if (param>=MeteoGrids::nrOfParameters) continue; //in order to be supported, a parameter must be declared in MeteoGrids
			addToVars(param); pushVar(nc_variables, param);
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
			
			const std::string units( vars[ param ].attributes.units );
			if (units=="h") {
				vars[param].attributes.units = "hours since " + date_str;
				vars[param].scale = 1./24.;
			} else if (units=="min") {
				vars[param].attributes.units = "minutes since " + date_str;
				vars[param].scale = 1./(24.*60);
			} else if (units=="s") {
				vars[param].attributes.units = "seconds since " + date_str;
				vars[param].scale = 1./(24.*3600.);
			} else 
				throw InvalidArgumentException("Unsupported time unit specified in schema: '"+units+"'", AT);
			vars[param].offset = ref_date_simplified.getJulian();
		}
		
		ncpp::create_variable(ncid, vars[ param ]);
		return true;
	}
	return false;
}

const std::vector<double> ncParameters::fillBufferForVar(const std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& station_idx, ncpp::nc_variable& var) const
{
	const size_t param = var.attributes.param;
	const size_t ref_station_idx = (station_idx==IOUtils::npos)? 0 : station_idx;
	const size_t st_start = (station_idx==IOUtils::npos)? 0 : station_idx;
	const size_t st_end = (station_idx==IOUtils::npos)? vecMeteo.size() : station_idx+1;
	const size_t nrStations = (station_idx==IOUtils::npos)? vecMeteo.size() : 1;
	
	const bool varIsLocation = (param==MeteoGrids::DEM || param==MeteoGrids::SLOPE || param==MeteoGrids::AZI || param==ncpp::ZREF || param==ncpp::UREF);
	if (param>=ncpp::firstdimension || varIsLocation) { //associated nc_variables
		if (param==ncpp::TIME) {
			const size_t nrTimeSteps = vecMeteo[ref_station_idx].size();
			std::vector<double> data(nrTimeSteps, var.nodata);
			if (var.attributes.type==NC_INT) { //in this case, we pre-round the data so when libnetcdf will cast, it will fall on what we want
				double prev = IOUtils::nodata;
				for (size_t ll=0; ll<nrTimeSteps; ll++) {
					data[ll] = static_cast<double>( Optim::round( (vecMeteo[ref_station_idx][ll].date.getJulian() - var.offset) / var.scale) );
					if (prev!=IOUtils::nodata && data[ll]==prev) throw InvalidArgumentException("When writing time as INT, some timesteps are rounded to identical values. Please change your sampling rate!", AT);
					prev = data[ll];
				}
			} else {
				for (size_t ll=0; ll<nrTimeSteps; ll++)
					data[ll] = (vecMeteo[ref_station_idx][ll].date.getJulian() - var.offset) / var.scale;
			}
			return data;
		} else {
			std::vector<double> data(nrStations, var.nodata);
			for (size_t jj=st_start; jj<st_end; jj++) {
				double value = IOUtils::nodata;
				if (param==MeteoGrids::DEM) {
					value = vecMeteo[jj].front().meta.position.getAltitude();
				} else if (param==MeteoGrids::SLOPE) {
					value = vecMeteo[jj].front().meta.getSlopeAngle();
					if (value==IOUtils::nodata) value = dflt_slope;
				} else if (param==MeteoGrids::AZI) {
					value = vecMeteo[jj].front().meta.getAzimuth();
					if (value==IOUtils::nodata) value = dflt_azi;
				} else if (param==ncpp::EASTING) {
					value = vecMeteo[jj].front().meta.position.getEasting();
				} else if (param==ncpp::NORTHING) {
					value = vecMeteo[jj].front().meta.position.getNorthing();
				} else if (param==ncpp::LATITUDE) {
					value = vecMeteo[jj].front().meta.position.getLat();
				} else if (param==ncpp::LONGITUDE) {
					value = vecMeteo[jj].front().meta.position.getLon();
				} else if (param==ncpp::ZREF) { //this two are required by Crocus and we don't have anything better for now...
					value = dflt_zref;
				} else if (param==ncpp::UREF) {
					value = dflt_uref;
				} else
					throw UnknownValueException("Unknown dimension found when trying to write out netcdf file", AT);
				
				if (value!=IOUtils::nodata) data[jj-st_start] = value;
			}
			return data;
		}
	} else { //normal nc_variables
		const MeteoData md( vecMeteo[ref_station_idx].front() );
		const size_t meteodata_param = md.getParameterIndex( MeteoGrids::getParameterName( param ) ); //retrieve the equivalent parameter in vecMeteo
		if (meteodata_param==IOUtils::npos) return std::vector<double>(); //this should not have happened...
		
		const size_t nrTimeSteps = vecMeteo[ref_station_idx].size();
		std::vector<bool> stationHasParameter(nrStations, false);
		for (size_t jj=st_start; jj<st_end; jj++) {
			if(!vecMeteo[jj].empty()) stationHasParameter[jj-st_start] = (vecMeteo[jj].front().getNrOfParameters()>meteodata_param);
		}
		
		std::vector<double> data(nrTimeSteps*nrStations, var.nodata);
		for (size_t ll=0; ll<nrTimeSteps; ll++) {
			for (size_t jj=st_start; jj<st_end; jj++) {
				if (stationHasParameter[jj-st_start] && vecMeteo[jj][ll]( meteodata_param )!=IOUtils::nodata) data[ll*nrStations + (jj-st_start)] = vecMeteo[jj][ll]( meteodata_param );
			}
		}

		//perform some units corrections, if necessary
		if (var.attributes.units=="%") for (size_t ii=0; ii<data.size(); ii++) data[ii] *= 100.;

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
		std::vector<double> data( nrPts, var.nodata );
		
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
		std::vector<double> data( grid.size(), var.nodata );
		const size_t nrows = grid.getNy(), ncols = grid.getNx();
		for (size_t kk=0; kk<nrows; ++kk) {
			for (size_t ll=0; ll<ncols; ++ll) {
				if (grid.grid2D(ll,kk)!=IOUtils::nodata) data[kk*ncols + ll] = grid.grid2D(ll,kk);
			}
		}
		return data;
	}
}

//bring back known units to MKSA 
void ncParameters::applyUnits(Grid2DObject& grid, const std::string& units, const size_t& time_pos, const bool& m2mm) const
{
	if (units.empty()) return;
	
	if (units=="m**2 s**-2") grid /= Cst::gravity;
	else if (units=="%") grid /= 100.;
	else if (units=="J m**-2") {
		if (vecTime.size()>1 && time_pos!=IOUtils::npos) {
			const Date integration_period = (time_pos>0)? (vecTime[time_pos] - vecTime[time_pos-1]) : (vecTime[time_pos+1] - vecTime[time_pos]);
			grid /= (integration_period.getJulian()*24.*3600.); //converting back to W/m2
		}
	}
	else if (m2mm && units=="m") grid *= 1000.;
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

		ncpp::nc_variable tmp_var( tmp_attr, schema_nodata );
		const bool readTimeTransform = (tmp_var.attributes.param==ncpp::TIME && tmp_var.attributes.type!=NC_CHAR);
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
	const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( ncpp::TIME );
	if (it==vars.end() || it->second.varid==-1) throw InvalidArgumentException("Could not find parameter \""+ncpp::getParameterName(ncpp::TIME)+"\" in file \""+file_and_path+"\"", AT);
	const bool timestamps_as_str = (it->second.attributes.type==NC_CHAR);
	
	if (!timestamps_as_str) {
		const std::vector<double> tmp_results( read_1Dvariable(ncid, ncpp::TIME) );
		std::vector<Date> results(tmp_results.size());
		for (size_t ii=0; ii<tmp_results.size(); ii++)
			results[ii].setDate(tmp_results[ii]*it->second.scale + it->second.offset, TZ);
		return results;
	} else {
		std::vector<std::string> tmp_results( read_1Dstringvariable(ncid, ncpp::TIME) );
		const size_t length = tmp_results.size();
		std::vector<Date> results( length );
		for(size_t ii=0; ii<length; ii++) {
			tmp_results[ii][ 10 ] = 'T';
			IOUtils::convertString(results[ii], tmp_results[ii], TZ);
		}
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
	ncpp::read_data(ncid, it->second, data);
	std::copy(data, data+length, results.begin());
	delete[] data;
	return results;
}

std::vector<std::string> ncParameters::read_1Dstringvariable(const int& ncid, const size_t& param) const
{
	const std::map<size_t, ncpp::nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end() || it->second.varid==-1) throw InvalidArgumentException("Could not find parameter \""+ncpp::getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	if (it->second.attributes.type!=NC_CHAR) throw InvalidArgumentException("Wrong data type for parameter \""+ncpp::getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	
	const size_t length = read_1DvariableLength(it->second); //this also checks that it depends on 2 dimensions for NC_CHAR
	const size_t strMaxLen = readDimension(it->second.dimids[1]); //string lenght is in second position
	char *data = (char*)calloc(length, sizeof(char)*strMaxLen);
	const int status = nc_get_var_text(ncid, it->second.varid, data);
	if (status != NC_NOERR) throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);
	
	std::vector<std::string> results(length);
	for(size_t ii=0; ii<length; ii++) {
		for(size_t jj=0; jj<strMaxLen; jj++) 
			results[ii].push_back( data[ii*strMaxLen+jj] );
	}
	free( data );
	return results;
}

size_t ncParameters::read_1DvariableLength(const ncpp::nc_variable& var) const
{
	const size_t ndims = var.dimids.size();
	if ((var.attributes.type==NC_CHAR && ndims!=2) || (var.attributes.type!=NC_CHAR && ndims!=1)) 
		throw InvalidArgumentException("Parameter \""+ncpp::getParameterName(var.attributes.param)+"\" in file \""+file_and_path+"\" is not a 1D variable", AT);
	
	return readDimension( var.dimids[0] ); //in the case of vector of strings, the first dimension is the vector size
}

size_t ncParameters::readDimension(const int& dimid) const
{
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
	
	return ncpp::var_attr(var, schema_dflt_type);
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
	
	return ncpp::var_attr(schema_dflt_type);
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
