/***********************************************************************************/
/*  Copyright 2009-2012 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   */
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
#include <meteoio/IOHandler.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/dataClasses/MeteoData.h> //needed for the merge strategies

#include <algorithm>
#include <fstream>

#cmakedefine PLUGIN_ALPUG
#cmakedefine PLUGIN_ARCIO
#cmakedefine PLUGIN_ARGOSIO
#cmakedefine PLUGIN_A3DIO
#cmakedefine PLUGIN_ARPSIO
#cmakedefine PLUGIN_CSVIO
#cmakedefine PLUGIN_DBO
#cmakedefine PLUGIN_GRASSIO
#cmakedefine PLUGIN_GEOTOPIO
#cmakedefine PLUGIN_SMETIO
#cmakedefine PLUGIN_SNIO
#cmakedefine PLUGIN_PGMIO
#cmakedefine PLUGIN_PMODIO
#cmakedefine PLUGIN_IMISIO
#cmakedefine PLUGIN_OSHDIO
#cmakedefine PLUGIN_GRIBIO
#cmakedefine PLUGIN_GOESIO
#cmakedefine PLUGIN_PNGIO
#cmakedefine PLUGIN_COSMOXMLIO
#cmakedefine PLUGIN_NETCDFIO
#cmakedefine PLUGIN_PSQLIO
#cmakedefine PLUGIN_SASEIO
#cmakedefine PLUGIN_ZRXPIO

#include <meteoio/plugins/ALPUG.h>
#include <meteoio/plugins/ARCIO.h>
#include <meteoio/plugins/Argos.h>
#include <meteoio/plugins/A3DIO.h>
#include <meteoio/plugins/ARPSIO.h>
#include <meteoio/plugins/CsvIO.h>
#include <meteoio/plugins/Goes.h>
#include <meteoio/plugins/GrassIO.h>
#include <meteoio/plugins/GeotopIO.h>
#include <meteoio/plugins/PGMIO.h>
#include <meteoio/plugins/SMETIO.h>
#include <meteoio/plugins/SNIO.h>

#ifdef PLUGIN_COSMOXMLIO
#include <meteoio/plugins/CosmoXMLIO.h>
#endif

#ifdef PLUGIN_DBO
#include <meteoio/plugins/DBO.h>
#endif

#ifdef PLUGIN_IMISIO
#include <meteoio/plugins/ImisIO.h>
#endif

#ifdef PLUGIN_OSHDIO
#include <meteoio/plugins/OshdIO.h>
#endif

#ifdef PLUGIN_GRIBIO
#include <meteoio/plugins/GRIBIO.h>
#endif

#ifdef PLUGIN_NETCDFIO
#include <meteoio/plugins/NetCDFIO.h>
#endif

#ifdef PLUGIN_PMODIO
#include <meteoio/plugins/PmodIO.h>
#endif

#ifdef PLUGIN_PNGIO
#include <meteoio/plugins/PNGIO.h>
#endif

#ifdef PLUGIN_PSQLIO
#include <meteoio/plugins/PSQLIO.h>
#endif

#ifdef PLUGIN_SASEIO
#include <meteoio/plugins/SASEIO.h>
#endif

#ifdef PLUGIN_ZRXPIO
#include <meteoio/plugins/ZRXPIO.h>
#endif

using namespace std;

namespace mio {
 /**
 * @page data_sources Data input/output overview
 * The data access is handled by a system of plugins. They all offer the same interface, meaning that a plugin can transparently be replaced by another one. Since they
 * might rely on third party libraries for accessing the data, they have been created as plugins, that is they are only compiled if requested when configuring the
 * compilation with cmake. A plugin can therefore fail to run if it has not been compiled.
 *
 * Please have a look at the support for \subpage coords "coordinate systems".
 *
 * @section available_categories Data sources categories
 * Several data sources categories have been defined that can be provided by a different plugin. Each data source category is defined by a specific key in the configuration file (usually, io.ini):
 * - METEO, for meteorological time series
 * - DEM, for Digital Elevation Maps
 * - LANDUSE, for land cover information
 * - GRID2D, for generic 2D grids (they can contain meteo fields and be recognized as such or arbitrary gridded data)
 * - POI, for a list of Points Of Interest that can be used for providing extra information at some specific location (extracting time series at a few selected points, etc)
 *
 * A plugin is "connected" to a given data source category simply by giving its keyword as value for the data source key:
 * @code
 * METEO = SMET
 * DEM = ARC
 * @endcode
 * Each plugin might have its own specific options, meaning that it might require its own keywords. Please check in each plugin documentation the supported options and keys (see links below).
 * Moreover, a given plugin might only support a given category for read or write (for example, PNG: there is no easy and safe way to interpret a given color as a given numeric value without knowing its color scale, so reading a png has been disabled).
 * Finally, the plugins usually don't implement all these categories (for example, ArcGIS file format only describes 2D grids, so the ARC plugin will only deal with 2D grids), so please check what a given plugin implements before connecting it to a specific data source category.
 *
 * @subsection available_plugins Available plugins
 * So far the following plugins have been implemented (by keyword for the io.ini key/value config file). Please read the documentation for each plugin in order to know the plugin-specific keywords:
 * <center><table border="1">
 * <tr><th>Plugin keyword</th><th colspan="2">Provides</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><th></th><th>Input</th><th>Output</th>		<th></th><th></th></tr>
 * <tr><td>\subpage a3d "A3D"</td><td>meteo, poi</td><td>meteo</td>		<td>original Alpine3D meteo files</td><td></td></tr>
 * <tr><td>\subpage alpug "ALPUG"</td><td>meteo</td><td></td>		<td>data files generated by the %ALPUG meteo stations</td><td></td></tr>
 * <tr><td>\subpage arc "ARC"</td><td>dem, landuse, grid2d</td><td>grid2d</td>		<td>ESRI/ARC ascii grid files</td><td></td></tr>
 * <tr><td>\subpage argosio "ARGOS"</td><td>meteo</td><td></td>		<td>ARGOS satellite raw files - <b>plugin not functional yet!</b></td><td></td></tr>
 * <tr><td>\subpage arps "ARPS"</td><td>dem, grid2d, grid3d</td><td></td>		<td>ARPS ascii formatted grids</td><td></td></tr>
 * <tr><td>\subpage cosmoxml "COSMOXML"</td><td>meteo</td><td></td>		<td>MeteoSwiss COSMO's postprocessing XML format</td><td><A HREF="http://xmlsoft.org/">libxml2</A></td></tr>
 * <tr><td>\subpage csvio "CSV"</td><td>meteo</td><td></td>		<td>flexible reading of CSV files</td><td></td></tr>
 * <tr><td>\subpage dbo "DBO"</td><td>meteo</td><td></td>		<td>connects to SLF's DBO web service interface</td><td><A HREF="http://curl.haxx.se/libcurl/">libcurl</A></td></tr>
 * <tr><td>\subpage geotop "GEOTOP"</td><td>meteo</td><td>meteo</td>		<td>GeoTop meteo files</td><td></td></tr>
 * <tr><td>\subpage goesio "GOES"</td><td>meteo</td><td></td>		<td>Meteo files transmitted by the GOES satellites</td><td></td></tr>
 * <tr><td>\subpage grass "GRASS"</td><td>dem, landuse, grid2d</td><td>grid2d</td>		<td>Grass grid files</td><td></td></tr>
 * <tr><td>\subpage gribio "GRIB"</td><td>meteo, dem, grid2d</td><td></td>		<td>GRIB meteo grid files</td><td><A HREF="http://www.ecmwf.int/products/data/software/grib_api.html">grib-api</A></td></tr>
 * <tr><td>\subpage imis "IMIS"</td><td>meteo</td><td></td>		<td>connects to the IMIS database</td><td><A HREF="http://docs.oracle.com/cd/B12037_01/appdev.101/b10778/introduction.htm">Oracle's OCCI library</A></td></tr>
 * <tr><td>\subpage netcdf "NETCDF"</td><td>meteo, dem, grid2d</td><td>meteo, grid2d</td>		<td>NetCDF grids and timeseries</td><td><A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF-C library</A></td></tr>
 * <tr><td>\subpage oshd "OSHD"</td><td>meteo</td><td></td>		<td>OSHD generated binary Matlab files</td><td><A HREF="https://sourceforge.net/projects/matio">libmatio</A></td></tr>
 * <tr><td>\subpage pgmio "PGM"</td><td>dem, grid2d</td><td>grid2d</td>		<td>PGM grid files</td><td></td></tr>
 * <tr><td>\subpage pmodio "PMOD"</td><td>meteo</td><td></td>		<td>Raw data files from Pmod/Wrc (experimental!)</td><td></td></tr>
 * <tr><td>\subpage pngio "PNG"</td><td></td><td>grid2d</td>		<td>PNG grid files</td><td><A HREF="http://www.libpng.org/pub/png/libpng.html">libpng</A></td></tr>
 * <tr><td>\subpage psqlio "PSQL"</td><td>meteo</td><td>meteo</td>		<td>connects to PostgreSQL database</td><td><A HREF="http://www.postgresql.org/">PostgreSQL</A>'s libpq</td></tr>
 * <tr><td>\subpage sase "SASE"</td><td>meteo</td><td></td>		<td>connects to the SASE database</td><td><A HREF="https://dev.mysql.com/doc/refman/5.0/en/c-api.html">MySQL's C API</A></td></tr>
 * <tr><td>\subpage smetio "SMET"</td><td>meteo, poi</td><td>meteo</td>		<td>SMET data files</td><td></td></tr>
 * <tr><td>\subpage snowpack "SNOWPACK"</td><td>meteo</td><td>meteo</td>		<td>original SNOWPACK meteo files</td><td></td></tr>
 * <tr><td>\subpage zrxpio "ZRXP"</td><td></td><td>meteo</td>		<td>WISKI database input files</td><td></td></tr>
 * </table></center>
 *
 * @note In order to optimize the data retrieval, the raw data is buffered. This means that up to \b BUFFER_SIZE days of data will be read at once by the plugin
 * so subsequent reads will not have to get back to the data source (this key is in the [General] section). It is usually a good idea to configure \b BUFFER_SIZE
 * to the intended duration of the simulation (in days).
 *
 * @subsection Multiple_input_plugins Multiple data sources
 * It is possible to use multiple plugins to read \b meteorological \b timeseries from multiple sources and combine them into one stream of data. This is
 * achieved by declaring as many \em DATASOURCExxx sections as necessary (where xxx represent any number, to make sure that not two
 * sections have the same name) and declaring \em METEO plugins in each of them. Please make sure that all required keys are defined within each
 * new such section (such as TIME_ZONE) because the plugins created this way won't have access to the original \em INPUT section. The other
 * plugins (such as for reading grids, dem, etc) as well as the raw data editing will only be read from the standard \em INPUT section and
 * performed as usual after the data has been provided by the plugins.
 *
 * If reading the same station from multiple sources (for example providing different time coverage), it might be useful to use
 * the EditingAutoMerge feature to merge all streams belonging to a station into one single stream. In this case, the data coming from [Input]
 * has priority over the various [DataSourcexxx] data sources.
 *
 */

IOInterface* IOHandler::getPlugin(std::string plugin_name, const Config& i_cfg) const
{
	IOUtils::toUpper( plugin_name );
#ifdef PLUGIN_ALPUG
	if (plugin_name == "ALPUG") return new ALPUG(i_cfg);
#endif
#ifdef PLUGIN_ARCIO
	if (plugin_name == "ARC") return new ARCIO(i_cfg);
#endif
#ifdef PLUGIN_ARGOSIO
	if (plugin_name == "ARGOS") return new ArgosIO(i_cfg);
#endif
#ifdef PLUGIN_A3DIO
	if (plugin_name == "A3D") return new A3DIO(i_cfg);
#endif
#ifdef PLUGIN_ARPSIO
	if (plugin_name == "ARPS") return new ARPSIO(i_cfg);
#endif
#ifdef PLUGIN_CSVIO
	if (plugin_name == "CSV") return new CsvIO(i_cfg);
#endif
#ifdef PLUGIN_GOESIO
	if (plugin_name == "GOES") return new GoesIO(i_cfg);
#endif
#ifdef PLUGIN_GRASSIO
	if (plugin_name == "GRASS") return new GrassIO(i_cfg);
#endif
#ifdef PLUGIN_GEOTOPIO
	if (plugin_name == "GEOTOP") return new GeotopIO(i_cfg);
#endif
#ifdef PLUGIN_SMETIO
	if (plugin_name == "SMET") return new SMETIO(i_cfg);
#endif
#ifdef PLUGIN_SNIO
	if (plugin_name == "SNOWPACK") return new SNIO(i_cfg);
#endif
#ifdef PLUGIN_PGMIO
	if (plugin_name == "PGM") return new PGMIO(i_cfg);
#endif
#ifdef PLUGIN_IMISIO
	if (plugin_name == "IMIS") return new ImisIO(i_cfg);
#endif
#ifdef PLUGIN_OSHDIO
	if (plugin_name == "OSHD") return new OshdIO(i_cfg);
#endif
#ifdef PLUGIN_GRIBIO
	if (plugin_name == "GRIB") return new GRIBIO(i_cfg);
#endif
#ifdef PLUGIN_PMODIO
	if (plugin_name == "PMOD") return new PmodIO(cfg);
#endif
#ifdef PLUGIN_PNGIO
	if (plugin_name == "PNG") return new PNGIO(i_cfg);
#endif
#ifdef PLUGIN_COSMOXMLIO
	if (plugin_name == "COSMOXML") return new CosmoXMLIO(i_cfg);
#endif
#ifdef PLUGIN_DBO
	if (plugin_name == "DBO") return new DBO(i_cfg);
#endif
#ifdef PLUGIN_NETCDFIO
	if (plugin_name == "NETCDF") return new NetCDFIO(i_cfg);
#endif
#ifdef PLUGIN_PSQLIO
	if (plugin_name == "PSQL") return new PSQLIO(i_cfg);
#endif
#ifdef PLUGIN_SASEIO
	if (plugin_name == "SASE") return new SASEIO(i_cfg);
#endif
#ifdef PLUGIN_ZRXPIO
	if (plugin_name == "ZRXP") return new ZRXPIO(i_cfg);
#endif

	return NULL; //no plugin found
}

/**
 * @brief Return a pointer to a plugin object matching the provided key.
 * @details If the plugin has not yet been constructed, it will be constructed (acting as an object factory) and then cached. It is cached per
 * section / plugin type / plugin combination, for example "INPUT::METEO::SMET" so calling the same kind of plugin
 * for the inputs and outputs will lead to different instances being cached.
 * @param[in] cfgkey Configuration key giving the plugin name (example: METEO);
 * @param[in] cfgsection Section where to find this configuration key (example: INPUT);
 * @param[in] sec_rename New section name if the section should be renamed before being passed to the plugin's constructor
 * (default: empty string, so no renaming)
 * @return Pointer to the constructed plugin
 *
 */
IOInterface* IOHandler::getPlugin(const std::string& cfgkey, const std::string& cfgsection, const std::string& sec_rename)
{
	const std::string op_src = cfg.get(cfgkey, cfgsection);
	const std::string plugin_key( cfgsection+"::"+cfgkey+"::"+op_src ); //otherwise, reading meteo+grids with the same plugin would rely on the same object

	if (mapPlugins.find(plugin_key) == mapPlugins.end()) { //the plugin has not already been constructed
		IOInterface *ioPtr = NULL;

		if (sec_rename.empty() || sec_rename==cfgsection) {
			ioPtr = getPlugin(op_src, cfg);
		} else {
			Config cfg2( cfg );
			cfg2.moveSection(cfgsection, sec_rename, true);
			ioPtr = getPlugin(op_src, cfg2);
		}

		if (ioPtr==NULL)
			throw IOException("Cannot find plugin " + op_src + " as requested in file " + cfg.getSourceName() + ". Has it been activated through ccmake? Is it declared in IOHandler::getPlugin?", AT);
		else
			mapPlugins[plugin_key] = ioPtr;
	}

	return mapPlugins[plugin_key];
}

//Copy constructor
IOHandler::IOHandler(const IOHandler& aio)
           : IOInterface(), cfg(aio.cfg), preProcessor(aio.cfg), mapPlugins(aio.mapPlugins)
{}

IOHandler::IOHandler(const Config& cfgreader)
           : IOInterface(), cfg(cfgreader), preProcessor(cfgreader), mapPlugins()
{}

IOHandler::~IOHandler() throw()
{
	// Get rid of the objects
	std::map<std::string, IOInterface*>::iterator mapit( mapPlugins.begin() );
	for (; mapit!=mapPlugins.end(); ++mapit) {
		delete mapit->second;
	}
}

IOHandler& IOHandler::operator=(const IOHandler& source) {
	if (this != &source) {
		preProcessor = source.preProcessor;
		mapPlugins = source.mapPlugins;
	}
	return *this;
}

//return the list of sections that declare a certain plugin key and matching a given pattern
//The [INPUT] section is currently always returned even if it does not match the pattern
std::vector<std::string> IOHandler::getListOfSources(const std::string& plugin_key, const std::string& sec_pattern) const
{
	const std::set<std::string> sections( cfg.getSections() );
	std::vector<std::string> results;

	//the [Input] section should always be returned if available and should come first (for priority in case of a later merge)
	if (cfg.keyExists(plugin_key, "INPUT")) results.push_back( "INPUT" );

	for (std::set<std::string>::const_iterator it = sections.begin(); it!=sections.end(); ++it) {
		if (*it=="INPUT") continue; //the {input] section has already been processed
		const size_t found_pos = it->find(sec_pattern, 0);
		if (found_pos==0 && cfg.keyExists(plugin_key, *it)) results.push_back( *it );
	}

	return results;
}

bool IOHandler::list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> > &list)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	return plugin->list2DGrids(start, end, list);
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const std::string& i_filename)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, i_filename);
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, parameter, date);
}

void IOHandler::read3DGrid(Grid3DObject& grid_out, const std::string& i_filename)
{
	IOInterface *plugin = getPlugin("GRID3D", "Input");
	plugin->read3DGrid(grid_out, i_filename);
}

void IOHandler::read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID3D", "Input");
	plugin->read3DGrid(grid_out, parameter, date);
}

void IOHandler::readDEM(DEMObject& dem_out)
{
	IOInterface *plugin = getPlugin("DEM", "Input");
	plugin->readDEM(dem_out);
	dem_out.update();
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	IOInterface *plugin = getPlugin("LANDUSE", "Input");
	plugin->readLanduse(landuse_out);
}

void IOHandler::readGlacier(Grid2DObject& glacier_out)
{
	IOInterface *plugin = getPlugin("GLACIER", "Input");
	plugin->readGlacier(glacier_out);
}

void IOHandler::readStationData(const Date& date, STATIONS_SET& vecStation)
{
	const std::vector<std::string> sources( getListOfSources("METEO", "DATASOURCE") ); //[INPUT] is included anyway
	if (sources.empty()) throw UnknownValueException("No plugin defined for METEO", AT);;

	for (size_t ii=0; ii<sources.size(); ii++) {
		IOInterface *plugin = getPlugin("METEO", sources[ii], "INPUT");

		if (ii==0) {
			plugin->readStationData(date, vecStation);
		} else  {
			STATIONS_SET vectmp;
			plugin->readStationData(date, vectmp);
			for (size_t jj=0; jj<vectmp.size(); jj++) vecStation.push_back( vectmp[jj] );
		}
	}
	
	preProcessor.editTimeSeries( vecStation );
}

void IOHandler::readMeteoData(const Date& dateStart, const Date& dateEnd,
                              std::vector<METEO_SET>& vecMeteo)
{
	const std::vector<std::string> sources( getListOfSources("METEO", "DATASOURCE") ); //[INPUT] is included anyway
	if (sources.empty()) throw UnknownValueException("No plugin defined for METEO", AT);

	//some time filters change the requested dates (for example, time loop)
	Date fakeStart( dateStart ),fakeEnd( dateEnd );
	preProcessor.timeproc.process(fakeStart, fakeEnd);

	for (size_t ii=0; ii<sources.size(); ii++) {
		IOInterface *plugin = getPlugin("METEO", sources[ii], "INPUT");

		if (ii==0) {
			plugin->readMeteoData(fakeStart, fakeEnd, vecMeteo);
		} else  {
			std::vector<METEO_SET> vectmp;
			plugin->readMeteoData(fakeStart, fakeEnd, vectmp);
			for (size_t jj=0; jj<vectmp.size(); jj++) vecMeteo.push_back( vectmp[jj] );
		}
	}

	preProcessor.editTimeSeries( vecMeteo );
}

void IOHandler::writeMeteoData(const std::vector<METEO_SET>& vecMeteo,
                               const std::string& name)
{
	IOInterface *plugin = getPlugin("METEO", "Output");
	plugin->writeMeteoData(vecMeteo, name);
}

void IOHandler::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	IOInterface *plugin = getPlugin("DA", "Input");
	plugin->readAssimilationData(date_in, da_out);
}

void IOHandler::readPOI(std::vector<Coords>& pts) {
	IOInterface *plugin = getPlugin("POI", "Input");
	plugin->readPOI(pts);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, name);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, parameter, date);
}

void IOHandler::write3DGrid(const Grid3DObject& grid_out, const std::string& options)
{
	IOInterface *plugin = getPlugin("GRID3D", "Output");
	plugin->write3DGrid(grid_out, options);
}

void IOHandler::write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID3D", "Output");
	plugin->write3DGrid(grid_out, parameter, date);
}

const std::string IOHandler::toString() const
{
	std::ostringstream os;
	os << "<IOHandler>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";

	os << "<mapPlugins>\n";
	std::map<std::string, IOInterface*>::const_iterator it1;
	for (it1=mapPlugins.begin(); it1 != mapPlugins.end(); ++it1) {
		os << setw(10) << it1->first << " = " << hex <<  it1->second << dec << "\n";
	}
	os << "</mapPlugins>\n";
	
	os << "<preProcessor>\n";
	os << preProcessor.toString();
	os << "</preProcessor>\n";

	os << "</IOHandler>\n";
	return os.str();
}

} //end namespace
