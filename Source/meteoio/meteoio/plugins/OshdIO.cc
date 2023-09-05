// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/OshdIO.h>
#include <meteoio/plugins/libMatioWrapper.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/FileUtils.h>

#include <matio.h>
#include <algorithm>
#include <regex>

using namespace std;

namespace mio {
/**
 * @page oshd OshdIO
 * This plugin reads the meteorological forecast data from COSMO either as gridded data or downscaled for each of the Swiss meteorological
 * networks IMIS/ANETZ stations as preprocessed by the
 * <A HREF="www.wsl.ch/fe/gebirgshydrologie/schnee_hydro/oshd/index_EN">Operational Snow-Hydrological Service</A>
 * of the <A HREF="www.wsl.ch">WSL/SLF</A>. The data is written as Matlab
 * <A HREF="http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf">binary files (.mat)</A>, one per timestep for all stations and meteorological parameters, 
 * available on an access-controlled server after each new <A HREF="www.cosmo-model.org/">COSMO</A> run. It therefore requires a third party 
 * library to read this file format: the Open Source <A HREF="https://sourceforge.net/projects/matio/">MatIO</A> library. This can be installed directly from
 * the default repositories under Linux or installed by downloading the proper package for Windows or OsX.
 * 
 * \note If non-ascii characters have been used and the file has been created under Windows, some of the strings might end up using the UTF-16 encoding.
 * This requires a recent version of libmatio (see <A HREF="https://github.com/tbeu/matio/issues/34">this issue</A>). Another option would be to 
 * add at the beginning of the Matlab routine a call to *feature('DefaultCharacterSet', 'UTF8')* in order to switch from the current default (which can be read by the same call, 
 * ommitting the 'UTF8' option) to the (<A HREF="http://blog.omega-prime.co.uk/?p=150">partial</A>) UTF-8  encoding of Matlab.
 * 
 * @section oshd_data_structure Data structure
 * The files are named with the following schema: <i>COSMODATA_{timestep}_{cosmo model version}{mode}_{runtime}.{meteo_ext}</i> with the following possible values:
 *     + *timestep* is written as purely numeric ISO with minute resolution;
 *     + *cosmo model version* is currently C1E;
 *     + *mode* is a one letter code (FC for Forecast);
 *     + *run time* is the purely numeric ISO date and time of when COSMO produced the dataset;
 *     + *meteo_ext* is currently .mat.
 * 
 * The station data files have the following internal data structure (represented as "name {data type}" and with meteo parameter as one of prcs, wnss, wnds, wnsc, lwrc, sdri, sdrd, sdfd, tais, taic, rhus, pais, pail):
 * @verbatim
        ├── acro {1x688 array of arrays of char variable}
        ├── NODATA_value {1x1 double variable}
        ├── meteo parameter {struct}
        │   ├── data {1x688 array of doubles variable}
        │   ├── name {array of char variable}
        │   ├── unit {array of char variable}
        │   └── source {array of char variable}
        ├── … more meteo parameters
        └── time {1x2 array of doubles variable}
  @endverbatim
 * 
 * The stations' acronyms follow a fixed order but their coordinates must be provided in a separate file, given as *METAFILE* key (see below). 
 * If the number of stations in this list does not match the number of stations declared in one of the meteorological file, an exception will be thrown. 
 * This file must have the following structure (the *x* and *y* coordinates being the CH1903 easting and northing, respectively): 
 * @verbatim
      statlist {1x1 struct}
        ├── acro {1x688 array of arrays of char variable}
        ├── name {1x688 array of arrays of char variable}
        ├── x {1x688 array of doubles variable}
        ├── y {1x688 array of doubles variable}
        ├── z {1x688 array of doubles variable}
        └── … various other arrays
  @endverbatim
 *
 * The gridded data have the following structure:
 * @verbatim
      grid {1x1 struct}
        ├── ncols {1x1 array of doubles}
        ├── nrows {1x1 array of doubles}
        ├── xllcorner {1x1 array of doubles}
        ├── yllcorner {1x1 array of doubles}
        ├── cellsize {1x1 array of doubles}
        ├── NODATA_value {1x1 array of doubles}
        ├── data {nrows x ncols array of double}
        └── desc {array of char}
  @endverbatim
 *
 *
 * @section oshd_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - METEOPATH: directory containing all the data files with the proper file naming schema; [Input] section
 * - METEOPATH_RECURSIVE: should *meteopath* be searched recursively for files? (default: false); [Input] section
 * - STATION#: input stations' IDs (in METEOPATH). As many stations' IDs as needed may be specified
 * - METAFILE: file containing the stations' IDs, names and location; [Input] section (either within METEOPATH if not path is 
 provided or within the provided path)
 * - DEMFILE: for reading the data as a DEMObject
 * - GRID2DPATH: meteo grids directory where to read/write the grids; [Input] and [Output] sections
 * - GRIDPATH_RECURSIVE: if set to true, grids will be searched recursively in GRID2DPATH (default: false)
 * - OSHD_DEBUG: write out extra information to better show what is in the files
 *
 * @section oshd_example Example use
 * @code
 * [Input]
 * METEO = OSHD
 * METEOPATH = /local/LATEST_03h_RUN
 * METEOPATH_RECURSIVE = true
 * METAFILE  = STAT_LIST.mat ;another possibility could be /local/metadata/STAT_LIST.mat
 * STATION1  = ATT2
 * STATION2  = WFJ2
 * @endcode
 *
 */

const char* OshdIO::meteo_ext = "mat";
const double OshdIO::in_dflt_TZ = 0.; //COSMO data is always GMT
std::map< std::string, MeteoGrids::Parameters > OshdIO::params_map;
std::map< MeteoGrids::Parameters, std::string > OshdIO::grids_map;
const bool OshdIO::__init = OshdIO::initStaticData();

OshdIO::OshdIO(const std::string& configfile) : cfg(configfile), cache_meteo_files(), cache_grid_files(), vecMeta(), statIDs(),
               coordin(), coordinparam(), grid2dpath_in(), in_meteopath(), in_metafile(), nrMetadata(0), debug(false)
{
	parseInputOutputSection();
}

OshdIO::OshdIO(const Config& cfgreader) : cfg(cfgreader), cache_meteo_files(), cache_grid_files(), vecMeta(), statIDs(),
               coordin(), coordinparam(), grid2dpath_in(), in_meteopath(), in_metafile(), nrMetadata(0), debug(false)
{
	parseInputOutputSection();
}

void OshdIO::parseInputOutputSection()
{
	cfg.getValue("OSHD_DEBUG", "INPUT", debug, IOUtils::nothrow);
	cfg.getValue("COORDSYS", "Input", coordin);
	cfg.getValue("COORDPARAM", "Input", coordinparam, IOUtils::nothrow);
	
	const std::string meteo_in = IOUtils::strToUpper( cfg.get("METEO", "Input", "") );
	if (meteo_in == "OSHD") {//keep it synchronized with IOHandler.cc for plugin mapping!!
		std::vector< std::string > vecIDs;
		cfg.getValues("STATION", "INPUT", vecIDs);
		for (auto& ID : vecIDs) statIDs.push_back( station_index( ID ) );
		
		cfg.getValue("METEOPATH", "Input", in_meteopath);
		const bool is_recursive = cfg.get("METEOPATH_RECURSIVE", "Input", false);
		cache_meteo_files = scanMeteoPath(in_meteopath, is_recursive);
		if (debug) {
			std::cout << "Meteo files cache content:\n";
			for(size_t ii=0; ii<cache_meteo_files.size(); ii++) std::cout << cache_meteo_files[ii].toString() << "\n";
		}

		cfg.getValue("METAFILE", "INPUT", in_metafile);
		if (FileUtils::getFilename(in_metafile) == in_metafile) { //ie there is no path in the provided filename
			in_metafile = in_meteopath + "/" + in_metafile;
		}
	}

	const std::string grid_in = IOUtils::strToUpper( cfg.get("GRID2D", "Input", "") );
	if (grid_in == "OSHD") {//keep it synchronized with IOHandler.cc for plugin mapping!!
		grid2dpath_in.clear();
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
		const bool is_recursive = cfg.get("GRIDPATH_RECURSIVE", "Input", false);
		cache_grid_files = scanMeteoPath(grid2dpath_in, is_recursive);
		if (debug) {
			std::cout << "Grid files cache content:\n";
			for(size_t ii=0; ii<cache_grid_files.size(); ii++) std::cout << cache_grid_files[ii].toString() << "\n";
		}
	}
}

bool OshdIO::initStaticData()
{
	params_map[ "lwrc" ] = MeteoGrids::ILWR;
	params_map[ "pail" ] = MeteoGrids::P;
	params_map[ "prcs" ] = MeteoGrids::PSUM;
	params_map[ "rhus" ] = MeteoGrids::RH;
	params_map[ "tais" ] = MeteoGrids::TA;
	params_map[ "wnsc" ] = MeteoGrids::VW;
	params_map[ "wnds" ] = MeteoGrids::DW;
	params_map[ "sdrd" ] = MeteoGrids::ISWR_DIR;
	params_map[ "sdfd" ] = MeteoGrids::ISWR_DIFF;

	grids_map[ MeteoGrids::ILWR ] = "ilwc";
	grids_map[ MeteoGrids::P ] = "pair";
	grids_map[ MeteoGrids::PSUM ] = "prec"; //in mm/ts
	grids_map[ MeteoGrids::RH ] = "rcor"; //old: rhum
	grids_map[ MeteoGrids::TA ] = "tcor"; //old:tair
	grids_map[ MeteoGrids::VW ] = "wcor"; //old: wind
	grids_map[ MeteoGrids::DW ] = "wdir";
	grids_map[ MeteoGrids::ISWR_DIFF ] = "idfc";
	grids_map[ MeteoGrids::ISWR_DIR ] = "idrc";
	grids_map[ MeteoGrids::ALB ] = "albd";

	return true;
}

//This builds an index of which timesteps are provided by which files, always keeping the most recent run when
//multiple files provide the same timesteps. The file names are read as: COSMODATA_{timestep}_C1EFC_{runtime}.{meteo_ext}
std::vector< struct OshdIO::file_index > OshdIO::scanMeteoPath(const std::string& meteopath_in, const bool& is_recursive)
{
	//matching file names such as COSMODATA_202211022300_C1EFC_202211010300.mat
	//dirty trick: make sure that {meteo_ext} does not contain unescaped regex special chars!!
	static const std::regex filename_regex("COSMODATA_([0-9]{12})_C1EFC_([0-9]{12})\\." + std::string(meteo_ext));
	std::smatch filename_matches;
	
	std::vector< struct OshdIO::file_index > data_files;
	const std::list<std::string> dirlist( FileUtils::readDirectory(meteopath_in, "COSMODATA", is_recursive) );

	std::map<std::string, size_t> mapIdx; //make sure each timestamp only appears once, ie remove duplicates
	for (const auto& file_and_path : dirlist) {
		const std::string filename( FileUtils::getFilename(file_and_path) );
		
		if (!std::regex_match(filename, filename_matches, filename_regex)) continue;
		const std::string date_str( filename_matches.str(1) );
		const std::string run_date( filename_matches.str(2) );

		//do we already have an entry for this date?
		size_t idx = IOUtils::npos;
		const std::map<std::string, size_t>::const_iterator it_map = mapIdx.find( date_str );
		if (it_map!=mapIdx.end()) {
			idx = it_map->second;
			if (data_files[idx].run_date>run_date) continue;
		}

		//we don't have an entry or it is too old -> create new entry / replace existing one
		const std::string path( FileUtils::getPath(file_and_path) );
		Date date;
		IOUtils::convertString(date, date_str, in_dflt_TZ);
		const file_index elem(date, path, filename, run_date);
		if (idx==IOUtils::npos) {
			data_files.push_back( elem );
			mapIdx[ date_str ] = data_files.size()-1;
		} else {
			data_files[ idx] = elem;
			mapIdx[ date_str ] = idx;
		}
	}

	std::sort(data_files.begin(), data_files.end());
	return data_files;
}

bool OshdIO::list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& results)
{
	results.clear();
	
	for (size_t ii=0; ii<cache_grid_files.size(); ii++) {
		const Date date( cache_grid_files[ii].date );
		if (date<start) continue;
		if (date>end) break;
		
		//we consider that all parameters are present at all timesteps
		for (auto const& grid : grids_map) {
			results[date].insert( grid.first );
		}
	}
	
	return true;
}

size_t OshdIO::getFileIdx(const std::vector< struct file_index >& cache, const Date& start_date)
{
	if (cache.empty()) throw InvalidArgumentException("No input files found or configured!", AT);
	if (cache.size()==1) return 0; //no other possibility

	//try to find a good match
	for (size_t idx=1; idx<cache.size(); idx++) {
		if (start_date>=cache[idx-1].date && start_date<cache[idx].date) {
			return --idx;
		}
	}

	//not found, we take the closest timestamp we have (ie very beginning or very end)
	if (start_date<cache.front().date) return 0;
	else return cache.size()-1;
}

void OshdIO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	if (vecMeta.empty()) fillStationMeta();
	vecStation = vecMeta;
}

void OshdIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	size_t file_idx = getFileIdx( cache_meteo_files, dateStart );
	Date station_date( cache_meteo_files[file_idx].date );
	if (station_date>dateEnd || cache_meteo_files.back().date<dateStart) return; //the requested period is NOT in the available files

	const size_t nr_files = cache_meteo_files.size();
	const size_t nrIDs = statIDs.size();
	
	if (vecMeta.empty()) fillStationMeta(); //this also fills vecIdx
	vecMeteo.resize( nrIDs );
	do {
		//create empty MeteoData for the current timestep
		for (size_t jj=0; jj<nrIDs; jj++) {
			MeteoData md( station_date, vecMeta[jj] );
			md.addParameter("ISWR_DIFF");
			md.addParameter("ISWR_DIR");
			vecMeteo[jj].push_back( md );
		}
		
		//read all stations and all parameters for the current time step
		const std::string path( in_meteopath + "/" + cache_meteo_files[ file_idx ].path );
		const std::string file_and_path( path + "/" + cache_meteo_files[ file_idx ].filename );
		readFromFile(file_and_path, station_date, vecMeteo);
		
		//for convenience, we directly compute ISWR global from DIR and DIFF
		for (size_t jj=0; jj<nrIDs; jj++) 
			vecMeteo[jj].back()( MeteoData::ISWR ) = vecMeteo[jj].back()( "ISWR_DIR" ) + vecMeteo[jj].back()( "ISWR_DIFF" );
		
		file_idx++;
		station_date = ((file_idx)<nr_files)? cache_meteo_files[file_idx].date : dateEnd+1.;
	} while (file_idx<nr_files && station_date<=dateEnd);
}

void OshdIO::readFromFile(const std::string& file_and_path, const Date& in_timestep, std::vector< std::vector<MeteoData> >& vecMeteo) const
{
	if (debug) matWrap::printFileStructure(file_and_path, in_dflt_TZ);
	mat_t *matfp = Mat_Open(file_and_path.c_str(), MAT_ACC_RDONLY);
	if ( nullptr == matfp ) throw AccessException(file_and_path, AT);

	//open the file 
	matvar_t *matvar = nullptr;
	const size_t nrIDs = statIDs.size();
	double nodata = -999.;
	std::vector<std::string> vecAcro;
	
	//extract each parameter one by one for the selected stations
	while ( (matvar = Mat_VarReadNextInfo(matfp)) != nullptr ) {
		const std::string varName( matvar->name );
		
		if (varName=="NODATA_value") {
			nodata = matWrap::readDouble(file_and_path, "NODATA_value", matfp, matvar);
			continue;
		}
		
		if (varName=="time") {
			const std::vector<double> vecTime( matWrap::readDoubleVector(file_and_path, "time", matfp, matvar) );
			//vecTime is a vector of two elements: begin date and end date of timestep
			if (vecTime.size()!=2) {
				throw InvalidFormatException("Two time step must be present in the 'time' vector: runtime and simulation time", AT);
			}
			Date timestep;
			timestep.setMatlabDate( vecTime[0], in_dflt_TZ ); //take the begin date
			if (in_timestep!=timestep) throw InvalidArgumentException("the in-file timestep and the filename time step don't match for for '"+file_and_path+"'", AT);
			continue;
		}
		
		if (varName=="acro") {
			vecAcro = matWrap::readStringVector(file_and_path, "acro", matfp, matvar);
			if (vecAcro.size()!=nrMetadata) {
				std::ostringstream ss;
				ss << "the number of stations changes between the metadata file (" << nrMetadata << ") and the data file '";
				ss << FileUtils::getFilename( file_and_path ) << "' (" << vecAcro.size() << ")\n";
				throw IndexOutOfBoundsException(ss.str(), AT);
			}
			for (size_t ii=0; ii<nrIDs; ii++) { //check that the IDs still match
				if (statIDs[ii].oshd_acro != vecAcro[ statIDs[ii].idx ]) {
					throw InvalidFormatException("station '"+statIDs[ii].id+"' is not listed in the same position as previously in file '"+file_and_path+"'", AT);
				}
			}
			continue;
		}
		
		//processing of meteo parameters
		const std::map< std::string, MeteoGrids::Parameters>::iterator it = params_map.find( varName );
		if (it==params_map.end()) continue;
		
		const std::string parname( MeteoGrids::getParameterName( it->second ) );
		const size_t parindex = vecMeteo.back().back().getParameterIndex( parname );
		const std::string units( matWrap::readString(file_and_path, "unit", matfp, matvar) );
		const std::vector<double> vecRaw( matWrap::readDoubleVector(file_and_path, "data", matfp, matvar) );
		if (vecAcro.size() != vecRaw.size()) throw InvalidFormatException("'acro' and 'data' arrays don't match in file '"+file_and_path+"'", AT);
		
		for (size_t ii=0; ii<nrIDs; ii++) {
			const double raw_value = vecRaw[ statIDs[ii].idx ];
			if (raw_value!=nodata) { //otherwise it keeps its initial value of IOUtils::nodata
				vecMeteo[ii].back()( parindex ) = convertUnits( raw_value, units, parindex, file_and_path);
			}
		}

		//cleanup pointers before looping to avoid memory leaks
		Mat_VarFree(matvar);
		matvar = nullptr;
	}
	
	Mat_Close(matfp);
}

//NOTE It seems that recent versions contain multibyte encoding and this is not supported by matio, leading to trucated units (at best)
double OshdIO::convertUnits(const double& val, const std::string& units, const size_t& param, const std::string& filename)
{
	if (units=="%") return val/100.;
	if (units=="m") return val;
	if (units=="cm") return val/100.;
	if (units=="mm") {
		if (param==MeteoData::PSUM) return val;
		else return val/1000.;
	}
	if (units=="\x3F\x43") return val+Cst::t_water_freezing_pt; //unknown encoding hex for '°C'
	if (units=="\xB0\x43") return val+Cst::t_water_freezing_pt; //ISO-8859-1 hex for '°C'
	if (units=="\xB0" && param==MeteoData::TA) return val+Cst::t_water_freezing_pt; //ISO-8859-1 hex for '°'
	if (units=="\x3F" && param==MeteoData::DW) return val; //unknown encoding hex for '°'
	if (units=="\xB0" && param==MeteoData::DW) return val; //ISO-8859-1 hex for '°'
	if (units=="deg") return val;
	if (units=="K") return val;
	if (units.empty()) return val;
	if (units=="Pa") return val;
	if (units=="W/m2") return val;
	if (units=="m/s") return val;
	else {
		std::ostringstream os;
		for(size_t ii=0; ii<units.size(); ++ii)
			os << " " << std::hex << static_cast<unsigned int>( units[ii] );
		throw IOException("Unknown units '"+units+"' (#"+os.str()+") for parameter "+MeteoData::getParameterName(param)+" in file "+filename, AT);
	}
	
	return val;
}

void OshdIO::fillStationMeta()
{
	//the STATION_LIST.mat file contains entries formated such as:
	//acro="SLF.WFJ2", name="WeissfluhjochSchneestation (ENET)"
	//regex for removing appended networks names
	static const std::regex name_regex("([^\\(\\)]+) (\\([a-zA-Z]+\\))(.*)");
	std::smatch name_matches;

	if (debug) matWrap::printFileStructure(in_metafile, in_dflt_TZ);
	vecMeta.resize( statIDs.size(), StationData() );
	mat_t *matfp = Mat_Open(in_metafile.c_str(), MAT_ACC_RDONLY);
	if ( nullptr == matfp ) throw AccessException(in_metafile, AT);

	matvar_t *matvar = Mat_VarReadInfo(matfp, "statlist");
	if (matvar==nullptr) throw NotFoundException("structure 'statlist' not found in file '"+in_metafile+"'", AT);
	
	const std::vector<std::string> vecAcro( matWrap::readStringVector(in_metafile, "acro", matfp, matvar) );
	const std::vector<std::string> vecNames( matWrap::readStringVector(in_metafile, "name", matfp, matvar) );
	const std::vector<double> easting( matWrap::readDoubleVector(in_metafile, "x", matfp, matvar) );
	const std::vector<double> northing( matWrap::readDoubleVector(in_metafile, "y", matfp, matvar) );
	const std::vector<double> altitude( matWrap::readDoubleVector(in_metafile, "z", matfp, matvar) );
	nrMetadata = vecAcro.size();
	Mat_VarFree(matvar);
	Mat_Close(matfp);
	
	if (debug) {
		for (size_t ii=0; ii<nrMetadata; ii++) 
			std::cout << std::setw(8) << vecAcro[ii] << std::setw(40) << vecNames[ii] << std::setw(8) << easting[ii] << std::setw(8) << northing[ii] << std::setw(8) << altitude[ii] << "\n";
		std::cout << endl;
	}
	
	buildVecIdx(vecAcro);
	for (size_t ii=0; ii<statIDs.size(); ii++) {
		const size_t idx = statIDs[ii].idx;
		Coords location(coordin, coordinparam);
		location.setXY(easting[idx], northing[idx], altitude[idx]);
		
		//if the network name has been appended, remove it. We also remove spaces, just in case
		std::string name( vecNames[idx] );
		if (std::regex_match(name, name_matches, name_regex)) {
			name = name_matches.str(1);
		}

		const StationData sd(location, statIDs[ii].id, name);
		vecMeta[ii] = sd;
	}
}

//Fill vecIdx so it contains for all IDs in the order of their appearance in the ini file, their index in the .mat files
//the STATION_LIST.mat file contains entries formated such as:
//acro="SLF.WFJ2", name="WeissfluhjochSchneestation (ENET)"
void OshdIO::buildVecIdx(std::vector<std::string> vecAcro)
{
	const size_t nrIDs = statIDs.size();
	if (nrIDs==0) throw InvalidArgumentException("Please provide at least one station ID to read!", AT);
	
	const size_t nrAcro = vecAcro.size();
	//regex for correcting the stations' Acro into the correct ones: matching things like 'SLF.WFJ2'
	static const std::regex acro_regex("([a-zA-Z]+)\\.([a-zA-Z]+)([0-9]*)");
	std::smatch acro_matches;
	
	//cleaning up all stations IDs and looking for the user requested station IDs
	for (size_t jj=0; jj<nrAcro; jj++) {
		std::string cleanAcro = vecAcro[jj];
		
		if (std::regex_match(vecAcro[jj], acro_matches, acro_regex)) {
			//rules applied by oshd: stations like '*WFJ' have been renamed as 'MCH.WFJ2'
			//stations  like '#DOL' have been renamed as 'MCH.DOL1', stations in AUT, DE, FR, IT have been added
			const std::string provider( acro_matches.str(1) );
			if (provider=="MCH") {
				const std::string st_nr( acro_matches.str(3) );
				if (st_nr=="2") cleanAcro = "*" + acro_matches.str(2);
				if (st_nr=="1") cleanAcro = "#" + acro_matches.str(2);
			} else {
				cleanAcro = acro_matches.str(2) + acro_matches.str(3);
			}
		}
		
		//does the cleaned station ID match a user-requested station ID?
		for (size_t ii=0; ii<nrIDs; ii++) {
			if (statIDs[ii].id==cleanAcro) {
				statIDs[ii].idx = jj;
				statIDs[ii].oshd_acro = vecAcro[jj];
				break;
			}
		}
	}
	
	//making sure all user requested station IDs have been found
	for (size_t ii=0; ii<nrIDs; ii++) {
		if (statIDs[ii].idx==IOUtils::npos)
			throw NotFoundException("station ID '"+statIDs[ii].id+"' could not be found in the provided metadata", AT);
	}
}

void OshdIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename)
{
	if (debug) matWrap::printFileStructure(filename, in_dflt_TZ);
	mat_t *matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
	if ( nullptr == matfp ) throw AccessException(filename, AT);

	matvar_t *matvar = Mat_VarReadInfo(matfp, "grid");
	if (matvar==nullptr) throw NotFoundException("structure 'grid' not found in file '"+filename+"'", AT);
	if (matvar->class_type!=MAT_C_STRUCT) throw InvalidFormatException("The matlab file should contain 1 structure", AT);

	const double xllcorner( matWrap::readDouble(filename, "xllcorner", matfp, matvar) );
	const double yllcorner( matWrap::readDouble(filename, "yllcorner", matfp, matvar) );
	Coords location(coordin, coordinparam);
	location.setXY(xllcorner, yllcorner, IOUtils::nodata);

	const double cellsize( matWrap::readDouble(filename, "cellsize", matfp, matvar) );
	const double ncols( matWrap::readDouble(filename, "ncols", matfp, matvar) );
	const double nrows( matWrap::readDouble(filename, "nrows", matfp, matvar) );

	//Initialize the 2D grid
	grid_out.set(static_cast<size_t>(ncols), static_cast<size_t>(nrows), cellsize, location);
	matWrap::readDoubleArray(filename, "data", matfp, matvar, grid_out.grid2D);

	Mat_Close(matfp);
}

void OshdIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const size_t file_idx = getFileIdx( cache_grid_files, date );
	Date grid_date( cache_grid_files[file_idx].date );
	if (grid_date!=date) return; //the requested date is NOT in the available files

	//build the proper file name
	const std::string file_suffix( cache_grid_files[ file_idx ].filename );
	const std::string path( grid2dpath_in + "/" + cache_grid_files[ file_idx ].path );

	if (grids_map.find(parameter)!=grids_map.end()) {
		const std::string filename( path + "/" + grids_map[parameter] +file_suffix );
		read2DGrid(grid_out, filename);
	} else
		throw NotFoundException("Parameter "+MeteoGrids::getParameterName(parameter)+" currently not supported", AT);

	//units corrections
	if (parameter==MeteoGrids::TA) grid_out += Cst::t_water_freezing_pt;
	else if (parameter==MeteoGrids::RH) grid_out /= 100.;
}

void OshdIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	read2DGrid(dem_out, filename);
}

} //namespace
