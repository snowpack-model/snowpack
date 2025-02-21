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
#include <meteoio/FileUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/plugins/GRIBIO.h>
#include <algorithm>


namespace mio {
using namespace codes;
/**
* @page gribio GRIBIO
* @section gribio_format Format
* 
* @note This plugin requires the <a href="https://confluence.ecmwf.int/display/ECC/ecCodes+Home">ecCodes library</a> to be installed (at least in version 2.27.0, lower versions won't compile).
* 
* This plugin reads Meteorological data from GRIB files, specified by the WMO in their Manual of Codes 306 (have a look at this <a href="https://en.wikipedia.org/wiki/GRIB">short introduction to GRIB</a>). 
* 
* In essence a Grib file contains an arbitrary number of messages. Each message is a self-explanatory record of some gridded data. Messages contain all the Metadata 
* and other relevant information, to correctly identify and interpret the data. The GRIB format is binary, and is therefore not human readable. There exists a variety 
* of software to read GRIB files, see Wiki entry. Here we use the <a href="https://confluence.ecmwf.int/display/ECC/ecCodes+installation">ecCodes library</a> from ECMWF to read GRIB files.
* 
* One big problem is, that the GRIB format is very flexible, and therefore the encoding of parameters... is specified in Tables provided by the data originating centers. 
* ecCodes uses many common tables to decode the data. However, (very unlikely) it might happen, that there is no table available.
* 
* The ECMWF provides a <a href="https://codes.ecmwf.int/grib/param-db/">parameter database</a> where all the information about parameter codes can be found. This is why we decided to 
* use their paramId as a key to identify the parameters. As there are many ways of specifing data, e.g. Air temperature at 2m as its own parameter 2mTempAir, the temperature read at 
* 2m: Temp at 2m, or the Atmospheric temperature at 2m... This is why we use a parameter table, that can be adjusted to the users needs. It can be found under doc/resources/GRIB_param.tbl. 
* It is a very basic text file, where empty lines, and lines starting wiht '#' are ignored. The default table is self explanatory, so read through the comments to see how to specify the 
* desired parameters.
* 
* Additionally, to read vertically spaced parameters the typeOfLevel and level values need to be specified. If there is no level needed, such as for a 
* surface level, just set the level to 0.
* 
* For now it is only possible to read parameters, that are predefined in MeteoIO (see MeteoGrids::Parameters, and MeteoData::Parameters). 
* 
* The GRIB_TABLE template can be found at <a href="https://gitlabext.wsl.ch/snow-models/meteoio/-/blob/master/doc/resources/GRIB_param.tbl?ref_type=heads">GRIB_param.tbl</a>
* If these are the correct Ids there is nothing to do, if you need to change something, just copy the file, and let meteoio know where to find it.
* 
* As it is incredibly flexible, we narrowed down the structure of GRIB files that MeteoIO is able to read. Each file can contain one or multiple timepoints, 
* but needs to contain all the parameters in one file. So having multiple files for the same timepoint is not supported.
* 
*
* @section paramId Finding the paramId
* If the default paramId is not the correct one for your GRIB file, you need to specify the correct ids in the parameer table. To find out, which 
* parameter maps to which paramId, we recommend using either the <a href="https://confluence.ecmwf.int/display/ECC/GRIB+tools">ecCodes grib tools </a> or the 
* <a href="https://metview.readthedocs.io/en/latest/">metview</a> software. Both tools can be used to inspect the GRIB file and find the correct paramId. 
* If you want to use the Grib Tools from ecCodes, you can use the following command:
* @code
* grib_ls -p paramId,typeOfLevel,name filename.gribext 
* @endcode
* This will print the paramId and the name of the parameter for each message in the file. If there are a lot of messages in the file, i.e. many timepoints, 
* you will have a large output, which will be repeated a lot.
* 
* Metview is the easier option, as it is installable with conda via:
* @code
* conda install metview -c conda-forge
* @endcode
* After installation, you can open the GRIB file in Metview, and inspect the parameters with:
* @code
* metview -e grib filename.gribext
* @endcode
* This will open the Metview GUI, where you can inspect the File. You will get a list of all messages on the left, and on the right you can find information about 
* the selected message. You will be able to easily see, which message contains which parameter. To find out the paramId, and typeOfLevel, 
* select the Namespaces in the toolbar on the right, and choose "parameters" under which you can find the paramId and "vertical" for the typeOfLevel.
*
* @section gribio_keywords Keywords
* This plugin uses the following keywords:
* When reading Grid2d files:
* - GRID2DPATH: Path to the repository containing the GRIB files
* - GRID2DEXT: Extension of the GRIB files. If not specified, all files in the directory will be read
* - GRID2DPATTERN: Pattern to filter the files in the directory. If not specified, all files in the directory will be read
* - GRIB_TABLE: Path to the parameter table. If not specified, the default table will be used
* - VERBOSE: If set to true, the plugin will print out information about the read files
* - RECURSIVE: If set to true, the plugin will search the directory recursively
* 
* When reading DEM files:
* - DEMFILE: Path to the DEM file
* - GRIB_DEM_UPDATE: If set to true, the plugin will update the DEM 
* - VERBOSE: If set to true, the plugin will print out information about the read files
* - GRIB_TABLE: Path to the parameter table. If not specified, the default table will be used
* 
* When reading virtual stations from the files:
* - METEOPATH: Path to the repository containing the GRIB files
* - METEOEXT: Extension of the GRIB files. If not specified, all files in the directory will be read
* - METEOPATTERN: Pattern to filter the files in the directory. If not specified, all files in the directory will be read
* - STATION#: The position of station # (latlon,xy...)
* - GRIB_TABLE: Path to the parameter table. If not specified, the default table will be used
* - VERBOSE: If set to true, the plugin will print out information about the read files
* - RECURSIVE: If set to true, the plugin will search the directory recursively
* 
* @note Grid2d and Meteo Files need to be seperated by either pattern or extension, if both are used.
* @note It is only possible to provide one GRIB table for all the different types of files.
* 
* @section example Example
* @code 
* [INPUT]
* METEO = GRIB
* METEOPATH = ../
* METEOEXT = .grib
* METEOPATTERN = meteo
* STATION1 = latlon(5,9,-1)
* 
* DEM = GRIB
* DEMFILE = ../dem.grib
* GRIB_DEM_UPDATE = true
* 
* GRID2D = GRIB
* GRID2DPATH = ../
* GRID2DEXT = .grib
* GRID2DPATTERN = grid
* 
* GRIB_TABLE = test.tbl
* VERBOSE = true
* RECURSIVE = false
* @endcode
* 
* 
* 
* @todo When add parameter is available for MeteoGrids, support adding unknown parameters via the parameter table for GRIB
*/

    static const double plugin_nodata = IOUtils::nodata; // plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
    static const std::string default_extension = ".grib";

    // ----------------------------- INITIALIZE -----------------------------
    GRIBIO::GRIBIO(const std::string &configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), meteopath_in(), grid2dpath_in(), table_path(), meteo_ext(default_extension),
                                                    meteo_pattern(), grid2d_ext(default_extension), grid_2d_pattern(), recursive_search(false), verbose(false), update_dem(false), 
                                                    bearing_offset(IOUtils::nodata), latitudeOfNorthernPole(), longitudeOfNorthernPole(), llcorner_initialized(false), llcorner(), 
                                                    cellsize(), factor_x(), factor_y(), grid_initialized(false), meteo_initialized(false), parameter_table(), 
                                                    cache_meteo(), cache_grid2d(), vecPts() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
        initialize();
    }

    GRIBIO::GRIBIO(const Config &cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), meteopath_in(), grid2dpath_in(), table_path(), meteo_ext(default_extension),
                                                    meteo_pattern(), grid2d_ext(default_extension), grid_2d_pattern(), recursive_search(false), verbose(false), update_dem(false), 
                                                    bearing_offset(IOUtils::nodata), latitudeOfNorthernPole(), longitudeOfNorthernPole(), llcorner_initialized(false), llcorner(), 
                                                    cellsize(), factor_x(), factor_y(), grid_initialized(false), meteo_initialized(false), parameter_table(), 
                                                    cache_meteo(), cache_grid2d(), vecPts() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
        initialize();
    }

    void GRIBIO::initialize() {
        setOptions();
        initTable();
    }

    void GRIBIO::setOptions() {
        const std::string in_2d = IOUtils::strToUpper(cfg.get("GRID2D", "Input", ""));
        const std::string in_meteo = IOUtils::strToUpper(cfg.get("METEO", "Input", ""));
        const std::string in_dem = IOUtils::strToUpper(cfg.get("DEM", "Input", ""));

        bool meteo_and_grid = in_2d == "GRIB" && in_meteo == "GRIB";

        if (in_2d == "GRIB" || in_meteo == "GRIB" || in_dem == "GRIB") {
            cfg.getValue("GRIB_TABLE", "Input", table_path, IOUtils::nothrow);
            cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);
            cfg.getValue("RECURSIVE", "Input", recursive_search, IOUtils::nothrow);
        }

        if (in_2d == "GRIB") { // keep it synchronized with IOHandler.cc for plugin mapping!!
            cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
            cfg.getValue("GRID2DEXT", "Input", grid2d_ext, IOUtils::nothrow);
            if (grid2d_ext == "none")
                grid2d_ext.clear();
            cfg.getValue("GRID2DPATTERN", "Input", grid_2d_pattern, IOUtils::nothrow);
            if (grid_2d_pattern == "none")
                grid_2d_pattern.clear();
        }

        if (in_meteo == "GRIB") {
            cfg.getValue("METEOEXT", "Input", meteo_ext, IOUtils::nothrow);
            if (meteo_ext == "none")
                meteo_ext.clear();
            cfg.getValue("METEOPATTERN", "Input", meteo_pattern, IOUtils::nothrow);
            if (meteo_pattern == "none")
                meteo_pattern.clear();
        }

        if (in_dem == "GRIB") {
            cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, IOUtils::nothrow);
        }

        if (meteo_and_grid && meteo_ext == grid2d_ext && meteo_pattern == grid_2d_pattern && meteopath_in == grid2dpath_in)
            throw InvalidArgumentException("Meteo and Grid2D files cannot have the same naming and be located in the same place.", AT);
    }

    void GRIBIO::initTable() {

        if (table_path.empty()) {
#ifdef DEBUG
            std::cerr << "No GRIB table specified, using default table" << std::endl;
#endif
            parameter_table = GRIBTable();
        } else {
            parameter_table = GRIBTable(table_path);
        }
#ifdef DEBUG
        parameter_table.printTable();
#endif
    }

    // scan the path for files with the correct extension and pattern
    void GRIBIO::scanPath(const std::string &in_path, const std::string &in_ext, const std::string &in_pattern, std::vector<GRIBFile> &cache) {
        std::list<std::string> dirlist;
        FileUtils::readDirectory(in_path, dirlist, in_ext, recursive_search);

        std::list<std::string>::const_iterator it = dirlist.begin();
        while ((it != dirlist.end())) {
            const std::string filename = *it;
            if (filename.find(in_pattern) != std::string::npos) {
                const std::string fullpath = in_path + "/" + filename;
                cache.push_back(GRIBFile(fullpath, parameter_table.getIndexes())); // Files will not be read twice, as we enforce a seperation between meteo and grid
            }
            it++;
        }
    }

    // read all the specified virtual stations
    void GRIBIO::readStations(std::vector<Coords> &vecPoints) {
        cfg.getValue("METEOPATH", "Input", meteopath_in);

        std::vector<std::string> vecStation;
        cfg.getValues("STATION", "INPUT", vecStation);
        for (size_t ii = 0; ii < vecStation.size(); ii++) {
            Coords tmp(coordin, coordinparam, vecStation[ii]);
            if (!tmp.isNodata())
                vecPoints.push_back(tmp);

            std::cout << "\tRead virtual station " << vecPoints.back().toString(Coords::LATLON) << "\n";
        }
    }

    // ---------------------- STATIC HELPERS -------------------------------------------
    // find the index of the file containing the date
    static size_t findDate(const std::vector<GRIBFile> &cache, const Date &date) {
        auto it = std::find_if(cache.begin(), cache.end(), [&date](const GRIBFile &file) { return file.isValidDate(date); });

        return (it != cache.end()) ? std::distance(cache.begin(), it) : IOUtils::npos;
    }

    // do some conversions for parameters we know are different in GRIB
    static void handleConversions(Grid2DObject& grid_out, const double& paramId) { 
        if (paramId == 163) { // slope
            grid_out.grid2D *= 90.;
        } else if (paramId == 162) { // azi
            grid_out.grid2D *= Cst::to_deg;
            grid_out.grid2D -= 90.;
        } else if (paramId == 129) { // geopotential
            grid_out.grid2D /= Cst::gravity;
        } 
    }

    // wrapper to find all the messages containing the parameter information, depending on the type of parameter id
    static std::vector<CodesHandlePtr> findMessages(GRIBFile &file, const std::string &param_key, const std::string &level_key, const std::string &level_type, const std::string &paramID_string,
                                                    double paramID_double, long paramID_long, const bool& verbose, const std::string& param_name) {
        double npos_double = static_cast<double>(IOUtils::npos);
        long npos_long = static_cast<long>(IOUtils::npos);

        if (!paramID_string.empty() && paramID_double == npos_double && paramID_long == npos_long) {
            return file.listParameterMessages(param_key, paramID_string, level_key, level_type);
        } else if (paramID_string.empty() && paramID_double != npos_double && paramID_long == npos_long) {
            return file.listParameterMessages(param_key, paramID_double, level_key, level_type);
        } else if (paramID_string.empty() && paramID_double == npos_double && paramID_long != npos_long) {
            return file.listParameterMessages(param_key, paramID_long, level_key, level_type);
        } else {
            if (verbose) {
                std::cout << "No parameter id provided for "+param_name << std::endl;
            }
            return {};
        }   
    }

    // get all the indexing information for the parameter and find the messages containing the parameter information
    static std::vector<CodesHandlePtr> extractParameterInfoAndFindMessages(GRIBFile &file, const std::string &param_name, const GRIBTable &parameter_table, long &level_no, const bool& verbose) {
        // get the paramID
        std::string paramID_string;
        double paramID_double;
        long paramID_long;
        parameter_table.getParamId(param_name, paramID_string, paramID_double, paramID_long);

        // get additional information from the parameter table
        std::string param_key = parameter_table.getParamKey();
        std::string level_key = parameter_table.getLevelKey();
        level_no = parameter_table.getLevelNo(param_name);
        std::string level_type = parameter_table.getLevelType(param_name);
        // check type of level and where it is lost
        return findMessages(file, param_key, level_key, level_type, paramID_string, paramID_double, paramID_long, verbose, param_name);
    }
    
    // wrapper to be able to use it with MeteoGrids::Parameters
    static std::vector<CodesHandlePtr> extractParameterInfoAndFindMessages(GRIBFile &file, const MeteoGrids::Parameters &parameter, const GRIBTable &parameter_table, long &level_no, const bool& verbose) {
        std::string param_name = MeteoGrids::getParameterName(parameter);
        return extractParameterInfoAndFindMessages(file, param_name, parameter_table, level_no, verbose);
    };


    // ----------------------------- GRIDDED DATA -----------------------------
    // legacy to support reading a single grid from a file
    void GRIBIO::read2DGrid(Grid2DObject &grid_out, const std::string &i_name) {
        const std::string filename(grid2dpath_in + "/" + i_name);
        if (!FileUtils::fileExists(filename))
            throw AccessException(filename, AT); // prevent invalid filenames

        std::vector<CodesHandlePtr> handles = getMessages(filename);

        // the file should contain only one grid
        if (handles.empty())
            throw IOException("No grid found in file \"" + filename + "\"", AT);
        if (handles.size() > 1)
            throw IOException("Multiple grids found in file \"" + filename + "\". Please specify which to load.", AT);
        
        // read the data in the file
        read2Dlevel(handles.front(), grid_out, getGridParameters(handles.front()));
    }

    // read the grid for a parameter and date from a file
    void GRIBIO::read2DGrid(Grid2DObject &grid_out, const MeteoGrids::Parameters &parameter, const Date &date) {
        if (!grid_initialized) {
            scanPath(grid2dpath_in, grid2d_ext, grid_2d_pattern, cache_grid2d);
            grid_initialized = true;
        }

        // find the file which actually contains the date
        size_t idx = findDate(cache_grid2d, date);
        if (idx == IOUtils::npos) {
            if (verbose) {
                std::cout << "No grid found for the specified date" << std::endl;
            }
            return;
        }
        if (verbose) {
            std::cout << "Reading grid for date " << date.toString() << std::endl;
        }

        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(cache_grid2d[idx], parameter, parameter_table, level_no, verbose);

        if (messages.empty()) {
            if (verbose) {
                std::cout << "No messages found for the specified parameter "+MeteoGrids::getParameterName(parameter) << std::endl;
            }
            return;
        }

        // loop through the messages for the correct parameter, and find the one with the correct level at the correct date
        for (auto &m : messages) {
            long level = 0;
            getParameter(m, "level", level);
            // if multiple timepoints in a grib file, we need to find the correct one:
            if (getMessageDateGrib(m, 0) != date) {
                continue;
            }
            if (level == level_no) {
                read2Dlevel(m, grid_out, cache_grid2d[idx].getGridParams());
                return;
            }
        }
        if (verbose)
            std::cout << "No messages found for the specified level and date" << std::endl;
        return;
    };

    void GRIBIO::read2Dlevel(CodesHandlePtr &h, Grid2DObject &grid_out, const std::map<std::string, double> &grid_params) {
        setMissingValue(h, plugin_nodata);
        std::vector<double> values;
        getGriddedValues(h, values);

        double Ni = grid_params.at("Ni");
        double Nj = grid_params.at("Nj");

        if (!llcorner_initialized) {
            initializeLLCorner(grid_params);
        }

        // set the values of the grid
        grid_out.set(static_cast<size_t>(Ni), static_cast<size_t>(Nj), cellsize, llcorner);
        size_t i = 0;
        for (size_t jj = 0; jj < (unsigned)Nj; jj++) {
            for (size_t ii = 0; ii < (unsigned)Ni; ii++)
                grid_out(ii, jj) = values[i++];
        }

        // cells were not square, we have to resample
        if (factor_x != IOUtils::nodata && factor_y != IOUtils::nodata) {
            grid_out.grid2D = LibResampling2D::Bilinear(grid_out.grid2D, factor_x, factor_y);
        }

        // finalizing llcorner, i.e. llcorner was set as center before
        if (!llcorner_initialized) { // take into account aspect ration conversion for computing true llcorner
            llcorner.moveByXY(-.5 * (double)grid_out.getNx() * cellsize, -.5 * (double)grid_out.getNy() * cellsize);
            llcorner_initialized = true;

            grid_out.llcorner = llcorner;
        }

        // handle some conversions for parameters we know are different in GRIB
        double paramId;
        getParameter(h, "paramId", paramId);
        handleConversions(grid_out, paramId);
        if (verbose) {
            std::cout << "Read " << values.size() << " values from GRIB file" << std::endl;
            std::cout << "Parameter " << paramId << std::endl;
        }
    }

    void GRIBIO::initializeLLCorner(const std::map<std::string, double> &grid_params) {
        // most important: get cellsize. llcorner will be finalized AFTER aspect ration correction
        double cellsize_x, cellsize_y;
        llcorner = getGeolocalization(cellsize_x, cellsize_y, grid_params); // this is the center cell

        cellsize = (double)Optim::round(std::min(cellsize_x, cellsize_y) * 100.) / 100.; // round to 1cm precision for numerical stability
        if (std::abs(cellsize_x - cellsize_y) / cellsize_x > 1. / 100.) {
            factor_x = cellsize_x / cellsize;
            factor_y = cellsize_y / cellsize;
        }
    
    }

    // handle shifted grids, by computing the center and offsets
    // returns the center point as reference
    Coords GRIBIO::getGeolocalization(double &cellsize_x, double &cellsize_y, const std::map<std::string, double> &grid_params) {
        latitudeOfNorthernPole = grid_params.at("latitudeOfNorthernPole");
        longitudeOfNorthernPole = grid_params.at("longitudeOfNorthernPole");
        double ll_latitude = grid_params.at("ll_latitude");
        double ll_longitude = grid_params.at("ll_longitude");
        double ur_latitude = grid_params.at("ur_latitude");
        double ur_longitude = grid_params.at("ur_longitude");
        double angleOfRotationInDegrees = grid_params.at("angleOfRotationInDegrees");
        double Ni = grid_params.at("Ni");
        double Nj = grid_params.at("Nj");

        if (angleOfRotationInDegrees != 0.) {
            throw InvalidArgumentException("Rotated grids not supported!", AT);
        }

        double ur_lat, ur_lon, ll_lat, ll_lon;
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, ur_latitude, ur_longitude, ur_lat, ur_lon);
        double cntr_lat, cntr_lon; // geographic coordinates
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, .5 * (ll_latitude + ur_latitude), .5 * (ll_longitude + ur_longitude), cntr_lat, cntr_lon);


        double bearing;
        cellsize_x = CoordsAlgorithms::VincentyDistance(cntr_lat, ll_lon, cntr_lat, ur_lon, bearing) / (double)Ni;
        cellsize_y = CoordsAlgorithms::VincentyDistance(ll_lat, cntr_lon, ur_lat, cntr_lon, bearing) / (double)Nj;

        // determining bearing offset
        double delta_lat, delta_lon; // geographic coordinates
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, .5 * (ll_latitude + ur_latitude) + 1., .5 * (ll_longitude + ur_longitude), delta_lat, delta_lon);
        CoordsAlgorithms::VincentyDistance(cntr_lat, cntr_lon, delta_lat, delta_lon, bearing_offset);
        bearing_offset = fmod(bearing_offset + 180., 360.) - 180.; // turn into [-180;180)

        // returning the center point as reference
        Coords cntr(coordin, coordinparam);
        cntr.setLatLon(cntr_lat, cntr_lon, IOUtils::nodata);

        return cntr;
    }
    

    // ---------------------------- DIGITAL ELEVATION MODEL -----------------------------
    // Read a DEM from a file
    // uses the same GRIBTable as METEO and GRID2D
    void GRIBIO::readDEM(DEMObject &dem_out) {
        const std::string filename = cfg.get("DEMFILE", "Input");

        GRIBTable dem_table(table_path);
        GRIBFile dem_file(filename, dem_table.getIndexes());

        processSingleMessage(dem_out, dem_file, dem_table, MeteoGrids::DEM);

        if (update_dem) {
            dem_out.update();
        } else {
            const int dem_ppt = dem_out.getUpdatePpt();
            if (dem_ppt & DEMObject::SLOPE) {
                Grid2DObject slope;
                processSingleMessage(slope, dem_file, dem_table, MeteoGrids::SLOPE);
                dem_out.slope = slope.grid2D;
                Grid2DObject azi;
                processSingleMessage(azi, dem_file, dem_table, MeteoGrids::AZI);
                dem_out.azi = azi.grid2D;
            }
            if (dem_ppt & DEMObject::NORMAL || dem_ppt & DEMObject::CURVATURE) {
                // we will only update the normals and/or curvatures, then revert update properties
                if (dem_ppt & DEMObject::NORMAL && dem_ppt & DEMObject::CURVATURE)
                    dem_out.setUpdatePpt((DEMObject::update_type)(DEMObject::NORMAL | DEMObject::CURVATURE));
                else if (dem_ppt & DEMObject::NORMAL)
                    dem_out.setUpdatePpt(DEMObject::NORMAL);
                else if (dem_ppt & DEMObject::CURVATURE)
                    dem_out.setUpdatePpt(DEMObject::CURVATURE);

                dem_out.update();
                dem_out.setUpdatePpt((DEMObject::update_type)dem_ppt);
            }

            dem_out.updateAllMinMax();
        }
    }


    void GRIBIO::processSingleMessage(Grid2DObject &dem_out, GRIBFile &dem_file, const GRIBTable &dem_table, const MeteoGrids::Parameters &parameter) {
        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(dem_file, parameter, dem_table, level_no, verbose);


        if (messages.empty()) {
            throw IOException("No messages containing DEM information found." AT);
        }
        if (messages.size() > 1) {
            throw IOException("Multiple messages containing DEM information found in file" + dem_file.getFilename(), AT);
        }

        read2Dlevel(messages.front(), dem_out, dem_file.getGridParams());
    }

    // ---------------------------- METEO DATA -----------------------------
    // ---------------------------- STATIC HELPERS for METEO DATA -----------------------------
    
    // compare files by their start date, to be able to order them
    static bool compareByDate(const GRIBFile &a, const GRIBFile &b) { return a.getStartDate() < b.getStartDate(); }

    // compare a file to a date, to be able to find the first file conatining the date
    static bool compareToDate(const GRIBFile &a, const Date &b) { return a.getStartDate() < b; } 

    // create a vector of MeteoData objects for each station
    static std::vector<MeteoData> createMeteoDataVector(const std::vector<StationData>& stations, const Date& date) {
        std::vector<MeteoData> vecMeteo;
        const size_t npoints = stations.size();

        for (size_t ii = 0; ii < npoints; ii++) {
            MeteoData md;
            md.meta = stations[ii];
            md.date = date;
            vecMeteo.push_back(md);
        }

        return vecMeteo;
    }

    // process all messages related to a parameter
    static void processMessages(std::vector<CodesHandlePtr>& messages, const long& level_no, std::vector<MeteoData>& vecMeteo, const std::vector<double>& lats, const std::vector<double>& lons, const size_t& npoints, const size_t& par_index) {
        Date curr_date = vecMeteo.front().date;
        // loop through all messages, to find the one with correct date and level
        for (auto &m : messages) {
            if (getMessageDateGrib(m, 0) != curr_date)
                continue;
            long level = 0;
            getParameter(m, "level", level);
            if (level == level_no) {
                // we need to set the missing value, as there is no default 
                setMissingValue(m, plugin_nodata);
                std::vector<double> outlats_vec(npoints), outlons_vec(npoints), distances_vec(npoints), values(npoints);
                std::vector<int> indexes_vec(npoints);
                // find the nearest points in the grid to lats, lons and get the values, as well as the locations and distances
                getNearestValues_grib(m, lats, lons, outlats_vec, outlons_vec, distances_vec, values, indexes_vec);

                // set the values for each station
                for (size_t ii = 0; ii < npoints; ii++) {
                    vecMeteo[ii](par_index) = values[ii];
                }
            }
        }
    }



    // ---------------------------- METEO READING -----------------------------
    void GRIBIO::readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecvecMeteo) {
        if (!meteo_initialized) {
            readStations(vecPts);
            scanPath(meteopath_in, meteo_ext, meteo_pattern, cache_meteo);
            meteo_initialized = true;
            std::sort(cache_meteo.begin(), cache_meteo.end(), compareByDate);
        }

        // initialize the vectors containing the meteo data
        vecvecMeteo.clear();
        std::vector<double> lats(vecPts.size());
        std::vector<double> lons(vecPts.size());
        std::vector<StationData> stations;
        
        // flag to see if we already read the metadata
        bool meta_ok = false;


        // find the first file containing the start date
        auto it = std::lower_bound(cache_meteo.begin(), cache_meteo.end(), dateStart, compareToDate);
        // set the start index to the first file containing the start date
        size_t start_idx = 0;
        if (it != cache_meteo.end())
            start_idx = std::distance(cache_meteo.begin(), it);

        // if the start date is not in the first file, we start before that, for interpolations...
        if (start_idx > 0)
            start_idx--;

        // loop through all files, and read the meteo data
        for (size_t i = start_idx; i < cache_meteo.size(); i++) {
            // go through all timepoints in a file
            for (const Date& current_date : cache_meteo[i].getDates()) {
                if (current_date > dateEnd) {
                    break;
                }
                // read metadata
                if (!meta_ok) {
                    if (!readMeteoMeta(cache_meteo[i], vecPts, stations, lats, lons)) {
                        // some points have been removed vecPts has been changed -> re-reading
                        lats.clear();
                        lons.clear();
                        readMeteoMeta(cache_meteo[i], vecPts, stations, lats, lons);
                    }
                    vecvecMeteo.insert(vecvecMeteo.begin(), vecPts.size(), std::vector<MeteoData>()); // allocation for the vectors now that we know how many true stations we have
                    meta_ok = true;
                }

                // create a vector of MeteoData objects for each station
                std::vector<MeteoData> vecMeteoStations = createMeteoDataVector(stations, current_date);
                const size_t npoints = vecMeteoStations.size();

                // try to read all known parameters from the file
                for (size_t par_index=0; par_index <= MeteoData::Parameters::lastparam; par_index++) {
                    std::string param_name = MeteoData::getParameterName(par_index);
                    long level_no;
                    std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(cache_meteo[i], param_name, parameter_table, level_no, verbose);

                    if (messages.empty()) {
                        if (verbose) {
                            std::cout << "No messages found for parameter " << param_name << " in file " << cache_meteo[i].getFilename() << std::endl;
                        }
                        continue; 
                    }

                    // read the data for a given parameter
                    processMessages(messages, level_no, vecMeteoStations, lats, lons, npoints, par_index);
                }
                // fill the vector containing all stations and dates
                for (size_t jj=0; jj<vecPts.size(); jj++)
                    vecvecMeteo[jj].push_back(vecMeteoStations[jj]);
            }
        }
    }


    // ---------------------------- METEO READING HELPERS -----------------------------
    bool GRIBIO::removeDuplicatePoints(std::vector<Coords> &vecPoints, std::vector<double> &lats, std::vector<double> &lons) { // remove potential duplicates. Returns true if some have been removed
        const size_t npoints = vecPoints.size();
        std::vector<size_t> deletions;
        deletions.reserve(npoints);
        for (size_t ii = 0; ii < npoints; ii++) {
            const double lat = lats[ii];
            const double lon = lons[ii];
            for (size_t jj = ii + 1; jj < npoints; jj++) {
                if (lat == lats[jj] && lon == lons[jj]) {
                    deletions.push_back(jj);
                }
            }
        }

        // we need to erase from the end in order to keep the index unchanged...
        for (size_t ii = deletions.size(); ii-- > 0;) {
            const size_t index = deletions[ii - 1];
            vecPoints.erase(vecPoints.begin() + index);
        }

        if (!deletions.empty())
            return true;
        return false;
    }

    bool GRIBIO::readMeteoMeta(GRIBFile& file ,std::vector<Coords>& vecPoints, std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons) {

        // level no is not needed, but we need to pass it
        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(file, MeteoGrids::DEM, parameter_table, level_no, verbose);

        if (messages.empty()) {
            throw IOException("No messages containing DEM information found." AT);
        }
        if (messages.size() > 1 && verbose) {
            std::cout << "Multiple messages containing DEM information found in file" << file.getFilename() << std::endl;
        }
        const size_t npoints = vecPoints.size();

        // we simply use the first dem message to get the grid parameters
        CodesHandlePtr& h = messages.front();
        // if something went weirdly wrongs
        if (h == nullptr) {
            throw IOException("Can not find DEM grid in GRIB file!", AT);
        }

        // get the coordinate parameters
        std::map<std::string, double> grid_params = getGridParameters(h);
        latitudeOfNorthernPole = grid_params.at("latitudeOfNorthernPole");
        longitudeOfNorthernPole = grid_params.at("longitudeOfNorthernPole");
        long Ni = static_cast<long>(grid_params.at("Ni"));

        // build GRIB local coordinates for the points
        for (size_t ii = 0; ii < npoints; ii++) {
            CoordsAlgorithms::trueLatLonToRotated(latitudeOfNorthernPole, longitudeOfNorthernPole, vecPoints[ii].getLat(), vecPoints[ii].getLon(), lats[ii], lons[ii]);
        }


        // retrieve nearest points
        size_t n_lats = lats.size();
        if (n_lats != lons.size())
            throw InvalidArgumentException("lats and lons vectors must have the same size. Something went seriously wront", AT);
        // initialize the output vectors
        std::vector<double> outlats_vec(n_lats), outlons_vec(n_lats), altitudes_vec(n_lats), distances_vec(n_lats);
        std::vector<int> indexes_vec(n_lats);
        // find the nearest points in the grid to lats, lons and get the values, as well as the locations and distances
        getNearestValues_grib(h, lats, lons, outlats_vec, outlons_vec, distances_vec, altitudes_vec, indexes_vec);

        // remove potential duplicates
        if (removeDuplicatePoints(vecPoints, outlats_vec, outlons_vec) == true)
            return false;

        // fill metadata
        for (size_t ii = 0; ii < npoints; ii++) {
            StationData sd;
            sd.position.setProj(coordin, coordinparam);
            double true_lat, true_lon;
            CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, outlats_vec[ii], outlons_vec[ii], true_lat, true_lon);
            sd.position.setLatLon(true_lat, true_lon, altitudes_vec[ii]);
            sd.stationID = "Point_" + IOUtils::toString(indexes_vec[ii]);
            std::ostringstream ss2;
            ss2 << "GRIB point (" << indexes_vec[ii] % Ni << "," << indexes_vec[ii] / Ni << ")";
            sd.stationName = ss2.str();
            stations.push_back(sd);
        }
        return true;
    }  
} // namespace
