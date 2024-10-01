// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2024 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
// #include <meteoio/plugins/bufrio.h>
#include <cassert>
#include <meteoio/plugins/BUFRIO.h>
#include <meteoio/plugins/plugin_utils.h>
#include <unordered_set>

using namespace std;

namespace mio {
using namespace PLUGIN;
    /**
     * @page bufrio BUFRIO
     * @section bufrio_format Format
     * BUFR is a binary data format, that was introduced by the <a href="https://wmo.int/">World Meteorological Organization (WMO)</a> to 
     * universally store meteorological observation data. It is very flexible and powerful, but can be quite complex to handle. 
     * You can read <a href="https://en.wikipedia.org/wiki/BUFR">this introduction</a> to get an idea of what it is.
     *
     * In short, a BUFR file can contain multiple messages, each of which can contain multiple subsets. Where the data is stored with 
     * its metadata and can be read by its respective key. Here we use keys from the 
     * <a href="https://confluence.ecmwf.int/display/ECC/WMO%3D40+element+table">WMO element</a> table 40 to read the data. To inspect 
     * your BUFR files we recommend using <a href="https://confluence.ecmwf.int/display/ECC/BUFR+tools">bufr tools</a> or 
     * the <a href="https://metview.readthedocs.io/en/latest/">metview software</a>. These are provided by ECMWF and therefore rely on
     * the same BUFR reading library (<a href="https://confluence.ecmwf.int/display/ECC/ecCodes+Home"> ecCodes</a>) as MeteoIO.
     *
     * Prior to reading the data a scan of the file is done, to find all stations in the file together with their metadata and to check if it is 
     * valid (we require a position (lat,lon,alt) and a station ID). If no station ID but a station Name is found, we use this as station ID. 
     * Additionally we require unique timepoints per station. So far only reading parameters known to MeteoIO is supported (if you need 
     * to read additional parameters, please contact the developers, as this is an option that could be supported).
     * 
     * @note This plugin requires the ecCodes library to be installed. It is available at <a href="https://confluence.ecmwf.int/display/ECC/ecCodes+Home"> ecCodes</a>
     *
     * @subsection write Write
     * It is also possible to write BUFR files. However, only writing known parameters is supported.
     * However, it is possible to write profiles (of heights) of some parameters, i.e. air temperature, pressure and relative humidity .
     * To do that the parameters at different heights need to be named as "TA@2", "TA@10"... for air temperature at 2m and 10m etc.
     * This is the responsibility on the input part.
     * The profile will be written with a repition of the parameter together with the height
     *
     * The filenames are generated on the fly, as meteoio_(timestamp).bufr, where the timestamp is the current date and time.
     * Each datapoint is its own message in the end, as it is not possible to create subsets with eccodes.
     *
     * @section bufrio_keywords Keywords
     * This plugin uses the following keywords:
     * INPUT:
     * - METEOPATH: Path to the BUFR files
     * - BUFREXT: Extension of the BUFR files (default: .bufr)
     * - STATION#: File name per station
     * - ADDITIONAL_PARAMS: Additional parameters to read from the BUFR file, as a comma separated list of WMO element table 40 keys
     * - FALLBACKTZ: Fallback timezone in hours, if no timezone information is found in the BUFR file (default: GMT+0)
     * - VERBOSE: Print additional information during reading (default: false)
     *OUTPUT:
     * - METEOPATH: Path to the BUFR files
     * - SEPARATESTATIONS: Write each station to a separate file (default: false)
     *
     * @section bufrio_example Example
     * @code
     * [INPUT]
     * METEO = BUFR
     * METEOPATH = path/to/bufr/files
     * BUFREXT = .bfr
     * ADDITIONAL_PARAMS = airTemperature, dewPointTemperature
     * STATION1 = example.bfr
     * VERBOSE = TRUE
     *
     * [OUTPUT]
     * METEO = BUFR
     * METEOPATH = path/to/output
     * SEPARATESTATIONS = TRUE
     * @endcode
     *
     * @section bufrio_cryo Cryo Station
     * The Cryo Staion is a specific BUFR Template for Cryospheric Stations containing snow/ice... profiles. It is possible to write those, by setting
     * - CRYOSTATION = TRUE 
     * in the output section. The following additional keywords are required:
     * - WIGOSIDSERIES: WIGOS ID Series (0-15)
     * - WIGOSISSUER: WIGOS Issuer
     * - WIGOSISSUENO: WIGOS Issue Number
     * - WIGOSLOCALID: WIGOS Local ID
     * - STATIONTYPE: Station Type (default: 3)
     * - SURFACETYPE: Surface Type (default: 255)
     * The station type and surface type are WMO codes, and can be found at <a href="https://confluence.ecmwf.int/display/ECC/WMO%3D40+code-flag+table#WMO=40codeflagtable-CF_008018">Flag Tables</a> 
     * 
     * @subsection allowed_params Allowed Parameters
     * The Cryo Station follows a specific template, where some parameters are not MeteoIO standard parameters. Additionally, 
     * some MeteoIO parameters are not supported by the Cryo Station (e.g. VW, ISWR...). Therefore we list the parameters that can be written, with their 
     * representation in MeteoIO (which should be the parameter name in the input file).
     * 
     * Single Parameters:
     * - Ice Thickness: ICE_THICKNESS
     * - Ground State: GROUNDSTATE (https://vocabulary-manager.eumetsat.int/vocabularies/BUFR/WMO/40/TABLE_CODE_FLAG/020062)
     * - Snow Depth: HS (MeteoIO standard)
     * - Surface Qualifier for Temperature Data: SURFACEQUALIFIER (What type of Surface: https://vocabulary-manager.eumetsat.int/vocabularies/BUFR/WMO/40/TABLE_CODE_FLAG/008010) 
     * - Skin Temperature: TSURFACE
     * Parameters with Height:
     * - Sensor Type: SENSORTYPE (https://confluence.ecmwf.int/display/ECC/WMO%3D40+code-flag+table#WMO=40codeflagtable-CF_002096)
     * - Air Temperature: TA (MeteoIO standard)
     * - Snow Temperature: TSNOW
     * - Soil Temperature: TSOIL (MeteoIO standard)
     * - Water Temperature: TWATER
     * - Ice Temperature: TICE
     * 
     * @subsection bufrio_cryo_surface Surface Type
     * Available surface types are:
     * 
     * - 0: OPEN OCEAN OR SEMI-ENCLOSED SEA
     * - 1: ENCLOSED SEA OR LAKE
     * - 2: CONTINENTAL ICE
     * - 3: LAND
     * - 4: LOW INLAND (BELOW SEA LEVEL)
     * - 5: MIX OF LAND AND WATER
     * - 6: MIX OF LAND AND LOW INLAND
     * - 11: RIVER
     * - 12: LAKE
     * - 13: SEA
     * - 14: GLACIER
     * - 15: URBAN LAND
     * - 16: RURAL LAND
     * - 17: SUBURBAN LAND
     * - 18: SEA ICE
     * - 255: MISSING VALUE
     * 
     * @subsection bufrio_cryo_station Station Type
     * Available station types are:
     * - 0: Automatic
     * - 1: Manned
     * - 2: Hybrid
     * - 3: Missing Value
     * 
     * @subsection bufrio_cryo_snow Snow Depth Method
     * Available snow depth methods can be found at https://confluence.ecmwf.int/display/ECC/WMO%3D40+code-flag+table#WMO=40codeflagtable-CF_002177
     * 
     * @note the given flags in these sub section might be outdated, as this is not regularly updated
     * 
     *
     *
     * @section dev_notes Developer Notes
     * If additional parameters need to be written to BUFR:
     * - Add the parameter to the constant maps in libcodes.cc
     * - Make Sure to handle specifics of parametrs in setMeteoData or setCryoData
     * - Add the parameter to the POSSIBLE_MULTIPLE_PARAMETERS vector if needed
     * The descriptor can be complete BUFR sequences, however, as soon as a parameter appears twice in two different sequences,
     * it will silentyly fail, and only the first descriptor is used. So be careful that each parameter only appears once.
     * - We might need the flexibility to add parameters outside of the cryo template, a way to flexibly add standard parameters
     *  @todo  The Cryo parameters which will be added to the WMO table need to be changed when available ( look for TODO in libcodes.cc)
     * 
     * For now Wigos id... metadata has to be set by ini, maybe make it available as StationData
     * Additionally, some parameters from the Cryo Station are not standard parameters, so they wont be handled-> needs a solution
     * 
     * 
     * @note If writing arbitrary parameters is required, we need descriptors and parameter keys, which are not hardcoded
     * @todo map precip type from code to...
     * 
     *
     */
    const double BUFRIO::plugin_nodata = IOUtils::nodata; // plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
    const std::string dflt_extension_BUFR = ".bufr";
    const std::string BUFRIO::template_filename = "MeteoIO.bufr";
    const std::set<std::string> ADDITIONAL_CRYO_PARAMS({"SURFACEQUALIFIER", "SENSORTYPE", "TSURFACE", "ICE_THICKNESS", "GROUNDSTATE", "SENSORTYPE","TICE", "TWATER"});
    const std::set<std::string> ALLOWED_CRYO_PARAMS({"ICE_THICKNESS", "GROUNDSTATE", "HS", "SURFACEQUALIFIER", "TSURFACE", "SENSORTYPE", "TA", "TSNOW", "TSOIL", "TWATER", "TICE"});
    const std::vector<MeteoParam> POSSIBLE_MULTIPLE_PARAMETERS = {
        MeteoParam::P, MeteoParam::TA, MeteoParam::RH, MeteoParam::TSOIL, MeteoParam::TSNOW}; // DEV note: the order of setting parameters in the BUFR file must follow the order here, otherwise heights get all mixed up

    BUFRIO::BUFRIO(const std::string &configfile)
        : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), station_files(), additional_params(), outpath(), separate_stations(false), verbose_out(false), write_cryo(false),
          wigos_id_series(IOUtils::nodata), wigos_issuer(IOUtils::nodata), wigos_issue_no(IOUtils::nodata), station_type(3), surface_type(255), snow_depth_method(15), wigos_local_id() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

        parseInputSection();
        parseOutputSection();
    }

    BUFRIO::BUFRIO(const Config &cfgreader)
        : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), station_files(), additional_params(), outpath(), separate_stations(false), verbose_out(false), write_cryo(false),
          wigos_id_series(IOUtils::nodata), wigos_issuer(IOUtils::nodata), wigos_issue_no(IOUtils::nodata), station_type(3), surface_type(255), snow_depth_method(15),   wigos_local_id() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
        parseInputSection();
        parseOutputSection();
    }

    void BUFRIO::parseInputSection() {
        const std::string in_meteo = IOUtils::strToUpper(cfg.get("METEO", "Input", ""));
        if (in_meteo == "BUFR") { // keep it synchronized with IOHandler.cc for plugin mapping!!
            const std::string inpath = cfg.get("METEOPATH", "Input");
            std::vector<std::string> vecFilenames;
            cfg.getValues("STATION", "INPUT", vecFilenames);

            std::string bufr_ext = dflt_extension_BUFR;
            cfg.getValue("BUFREXT", "Input", bufr_ext, IOUtils::nothrow);
            if (bufr_ext == "none")
                bufr_ext.clear();

            cfg.getValue("ADDITIONAL_PARAMS", "INPUT", additional_params, IOUtils::nothrow);

            double fallback_tz = 0;
            cfg.getValue("FALLBACKTZ", "Input", fallback_tz, IOUtils::nothrow);

            bool verbose = false;
            cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);

            if (vecFilenames.empty())
                scanMeteoPath(cfg, inpath, vecFilenames, bufr_ext);

            const std::vector<std::string> all_files_and_paths = getFilesWithPaths(vecFilenames, inpath, dflt_extension_BUFR);
            for (const auto &filename : all_files_and_paths) {
                station_files.push_back(BUFRFile(filename, coordin, verbose, fallback_tz));
            }
        }
    }

    void BUFRIO::parseOutputSection() {
        const std::string out_meteo = IOUtils::strToUpper(cfg.get("METEO", "Output", ""));
        if (out_meteo == "BUFR") {
            outpath = cfg.get("METEOPATH", "Output", "");
            if (outpath.empty()) {
                throw IOException("No output path specified for BUFR output", AT);
            }
            separate_stations = cfg.get("SEPARATESTATIONS", "Output", false);
            verbose_out = cfg.get("VERBOSE", "Output", false);
            write_cryo = cfg.get("CRYOSTATION", "Output", false);
            if (write_cryo) {
                cfg.getValue("WIGOSIDSERIES", "Output", wigos_id_series); // we might need to have this as metadata
                cfg.getValue("WIGOSISSUER", "Output", wigos_issuer);      // we might need to have this as metadata
                cfg.getValue("WIGOSISSUENO", "Output", wigos_issue_no);   // we might need to have this as metadata
                cfg.getValue("WIGOSLOCALID", "Output", wigos_local_id);   // we might need to have this as metadata
                station_type = static_cast<long>(cfg.get("STATIONTYPE", "Output", 3.));      // missing val as default // we might need to have this as metadata
                surface_type = static_cast<long>(cfg.get("SURFACETYPE", "Output", 255.));     // missing val as default  // we might need to have this as metadata
                snow_depth_method = static_cast<long>(cfg.get("SNOWDEPTHMETHOD", "Output", 15.)); // missing val as default // we might need to have this as metadata
            }
        }
    }

    // ------------------------- READ -------------------------

    void BUFRIO::readMeteoData(const Date & /*dateStart*/, const Date & /*dateEnd*/, std::vector<std::vector<MeteoData>> &vecvecMeteo) {
        vecvecMeteo.clear();
        std::map<std::string, size_t> station_ids;
        for (auto &station_file : station_files) {
            station_file.readData(vecvecMeteo, station_ids, additional_params);
        }
        // sort the data, as multiple bufr files can have information about the same station, but order is not guaranteed
        for (auto &vecMeteo : vecvecMeteo) {
            std::sort(vecMeteo.begin(), vecMeteo.end());
        }
    }

    void BUFRIO::readStationData(const Date & /* date */, std::vector<StationData> &vecStation) {
        vecStation.clear();
        std::unordered_set<std::string> station_id_set;
        for (const auto &station_file : station_files) {
            for (const auto &station : station_file.getMetadata()) {
                auto success = station_id_set.insert(station.second.getStationID());
                if (success.second) {
                    vecStation.push_back(station.second);
                }
            }
        }
    }

    // ------------------------- WRITE -------------------------

    // STATIC WRITE HELPERS
    static bool isNumber(const std::string &s) {
        try {
            std::stod(s);
            return true;
        } catch (const std::invalid_argument &) {
            return false;
        }
    }

    static void setStationId(CodesHandlePtr &message, const StationData &station, const std::string &subset_prefix) {
        std::string station_id = station.getStationID();
        bool set_id = false;
        if (isNumber(station_id)) {
            set_id = setParameter(message, subset_prefix + "stationNumber", station_id);
        } else {
            set_id = setParameter(message, subset_prefix + "shortStationName", station_id);
        }
        if (!set_id) {
            throw IOException("Station ID could not be set", AT);
        }
    };

    void BUFRIO::setWIGOSId(CodesHandlePtr &message, const std::string &subset_prefix) {
        if (!setParameter(message, subset_prefix + "wigosIdentifierSeries", wigos_id_series))
            throw IOException("WIGOS ID Series could not be set", AT);
        if (!setParameter(message, subset_prefix + "wigosIssuerOfIdentifier", wigos_issuer))
            throw IOException("WIGOS Issuer could not be set", AT);
        if (!setParameter(message, subset_prefix + "wigosIssueNumber", wigos_issue_no))
            throw IOException("WIGOS Issue Number could not be set", AT);
        if (!setParameter(message, subset_prefix + "wigosLocalIdentifierCharacter", wigos_local_id))
            throw IOException("WIGOS Local ID could not be set", AT);
    }

    void BUFRIO::setStationData(CodesHandlePtr &message, const StationData &station, const Coords &position, const std::string &subset_prefix) {
        if (write_cryo) {
            setWIGOSId(message, subset_prefix);
            if (!setParameter(message, subset_prefix + "stationType", station_type))
                throw IOException("Station Type could not be set", AT);
            if (!setParameter(message, subset_prefix + "surfaceType", surface_type))
                throw IOException("Surface Type could not be set", AT);
            if (!setParameter(message, subset_prefix + "methodOfSnowDepthMeasurement", snow_depth_method))
                throw IOException("Snow Depth Method could not be set", AT);
        } else {
            setStationId(message, station, subset_prefix);
        }

        std::string station_name_key = write_cryo ? "longStationName" : "stationOrSiteName";
        if (!setParameter(message, subset_prefix + station_name_key, station.getStationName()))
            throw IOException("Station Name could not be set", AT);

        if (!setParameter(message, subset_prefix + "latitude", position.getLat()))
            throw IOException("latitude could not be set", AT);
        if (!setParameter(message, subset_prefix + "longitude", position.getLon()))
            throw IOException("longitude could not be set", AT);
        if (!setParameter(message, subset_prefix + "heightOfStationGroundAboveMeanSeaLevel", position.getAltitude()))
            throw IOException("station height could not be set", AT);
    }

    // setMeteoData and setCryoData essentially do the same, however the order in which setting height and parameters is done is different
    // Therefore the functions are very simliar but couldnt be refactored

    static void setMeteoData(CodesHandlePtr &message, const MeteoData &meteo, const std::vector<std::string> &sorted_params, std::map<std::string, int> &parameter_occurences, const bool &verbose) {
        for (const auto &param : sorted_params) {
            if (param == "RSWR")
                continue; // skip RSWR, as it is not in the BUFR template

            size_t param_id = meteo.getParameterIndex(param);

            std::string parsed_param_name;
            double parsed_param_number;

            bool success = false;
            // this is a preliminary way to write a type of profiles
            if (MeteoData::getTypeAndNo(param, parsed_param_name, parsed_param_number)) {
                std::string prefix = "#" + std::to_string(parameter_occurences[parsed_param_name]) + "#";

                // handle all the rest
                if (parsed_param_number == IOUtils::nodata) { // is a standard parameter, so no need to worry about height
                    success = setParameter(message, prefix + BUFR_PARAMETER.at(parsed_param_name), meteo(param_id));
                    if (parsed_param_name == "TSS") // set the surface qualifier, as it is a special case
                        success = setParameter(message, prefix + BUFR_PARAMETER.at("SURFACEQUALIFIER"), SNOW_SURFACE_QUALIFIER); 
                } else { // set the height together with the parameter
                    if (std::find(POSSIBLE_MULTIPLE_PARAMETERS.begin(), POSSIBLE_MULTIPLE_PARAMETERS.end(), MeteoData::toParameter(parsed_param_name)) == POSSIBLE_MULTIPLE_PARAMETERS.end()) {
                        if (verbose)
                            std::cout << "Parameter " << param << " is not suppoprted for profiles, and will not be written to BUFR" << std::endl;
                        continue; // we cant handle those height dependent parameters yet                        
                    }
                    // parameter
                    std::string param_key = parsed_param_name == "TA" ? BUFR_PARAMETER_ALT.at(parsed_param_name) : BUFR_PARAMETER.at(parsed_param_name);
                    success = setParameter(message, prefix + param_key, meteo(param_id));

                    // height
                    std::string height_prefix = "#" + std::to_string(parameter_occurences[BUFR_HEIGHT_KEY]) + "#";
                    success = setParameter(message, height_prefix + BUFR_HEIGHT_KEY, parsed_param_number);
                    parameter_occurences[BUFR_HEIGHT_KEY] += 1; // increment occurences
                }
                parameter_occurences[parsed_param_name] += 1; // increment occurences

                if (!success) // only throw error if it is a known parameter
                    throw IOException("Parameter " + param + " could not be set", AT);
            } else {
                if (verbose)
                    std::cout << "Parameter " << param << " is not a known parameter, and will not be written to BUFR" << std::endl;
            }
        }
    }

    static void setCryoData(CodesHandlePtr &message, const MeteoData &meteo, const std::vector<std::pair<double, std::vector<std::string>>> &params_at_heights,
                            std::map<std::string, int> &parameter_occurences, const bool &/* verbose */) {
        // we need to set all parameters for one height that we have, and then move to the next height
        for (const auto &height_and_params : params_at_heights) {
            const double height = height_and_params.first;
            const std::vector<std::string> &params = height_and_params.second;

            bool success = false;
            // we have every parameter in each height layer, even when they are missing, so we need this prefix for all repeated params
            std::string height_prefix = "#" + std::to_string(parameter_occurences[BUFR_HEIGHT_KEY]) + "#";

            for (const auto &param : params) {
                size_t param_id = meteo.getParameterIndex(param);

                std::string parsed_param_name;
                double parsed_param_number;

                if (MeteoData::getTypeAndNo(param, parsed_param_name, parsed_param_number, ADDITIONAL_CRYO_PARAMS)) {

                    if (parsed_param_number == IOUtils::nodata) { // is a standard parameter, so no need to worry about height
                        if (ALLOWED_CRYO_PARAMS.find(parsed_param_name) == ALLOWED_CRYO_PARAMS.end()) {
                            throw IOException("Parameter " + param + " is not allowed for Cryo Station", AT);                    
                        }
                        std::string prefix = "#" + std::to_string(parameter_occurences[parsed_param_name]) + "#";
                        success = setParameter(message, prefix + BUFR_PARAMETER.at(parsed_param_name), meteo(param_id));
                    } else { // set the height together with the parameter

                        assert(parsed_param_number == height); // height should be the same as the one we are setting

                        std::string param_key =
                            parsed_param_name == "TA" ? BUFR_PARAMETER_ALT.at(parsed_param_name) : BUFR_PARAMETER.at(parsed_param_name); // TA standard is has a different key then when settin repeated
                        success = setParameter(message, height_prefix + param_key, meteo(param_id));
                    }
                    parameter_occurences[parsed_param_name] += 1; // increment occurences

                    if (!success) // only throw error if it is a known parameter
                        throw IOException("Parameter " + param + " could not be set", AT);
                } else {
                    throw IOException("Parameter " + param + " is not supported for CRYO STATIONS", AT);
                }
            }

            if (height != IOUtils::nodata) { // set height and increment at each level, except nodata(==single occurences)
                success = setParameter(message, height_prefix + BUFR_HEIGHT_KEY, height);
                parameter_occurences[BUFR_HEIGHT_KEY] += 1;
            }
        }
    }

    static std::set<double> getHeights(const std::set<std::string> &available_params) {
        std::set<double> heights = MeteoData::retrieveAllHeights(available_params, ADDITIONAL_CRYO_PARAMS); // params without heights, have IOUtils::nodata as height
        if (heights.empty()) {
            std::cout << "No profile parameters found" << std::endl;
        }
        return heights;
    }

    static std::map<MeteoParam, size_t> getOccurencesOfParameters(const std::vector<MeteoData> &vecMeteo) {
        std::map<MeteoParam, size_t> occurences;
        for (const auto &param : POSSIBLE_MULTIPLE_PARAMETERS) {
            occurences[param] = vecMeteo.front().getOccurencesOfParameter(param, ADDITIONAL_CRYO_PARAMS); // all meteodata will have the same number of parameters
        }
        return occurences;
    }

    static std::map<std::string, int> initParameterOccurences(const std::set<std::string> &available_params) {
        std::map<std::string, int> parameter_occurences;
        for (const std::string &param : available_params) {
            parameter_occurences[param] = 1;
        }
        parameter_occurences[BUFR_HEIGHT_KEY] = 1;
        return parameter_occurences;
    }

    // we need the method overload to handle setting cryo station or meteoio template
    static void populateSubset(CodesHandlePtr &message, const std::string &subset_prefix, const MeteoData &meteo, const std::vector<std::string> &sorted_params,
                               std::map<std::string, int> &parameter_occurences, const bool &verbose) {
        if (sorted_params.empty()) {
            throw IOException("No parameters to write to BUFR", AT);
        }
        setTime(message, meteo.date, subset_prefix);
        setMeteoData(message, meteo, sorted_params, parameter_occurences, verbose);
    }
    static void populateSubset(CodesHandlePtr &message, const std::string &subset_prefix, const MeteoData &meteo, const std::vector<std::pair<double, std::vector<std::string>>> &params_at_heights,
                               std::map<std::string, int> &parameter_occurences, const bool &verbose) {
        if (params_at_heights.empty()) {
            throw IOException("No parameters to write to BUFR", AT);
        }
        setTime(message, meteo.date, subset_prefix);
        setCryoData(message, meteo, params_at_heights, parameter_occurences, verbose);
    }

    static std::vector<std::string> sortAvailableParameters(const std::set<std::string> &available_params, const bool &verbose_out) {
        std::set<std::string> aggregated_params;
        std::vector<std::string> sorted_params;
        // We need to enforce this order, of paramters occuring multiple times, otherwise tracking the height bufr paramter is a nightmare
        // Sort, so that repeated parameters occur right after each other
        for (const auto &param_type : POSSIBLE_MULTIPLE_PARAMETERS) {
            if (!aggregated_params.insert(MeteoData::parToString(param_type)).second)
                throw IOException("Parameter " + MeteoData::parToString(param_type) + " already in aggregated_params. Probably same Multiple Possible Parameters", AT);
            std::string param_name(MeteoData::parToString(param_type));
            std::vector<std::string> all_params_of_type_sorted(MeteoData::retrieveAllHeightsForParam(available_params, param_name));
            sorted_params.insert(sorted_params.end(), all_params_of_type_sorted.begin(), all_params_of_type_sorted.end());
        }
        for (const auto &param : available_params) {
            std::string parsed_param_name;
            double parsed_param_number;
            MeteoData::getTypeAndNo(param, parsed_param_name, parsed_param_number);
            if (aggregated_params.insert(parsed_param_name).second) {
                sorted_params.push_back(param); // these will not occur multiple times, and if they do, it will not work
            }
        }
        if (verbose_out) {
            std::cout << "Sorted parameters: " << std::endl;
            for (const auto &param : sorted_params) {
                std::cout << param << std::endl;
            }
        }

        return sorted_params;
    }

    // Sorting for Cryo Station
    static std::vector<std::pair<double, std::vector<std::string>>> sortAvailableParametersByHeight(const std::set<std::string> &available_params, const bool &/* verbose_out */, const std::set<double>& heights) {
        std::vector<std::pair<double, std::vector<std::string>>> params_at_heights;
        for (const auto &height : heights) {
            std::vector<std::string> params_at_h(MeteoData::retrieveAllParametersAtHeight(available_params, height, ADDITIONAL_CRYO_PARAMS));
            params_at_heights.push_back(std::make_pair(height, params_at_h));
        }
        return params_at_heights;
    }

    // WRITE METEO

    void BUFRIO::writeMeteoData(const std::vector<std::vector<MeteoData>> &vecStations, const std::string & /* name */) {
        const std::string file_prefix = write_cryo ? "/meteoio_cryo_" : "/meteoio_";
        const std::string extension = FileUtils::getDateTime() + ".bufr";

        // create and write a message for each station
        for (const auto &vecMeteo : vecStations) {
            // get header and other constant information
            const StationData station = vecMeteo.front().meta;
            const std::set<std::string> available_params(MeteoData::listAvailableParameters(vecMeteo));
            
            if (vecMeteo.front().listUnknownParameters(ADDITIONAL_CRYO_PARAMS) > 0) {
                std::cout << "Warning: Unknown parameters in station " << station.getStationID() << " will not be written to BUFR" << std::endl;
            }

            // we need different sorting of parameters for cryo station, and meteoio template
            std::vector<std::string> sorted_params;
            std::vector<std::pair<double, std::vector<std::string>>> params_at_heights;
            std::set<double> heights;
            long num_heights = 0;
            if (!write_cryo) {
                sorted_params =
                    sortAvailableParameters(available_params, verbose_out); // need to set the parameters sorted, according to POSSIBLE_MULIT... otherwise tracking the height occurences is a nightmare
            } else {
                heights = getHeights(available_params);
                num_heights = heights.find(IOUtils::nodata) == heights.end() ? heights.size() : heights.size() - 1; // if we have a nodata height, we dont count it
                params_at_heights = sortAvailableParametersByHeight(available_params, verbose_out, heights);
            }

            // the number of repitions of a parameter
            const std::map<MeteoParam, size_t> total_occurences(getOccurencesOfParameters(vecMeteo));

            CodesHandlePtr message = createBUFRMessageFromSample(vecMeteo.size(), total_occurences, available_params, POSSIBLE_MULTIPLE_PARAMETERS, write_cryo, num_heights);

            setMissingValue(message, plugin_nodata);
            // typical time is the starting date
            setTypicalTime(message, vecMeteo.front().date);

            // map all the occurences of parameters, so the right one is set (incluing height above ground, as it is used often)
            std::map<std::string, int> parameter_occurences_in_file(initParameterOccurences(MeteoData::retriveUniqueParameters(available_params, ADDITIONAL_CRYO_PARAMS)));

            const Coords position = station.getPosition();
            // set subset specific information, i.e. each time point
            for (size_t i = 1; i <= vecMeteo.size(); ++i) {
                const MeteoData meteo = vecMeteo[i - 1];
                const std::string subset_prefix = "#" + std::to_string(i) + "#"; // we use occurence number as prefix
                setStationData(message, station, position, subset_prefix);       // is the same for each subset, but needs to be set, otherwise missing
                if (write_cryo) {
                    populateSubset(message, subset_prefix, meteo, params_at_heights, parameter_occurences_in_file, verbose_out); // overload for cryo station
                } else {
                    populateSubset(message, subset_prefix, meteo, sorted_params, parameter_occurences_in_file, verbose_out);
                }
            }
            // pack and write message, is always appended if the file already exists
            packMessage(message);

            const std::string outfile = outpath + file_prefix + (separate_stations ? station.getStationID() : "") + extension;
            writeToFile(message, outfile);
        }
    }
} // namespace
