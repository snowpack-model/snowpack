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
#include <fstream>
#include <iostream>
#include <meteoio/FStream.h>
#include <meteoio/plugins/plugin_utils.h>
#include <meteoio/plugins/iCSVIO.h>
#include <regex>

// using namespace std;

namespace mio {
/**
* @page icsvio iCSV
* @section template_format Format
*
* @brief The <b>interoperable CSV</b> format (iCSV) is a versatile and intuitive format that merges the self-documenting
* capabilities of NetCDF with the human-friendly readability and writeability of CSV.
*
* @details The iCSV format has been primarily designed for the exchange and preservation of time series data in environmental
* data repositories but can also support a much broader use as very few metadata are mandatory (but many are available). Being
* a text format, it should remain easily readable in decades to come, thus supporting long term data preservation.
* For comprehensive details about the format, refer to its
* <a href="http://envidat.gitlab-pages.wsl.ch/icsv/">official format documentation</a>.
*
* When using iCSV with MeteoIO, the following additional MetaData is required:
* - Timestamp or Julian field: This represents the timestamps of the measurement in ISO format.
* - Location: This should be specified either in the header or in the appropriate column as a location for a station. Only POINT(x y)
* or
* POINTZ(x y z) formats are supported.
* - Timezone: This should be specified in the METADATA section, or as part of the ISO timestamp string.
*
* Please note that MeteoIO currently does not support handling additional information on the fields, such as long_name. Such
* information will be ignored.
*
* @section template_units Units
* All units are <a href="https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.330-2019.pdf">coherent derived SI units</a>
* (section 2.3.4 in the SI-Brochure). It is however
* possible to use  multipliers and offsets (but they must be specified in the FIELDS section).
*
* @section template_keywords Keywords
* This plugin uses the following keywords:
* - METEOPATH: meteo files directory where to read/write the meteofiles; [Input] and [Output] sections
* - METEOFILE#: input filename (in METEOPATH). As many meteofiles as needed may be specified. If nothing is specified, the METEOPATH
* directory
* will be scanned for files ending in ".icsv" and sorted in ascending order;
* - METEOPATH_RECURSIVE: if set to true, the scanning of METEOPATH is performed recursively (default: false); [Input] section;
* - SNOWPACK_SLOPES: if set to true and no slope information is found in the input files,
* the <a
* href="https://www.slf.ch/en/avalanche-bulletin-and-snow-situation/measured-values/description-of-automated-stations.html">IMIS/Snowpack</a>
* naming scheme will be used to derive the slope information (default: false, [Input] section).     *
* - iCSV_APPEND: when an output file already exists, should the plugin try to append data (default: false); [Output] section
* - iCSV_OVERWRITE: when an output file already exists, should the plugin overwrite it (default: true)? [Output] section
* - ACDD_WRITE: add the Attribute Conventions Dataset Discovery <A href="http://wiki.esipfed.org/index.php?title=Category:Attribute_Conventions_Dataset_Discovery">(ACDD)</A> 
* metadata to the headers (then the individual keys are provided according to the ACDD class documentation) (default: false, [Output] section)
* - iCSV_SEPARATOR: choice of field delimiter, options are: [,;:|/\]; [Output] section
* 
* @note There is a python module available to read iCSV files, see <a href="http://patrick.leibersperger.gitlab-pages.wsl.ch/snowpat/">snowpat</a>
* 
*/
using namespace PLUGIN;
// clang-format off
static const std::string dflt_extension_iCSV = ".icsv";
const double iCSVIO::snVirtualSlopeAngle = 38.; //in Snowpack, virtual slopes are 38 degrees
static const std::streamsize MAXMEMORY = static_cast<std::streamsize>(20) * 1024 * 1024 * 1024; // if i use more than 20GB of memory for data read it sequentially
// clang-format on

// ------------------ IO Helpers ------------------
/**
* Checks if the total size of the files specified by the given filenames
* is within the memory limit.
*
* @param filenames The vector of filenames to check.
* @return True if the total size is within the memory limit, false otherwise.
*/
static bool areFilesWithinLimit(const std::vector<std::string> &filenames) {
    std::streamsize totalSize = 0;
    for (const auto &filename : filenames) {
        std::ifstream file(filename, std::ios::binary | std::ios::ate);

        const std::streamsize size = file.tellg();
        file.close();

        totalSize += size;
        if (totalSize > MAXMEMORY) {
            return false;
        }
    }
    return true;
}

static std::string joinVector(const std::vector<std::string> &vec, const char &delimiter) {
    std::string result;
    for (size_t i = 0; i < vec.size(); i++) {
        result += vec[i];
        if (i != vec.size() - 1) {
            result += delimiter;
        }
    }
    return result;
}

static std::vector<Coords> getUniqueLocations(iCSVFile &file) {
    std::vector<Coords> locations;
    if (file.location_in_header) {
        locations.push_back(file.station_location.toCoords(file.METADATA.epsg));
    } else {
        for (size_t ii = 0; ii < file.getAllLocationsInData().size(); ii++) {
            const Coords latest = file.getLocationAt(ii).toCoords(file.METADATA.epsg);
            if (latest != locations.back()) {
                locations.push_back(latest);
            }
        }
    }
    return locations;
}

using namespace iCSV;

// ----------------- iCSVIO -----------------
iCSVIO::iCSVIO(const std::string &configfile)
    : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), snowpack_slopes(false), read_sequential(false),
        stations_files(), acdd_metadata(false), TZ_out(0), outpath(""), allow_overwrite(false), allow_append(false), out_delimiter(','), file_extension_out(dflt_extension_iCSV) {
    parseInputSection();
    parseOutputSection();
}

iCSVIO::iCSVIO(const Config &cfgreader)
    : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), snowpack_slopes(false), read_sequential(false),
        stations_files(), acdd_metadata(false), TZ_out(0), outpath(""), allow_overwrite(false), allow_append(false), out_delimiter(','), file_extension_out(dflt_extension_iCSV) {
    parseInputSection();
    parseOutputSection();
}

// ----------------- iCSVIO parse -----------------
void iCSVIO::parseInputSection() {
    // default timezones

    // Parse the [Input] and [Output] sections within Config object cfg
    IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

    // Parse input section: extract number of files to read and store filenames in vecFiles
    const std::string in_meteo = IOUtils::strToUpper(cfg.get("METEO", "Input", ""));
    if (in_meteo == "ICSV") { // keep it synchronized with IOHandler.cc for plugin mapping!!
        cfg.getValue("SNOWPACK_SLOPES", "Input", snowpack_slopes, IOUtils::nothrow);
        const std::string inpath = cfg.get("METEOPATH", "Input");

        //handle the deprecated STATION# syntax
        std::string meteofiles_key( "METEOFILE" );
        //const std::vector< std::pair<std::string, std::string> > vecDeprecated( cfg.getValues("STATION", "INPUT") );
       const std::vector< std::string > vecDeprecated( cfg.getKeys("STATION", "INPUT") );
		if (!vecDeprecated.empty()) {
			meteofiles_key = "STATION";
			std::cerr << "[W] The STATION# syntax has been deprecated for the iCSV input plugin, please rename these keys as METEOFILE#!\n";
			//throw InvalidArgumentException("The STATION# syntax has been deprecated for the SMET plugin, please rename these keys as METEOFILE#!", AT);
		}
        
        std::vector<std::string> vecFilenames;
        cfg.getValues(meteofiles_key, "INPUT", vecFilenames);

        if (vecFilenames.empty())
            scanMeteoPath(cfg, inpath, vecFilenames, dflt_extension_iCSV);

        const std::vector<std::string> all_files_and_paths = getFilesWithPaths(vecFilenames, inpath, dflt_extension_iCSV);
        read_sequential = !areFilesWithinLimit(all_files_and_paths);
        for (auto &filename : all_files_and_paths)
            stations_files.push_back(iCSVFile(filename, read_sequential));
    } 
}

void iCSVIO::parseOutputSection() {
    // Parse output section: extract info on whether to write ASCII or BINARY, gzipped or not, acdd...
    cfg.getValue("TIME_ZONE", "Output", TZ_out, IOUtils::nothrow);

    outpath.clear();
    const std::string out_meteo = IOUtils::strToUpper(cfg.get("METEO", "Output", ""));
    if (out_meteo == "ICSV") { // keep it synchronized with IOHandler.cc for plugin mapping!!

        cfg.getValue("EXTENSION_OUT", "Output", file_extension_out, IOUtils::nothrow);
        static const std::regex valid_extension("^[.][a-z0-9]+$", std::regex::icase | std::regex::optimize);
        if (!std::regex_match(file_extension_out, valid_extension))
            throw InvalidFormatException("Invalid extension format (valid: \".abc123\") : " + file_extension_out, AT);

        bool write_acdd = false;
        cfg.getValue("ACDD_WRITE", "Output", write_acdd, IOUtils::nothrow);
        if (write_acdd) {
            acdd_metadata.setEnabled(true);
            acdd_metadata.setUserConfig(cfg, "Output", false); // do not allow multi-line keys
        }

        cfg.getValue("METEOPATH", "Output", outpath, IOUtils::nothrow);
        cfg.getValue("iCSV_APPEND", "Output", allow_append, IOUtils::nothrow);
        cfg.getValue("iCSV_OVERWRITE", "Output", allow_overwrite, IOUtils::nothrow);
        cfg.getValue("iCSV_SEPARATOR", "Output", out_delimiter,
                        IOUtils::nothrow); // allow specifying a different field separator as required by some import programs
        std::vector<std::string> vecArgs;

        if (allow_overwrite && allow_append)
            throw InvalidFormatException("Cannot allow both iCSV_APPEND and iCSV_OVERWRITE", AT);
    } 
}

// ------------------------------------- iCSVIO read -----------------------------------------

void iCSVIO::readStationData(const Date & /*date*/, std::vector<StationData> &vecStation) {
    vecStation.clear();
    vecStation.resize(stations_files.size());

    // Now loop through all requested stations, open the respective files and parse them
    for (size_t ii = 0; ii < stations_files.size(); ii++) {
        StationData sd;
        iCSVFile &current_file = stations_files[ii];

        read_meta_data(current_file, sd);
        vecStation[ii] = sd;
    }
}

void iCSVIO::readMeteoData(const Date &start_date, const Date &end_date, std::vector<std::vector<MeteoData>> &vecvecMeteo) {
    vecvecMeteo.clear();
    vecvecMeteo.resize(stations_files.size());

    for (size_t ii = 0; ii < stations_files.size(); ii++) {
        iCSVFile current_file = stations_files[ii];
        if (read_sequential) {
            readDataSequential(current_file);
        }

        std::vector<Date> date_vec = current_file.getDatesInFile(start_date, end_date);
        std::vector<geoLocation> location_vec = current_file.getLocationsInData(start_date, end_date);

        std::vector<MeteoData> vecMeteo = createMeteoDataVector(current_file, date_vec, location_vec);
        vecvecMeteo.push_back(vecMeteo);
    }
}

// --------------------------- iCSVIO read helper functions ---------------------------------------

/**
* @brief Reads the meta data from the iCSV file and populates the StationData object.
*
* This function reads the meta data from the iCSV file specified by current_file and populates
* the StationData object specified by meta. It sets the meta data position, station ID, station name,
* meta data slope, and extra information. If the location is not specified in the header of the iCSV file,
* it sets the EPSG code for the position. The nodata value is used to handle missing data.
*
* @param current_file The iCSV file from which to read the meta data.
* @param meta The StationData object to populate with the meta data.
*/
void iCSVIO::read_meta_data(const iCSVFile &current_file, StationData &meta) {
    const double nodata_value = current_file.getNoData();

    setMetaDataPosition(current_file, meta, nodata_value);

    meta.stationID = current_file.METADATA.station_id;
    meta.stationName = current_file.METADATA.getOptionalMetaData("station_name");

    setMetaDataSlope(current_file, meta, nodata_value);

    meta.extra = current_file.METADATA.toMetaMap();
}

void iCSVIO::setMetaDataPosition(const iCSVFile &current_file, StationData &meta, const double &nodata_value) {
    meta.position.setProj(coordin, coordinparam);
    if (current_file.location_in_header) {
        const double east = IOUtils::standardizeNodata(current_file.station_location.x, nodata_value);
        const double north = IOUtils::standardizeNodata(current_file.station_location.y, nodata_value);
        const double alt = IOUtils::standardizeNodata(current_file.station_location.z, nodata_value);
        meta.position.setPoint(east, north, alt, current_file.METADATA.epsg);
    }
    if (!meta.position.isNodata())
        meta.position.check("Inconsistent geographic coordinates in file \"" + current_file.filename + "\": ");
}

void iCSVIO::setMetaDataSlope(const iCSVFile &current_file, StationData &meta, const double &nodata_value) {
    const double slope_angle = IOUtils::standardizeNodata(current_file.station_location.slope_angle, nodata_value);
    const double slope_azi = IOUtils::standardizeNodata(current_file.station_location.slope_azi, nodata_value);
    if (slope_angle != IOUtils::nodata && slope_azi != IOUtils::nodata) {
        meta.setSlope(slope_angle, slope_azi);
    } else if (slope_angle == 0.) {
        meta.setSlope(slope_angle, 0.);
    } else if (snowpack_slopes) {
        const double exposition = getSnowpackSlope(meta.stationID);
        if (exposition == 0.)
            meta.setSlope(0., 0.);
        else if (exposition != IOUtils::nodata)
            meta.setSlope(snVirtualSlopeAngle, (exposition - 1.) * 90.);
    }
}

/**
* Reads data sequentially from a file.
*
* @param current_file The iCSVFile object representing the file to read from.
* @throws IOException if the file cannot be opened.
*/
void iCSVIO::readDataSequential(iCSVFile &current_file) {
    std::ifstream file(current_file.filename);
    if (!file.is_open()) {
        throw IOException("Unable to open file " + current_file.filename, AT);
    }

    std::string line;
    size_t line_count = 0;
    while (std::getline(file, line)) {
        if (line_count < current_file.skip_lines_to_data) {
            line_count++;
            continue;
        }
        IOUtils::trim(line);
        if (line.empty()) {
            continue;
        }
        current_file.processData(line);
    }
}

/**
* @brief Creates a vector of MeteoData objects, i.e. a time series.
*
* The MeteoData objects are created based on the
* provided iCSVFile, Date, and geoLocation information.
*
* @param current_file The iCSVFile object containing the necessary data.
* @param date_vec The vector of Dates
* @param location_vec The vector of glocations
* @return std::vector<MeteoData> The vector of MeteoData objects created.
*/
std::vector<MeteoData> iCSVIO::createMeteoDataVector(iCSVFile &current_file, std::vector<Date> &date_vec,
                                                        std::vector<geoLocation> &location_vec) {
    std::vector<MeteoData> vecMeteo;
    vecMeteo.reserve(date_vec.size());

    std::vector<size_t> indexes;
    MeteoData md;
    identify_fields(current_file.FIELDS.fields, indexes, md, current_file.METADATA.geometry);

    double nodata = current_file.getNoData();

    MeteoData tmp_md(md);
    for (size_t d_idx = 0; d_idx < date_vec.size(); d_idx++) {
        tmp_md.reset();
        Date &date = date_vec[d_idx];

        tmp_md.setDate(date);

        read_meta_data(current_file, tmp_md.meta);
        if (location_vec.size() == date_vec.size()) {
            setMeteoDataLocation(tmp_md, location_vec[d_idx], current_file, nodata);
        }

        setMeteoDataFields(tmp_md, current_file, date, indexes, nodata);
        vecMeteo.push_back(tmp_md);
    }
    return vecMeteo;
}

void iCSVIO::setMeteoDataLocation(MeteoData &tmp_md, geoLocation &loc, iCSVFile &current_file, double nodata) {
    tmp_md.meta.position.setPoint(IOUtils::standardizeNodata(loc.x, nodata), IOUtils::standardizeNodata(loc.y, nodata),
                                    IOUtils::standardizeNodata(loc.z, nodata), current_file.METADATA.epsg); // TODO: what happens if alt=nodata?
    tmp_md.meta.position.check("Inconsistent geographic coordinates in file \"" + current_file.filename + "\": ");
}

void iCSVIO::setMeteoDataFields(MeteoData &tmp_md, iCSVFile &current_file, Date &date, std::vector<size_t> &indexes, double nodata) {
    double offset = 0;
    double multiplier = 1;
    for (size_t field_idx = 0; field_idx < current_file.FIELDS.fields.size(); field_idx++) {
        if (!current_file.FIELDS.units_offsets.empty())
            offset = current_file.FIELDS.units_offsets[field_idx];
        if (!current_file.FIELDS.units_multipliers.empty())
            multiplier = current_file.FIELDS.units_multipliers[field_idx];
        const std::string &fieldname = current_file.FIELDS.fields[field_idx];
        const double value = current_file.readData(date, fieldname);

        if (fieldname == "timestamp" || fieldname == "julian")
            continue;
        // TODO: handle geometry = multiple columns
        if (fieldname == current_file.METADATA.geometry)
            continue;

        if (value == nodata) {
            tmp_md(indexes[field_idx]) = IOUtils::nodata;
        } else {
            tmp_md(indexes[field_idx]) = (value + offset) * multiplier;
        }
    }
}

// assume an operational snowpack virtual slopes naming: the station ID is made of letters followed by a station
// number (1 digit) and an optional virtual slope (1 digit, between 1 and 4)
double iCSVIO::getSnowpackSlope(const std::string &id) {
    static const char ALPHA[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

    const std::size_t end_name_pos = id.find_first_not_of(ALPHA);
    if (end_name_pos == std::string::npos)
        return IOUtils::nodata; // this is not a Snowpack virtual station naming
    if (id[end_name_pos] < '0' || id[end_name_pos] > '9')
        return IOUtils::nodata;

    if (id.size() - end_name_pos == 1)
        return 0.;                                // this is the flat field station
    const char slope_code = id[end_name_pos + 1]; // we are now sure that there is another char at this position

    if (slope_code == '1')
        return 1.;
    else if (slope_code == '2')
        return 2.;
    else if (slope_code == '3')
        return 3.;
    else if (slope_code == '4')
        return 4.;
    else
        return IOUtils::nodata; // this is not a digit or not in the [1-4] range
}

/**
* @brief Associate MeteoData parameter index with each iCSV field.
* @details This function associates a parameter index for MeteoData objects with the
* lineup of field types in iCSV FIELDS
*
* If a paramter is unknown in the fields section, then it is added as separate field to MeteoData
*
* @param[in] fields the fields coming from the SMET file
* @param[out] indexes the matching parameter indexes
* @param[out] julian_present set to true if a column contains a julian date
* @param[out] md a MeteoData object where extra parameters would be added
*/
void iCSVIO::identify_fields(const std::vector<std::string> &fields, std::vector<size_t> &indexes, MeteoData &md,
                                const std::string &geometry_field) {
    for (size_t ii = 0; ii < fields.size(); ii++) {
        const std::string &key = fields[ii];
        if (key == "julian" || key == "timestamp") {
            indexes.push_back(IOUtils::npos);
            continue;
        }

        // TODO: handle geometry = multiple columns
        if (key == geometry_field) {
            indexes.push_back(IOUtils::npos);
            continue;
        }

        if (md.param_exists(key)) {
            indexes.push_back(md.getParameterIndex(key));
            continue;
        }

        // specific key mapping
        if (key == "OSWR") {
            std::cerr << "The OSWR field name has been deprecated, it should be renamed into RSWR. Please update your files!!\n";
            indexes.push_back(md.getParameterIndex("RSWR"));
        } else if (key == "OLWR") {
            md.addParameter("OLWR");
            indexes.push_back(md.getParameterIndex("OLWR"));
        } else if (key == "PINT") { // in mm/h
            md.addParameter("PINT");
            indexes.push_back(md.getParameterIndex("PINT"));
        } else {
            // this is an extra parameter, we convert to uppercase
            const std::string extra_param(IOUtils::strToUpper(key));
            md.addParameter(extra_param);
            indexes.push_back(md.getParameterIndex(extra_param));
        }
    }
}

// ---------------------------- iCSVIO write -----------------------------------

void iCSVIO::writeMeteoData(const std::vector<std::vector<MeteoData>> &vecvecMeteo, const std::string &) {
    for (size_t ii = 0; ii < vecvecMeteo.size(); ii++) {
        if (vecvecMeteo[ii].empty())
            continue;

        std::vector<MeteoData> vecMeteo = vecvecMeteo[ii];
        iCSVFile outfile;
        bool file_exists = createFilename(outfile, vecMeteo[0].meta, ii);

        prepareOutfile(outfile, vecMeteo, file_exists);
        outfile.aggregateData(vecMeteo);

        if (acdd_metadata.isEnabled()) {
            acdd_metadata.setTimeCoverage(vecMeteo);
            acdd_metadata.setGeometry(getUniqueLocations(outfile), true);
        }

        writeToFile(outfile);
    }
}

// ----------------- iCSVIO write helper functions -----------------
/**
* @brief Creates a filename for the iCSV file based on the given station data.
*
* Creates a filename based on the station id, and if none is provided, the internal index
* of the station is used
*
* @param[out] outfile The iCSVFile object to store the generated filename.
* @param[in] station The StationData object containing the station information.
* @param[in] station_number The index used to generate the default station ID.
* @return True if the generated filename exists, false otherwise.
* @throws InvalidNameException if the generated filename or path is invalid.
*/
bool iCSVIO::createFilename(iCSVFile &outfile, const StationData &station, size_t ii) {
    outfile.METADATA.station_id = station.getStationID().empty() ? "Station" + IOUtils::toString(ii + 1) : station.getStationID();
    outfile.filename = outpath + "/" + outfile.METADATA.station_id + file_extension_out;
    if (!FileUtils::validFileAndPath(outfile.filename)) {
        throw InvalidNameException(outfile.filename, AT);
    }
    return FileUtils::fileExists(outfile.filename);
}

/**
* Prepares the output file for writing data.
* If the file exists and appending is allowed, it reads and parses the file
* If the file does not exist or appending is not allowed, it handles creating a new file
*
* checks the format validity and MeteoIO compatibility, and appends data to the file if needed.
*
* @param outfile The output file to be prepared.
* @param vecMeteo The vector of MeteoData containing the data to be written.
* @param file_exists A flag indicating whether the file already exists.
*/
void iCSVIO::prepareOutfile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists) {
    if (file_exists && allow_append) {
        outfile.readFile(outfile.filename, false);
        outfile.parseGeometry();
    } else {
        handleNewFile(outfile, vecMeteo, file_exists);
    }

    outfile.firstline = iCSV_firstline;
    outfile.checkFormatValidity();
    outfile.checkMeteoIOCompatibility();

    if (file_exists && allow_append) {
        handleFileAppend(outfile, vecMeteo);
    }
}

void iCSVIO::handleNewFile(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists) {
    if (file_exists && !allow_overwrite) {
        throw IOException("File " + outfile.filename + " already exists and iCSV_OVERWRITE is not allowed", AT);
    }
    createMetaDataSection(outfile, vecMeteo);
    createFieldsSection(outfile, vecMeteo);
}

void iCSVIO::handleFileAppend(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo) {
    if (outfile.location_in_header && outfile.station_location != toiCSVLocation(vecMeteo[0].meta.position, outfile.METADATA.epsg)) {
        throw IOException("Inconsistent geographic coordinates between header and data in file \"" + outfile.filename + "\": " +
                                outfile.station_location.toString() + " != " + toiCSVLocation(vecMeteo[0].meta.position, outfile.METADATA.epsg).toString(),
                            AT);
    }
    std::vector<std::string> columns_to_append = outfile.columnsToAppend(vecMeteo);
    if (!columns_to_append.empty()) {
        std::cerr << "Will append the following columns to file " << outfile.filename << ": "
                    << joinVector(columns_to_append, outfile.METADATA.field_delimiter) << "\n";
        outfile.FIELDS.fields.insert(outfile.FIELDS.fields.end(), columns_to_append.begin(), columns_to_append.end());
    }
}

/**
* @brief Creates the metadata section for the file to write
*
* This function adds metadata information to the output file, such as the timestamp, timezone, nodata value,
* field delimiter, EPSG code, location geometry, and optional metadata. It also checks the consistency of the
* location information and updates the location header accordingly.
*
* @param outfile The iCSVFile object representing the output file.
* @param vecMeteo A vector of MeteoData objects containing the meteorological data.
*/
void iCSVIO::createMetaDataSection(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo) {
    outfile.FIELDS.fields.push_back("timestamp"); // force time stamp to be the first field

    outfile.METADATA.timezone = TZ_out;
    outfile.METADATA.nodata = IOUtils::nodata;
    outfile.METADATA.field_delimiter = out_delimiter;

    Coords loc = vecMeteo[0].meta.position;
    loc.setProj(coordout, coordoutparam);
    int epsg = loc.getEPSG();
    outfile.METADATA.setEPSG(epsg);

    outfile.location_in_header = checkLocationConsistency(vecMeteo);
    if (outfile.location_in_header) {
        double x = loc.getEasting();  // is this the correct way to get xy according to epsg?
        double y = loc.getNorthing(); // is this the correct way to get xy according to epsg?
        double z = loc.getAltitude(); // is this the correct way to get xy according to epsg?
        outfile.station_location = geoLocation(x, y, z);
        outfile.METADATA.geometry = getGeometry(outfile.station_location);
    } else {
        outfile.METADATA.geometry = "geometry";
        outfile.FIELDS.fields.push_back("geometry");
    }
    outfile.METADATA.optional_metadata = vecMeteo[0].meta.extra;
    outfile.findLocation();
}

// TODO: somehow get other information like longname etc.
void iCSVIO::createFieldsSection(iCSVFile &outfile, const std::vector<MeteoData> &vecMeteo) {
    std::set<std::string> available_params = MeteoData::listAvailableParameters(vecMeteo);
    for (const auto &param : available_params) {
        outfile.FIELDS.fields.push_back(param);
    }
}

bool iCSVIO::checkLocationConsistency(const std::vector<MeteoData> &vecMeteo) {
    for (size_t ii = 1; ii < vecMeteo.size(); ii++) {
        const Coords &p1 = vecMeteo[ii - 1].meta.position;
        const Coords &p2 = vecMeteo[ii].meta.position;
        if (p1 != p2) {
            // we don't mind if p1==nodata or p2==nodata
            if (p1.isNodata() == false && p2.isNodata() == false)
                return false;
        }
    }

    return true;
}

std::string iCSVIO::getGeometry(const geoLocation &loc) {
    bool dim_2 = loc.z == IOUtils::nodata;
    std::string geometry = dim_2 ? "POINT(" : "POINTZ(";
    geometry += std::to_string(loc.x) + " " + std::to_string(loc.y);
    if (!dim_2) {
        geometry += " " + std::to_string(loc.z);
    }
    geometry += ")";
    return geometry;
}

void iCSVIO::writeToFile(const iCSVFile &outfile) {
    ofilestream file(outfile.filename);
    if (!file.is_open()) {
        throw IOException("Unable to open file " + outfile.filename, AT);
    }
    file << outfile.firstline << "\n";
    file << "# [METADATA]\n";

    std::map<std::string, std::string> metadata = outfile.METADATA.toOutputMap();
    for (const auto &meta : metadata) {
        file << "# " << meta.first << " = " << meta.second << "\n";
    }

    // print acdd metadata
    if (acdd_metadata.isEnabled()) {
        // print ACDD headers
        for (auto it = acdd_metadata.cbegin(); it != acdd_metadata.cend(); ++it) {
            const std::string header_field(it->second.getName());
            const std::string value(it->second.getValue());
            if (header_field.empty() || value.empty())
                continue;

            file << "# " << header_field << " = " << value << "\n";
        }
    }

    file << "# [FIELDS]\n";
    std::map<std::string, std::vector<std::string>> fields = outfile.FIELDS.toMap();
    for (const auto &field : fields) {
        file << "# " << field.first << " = " << joinVector(field.second, outfile.METADATA.field_delimiter) << "\n";
    }
    file << "# [DATA]\n";
    // timestamp and geometry will not be in row data
    size_t num_data_fields = outfile.FIELDS.fields.size() - 1;
    if (!outfile.location_in_header) {
        num_data_fields--;
    }
    const auto& out_data = outfile.getRowData();
    const auto& out_dates = outfile.getAllDatesInFile();
    const auto& out_locations = outfile.getAllLocationsInData();
    const double nodata = outfile.getNoData();
    for (size_t ii = 0; ii < outfile.getRowData().size(); ii++) {
        size_t data_idx = 0;
        for (size_t jj = 0; jj < outfile.FIELDS.fields.size(); jj++) {
            if (outfile.FIELDS.fields[jj] == "timestamp") {
                file << out_dates[ii].toString(Date::ISO) << outfile.METADATA.field_delimiter;
            } else if (outfile.FIELDS.fields[jj] == outfile.METADATA.geometry) {
                file << getGeometry(out_locations[ii]) << outfile.METADATA.field_delimiter;
            } else {
                if (data_idx < num_data_fields) {
                    file << IOUtils::standardizeNodata(out_data[ii][data_idx], nodata);
                } else {
                    file << IOUtils::nodata;
                }
                data_idx++;
                if (jj != outfile.FIELDS.fields.size() - 1) {
                    file << outfile.METADATA.field_delimiter;
                }
            }
        }
        if (ii != out_data.size() - 1) {
            file << "\n";
        }
    }
    file.close();
}

} // namespace
