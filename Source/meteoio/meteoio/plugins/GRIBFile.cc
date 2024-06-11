// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/GRIBFile.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOUtils.h>

#include <fstream>


namespace mio {
    using namespace codes;

    // ------------------------- GRIBTable -------------------------

    GRIBTable::GRIBTable()
        : filename(), param_indexing("paramId"), level_indexing("typeOfLevel"), param_table(), param_table_double(), param_table_long(GRIB_DEFAULT_PARAM_TABLE), level_type_table(GRIB_DEFAULT_LEVELTYPE_TABLE), level_no_table(GRIB_DEFAULT_LEVELNO_TABLE),
          parameter_id_type(PARAM_TYPE::LONG), known_params() {
        
        // add the default parameters to the known parameters
        for (auto &p : GRIB_DEFAULT_PARAM_TABLE) {
            known_params.insert(p.first);
        }

    }



    GRIBTable::GRIBTable(const std::string &in_filename)
        : filename(in_filename), param_indexing(), level_indexing(), param_table(), param_table_double(), param_table_long(), level_type_table(), level_no_table(),
          parameter_id_type(PARAM_TYPE::LONG), known_params() {
        // Open the file
        init_known_params();
        readTable();
    }

    void GRIBTable::init_known_params() {
        for (size_t i=0; i<MeteoData::nrOfParameters; i++) {
            std::string param_name = MeteoData::getParameterName(i);
            known_params.insert(param_name);
            param_table[param_name] = "";
            param_table_double[param_name] = static_cast<double>(IOUtils::npos);;
            param_table_long[param_name] = static_cast<long>(IOUtils::npos);;
            level_type_table[param_name] = "";
            level_no_table[param_name] = static_cast<long>(IOUtils::npos);
        }
        for (size_t i=0; i<MeteoGrids::nrOfParameters; i++) {
            std::string param_name = MeteoGrids::getParameterName(i);
            known_params.insert(param_name);
            param_table[param_name] = "";
            param_table_double[param_name] = static_cast<double>(IOUtils::npos);;
            param_table_long[param_name] = static_cast<long>(IOUtils::npos);;
            level_type_table[param_name] = "";
            level_no_table[param_name] = static_cast<long>(IOUtils::npos);
        }
    };

    // ------------------------- reading the GRIB table -------------------------
    void GRIBTable::readTable() {
        std::ifstream file(filename);

        if (!file.is_open()) {
            throw AccessException("Could not open GRIB table: " + filename, AT);
        }

        // Read the file line by line and parse it
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#' || line[0] == ' ') {
                continue;
            }

            std::vector<std::string> values = IOUtils::split(line, ':');

            if (parseIndexing(values))
                continue;
            if (parseParamType(values))
                continue;

            // only comments or parameters left so fill the parameter table
            fillParameterTables(values);
        }
    }

    bool GRIBTable::parseIndexing(const std::vector<std::string> &line_vals) {
        if (line_vals[0] == "index") {
            if (line_vals.size() != 3) {
                std::stringstream ss;
                ss << "Invalid number of indexing values to retrieve data, got: ";
                for (const auto &v : line_vals) {
                    ss << v << " ";
                }
                ss << std::endl;
                ss << "Expected: index:param:levelType" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            param_indexing = line_vals[1];
            level_indexing = line_vals[2];
            return true;
        }
        return false;
    }

    bool GRIBTable::parseParamType(const std::vector<std::string> &line_vals) {
        if (line_vals[0] == "paramIdType") {
            if (line_vals.size() != 2) {
                std::stringstream ss;
                ss << "Invalid parameter type specification, got: ";
                for (const auto &v : line_vals) {
                    ss << v << " ";
                }
                ss << std::endl;
                ss << "Expected: paramIdType:[doubele|str|long]" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            if (line_vals[1] == "double") {
                parameter_id_type = PARAM_TYPE::DOUBLE;
            } else if (line_vals[1] == "string") {
                parameter_id_type = PARAM_TYPE::STRING;
            } else if (line_vals[1] == "long") {
                parameter_id_type = PARAM_TYPE::LONG;
            } else {
                std::stringstream ss;
                ss << "Invalid parameter type specification, got: " << line_vals[1] << std::endl;
                ss << "Expected: paramIdType:[doubele|str]" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            return true;
        }
        return false;
    }

    void GRIBTable::fillParameterTables(const std::vector<std::string> &line_vals) {
        std::string key, paramID, levelType;
        long levelNo;
        if (line_vals.size() != 4) {
            std::stringstream ss;
            ss << "Invalid number of values to retrieve data, got: ";
            for (const auto &v : line_vals) {
                ss << v << " ";
            }
            ss << std::endl;
            ss << "Expected: paramName:paramID:levelType:levelNo" << std::endl;
            throw InvalidFormatException(ss.str(), AT);
        }

        key = line_vals[0];
        paramID = line_vals[1];
        levelType = line_vals[2];
        if (!IOUtils::convertString(levelNo, line_vals[3])) {
            throw InvalidArgumentException("Invalid level number: " + line_vals[3], AT);
        };

        if (known_params.find(key) == known_params.end()) {
            throw InvalidArgumentException("Unknown parameter: " + key, AT);
        }

        if (parameter_id_type == PARAM_TYPE::DOUBLE) {
            double paramID_double;
            if (!IOUtils::convertString(paramID_double, paramID)) {
                throw InvalidArgumentException("Could not convert parameter ID to number (as specified in paramType): " + paramID, AT);
            }
            param_table_double[key] = paramID_double;
        } else if (parameter_id_type == PARAM_TYPE::STRING) {
            param_table[key] = paramID;
        } else if (parameter_id_type == PARAM_TYPE::LONG) {
            long paramID_long;
            if (!IOUtils::convertString(paramID_long, paramID)) {
                throw InvalidArgumentException("Could not convert parameter ID to number (as specified in paramType): " + paramID, AT);
            }
            param_table_long[key] = paramID_long;
        } else {
            throw InvalidArgumentException("Invalid parameter type, if this appears, please contact the developers", AT);
        }
        
        level_type_table[key] = levelType;
        level_no_table[key] = levelNo;
    }

    // ------------------------- GETTERS -------------------------
    std::vector<std::string> GRIBTable::getIndexes() const {
        return {param_indexing, level_indexing};
    }

    void GRIBTable::getParamId(const std::string &param_name, std::string &paramId, double &paramId_num, long &paramId_long) const {
        if (parameter_id_type == PARAM_TYPE::DOUBLE) {
            paramId_num = param_table_double.at(param_name);
            paramId = "";
            paramId_long = static_cast<long> (IOUtils::npos);
        } else if (parameter_id_type == PARAM_TYPE::STRING) {
            paramId = param_table.at(param_name);
            paramId_num = static_cast<double> (IOUtils::npos);
            paramId_long = static_cast<long> (IOUtils::npos);
        } else if (parameter_id_type == PARAM_TYPE::LONG) {
            paramId_long = param_table_long.at(param_name);
            paramId = "";
            paramId_num = static_cast<double> (IOUtils::npos);
        } else {
            throw InvalidArgumentException("Invalid parameter type, if this appears, please contact the developers", AT);
        }   
    }

    std::string GRIBTable::getLevelType(const std::string &param_name) const { return level_type_table.at(param_name); }

    long GRIBTable::getLevelNo(const std::string &param_name) const { return level_no_table.at(param_name); }

#ifdef DEBUG
    void GRIBTable::printTable() const {
        std::cerr << "GRIB table: " << filename << std::endl;
        std::cerr << "Indexing: " << param_indexing << " " << level_indexing << std::endl;
        std::cerr << "ParamIdType: ";
        if (parameter_id_type == PARAM_TYPE::DOUBLE) {
            std::cerr << "double" << std::endl;
        } else if (parameter_id_type == PARAM_TYPE::STRING) {
            std::cerr << "string" << std::endl;
        } else if (parameter_id_type == PARAM_TYPE::LONG) {
            std::cerr << "long" << std::endl;
        } else {
            std::cerr << "unknown" << std::endl;
        }
        for (const auto &p : param_table) {
            std::cerr << p.first << " " << p.second << " " << param_table_double.at(p.first) << " " << param_table_long.at(p.first) << " " << level_type_table.at(p.first) << " " << level_no_table.at(p.first) << std::endl;
        }
        for (const auto &p : param_table_double) {
            std::cerr << p.first << " " << p.second << " " << param_table.at(p.first) << " " << param_table_long.at(p.first) << " " << level_type_table.at(p.first) << " " << level_no_table.at(p.first) << std::endl;
        }
        for (const auto& p : param_table_long) {
            std::cerr << p.first << " " << p.second << " " << param_table.at(p.first) << " " << param_table_double.at(p.first) << " " << level_type_table.at(p.first) << " " << level_no_table.at(p.first) << std::endl;
        }
    }
#endif


    // ------------------------- GRIBFile -------------------------
    GRIBFile::GRIBFile(const std::string &in_filename, const std::vector<std::string> &indexes)
        : filename(in_filename), file(), grid_params(), timepoints() {

            file = indexFile(filename, indexes, false); // true = verbose output
            checkValidity();
        }
    
    void GRIBFile::checkValidity() {
        std::set<std::pair<Date, long>> timepoints_for_parameter;
        std::vector<CodesHandlePtr> messages = getMessages(filename, PRODUCT_GRIB);
        if (messages.empty()) {
            throw AccessException("No messages found in GRIB file: " + filename);
        }

        long edition = -1;
        for (auto &m : messages) {
            long curr_edition;
            getParameter(m, "editionNumber", curr_edition);
            if (edition == -1)
                edition = curr_edition;
            else if (edition != curr_edition) {
                // throw InvalidFormatException("Different GRIB editions found in GRIB file: " + filename, AT);
            }
            std::map<std::string, double> new_grid_params = getGridParameters(m);
            if (grid_params.empty()) {
                grid_params = new_grid_params;
            } else if (!new_grid_params.empty() && grid_params != new_grid_params) {
                std::cerr << "Grid parameters do not match in GRIB file: " << filename << std::endl; // TODO: change back
                // throw InvalidFormatException("Grid parameters do not match in GRIB file: " + filename,AT);
            }

            long paramId;
            Date curr_date = getMessageDateGrib(m, 0);
            getParameter(m, "paramId", paramId);
            std::pair<std::set<std::pair<Date,long>>::iterator,bool> res = timepoints_for_parameter.insert(std::make_pair(curr_date, paramId));
            timepoints.insert(curr_date);
            if (!res.second) {
                throw InvalidFormatException("Duplicate timepoints (" + curr_date.toString(Date::ISO) + ") found in GRIB file: " + filename, AT);
            }
        }
        if (grid_params.empty()) {
            throw InvalidFormatException("No grid parameters found in GRIB file: " + filename);
        }
        if (timepoints.empty()) throw InvalidFormatException("No timepoint found in GRIB file: " + filename,AT);
    }  
} // namespace mio