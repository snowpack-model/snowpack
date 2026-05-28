// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*                   Copyright GridGroup, EIA-FR 2010                              */
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef DATAQA_H
#define DATAQA_H

#include <meteoio/Config.h>
#include <meteoio/dataClasses/Date.h>
#include <string>

namespace mio {

/**
 * @class DataQA
 * @brief A class for handling data quality assurance logging.
 *
 * This class centralizes all data QA logging functionality using a generic
 * approach with operation types specified via enum.
 */
class DataQA {
public:
	/**
	 * @brief Enum defining the types of QA operations
	 */
	enum OperationType {
		GENERATION,    ///< Data generation operation
		RESAMPLING,    ///< Data resampling operation
		FILTERING,     ///< Data filtering operation
		MISSING,       ///< Missing data operation
		CUSTOM         ///< Custom operation type
	};

	/**
	 * @brief Constructor
	 * @param cfg Configuration object to read DATA_QA_LOGS setting
	 */
	explicit DataQA(const Config& cfg);

	/**
	 * @brief Constructor with explicit enable flag
	 * @param enabled Whether data QA logging is enabled
	 */
	explicit DataQA(bool i_enabled = false);

	/**
	 * @brief Print a QA log message with specific operation type
	 * @param operationType Type of operation (GENERATION, RESAMPLING, FILTERING, etc.)
	 * @param stat Station ID
	 * @param parname Parameter name
	 * @param processingName Processing algorithm name
	 * @param date Date of the operation
	 */
	void printQA(OperationType operationType, const std::string& stat, const std::string& parname, const std::string& processingName, const Date& date) const;

	/**
	 * @brief Print a QA log message with custom operation type string
	 * @param operationTypeString Custom operation type string
	 * @param stat Station ID
	 * @param parname Parameter name
	 * @param processingName Processing algorithm name
	 * @param date Date of the operation
	 */
	void printQA(const std::string& operationTypeString, const std::string& stat, const std::string& parname, const std::string& processingName, const Date& date) const;

	/**
	 * @brief Check if data QA logging is enabled
	 * @return true if enabled, false otherwise
	 */
	bool isEnabled() const { return enabled; }

	/**
	 * @brief Enable or disable data QA logging
	 * @param enable true to enable, false to disable
	 */
	void setEnabled(bool enable) { enabled = enable; }

	/**
	 * @brief Convert operation type enum to string
	 * @param operationType The operation type enum value
	 * @return String representation of the operation type
	 */
	static std::string operationTypeToString(OperationType operationType);

private:
	bool enabled; ///< Whether data QA logging is enabled
};

} //end namespace mio
#endif
