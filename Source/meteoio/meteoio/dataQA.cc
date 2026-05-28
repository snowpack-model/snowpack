// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*                   Copyright GridGroup, EIA-FR 2010                              */
/*  Copyright 2010-2013 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   */
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

#include <iostream>
#include <meteoio/dataQA.h>

namespace mio {

DataQA::DataQA(const Config& cfg) : enabled(false)
{
	cfg.getValue("DATA_QA_LOGS", "GENERAL", enabled, IOUtils::nothrow);
}

DataQA::DataQA(bool i_enabled) : enabled(i_enabled) {}

void DataQA::printQA(OperationType operationType, const std::string& stat, const std::string& parname, const std::string& processingName, const Date& date) const
{
	printQA(operationTypeToString(operationType), stat, parname, processingName, date);
}

void DataQA::printQA(const std::string& operationTypeString, const std::string& stat, const std::string& parname, const std::string& processingName, const Date& date) const
{
	std::cout << "[DATA_QA] " << operationTypeString << " " << stat << "::" << parname;
	if (!processingName.empty()) {
		std::cout << "::" << processingName;
	}
	std::cout << " " << date.toString(Date::ISO_TZ) << " [" << date.toString(Date::ISO_WEEK) << "]\n";
}

std::string DataQA::operationTypeToString(OperationType operationType)
{
	switch (operationType) {
		case GENERATION: return "Generating";
		case RESAMPLING: return "Resampling";
		case FILTERING: return "Filtering";
		case MISSING: return "Missing";
		case CUSTOM: return "Custom";
		default: return "Unknown";
	}
}

} //namespace

