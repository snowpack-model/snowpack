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
#include <meteoio/IOUtils.h>
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/dataClasses/StationData.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>

#include <regex>

using namespace std;

namespace mio {

const std::string Meteo1DInterpolator::interpol_section("Interpolations1D");
const std::string Meteo1DInterpolator::interpol_pattern("::RESAMPLE");


static void eraseArg(std::vector<std::pair<std::string, std::string>> &vecArgs, const std::string &argname)
{
	std::vector<std::pair<std::string, std::string>>::iterator it( vecArgs.begin() );
	while (it != vecArgs.end()) {
		if (it->first == argname) {
			it = vecArgs.erase(it);
		} else {
			++it;
		}
	}
}


// --------------------- Meteo1DInterpolator Implementation --------------------------------------

static bool checkDeprecatedWindowSize(const Config &cfg)
{
	if (cfg.keyExistsRegex(".*WINDOW_SIZE.*", Meteo1DInterpolator::interpol_section)) {
		std::cout << "[W] Using deprecated key: WINDOW_SIZE in section " + Meteo1DInterpolator::interpol_section + ".\n Please use MAX_GAP_SIZE instead.\n";
		return true;
	}
	return false;
}

Meteo1DInterpolator::Meteo1DInterpolator(const Config &in_cfg, const char &rank, const IOUtils::OperationMode &mode)
    : mapAlgorithms(), cfg(in_cfg), max_gap_size(86400.), enable_resampling(true), data_qa_logs(false), gap_size_key("MAX_GAP_SIZE")
{

	if (checkDeprecatedWindowSize(cfg)) gap_size_key = "WINDOW_SIZE";

	cfg.getValue("DATA_QA_LOGS", "GENERAL", data_qa_logs, IOUtils::nothrow);
	cfg.getValue("ENABLE_RESAMPLING", "Interpolations1D", enable_resampling, IOUtils::nothrow);

	// default max_gap_size is 2 julian days
	cfg.getValue("MAX_GAP_SIZE", "Interpolations1D", max_gap_size, IOUtils::nothrow);
	if (max_gap_size < 0.)
		throw IOException("MAX_GAP_SIZE not valid, it should be a duration in seconds at least greater than 0", AT);
	max_gap_size /= 86400.; // user uses seconds, internally julian day is used
	if (max_gap_size == 0.)
		enable_resampling = false;

	createResamplingStacks(mode, rank);
}

static int findIndex(const std::string& key)
{
	// Define the regular expression pattern
	static const std::regex pattern("^.+::resample(\\d+)$", std::regex::icase);

	std::smatch match;
	// Search for the pattern in the key
	if (!std::regex_search(key, match, pattern))
		throw IOException("Resampling Key " + key + " does not contain a valid index; I.e. it does not match the pattern PARAM::resample#",AT);

	// A match has been found, return the number
	return std::stoi(match[1]);
}

static std::vector<std::pair<int, std::string>> orderAlgoStack(const std::vector<std::pair<std::string, std::string>> &vecAlgos)
{
	std::vector<std::pair<int, std::string>> orderedAlgos;
	for (const auto &key_algo : vecAlgos) {
		const std::string key( key_algo.first );
		const int index = findIndex(key);
		orderedAlgos.push_back(std::make_pair(index, key_algo.second));
	}

	auto orderByIndex = [](const std::pair<int, std::string> &a, const std::pair<int, std::string> &b) {
		return a.first < b.first;
	};

	std::sort(orderedAlgos.begin(), orderedAlgos.end(), orderByIndex);
	return orderedAlgos;
}

void Meteo1DInterpolator::createResamplingStacks(const IOUtils::OperationMode &mode, const char &rank)
{
	for (size_t ii = 0; ii < MeteoData::nrOfParameters; ii++) {     // loop over all MeteoData member variables
		const std::string parname( MeteoData::getParameterName(ii) ); // Current parameter name

		// extract each interpolation algorithm and its arguments, then build the stack
		const std::vector<std::pair<std::string, std::string>> vecAlgos( cfg.getValues(parname + interpol_pattern, interpol_section) );
		const std::vector<std::pair<int, std::string>> orderedAlgos( orderAlgoStack(vecAlgos) );
		mapAlgorithms[parname] = ResamplingStack();
		processAlgorithms(parname, orderedAlgos, "", mode, rank);
	}
}

static void checkDeprecatedKeySyntax(const mio::Config& cfg, const std::string& parname, const std::string& algo_name, const std::string& interpol_section)
{
	if (cfg.keyExistsRegex(parname + "::" + algo_name +".*", interpol_section)) {
		std::stringstream ss;
		ss << "Using deprecated interpolation algorithm syntax " << parname << "::" << algo_name << " in section " << interpol_section << std::endl;
		ss << "Please use the new syntax " << parname << "::" << "ARG#" << std::endl;
		ss << "For detailed information, see the documentation: https://meteoio.slf.ch/doc-release/html/resampling.html" << std::endl;
		throw IOException(ss.str(), AT);
	}
}


void Meteo1DInterpolator::addAlgorithmToStack(const std::string& parname,const std::string& algo_name, const std::vector<std::pair<std::string, std::string>>& vecArgs, const double& i_max_gap_size)
{
	const std::shared_ptr<ResamplingAlgorithms> algo_ptr( ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, i_max_gap_size, vecArgs, cfg) );
	mapAlgorithms[parname].addAlgorithm(algo_ptr, i_max_gap_size);
}

void Meteo1DInterpolator::createDefaultAlgorithm(const std::string &parname)
{
	std::vector<std::pair<std::string, std::string>> vecArgs; // is empty anyways as we already iterated through all specified algorithms
	addAlgorithmToStack(parname, "LINEAR", vecArgs, max_gap_size);
}

void Meteo1DInterpolator::processAlgorithms(const std::string &parname, const std::vector<std::pair<int, std::string>> &vecAlgos, std::string base_parname,
											const IOUtils::OperationMode &mode, const char &rank)
{
	// ARG# is the syntax instad of ALGO# for the interpolation algorithms (to support repeated algorithms with different arguments)
	static const std::string arguments_ini_key( "ARG" );

	// in case we dont even care about the height
	if (base_parname.empty()) base_parname = parname;

	// process each algorithm, and possibly add it to the stack
	for (const auto &index_algo : vecAlgos) {
		const std::string algo_name( IOUtils::strToUpper(index_algo.second) );
		const int algo_index = index_algo.first;
		checkDeprecatedKeySyntax(cfg, base_parname, algo_name, interpol_section);

		double max_gap_size_override( cfg.get(base_parname + "::" + arguments_ini_key + std::to_string(algo_index) + "::" + gap_size_key, interpol_section, max_gap_size) );
		if (max_gap_size_override != max_gap_size) {
			if (max_gap_size_override>0)
				max_gap_size_override /= 86400.; // user uses seconds, internally julian day is used
			else
				throw IOException("Invalid MAX_GAP_SIZE, it should be a duration in seconds at least greater than 0", AT);
		}

		if (algo_name == "ACCUMULATE" && mode == IOUtils::VSTATIONS && rank == 1) {
			// handle accumulation algorithm
			const std::string vstations_refresh_rate = cfg.get("VSTATIONS_REFRESH_RATE", "InputEditing");
			const std::vector<std::pair<std::string, std::string>> vecArgs(1, std::make_pair("PERIOD", vstations_refresh_rate));
			addAlgorithmToStack(parname, algo_name, vecArgs, max_gap_size_override);
		} else {
			std::vector<std::pair<std::string, std::string>> vecArgs(cfg.getArgumentsForAlgorithm(base_parname, arguments_ini_key, algo_index));
			eraseArg(vecArgs, gap_size_key);
			addAlgorithmToStack(parname, algo_name, vecArgs, max_gap_size_override);
		}
	}
	if (mapAlgorithms[parname].empty()) { // create a default linear algorithm
		createDefaultAlgorithm(parname);
	}
}

void Meteo1DInterpolator::getWindowSize(ProcessingProperties &o_properties) const
{
	o_properties.points_before = 1;
	o_properties.points_after = 1;
	o_properties.time_before = Duration(max_gap_size, 0.); // we will need to cut a window 2x larger so we can interpolate each point in the window
	o_properties.time_after = Duration(max_gap_size, 0.);
}

bool Meteo1DInterpolator::resampleData(const Date &date, const std::string &stationHash, const std::vector<MeteoData> &vecM, MeteoData &md)
{
	if (vecM.empty()) // Deal with case of the empty vector
		return false; // nothing left to do

	// Find element in the vector or the next index
	size_t index = IOUtils::seek(date, vecM, false);

	// Three cases
	bool isResampled = true;
	ResamplingAlgorithms::ResamplingPosition elementpos = ResamplingAlgorithms::exact_match;
	if (index == IOUtils::npos) { // nothing found append new element at the left or right
		if (date < vecM.front().date) {
			elementpos = ResamplingAlgorithms::begin;
			index = 0;
		} else if (date >= vecM.back().date) {
			elementpos = ResamplingAlgorithms::end;
			index = vecM.size() - 1;
		}
	} else if ((index != IOUtils::npos) && (vecM[index].date != date)) { // element found nearby
		elementpos = ResamplingAlgorithms::before;
	} else { // element found at the right time
		isResampled = false;
	}
	md = vecM[index]; // create a clone of the found element

	if (!enable_resampling) {
		if (isResampled == false)
			return true; // the element was found at the right time
		else {           // not found or wrong time: return a nodata element
			md.reset();
			md.setDate(date);
			return true;
		}
	}

	md.reset(); // set all values to IOUtils::nodata
	md.setDate(date);
	md.setResampled(isResampled);

	// now, perform the resampling
	for (size_t ii = 0; ii < md.getNrOfParameters(); ii++) {
		const std::string parname( md.getNameForParameter(ii) ); // Current parameter name

		const std::map<std::string, ResamplingStack>::const_iterator it = mapAlgorithms.find(parname);
		if (it != mapAlgorithms.end()) { // the parameter has been found
			it->second.resample(stationHash, index, elementpos, ii, vecM, md, max_gap_size);

		} else { // we are dealing with an extra parameter, we need to add it to the map first, so it will exist next time...
			// all parameters with height will end up here.
			double height; std::string base_parname;
			if (!MeteoData::getTypeAndNo(parname, base_parname, height)) {
				base_parname = parname;
			}

			const std::vector<std::pair<std::string, std::string>> vecAlgos(cfg.getValues(base_parname + interpol_pattern, interpol_section));
			const std::vector<std::pair<int, std::string>> orderedAlgos(orderAlgoStack(vecAlgos));
			mapAlgorithms[parname] = ResamplingStack();

			// TODO: we need base name here
			processAlgorithms(parname, orderedAlgos, base_parname);
			mapAlgorithms[parname].resample(stationHash, index, elementpos, ii, vecM, md, max_gap_size);
		}

		if ((index != IOUtils::npos) && vecM[index](ii) != md(ii)) {
			md.setResampledParam(ii);
			if (data_qa_logs) {
				const std::map<std::string, ResamplingStack>::const_iterator it2 = mapAlgorithms.find(parname); // we have to re-find it in order to handle extra parameters
				const std::string statName(md.meta.getStationName());
				const std::string statID(md.meta.getStationID());
				const std::string stat( (!statID.empty()) ? statID : statName );
				const std::string algo_name(it2->second.getStackStr());
				std::cout << "[DATA_QA] Resampling " << stat << "::" << parname << "::" << algo_name << " " << md.date.toString(Date::ISO_TZ) << " [" << md.date.toString(Date::ISO_WEEK) << "]\n";
			}
		}
	} // endfor ii

	return true; // successfull resampling
}

void Meteo1DInterpolator::resetResampling()
{
	//std::map<std::string, ResamplingStack>::iterator it;
	for (auto it = mapAlgorithms.begin(); it != mapAlgorithms.end(); ++it) {
		it->second.resetResampling();
	}
}

std::string Meteo1DInterpolator::getAlgorithmsForParameter(const std::string &parname) const
{
	std::string algo("linear"); // default value
	cfg.getValue(parname + "::resample", "Interpolations1D", algo, IOUtils::nothrow);
	return algo;
}

Meteo1DInterpolator &Meteo1DInterpolator::operator=(const Meteo1DInterpolator &source)
{
	if (this != &source) {
		mapAlgorithms = source.mapAlgorithms;
		max_gap_size = source.max_gap_size;
		enable_resampling = source.enable_resampling;
		data_qa_logs = source.data_qa_logs;
	}
	return *this;
}

const std::string Meteo1DInterpolator::toString() const
{
	ostringstream os;
	os << "<Meteo1DInterpolator>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	if (enable_resampling) {
		os << "Resampling algorithms:\n";
		map<string, ResamplingStack>::const_iterator it;
		for (it = mapAlgorithms.begin(); it != mapAlgorithms.end(); ++it) {
			// os << setw(10) << it->first << "::" << it->second->getAlgo() << "\n";
			os << it->second.getStackStr() << "\n";
		}
	} else {
		os << "Resampling disabled\n";
	}
	os << "</Meteo1DInterpolator>\n";

	return os.str();
}

// --------------------- Resampling Stack Implementation --------------------------------------
ResamplingStack::ResamplingStack() : max_gap_sizes(), stack() {}

void ResamplingStack::addAlgorithm(std::shared_ptr<ResamplingAlgorithms> algo, const double &max_gap_size) {
	stack.push_back(algo);
	max_gap_sizes.push_back(max_gap_size);
}

std::vector<std::shared_ptr<ResamplingAlgorithms>> ResamplingStack::buildStack(const ResamplingAlgorithms::gap_info &gap) const
{
	std::vector<std::shared_ptr<ResamplingAlgorithms>> res;
	for (size_t ii = 0; ii < stack.size(); ii++) {
		if (max_gap_sizes[ii] == IOUtils::nodata)
			res.push_back(stack[ii]);
		else if (gap.size() <= max_gap_sizes[ii]) {
			res.push_back(stack[ii]);
		}
	}
	return res;
}

void ResamplingStack::resetResampling()
{
	for (size_t ii = 0; ii < stack.size(); ii++) {
		stack[ii]->resetResampling();
	}
}

void ResamplingStack::resample(const std::string &stationHash, const size_t &index, const ResamplingAlgorithms::ResamplingPosition elementpos, const size_t &ii, const std::vector<MeteoData> &vecM,
                                MeteoData &md, const double &i_max_gap_size) const
{
	const ResamplingAlgorithms::gap_info gap( ResamplingAlgorithms::findGap(index, ii, vecM, md.date, i_max_gap_size) );
	const std::vector<std::shared_ptr<ResamplingAlgorithms>> resampling_stack = buildStack(gap);

	for (size_t jj = 0; jj < resampling_stack.size(); jj++) {
		resampling_stack[jj]->resample(stationHash, index, elementpos, ii, vecM, md);
		if (jj > 0 && (index != IOUtils::npos) && vecM[index](ii) != md(ii)) {
			break;
		} else if (ResamplingAlgorithms::exact_match == elementpos && vecM[index](ii) == md(ii)) {
			break;
		}
	}
}

bool ResamplingStack::empty() const { return stack.empty(); }

std::string ResamplingStack::getStackStr() const
{
	ostringstream os;
	os << "[";
	for (size_t ii = 0; ii < stack.size(); ii++) {
		os << stack[ii]->getAlgo();
		if (ii < stack.size() - 1)
			os << ", ";
	}
	os << "]";
	return os.str();
}

} // namespace
