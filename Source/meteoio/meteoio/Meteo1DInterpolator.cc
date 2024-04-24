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
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/dataClasses/StationData.h>

#include <iostream>
#include <utility>

using namespace std;

namespace mio {

    const std::string Meteo1DInterpolator::interpol_section("Interpolations1D");
    const std::string Meteo1DInterpolator::interpol_pattern("::RESAMPLE");

    static void eraseArg(std::vector<std::pair<std::string, std::string>> &vecArgs, const std::string &argname) {
        std::vector<std::pair<std::string, std::string>>::iterator it = vecArgs.begin();
        while (it != vecArgs.end()) {
            if (it->first == argname) {
                it = vecArgs.erase(it);
            } else {
                ++it;
            }
        }
    }

    // --------------------- Meteo1DInterpolator Implementation --------------------------------------

    Meteo1DInterpolator::Meteo1DInterpolator(const Config &in_cfg, const char &rank, const IOUtils::OperationMode &mode)
        : mapAlgorithms(), cfg(in_cfg), window_size(86400.), enable_resampling(true), data_qa_logs(false) {
        cfg.getValue("DATA_QA_LOGS", "GENERAL", data_qa_logs, IOUtils::nothrow);
        cfg.getValue("ENABLE_RESAMPLING", "Interpolations1D", enable_resampling, IOUtils::nothrow);

        // default window_size is 2 julian days
        cfg.getValue("WINDOW_SIZE", "Interpolations1D", window_size, IOUtils::nothrow);
        if (window_size < 0.)
            throw IOException("WINDOW_SIZE not valid, it should be a duration in seconds at least greater than 0", AT);
        window_size /= 86400.; // user uses seconds, internally julian day is used
        if (window_size == 0.)
            enable_resampling = false;

        createResamplingStacks(mode, rank);
    }

    void Meteo1DInterpolator::createResamplingStacks(const IOUtils::OperationMode &mode, const char &rank) {
        for (size_t ii = 0; ii < MeteoData::nrOfParameters; ii++) {     // loop over all MeteoData member variables
            const std::string parname(MeteoData::getParameterName(ii)); // Current parameter name

            // extract each interpolation algorithm and its arguments, then build the stack
            const std::vector<std::pair<std::string, std::string>> vecAlgos(cfg.getValues(parname + interpol_pattern, interpol_section));
            mapAlgorithms[parname] = ResamplingStack();
            processAlgorithms(true, parname, vecAlgos, mode, rank);
        }
    }

    void Meteo1DInterpolator::processAlgorithms(bool first_time, const std::string &parname, const std::vector<std::pair<std::string, std::string>> &vecAlgos, const IOUtils::OperationMode &mode,
                                                const char &rank) {
        for (size_t ii = 0; ii < vecAlgos.size(); ii++) {
            std::string algo_name(IOUtils::strToUpper(vecAlgos[ii].second));
            if (algo_name == "NONE")
                algo_name = "LINEAR"; // default value
            double max_gap_size = cfg.get(parname + "::" + algo_name + "::MAX_GAP_SIZE", interpol_section, IOUtils::nodata);
            if (max_gap_size != IOUtils::nodata) {
                if (max_gap_size < 0.)
                    throw IOException("MAX_GAP_SIZE not valid, it should be a duration in seconds at least greater than 0", AT);
                else
                    max_gap_size /= 86400.; // user uses seconds, internally julian day is used
            }

            if (first_time && algo_name == "ACCUMULATE" && mode == IOUtils::VSTATIONS && rank == 1) {
                const std::string vstations_refresh_rate = cfg.get("VSTATIONS_REFRESH_RATE", "InputEditing");
                const std::vector<std::pair<std::string, std::string>> vecArgs(1, std::make_pair("PERIOD", vstations_refresh_rate));
                std::shared_ptr<ResamplingAlgorithms> algo_ptr(ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, window_size, vecArgs));
                mapAlgorithms[parname].addAlgorithm(algo_ptr, max_gap_size);
            } else {
                std::vector<std::pair<std::string, std::string>> vecArgs(cfg.getArgumentsForAlgorithm(parname, algo_name));
                eraseArg(vecArgs, "MAX_GAP_SIZE");
                std::shared_ptr<ResamplingAlgorithms> algo_ptr(ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, window_size, vecArgs));
                mapAlgorithms[parname].addAlgorithm(algo_ptr, max_gap_size);
            }
        }
        if (mapAlgorithms[parname].empty()) { // create a default linear algorithm
            std::vector<std::pair<std::string, std::string>> vecArgs(cfg.getArgumentsForAlgorithm(parname, "LINEAR"));
            std::shared_ptr<ResamplingAlgorithms> algo_ptr(ResamplingAlgorithmsFactory::getAlgorithm("LINEAR", parname, window_size, vecArgs));
            mapAlgorithms[parname].addAlgorithm(algo_ptr, IOUtils::nodata);        
        }
    }

    void Meteo1DInterpolator::getWindowSize(ProcessingProperties &o_properties) const {
        o_properties.points_before = 1;
        o_properties.points_after = 1;
        o_properties.time_before = Duration(window_size, 0.); // we will need to cut a window 2x larger so we can interpolate each point in the window
        o_properties.time_after = Duration(window_size, 0.);
    }

    bool Meteo1DInterpolator::resampleData(const Date &date, const std::string &stationHash, const std::vector<MeteoData> &vecM, MeteoData &md) {
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
            const std::string parname(md.getNameForParameter(ii)); // Current parameter name
            const std::map<std::string, ResamplingStack>::const_iterator it = mapAlgorithms.find(parname);
            if (it != mapAlgorithms.end()) { // the parameter has been found
                it->second.resample(stationHash, index, elementpos, ii, vecM, md, window_size);

            } else { // we are dealing with an extra parameter, we need to add it to the map first, so it will exist next time...
                const std::vector<std::pair<std::string, std::string>> vecAlgos(cfg.getValues(parname + interpol_pattern, interpol_section));
                mapAlgorithms[parname] = ResamplingStack();

                processAlgorithms(false, parname, vecAlgos);
                mapAlgorithms[parname].resample(stationHash, index, elementpos, ii, vecM, md, window_size);
            }

            if ((index != IOUtils::npos) && vecM[index](ii) != md(ii)) {
                md.setResampledParam(ii);
                if (data_qa_logs) {
                    const std::map<std::string, ResamplingStack>::const_iterator it2 = mapAlgorithms.find(parname); // we have to re-find it in order to handle extra parameters
                    const std::string statName(md.meta.getStationName());
                    const std::string statID(md.meta.getStationID());
                    const std::string stat = (!statID.empty()) ? statID : statName;
                    const std::string algo_name(it2->second.getStackStr());
                    cout << "[DATA_QA] Resampling " << stat << "::" << parname << "::" << algo_name << " " << md.date.toString(Date::ISO_TZ) << " [" << md.date.toString(Date::ISO_WEEK) << "]\n";
                }
            }
        } // endfor ii

        return true; // successfull resampling
    }

    void Meteo1DInterpolator::resetResampling() {
        std::map<std::string, ResamplingStack>::iterator it;
        for (it = mapAlgorithms.begin(); it != mapAlgorithms.end(); ++it) {
            it->second.resetResampling();
        }
    }

    std::string Meteo1DInterpolator::getAlgorithmsForParameter(const std::string &parname) const {
        std::string algo("linear"); // default value
        cfg.getValue(parname + "::resample", "Interpolations1D", algo, IOUtils::nothrow);
        return algo;
    }

    Meteo1DInterpolator &Meteo1DInterpolator::operator=(const Meteo1DInterpolator &source) {
        if (this != &source) {
            mapAlgorithms = source.mapAlgorithms;
            window_size = source.window_size;
            enable_resampling = source.enable_resampling;
            data_qa_logs = source.data_qa_logs;
        }
        return *this;
    }

    const std::string Meteo1DInterpolator::toString() const {
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

    std::vector<std::shared_ptr<ResamplingAlgorithms>> ResamplingStack::buildStack(const ResamplingAlgorithms::gap_info &gap) const {
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

    void ResamplingStack::resetResampling() {
        for (size_t ii = 0; ii < stack.size(); ii++) {
            stack[ii]->resetResampling();
        }
    }

    void ResamplingStack::resample(const std::string &stationHash, const size_t &index, const ResamplingAlgorithms::ResamplingPosition elementpos, const size_t &ii, const std::vector<MeteoData> &vecM,
                                   MeteoData &md, const double &i_window_size) const {
        const ResamplingAlgorithms::gap_info gap = ResamplingAlgorithms::findGap(index, ii, vecM, md.date, i_window_size);
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

    bool ResamplingStack::empty() const {
        return stack.empty();
    }

    std::string ResamplingStack::getStackStr() const {
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
