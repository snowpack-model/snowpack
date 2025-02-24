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
#include <fstream>
#include <limits>
#include <sstream>
#include <cerrno>
#include <cstring>

#include <meteoio/MeteoProcessor.h>		//required for some generic initializations static methods
#include <meteoio/FileUtils.h>
#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/meteoFilters/FilterMin.h>
#include <meteoio/meteoFilters/FilterMax.h>
#include <meteoio/meteoFilters/FilterMinMax.h>
#include <meteoio/meteoFilters/FilterMinMaxConditional.h>
#include <meteoio/meteoFilters/FilterPotentialSW.h>
#include <meteoio/meteoFilters/FilterStdDev.h>
#include <meteoio/meteoFilters/FilterRate.h>
#include <meteoio/meteoFilters/FilterUnheatedPSUM.h>
#include <meteoio/meteoFilters/FilterTukey.h>
#include <meteoio/meteoFilters/FilterMAD.h>
#include <meteoio/meteoFilters/ProcAggregate.h>
#include <meteoio/meteoFilters/ProcIIR.h>
#include <meteoio/meteoFilters/ProcDeAccumulate.h>
#include <meteoio/meteoFilters/ProcUndercatch_WMO.h>
#include <meteoio/meteoFilters/ProcUndercatch_Forland.h>
#include <meteoio/meteoFilters/ProcUndercatch_Hamon.h>
#include <meteoio/meteoFilters/ProcPSUMDistribute.h>
#include <meteoio/meteoFilters/ProcUnventilatedT.h>
#include <meteoio/meteoFilters/ProcQuantileMapping.h>
#include <meteoio/meteoFilters/ProcShade.h>
#include <meteoio/meteoFilters/ProcAdd.h>
#include <meteoio/meteoFilters/ProcMult.h>
#include <meteoio/meteoFilters/ProcExpSmoothing.h>
#include <meteoio/meteoFilters/ProcWMASmoothing.h>
#include <meteoio/meteoFilters/ProcRadProj.h>
#include <meteoio/meteoFilters/ProcReducePressure.h>
#include <meteoio/meteoFilters/ProcRHWaterToIce.h>
#include <meteoio/meteoFilters/ProcTransformWindVector.h>
#include <meteoio/meteoFilters/ProcShift.h>
#include <meteoio/meteoFilters/FilterNoChange.h>
#include <meteoio/meteoFilters/FilterTimeconsistency.h>
#include <meteoio/meteoFilters/FilterDeGrass.h>
#include <meteoio/meteoFilters/TimeFilters.h>
#include <meteoio/meteoFilters/FilterDespikingPS.h>
#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/meteoFilters/FilterKalman.h>
#include <meteoio/meteoFilters/FilterMaths.h>

namespace mio {
/**
 * @page processing Processing overview
 * The pre-processing infrastructure is described in ProcessingBlock (for its API). The goal of this page is to give an overview of the available filters 
 * and processing elements and their usage. Moreover, there is a special mode of operation where MeteoIO writes on the screen a line for each data point that gets modified 
 * (either filtered, resampled or generated). This is enabled by setting the DATA_QA_LOGS key to *true* in the [General] section. The outputs then look like the following:
 * @code
 * [DATA_QA] Filtering WFJ2::ILWR::MIN_MAX 2016-01-18T09:00:00+01:00 [2016-W03-01]
 * [DATA_QA] Resampling WFJ1::VW_MAX::LINEAR 2015-10-12T13:00:00+01:00 [2015-W42-01]
 * @endcode
 *
 * @section processing_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * It is possible to disable a given Processing Element for specific stations, using the **exclude** or **only**
 * options followed by a list of station IDs (see example below). This is supported automatically by all Processing Elements. Several Processing Elements
 * take arguments describing a processing window (for example, FilterStdDev). In such a case, they take the window parameters arguments as
 * defined in WindowedFilter::setWindowFParams().
 * 
 * It is also possible to rectrict any filter to a specific set of time ranges, using the **when** option followed by a comma delimited list of
 * date intervals (represented by two ISO formatted dates seperated by ' - ') or individual dates (also ISO formatted).
 * 
 * The processing is normally defined per meteorological parameter (either from the \ref MeteoData::Parameters "list of standard parameter names" or 
 * any other name of your choice). A special kind of processing is available on the timestamps themselves and takes place before any 
 * other processing (see below in \ref processing_available "Available processing elements").
 *
 * It is possible to apply filters to only specific heights/numbers of a parameter, when they are given as TA@15..., either 
 * by specifying the ids for which the filter should be used (AT_HEIGHTS), or for which it should not be used (EXCLUDE_HEIGHTS).
 * Per default a filter will run on all available heights, so every TA will be applied the same filter to.
 * 
 * The heights are given as a list of heights similar to the station IDs, if you have default parameters and 
 * some that use a height, i.e. TA, TA@1, TA@10, the default parameter has no height, i.e. "NO" in the list to 
 * include or exclude.
 * 
 * @section processing_section Filtering section
 * The filters are specified for each parameter in the [Filters] section. This section contains
 * a list of the various meteo parameters (see MeteoData) with their associated choice of filtering algorithms and
 * optional parameters.The filters are applied serially, in the order they are given in. 
 * 
 * 
 * An example of such section is given below:
 * @code
 * [Filters]
 * TA::filter1   = min_max
 * TA::arg1::AT_HEIGHTS = NO 1 
 * TA::arg1::min = 230
 * TA::arg1::max = 330
 *
 * RH::filter1   = min_max
 * RH::arg1::EXCLUDE_HEIGHTS = 10
 * RH::arg1::min = -0.2
 * RH::arg1::max = 1.2
 * RH::filter2    = min_max
 * RH::arg2::soft = true
 * RH::arg2::min  = 0.0
 * RH::arg2::max  = 1.0
 *
 * PSUM::filter1    = min
 * PSUM::arg1::min  = -0.1
 * PSUM::filter2    = min
 * PSUM::arg2::soft = true
 * PSUM::arg2::min  = 0.
 * PSUM::filter3    = undercatch_wmo
 * PSUM::arg3::type = Hellmannsh
 * PSUM::arg3::exclude = DAV3 WFJ2
 * 
 * #Correct a wrongly mounted wind sensor for two time periods
 * DW::filter1   = add
 * DW::arg1::type = Cst
 * DW::arg1::cst = 45
 * DW::arg1::when = 2020-07-01 - 2020-07-10 , 2020-07-20T12:00 - 2020-08-01
 * @endcode
 * 
 * @note It is possible to turn off all meteo filtering by setting the *Enable_Meteo_Filtering* key to false in the [Filters] section; 
 * the same can be done for timestamps filtering with the *Enable_Time_Filtering* key.
 *
 * @section processing_available Available processing elements
 * New filters can easily be developed. The filters that are currently available are the following:
 * - NONE: this does nothing (this is useful in an \ref config_import "IMPORT" to overwrite previous filters);
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - MIN_MAX: range check filter, see FilterMinMax
 * - MIN_MAX_CONDITIONAL: range check only if a different parameter holds true to a comparison, see FilterMinMaxConditional
 * - RATE: rate of change filter, see FilterRate
 * - UNHEATED_RAINGAUGE: detection of snow melting in a rain gauge, see FilterUnheatedPSUM
 * - DETECT_GRASS: detection of grass growing under the snow height sensor, see FilterDeGrass
 * - POTENTIALSW: ensuring physically realistic incoming short wave radiation, see FilterPotentialSW
 * - MATHS: evaluating arithmetic expressions with access to meteo data, see FilterMaths
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterStdDev
 * - MAD: median absolute deviation, see FilterMAD
 * - TUKEY: Tukey53H spike detection, based on median, see FilterTukey
 * - DESPIKING: despiking in phase space according to Goring and Nikora (2002), see FilterDespikingPS
 * - NO_CHANGE: reject data that changes too little (low variance), see FilterNoChange
 * - TIME_CONSISTENCY: reject data that changes too much, see FilterTimeconsistency
 * - KALMAN: dynamic state likelihood estimation via Bayesian statistics, see FilterKalman
 * - PARTICLE: Monte Carlo sampling method for dynamic state estimation, see FilterParticle
 *
 * Some data transformations are also supported besides filtering, both very basic and generic data transformations:
 * - SUPPR: delete all or some data, see FilterSuppr
 * - ADD: adds a given offset to the data, see ProcAdd
 * - MULT: multiply the data by a given factor, see ProcMult
 * - SHIFT: shift a specific meteo parameter in time, see ProcShift
 * - QM: quantile mapping, see ProcQuantileMapping
 *
 * As well as more specific data transformations:
 * - AGGREGATE: various data aggregation algorithms, see ProcAggregate
 * - DEACCUMULATE: recompute instantaneous values from accumulated values, see ProcDeAccumulate
 * - EXP_SMOOTHING: exponential smoothing of data, see ProcExpSmoothing
 * - WMA_SMOOTHING: weighted moving average smoothing of data, see ProcWMASmoothing
 * - IIR: Low Pass or High Pass critically damped filter, see ProcIIR
 * - UNDERCATCH_WMO: WMO rain gauge correction for undercatch, using various correction models, see ProcUndercatch_WMO
 * - UNDERCATCH_FORLAND: Forland1996 rain gauge correction for solid and liquid undercatch, using various correction models, see ProcUndercatch_Forland
 * - UNDERCATCH_HAMON: Hamon1973 rain gauge correction for undercatch, see ProcUndercatch_Hamon
 * - UNVENTILATED_T: unventilated temperature sensor correction, see ProcUnventilatedT
 * - PSUM_DISTRIBUTE: distribute accumulated precipitation over preceeding timesteps, see ProcPSUMDistribute
 * - SHADE: apply a shading mask to the Incoming or Reflected Short Wave Radiation, see ProcShade
 * - REDUCE_PRESSURE: reduce local pressure to sea level pressure, see ProcReducePressure
 * - RHWATERTOICE: correct relative humidity over water to over ice in case temperature is below freezing, see ProcRHWaterToIce
 * - TRANSFORMWINDVECTOR: transform wind direction and/or wind speed components, see ProcTransformWindVector
 * - PROJECT_RADIATION: project short wave radiation from one inclination to another one, see RadProj
 *
 * A few filters can be applied to the timestamps themselves:
 * - SUPPR: delete whole timesteps (based on a list or other criteria such as removing duplicates, etc), see TimeSuppr
 * - SHIFT: add offsets at specific times, for example to bring timestamps that contain Daylight Saving Time back to Winter time, see TimeShift
 * - SORT: sort the timestamps in increasing order, see TimeSort
 * - TIMELOOP: loop over a specific time period (for example for model spin-ups), see TimeLoop
 */

const double ProcessingBlock::default_height(IOUtils::nodata); // TODO: change if negative heights are supported

ProcessingBlock* BlockFactory::getBlock(const std::string& blockname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	//the indenting is a little weird, this is in order to show the same groups as in the documentation above

	//normal filters
	if (blockname == "MIN"){
		return new FilterMin(vecArgs, blockname, cfg);
	} else if (blockname == "MAX"){
		return new FilterMax(vecArgs, blockname, cfg);
	} else if (blockname == "MIN_MAX"){
		return new FilterMinMax(vecArgs, blockname, cfg);
	} else if (blockname == "MIN_MAX_CONDITIONAL"){
		return new FilterMinMaxConditional(vecArgs, blockname, cfg);
	} else if (blockname == "RATE"){
		return new FilterRate(vecArgs, blockname, cfg);
	} else if (blockname == "STD_DEV"){
		return new FilterStdDev(vecArgs, blockname, cfg);
	} else if (blockname == "MAD"){
		return new FilterMAD(vecArgs, blockname, cfg);
	} else if (blockname == "TUKEY"){
		return new FilterTukey(vecArgs, blockname, cfg);
	} else if (blockname == "UNHEATED_RAINGAUGE"){
		return new FilterUnheatedPSUM(vecArgs, blockname, cfg);
	} else if (blockname == "NO_CHANGE"){
		return new FilterNoChange(vecArgs, blockname, cfg);
	} else if (blockname == "TIME_CONSISTENCY"){
		return new FilterTimeconsistency(vecArgs, blockname, cfg);
	} else if (blockname == "DETECT_GRASS"){
		return new FilterDeGrass(vecArgs, blockname, cfg);
	} else if (blockname == "POTENTIALSW"){
		return new FilterPotentialSW(vecArgs, blockname, cfg);
	} else if (blockname == "DESPIKING"){
		return new FilterDespikingPS(vecArgs, blockname, cfg);
	} else if (blockname == "PARTICLE"){
		return new FilterParticle(vecArgs, blockname, cfg);
	} else if (blockname == "KALMAN"){
		return new FilterKalman(vecArgs, blockname, cfg);
	} else if (blockname == "MATHS"){
		return new FilterMaths(vecArgs, blockname, cfg);
	}

	//general data transformations
	else if (blockname == "SUPPR"){
		return new FilterSuppr(vecArgs, blockname, cfg);
	} else if (blockname == "ADD"){
		return new ProcAdd(vecArgs, blockname, cfg);
	} else if (blockname == "MULT"){
		return new ProcMult(vecArgs, blockname, cfg);
	} else if (blockname == "QM"){
		return new ProcQuantileMapping(vecArgs, blockname, cfg);
	} else if (blockname == "SHIFT"){
		return new ProcShift(vecArgs, blockname, cfg);
	} 

	//more specific data transformations
	else if (blockname == "EXP_SMOOTHING"){
		return new ProcExpSmoothing(vecArgs, blockname, cfg);
	} else if (blockname == "WMA_SMOOTHING"){
		return new ProcWMASmoothing(vecArgs, blockname, cfg);
	} else if (blockname == "IIR"){
		return new ProcIIR(vecArgs, blockname, cfg);
	} else if (blockname == "AGGREGATE"){
		return new ProcAggregate(vecArgs, blockname, cfg);
	} else if (blockname == "DEACCUMULATE"){
		return new ProcDeAccumulate(vecArgs, blockname, cfg);
	} else if (blockname == "UNDERCATCH_WMO"){
		return new ProcUndercatch_WMO(vecArgs, blockname, cfg);
	} else if (blockname == "UNDERCATCH_FORLAND"){
		return new ProcUndercatch_Forland(vecArgs, blockname, cfg);
	} else if (blockname == "UNDERCATCH_HAMON"){
		return new ProcUndercatch_Hamon(vecArgs, blockname, cfg);
	} else if (blockname == "UNVENTILATED_T"){
		return new ProcUnventilatedT(vecArgs, blockname, cfg);
	} else if (blockname == "PROJECT_RADIATION"){
		return new RadProj(vecArgs, blockname, cfg);
	} else if (blockname == "PSUM_DISTRIBUTE"){
		return new ProcPSUMDistribute(vecArgs, blockname, cfg);
	} else if (blockname == "SHADE"){
		return new ProcShade(vecArgs, blockname, cfg);
	} else if (blockname == "REDUCE_PRESSURE"){
		return new ProcReducePressure(vecArgs, blockname, cfg);
	} else if (blockname == "RHWATERTOICE"){
		return new ProcRHWaterToIce(vecArgs, blockname, cfg);
	} else if (blockname == "TRANSFORMWINDVECTOR"){
		return new ProcTransformWindVector(vecArgs, blockname, cfg);
	} else {
		throw IOException("The processing block '"+blockname+"' does not exist! " , AT);
	}
}

ProcessingBlock* BlockFactory::getTimeBlock(const std::string& blockname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	if (blockname == "SUPPR"){
		return new TimeSuppr(vecArgs, blockname, cfg);
	} else if (blockname == "SHIFT"){
		return new TimeShift(vecArgs, blockname, cfg);
	} else if (blockname == "SORT"){
		return new TimeSort(vecArgs, blockname, cfg);
	} else if (blockname == "TIMELOOP"){
		return new TimeLoop(vecArgs, blockname, cfg);
	} else {
		if (blockname == "UNDST") std::cerr << "[E] time filter UnDST has been renamed into Shift\n";
		throw IOException("The processing block '"+blockname+"' does not exist for the TIME parameter! " , AT);
	}
}

ProcessingBlock::ProcessingBlock(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
                            : excluded_stations( MeteoProcessor::initStationSet(vecArgs, "EXCLUDE") ), kept_stations( MeteoProcessor::initStationSet(vecArgs, "ONLY") ), 
                              time_restrictions( MeteoProcessor::initTimeRestrictions(vecArgs, "WHEN", "Filters::"+name, cfg.get("TIME_ZONE", "Input")) ), included_heights(), excluded_heights(), all_heights(true), properties(), block_name(name) {
	initHeightRestrictions(vecArgs);
}

void ProcessingBlock::initHeightRestrictions(const std::vector<std::pair<std::string, std::string>> vecArgs) {
	std::set<std::string> results_include;
	std::set<std::string> results_exclude;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first == "AT_HEIGHTS") {
			std::istringstream iss(vecArgs[ii].second);
			std::string word;
			while (iss >> word){
				results_include.insert(word);
			}
		}
		else if (vecArgs[ii].first == "EXCLUDE_HEIGHTS") {
			std::istringstream iss(vecArgs[ii].second);
			std::string word;
			while (iss >> word){
				results_exclude.insert(word);
			}
		}
	}

	if (!results_exclude.empty() || !results_include.empty()) {
		all_heights = false;
	}

	for (const auto& h_str : results_include) {
		double height;
		if (h_str == "NO") {
			height = default_height;
		} else {
			height = std::stod(h_str);
		}

		if (!included_heights.insert(height).second) { 
			throw IOException("Could not parse height of " + h_str + " for filter " + block_name);
		}
	}

	for (const auto& h_str : results_exclude) {
		double height;
		if (h_str == "NO") {
			height = default_height;
		} else {
			height = std::stod(h_str);
		}

		if (!excluded_heights.insert(height).second) { 
			throw IOException("Could not parse height of " + h_str + " for filter " + block_name);
		}
	}
}

/**
 * @brief Should the provided station be skipped in the processing?
 * @param[in] station_id stationID to test
 * @return true if the startion should be skipped, false otherwise
 */
bool ProcessingBlock::skipStation(const std::string& station_id) const
{
	if (excluded_stations.count(station_id)!=0) return true; //the station is in the excluded list -> skip
	if (kept_stations.empty()) return false; //there are no kept stations -> do not skip

	return (kept_stations.count(station_id)==0);
}

/**
 * @brief Should the provided height be skipped in the processing?
 * @param[in] height height to test
 * @return true if the startion should be skipped, false otherwise
 */
bool ProcessingBlock::skipHeight(const double& height) const {

	if (all_heights) return false;

	if (excluded_heights.count(height) != 0) return true;
	if (included_heights.empty()) return false;

	return (included_heights.count(height) == 0);
}

/**
 * @brief Read a data file structured as X Y value on each lines
 * @param[in] filter Calling filter name for error reporting
 * @param[in] filename file and path to open and read the data from
 * @param[out] X vector of X values
 * @param[out] Y vector of Y values
 */
void ProcessingBlock::readCorrections(const std::string& filter, const std::string& filename, std::vector<double> &X, std::vector<double> &Y)
{
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	try {
		size_t lcount=0;
		double xvalue, yvalue;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			iss >> std::skipws >> xvalue;
			if (!iss) {
				std::ostringstream ss;
				ss << "Invalid index in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			iss >> std::skipws >> yvalue;
			if ( iss.fail() ) {
				std::ostringstream ss;
				ss << "Invalid value in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			
			X.push_back( xvalue );
			Y.push_back( yvalue );
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}
}

/**
 * @brief Read a data file structured as X Y1 Y2 value on each lines
 * @param[in] filter Calling filter name for error reporting
 * @param[in] filename file and path to open and read the data from
 * @param[out] X vector of X values
 * @param[out] Y1 vector of Y1 values
 * @param[out] Y2 vector of Y2 values
 */
void ProcessingBlock::readCorrections(const std::string& filter, const std::string& filename, std::vector<double> &X, std::vector<double> &Y1, std::vector<double> &Y2)
{
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	try {
		size_t lcount=0;
		double xvalue, y1value, y2value;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			iss >> std::skipws >> xvalue;
			if (!iss) {
				std::ostringstream ss;
				ss << "Invalid index in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			iss >> std::skipws >> y1value;
			if ( iss.fail() ) {
				std::ostringstream ss;
				ss << "Invalid value in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			
			iss >> std::skipws >> y2value;
			if ( iss.fail() ) {
				std::ostringstream ss;
				ss << "Invalid value in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			
			X.push_back( xvalue );
			Y1.push_back( y1value );
			Y2.push_back( y2value );
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}
}

/**
 * @brief Read a correction file applicable to repeating time period.
 * @details This reads corrections to apply to repeating time periods, such as hours, dyas, months or years. 
 * Not all values must be provided as a default initial value is set. Depending on the time period, a check on
 * the index range is performed (hours must be <=24, days <=366, months <=12, years have no maximum).
 * @param[in] filter Calling filter name for error reporting
 * @param[in] filename file and path to open and read the data from
 * @param[in] col_idx column to read the data from (column 1 contains the timestamps)
 * @param[in] c_type expected time period
 * @param[in] init default value to initalize the results
 * @return a vector of corrections as read from the provided file, the indices matching the indices of the chosen time period
 */
std::vector<double> ProcessingBlock::readCorrections(const std::string& filter, const std::string& filename, const size_t& col_idx, const char& c_type, const double& init)
{
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	size_t maxIndex = 0;
	const size_t minIndex = (c_type=='h')? 0 : 1;
	if (c_type=='y') maxIndex = IOUtils::npos;
	else if (c_type=='m') maxIndex = 12;
	else if (c_type=='d') maxIndex = 366;
	else if (c_type=='h') maxIndex = 24;
	std::vector<double> corrections(maxIndex, init);

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

	try {
		size_t index, lcount=0;
		double value;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			iss >> std::skipws >> index;
			if ( !iss || index<minIndex || (index-minIndex)>=maxIndex) {
				std::ostringstream ss;
				ss << "Invalid index in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			size_t ii=2;
			do {
				iss >> std::skipws >> value;
				if ( iss.fail() ){
					std::ostringstream ss;
					ss << "In file " << filename << " at line " << lcount;
					if (!iss.eof())
						ss << ": invalid value";
					else
						ss << ": trying to read column " << col_idx << " of " << ii-1 << " columns";
					throw InvalidArgumentException(ss.str(), AT);
				}
			} while ((ii++) < col_idx);
			corrections.at( index-minIndex ) = value;
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}

	return corrections;
}

/**
 * @brief Read a correction file, ie a file structured as timestamps followed by values on each lines
 * @param[in] filter Calling filter name for error reporting
 * @param[in] filename file and path to open and read the data from
 * @param[in] TZ default timezone for the timestamps
 * @param[in] col_idx column to read the data from (column 1 contains the timestamps)
 * @return a vector of <Date, value> as read from the provided file
 */
std::vector<ProcessingBlock::offset_spec> ProcessingBlock::readCorrections(const std::string& filter, const std::string& filename, const double& TZ, const size_t& col_idx)
{
	if (col_idx<2)
		throw InvalidArgumentException("Filter "+filter+": the column index must be greater than 1!", AT);

	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	std::vector<offset_spec> corrections;

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

	try {
		size_t lcount=0;
		double value;
		Date date;
		std::string tmp;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss >> std::skipws >> tmp;
			const bool status = IOUtils::convertString(date, tmp, TZ);
			if ( !iss || !status) {
				std::ostringstream ss;
				ss << "Invalid date in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			size_t ii=2;
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			do {
				iss >> std::skipws >> value;
				if ( iss.fail() ){
					std::ostringstream ss;
					ss << "In file " << filename << " at line " << lcount;
					if (!iss.eof())
						ss << ": invalid value";
					else
						ss << ": trying to read column " << col_idx << " of " << ii-1 << " columns";
					throw InvalidArgumentException(ss.str(), AT);
				}
			} while ((ii++) < col_idx);
			corrections.push_back( offset_spec(date, value) );
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}

	std::sort(corrections.begin(), corrections.end());
	return corrections;
}

/**
 * @brief Read a list of date ranges by stationIDs from a file
 * @details Each station ID read in the provided file is attributed a list fo date ranges. For example, the file
 * to read the data from can contain:
 * @code
 * *WFJ 2015-11-10T06:00
 * *WFJ 2015-12-25T01:00 2015-12-27T13:30
 * @endcode
 * All these time ranges will populate a vector that will be mapped to the "*WFJ" station ID.
 * 
 * @param[in] filter Calling filter name for error reporting
 * @param[in] filename file and path to open and read the data from
 * @param[in] TZ default timezone for the timestamps
 * @return a map of <stationID, std::vector<DateRange>> as read from the provided file
 */
std::map< std::string, std::vector<DateRange> > ProcessingBlock::readDates(const std::string& filter, const std::string& filename, const double& TZ)
{
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename, AT);
	if (!FileUtils::fileExists(filename)) throw NotFoundException(filename, AT);

	std::ifstream fin(filename.c_str());
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	std::map< std::string, std::vector<DateRange> > dates_specs;

	Date d1, d2;
	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::vector<std::string> vecString;
			const size_t nrElems = IOUtils::readLineToVec(line, vecString);
			if (nrElems<2)
				throw InvalidFormatException("Invalid syntax for filter " + filter + " in file \"" + filename + "\": expecting at least 2 arguments", AT);

			const std::string station_ID( vecString[0] );
			if (nrElems==2) {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				const DateRange range(d1, d1);
				dates_specs[ station_ID ].push_back( range );
			} else if (nrElems==3) {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[2], TZ))
					throw InvalidFormatException("Could not process date "+vecString[2]+" in file \""+filename+"\"", AT);
				const DateRange range(d1, d2);
				dates_specs[ station_ID ].push_back( range );
			} else if (nrElems==4 && vecString[2]=="-") {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[3], TZ))
					throw InvalidFormatException("Could not process date "+vecString[3]+" in file \""+filename+"\"", AT);
				const DateRange range(d1, d2);
				dates_specs[ station_ID ].push_back( range );
			} else
				throw InvalidFormatException("Unrecognized syntax in file \""+filename+"\": '"+line+"'\n", AT);
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}

	//sort all the suppr_specs
	std::map< std::string, std::vector<DateRange> >::iterator station_it( dates_specs.begin() );
	for (; station_it!=dates_specs.end(); ++station_it) {
		std::sort(station_it->second.begin(), station_it->second.end());
	}
	return dates_specs;
}


void ProcessingBlock::extract_dbl_vector(const unsigned int& param, const std::vector<MeteoData>& ivec,
                                     std::vector<double>& ovec)
{
	ovec.resize( ivec.size() );
	for (size_t ii=0; ii<ivec.size(); ii++) {
		ovec[ii] = ivec[ii](param);
	}
}

void ProcessingBlock::extract_dbl_vector(const unsigned int& param, const std::vector<const MeteoData*>& ivec,
                                     std::vector<double>& ovec)
{
	ovec.resize( ivec.size() );
	for (size_t ii=0; ii<ivec.size(); ii++) {
		ovec[ii] =  (*ivec[ii])(param);
	}
}

const std::string ProcessingBlock::toString() const {
	std::ostringstream os;
	os << "[" << block_name << " ";
	os << properties.toString();
	os << "]";
	return os.str();
}

const std::string ProcessingProperties::toString() const
{
	std::ostringstream os;
	const double h_before = time_before.getJulian()*24.;
	const double h_after = time_after.getJulian()*24.;
	const size_t p_before = points_before;
	const size_t p_after = points_after;

	os << "{";
	if(h_before>0. || h_after>0.) os << "-" << h_before << " +" << h_after << " h; ";
	if(p_before>0 || p_after>0) os << "-" << p_before << " +" << p_after << " pts; ";
	if(stage==ProcessingProperties::first)
		os << "p¹";
	if(stage==ProcessingProperties::second)
		os << "p²";
	if(stage==ProcessingProperties::both)
		os << "p½";
	os << "}";
	return os.str();
}

} //end namespace
