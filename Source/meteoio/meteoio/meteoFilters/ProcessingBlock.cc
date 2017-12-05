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
#include <sstream>
#include <errno.h>
#include <cstring>

#include <meteoio/FileUtils.h>
#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/meteoFilters/FilterMin.h>
#include <meteoio/meteoFilters/FilterMax.h>
#include <meteoio/meteoFilters/FilterMinMax.h>
#include <meteoio/meteoFilters/FilterPotentialSW.h>
#include <meteoio/meteoFilters/FilterStdDev.h>
#include <meteoio/meteoFilters/FilterRate.h>
#include <meteoio/meteoFilters/FilterUnheatedPSUM.h>
#include <meteoio/meteoFilters/FilterTukey.h>
#include <meteoio/meteoFilters/FilterMAD.h>
#include <meteoio/meteoFilters/ProcAggregate.h>
#include <meteoio/meteoFilters/ProcIIR.h>
#include <meteoio/meteoFilters/ProcUndercatch_WMO.h>
#include <meteoio/meteoFilters/ProcUndercatch_Forland.h>
#include <meteoio/meteoFilters/ProcUndercatch_Hamon.h>
#include <meteoio/meteoFilters/ProcPSUMDistribute.h>
#include <meteoio/meteoFilters/ProcUnventilatedT.h>
#include <meteoio/meteoFilters/ProcShade.h>
#include <meteoio/meteoFilters/ProcAdd.h>
#include <meteoio/meteoFilters/ProcMult.h>
#include <meteoio/meteoFilters/ProcExpSmoothing.h>
#include <meteoio/meteoFilters/ProcWMASmoothing.h>
#include <meteoio/meteoFilters/FilterNoChange.h>
#include <meteoio/meteoFilters/FilterTimeconsistency.h>
#include <meteoio/meteoFilters/FilterDeGrass.h>
#include <meteoio/meteoFilters/TimeFilters.h>

namespace mio {
/**
 * @page processing Processing overview
 * The pre-processing infrastructure is described in ProcessingBlock (for its API). The goal of this page is to give an overview of the available filters and processing elements and their usage.
 *
 * @section processing_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * It is possible to disable a given Processing Element for specific stations, using the **exclude** or **only**
 * options followed by a list of station IDs (see example below). This is supported automatically by all Processing Elements. Several Processing Elements
 * take arguments describing a processing window (for example, FilterStdDev). In such a case, they take the window parameters arguments as
 * defined in WindowedFilter::setWindowFParams().
 *
 * @section processing_section Filtering section
 * The filters are specified for each parameter in the [Filters] section. This section contains
 * a list of the various meteo parameters (see MeteoData) with their associated choice of filtering algorithms and
 * optional parameters.The filters are applied serialy, in the order they are given in. An example of such section is given below:
 * @code
 * [Filters]
 * TA::filter1   = min_max
 * TA::arg1::min = 230
 * TA::arg1::max = 330
 *
 * RH::filter1   = min_max
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
 * @endcode
 *
 * @section processing_available Available processing elements
 * New filters can easily be developed. The filters that are currently available are the following:
 * - NONE: this does nothing (this is useful in an \ref config_import "IMPORT" to overwritte previous filters);
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - MIN_MAX: range check filter, see FilterMinMax
 * - RATE: rate of change filter, see FilterRate
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterStdDev
 * - MAD: median absolute deviation, see FilterMAD
 * - TUKEY: Tukey53H spike detection, based on median, see FilterTukey
 * - UNHEATED_RAINGAUGE: detection of snow melting in a rain gauge, see FilterUnheatedPSUM
 * - NO_CHANGE: reject data that changes too little (low variance), see FilterNoChange
 * - TIME_CONSISTENCY: reject data that changes too much , see FilterTimeconsistency
 * - DETECT_GRASS: detection of grass growing under the snow height sensor, see FilterDeGrass
 * - POTENTIALSW: ensuring physically realistic incoming short wave radiation, see FilterPotentialSW
 *
 * Some data transformations are also supported besides filtering, both very basic and generic data transformations:
 * - SUPPR: delete all or some data, see FilterSuppr
 * - ADD: adds a given offset to the data, see ProcAdd
 * - MULT: multiply the data by a given factor, see ProcMult
 *
 * As well as more specific data transformations:
 * - EXP_SMOOTHING: exponential smoothing of data, see ProcExpSmoothing
 * - WMA_SMOOTHING: weighted moving average smoothing of data, see ProcWMASmoothing
 * - IIR: Low Pass or High Pass critically damped filter, see ProcIIR
 * - AGGREGATE: various data aggregation algorithms, see ProcAggregate
 * - UNDERCATCH_WMO: WMO rain gauge correction for undercatch, using various correction models, see ProcUndercatch_WMO
 * - UNDERCATCH_FORLAND: Forland1996 rain gauge correction for solid and liquid undercatch, using various correction models, see ProcUndercatch_Forland
 * - UNDERCATCH_HAMON: Hamon1973 rain gauge correction for undercatch, see ProcUndercatch_Hamon
 * - UNVENTILATED_T: unventilated temperature sensor correction, see ProcUnventilatedT
 * - PSUM_DISTRIBUTE: distribute accumulated precipitation over preceeding timesteps, see ProcPSUMDistribute
 * - SHADE: apply a shading mask to the Incoming or Reflected Short Wave Radiation, see ProcShade
 * 
 * A few filters can be applied to the timestamps themselves:
 * - SUPPR: delete whole timesteps, see TimeSuppr
 * - UNDST: correct timestamps that contain Daylight Saving Time back to Winter time, see TimeUnDST
 */

ProcessingBlock* BlockFactory::getBlock(const std::string& blockname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	//the indenting is a little weird, this is in order to show the same groups as in the documentation above
	
	//normal filters
	if (blockname == "MIN"){
		return new FilterMin(vecArgs, blockname);
	} else if (blockname == "MAX"){
		return new FilterMax(vecArgs, blockname);
	} else if (blockname == "MIN_MAX"){
		return new FilterMinMax(vecArgs, blockname);
	} else if (blockname == "RATE"){
		return new FilterRate(vecArgs, blockname);
	} else if (blockname == "STD_DEV"){
		return new FilterStdDev(vecArgs, blockname);
	} else if (blockname == "MAD"){
		return new FilterMAD(vecArgs, blockname);
	} else if (blockname == "TUKEY"){
		return new FilterTukey(vecArgs, blockname);
	} else if (blockname == "UNHEATED_RAINGAUGE"){
		return new FilterUnheatedPSUM(vecArgs, blockname);
	} else if (blockname == "NO_CHANGE"){
		return new FilterNoChange(vecArgs, blockname);
	} else if (blockname == "TIME_CONSISTENCY"){
		return new FilterTimeconsistency(vecArgs, blockname);
	} else if (blockname == "DETECT_GRASS"){
		return new FilterDeGrass(vecArgs, blockname);
	} else if (blockname == "POTENTIALSW"){
		return new FilterPotentialSW(vecArgs, blockname);
	}
	
	//general data transformations
	else if (blockname == "SUPPR"){
		return new FilterSuppr(vecArgs, blockname, cfg.getConfigRootDir(), cfg.get("TIME_ZONE", "Input"));
	} else if (blockname == "ADD"){
		return new ProcAdd(vecArgs, blockname, cfg.getConfigRootDir());
	} else if (blockname == "MULT"){
		return new ProcMult(vecArgs, blockname, cfg.getConfigRootDir());
	}
	
	//more specific data transformations
	else if (blockname == "EXP_SMOOTHING"){
		return new ProcExpSmoothing(vecArgs, blockname);
	} else if (blockname == "WMA_SMOOTHING"){
		return new ProcWMASmoothing(vecArgs, blockname);
	} else if (blockname == "IIR"){
		return new ProcIIR(vecArgs, blockname);
	} else if (blockname == "AGGREGATE"){
		return new ProcAggregate(vecArgs, blockname);
	} else if (blockname == "UNDERCATCH_WMO"){
		return new ProcUndercatch_WMO(vecArgs, blockname);
	} else if (blockname == "UNDERCATCH_FORLAND"){
		return new ProcUndercatch_Forland(vecArgs, blockname);
	} else if (blockname == "UNDERCATCH_HAMON"){
		return new ProcUndercatch_Hamon(vecArgs, blockname);
	} else if (blockname == "UNVENTILATED_T"){
		return new ProcUnventilatedT(vecArgs, blockname);
	} else if (blockname == "PSUM_DISTRIBUTE"){
		return new ProcPSUMDistribute(vecArgs, blockname);
	} else if (blockname == "SHADE"){
		return new ProcShade(vecArgs, blockname, cfg);
	} else {
		throw IOException("The processing block '"+blockname+"' does not exist! " , AT);
	}
}

ProcessingBlock* BlockFactory::getTimeBlock(const std::string& blockname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	if (blockname == "SUPPR"){
		return new TimeSuppr(vecArgs, blockname, cfg.getConfigRootDir(), cfg.get("TIME_ZONE", "Input"));
	} else if (blockname == "UNDST"){
		return new TimeUnDST(vecArgs, blockname, cfg.getConfigRootDir(), cfg.get("TIME_ZONE", "Input"));
	} else {
		throw IOException("The processing block '"+blockname+"' does not exist for the TIME parameter! " , AT);
	}
}

const double ProcessingBlock::soil_albedo = .23; //grass
const double ProcessingBlock::snow_albedo = .85; //snow
const double ProcessingBlock::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

ProcessingBlock::ProcessingBlock(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
                            : excluded_stations( initStationSet(vecArgs, "EXCLUDE") ), kept_stations( initStationSet(vecArgs, "ONLY") ), properties(), block_name(name) {}

std::set<std::string> ProcessingBlock::initStationSet(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& keyword)
{
	std::set<std::string> results;
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first==keyword) {
			std::istringstream iss(vecArgs[ii].second);
			std::string word;
			while (iss >> word){
				results.insert(word);
			}
		}
	}

	return results;
}

bool ProcessingBlock::skipStation(const std::string& station_id) const
{
	if (excluded_stations.count(station_id)!=0) return true; //the station is in the excluded list -> skip
	if (kept_stations.empty()) return false; //there are no kept stations -> do not skip

	return (kept_stations.count(station_id)==0);
}

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
	if (c_type=='m') maxIndex = 12;
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
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
	
	return corrections;
}

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
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
	
	std::sort(corrections.begin(), corrections.end());
	return corrections;
}

std::map< std::string, std::vector<ProcessingBlock::dates_range> > ProcessingBlock::readDates(const std::string& filter, const std::string& filename, const double& TZ)
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
	std::map< std::string, std::vector<dates_range> > dates_specs;
	
	Date d1, d2;
	try {
		size_t lcount=0;
		do {
			lcount++;
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
				const dates_range range(d1, d1);
				dates_specs[ station_ID ].push_back( range );
			} else if (nrElems==3) {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[2], TZ))
					throw InvalidFormatException("Could not process date "+vecString[2]+" in file \""+filename+"\"", AT);
				const dates_range range(d1, d2);
				dates_specs[ station_ID ].push_back( range );
			} else if (nrElems==4 && vecString[2]=="-") {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[3], TZ))
					throw InvalidFormatException("Could not process date "+vecString[3]+" in file \""+filename+"\"", AT);
				const dates_range range(d1, d2);
				dates_specs[ station_ID ].push_back( range );
			} else
				throw InvalidFormatException("Unrecognized syntax in file \""+filename+"\": '"+line+"'\n", AT);
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}

	//sort all the suppr_specs
	std::map< std::string, std::vector<dates_range> >::iterator station_it( dates_specs.begin() );
	for (; station_it!=dates_specs.end(); ++station_it) {
		std::sort(station_it->second.begin(), station_it->second.end());
	}
	return dates_specs;
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
