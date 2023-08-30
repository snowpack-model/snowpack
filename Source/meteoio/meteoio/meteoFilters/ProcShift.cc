// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcShift.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/FileUtils.h>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <algorithm>

using namespace std;

namespace mio {

ProcShift::ProcShift(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
          : ProcessingBlock(vecArgs, name, cfg), corrections(), ref_param(), root_path(cfg.getConfigRootDir()), offsets_file(), 
            cst_offset(IOUtils::nodata), sampling_rate(IOUtils::nodata), offset_range(1.), width_d(2.), width_idx(0), 
            offsets_interp(cst), extract_offsets(false)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void ProcShift::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ivec.empty()) return;
	
	if (extract_offsets) { //compute the offsets between two parameters
		writeOffsets(param, ivec);
	} else { //apply a correction file or a fixed correction
		correctOffsets(param, ovec);
	}
}

void ProcShift::writeOffsets(const unsigned int& param, const std::vector<MeteoData>& ivec)
{
	//prepare the parameters that we need: reference parameter index, clean sampling rate and vectors with the two parameters
	const size_t param_ref = ivec[0].getParameterIndex(ref_param);
	if (sampling_rate==IOUtils::nodata) {
		sampling_rate = getMedianSampling(param_ref, ivec);
		//rounding the sampling rate to 1‰ in seconds so the resampled timestamps don't accumulate rounding errors
		const int sampling_decimal_range = static_cast<int>(floor( log10(sampling_rate*24.*3600.) ));
		const double rounding_precision = pow(10, sampling_decimal_range-3);
		sampling_rate = round(sampling_rate*24.*3600./rounding_precision) * rounding_precision / (3600.*24.);
	}
	width_idx = static_cast<size_t>( round( width_d / sampling_rate ) );
	const std::vector<offset_spec> vecX( resampleVector(ivec, param_ref) );
	const std::vector<offset_spec> vecY( resampleVector(ivec, param) );
	
	//check and open the output file that will contain the extracted offsets
	if (!FileUtils::validFileAndPath(offsets_file)) throw AccessException("Invalid file name \""+offsets_file+"\"", AT);
	errno = 0;
	
	//if the filter runs multiple times (because of rebuffering for example), we need to append the results
	//so the output file has been deleted in the constructor
	std::ofstream fout(offsets_file.c_str(), std::ofstream::out | std::ios_base::app | ios::binary);
	if (fout.fail())
		throw AccessException("Error opening file \"" + offsets_file + "\" for writing, possible reason: " + std::string(std::strerror(errno)), AT);
	
	//the offsets are written in seconds as corrections that should be applied to re-synchronize the data
	for (size_t ii=0; ii<vecX.size(); ii++) {
		const int offset = getOffset(vecX, vecY, ii);
		if (offset!=IOUtils::inodata)
			fout << vecX[ii].date.toString(Date::ISO) << " " << -offset*sampling_rate*24.*3600. << "\n";
	}
	
	fout.close();
}

//the strategy is: we created a new vector<MeteoData> that only contains the shifted parameter
//and then we merge it with ovec
void ProcShift::correctOffsets(const unsigned int& param, std::vector<MeteoData>& ovec)
{
	const std::string paramname( ovec.front().getNameForParameter(param) );
	std::vector<MeteoData> shifted_param;
	shifted_param.reserve( ovec.size() );
	
	if (offsets_interp==cst) {
		for (size_t ii=0; ii<ovec.size(); ii++) {
			MeteoData md( ovec[ii].date + cst_offset );
			md.addParameter(paramname);
			md(paramname) = ovec[ii](param);
			shifted_param.push_back( md );
			ovec[ii](param) = IOUtils::nodata;
		}
	} else {
		const size_t nCorr = corrections.size();
		if (nCorr==0) return;
		
		if (offsets_interp==stepwise) {
			size_t jj=0; //index for the corrections
			double offset = 0.;
			
			for (size_t ii=0; ii<ovec.size(); ii++) {
				const Date curr_date(ovec[ii].date);
				
				//update the offset?
				if (jj<nCorr && corrections[jj].date <= curr_date) {
					//find the next user provided offset
					while (jj<nCorr && corrections[jj].date <= curr_date) jj++;
					const double last_offset = offset;
					if (jj>0) offset = corrections[jj-1].value / (3600.*24.);
					
					//remove the points that would otherwise be duplicates
					//this happens for example if first moving data forward by 1 hour
					//and then moving data backward by 2 hours -> there is an overlap
					if (offset<last_offset) {
						const Date last_valid( curr_date + offset );
						size_t kk = shifted_param.size() - 1;
						while (kk>0 && shifted_param[kk].date>last_valid) kk--;
						shifted_param.erase( shifted_param.begin()+kk, shifted_param.end() );
					}
				}
				
				MeteoData md( curr_date + offset );
				md.addParameter(paramname);
				md(paramname) = ovec[ii](param);
				shifted_param.push_back( md );
				ovec[ii](param) = IOUtils::nodata;
			}
		}
	}
	
	//now merge back the results
	MeteoData::mergeTimeSeries(ovec, shifted_param, MeteoData::WINDOW_MERGE, MeteoData::CONFLICTS_PRIORITY_LAST);
}

bool ProcShift::isAllNodata(const std::vector<offset_spec>& vecX, const size_t& startIdx, const size_t& endIdx)
{
	for (size_t ii=startIdx; ii<endIdx; ii++) {
		if (vecX[ii].value!=IOUtils::nodata) return false;
	}
	
	return true;
}

/**
 * @brief Get the median sampling rate for a given parameter with a meteo timeseries
 * @param param Meteorological parameter to consider
 * @param ivec meteorological timeseries
 * @return median sampling rate in days
 */
double ProcShift::getMedianSampling(const size_t& param, const std::vector<MeteoData>& ivec)
{
	std::vector<double> vecSampling;
	
	for (size_t ii=0; ii<(ivec.size()-1); ii++){
		if (ivec[ii](param)!=IOUtils::nodata && ivec[ii+1](param)!=IOUtils::nodata)
			vecSampling.push_back( ivec[ii+1].date.getJulian(true) - ivec[ii].date.getJulian(true) );
	}
	
	return Interpol1D::getMedian(vecSampling, false);
}

std::vector<ProcessingBlock::offset_spec> ProcShift::resampleVector(const std::vector<MeteoData>& ivec, const size_t& param) const
{
	const size_t n_ivec = ivec.size();
	const Date dt_start( ivec.front().date );
	const size_t nrSteps = static_cast<size_t>(round( (ivec.back().date.getJulian(true) - dt_start.getJulian(true)) / sampling_rate ));
	
	std::vector<ProcessingBlock::offset_spec> vecResults;
	size_t jj=0; //position within ivec
	
	for (size_t ii=0; ii<nrSteps; ii++) {
		const Date dt( dt_start+(static_cast<double>(ii)*sampling_rate) );
		
		if (ivec[jj].date==dt) { //do we have the exact element that we are looking for?
			vecResults.push_back( ProcessingBlock::offset_spec(dt, ivec[jj](param)) );
			jj++;
			continue;
		}
		
		while (ivec[jj].date<dt) { //find the first element >= dt
			if (jj==(n_ivec-1)) {
				return vecResults;
			}
			jj++;
		}
		
		if (jj>0) { //interpolate the value based on the left and right neighbors
			if (jj==n_ivec) return vecResults;
			
			const double x1 = ivec[jj-1].date.getJulian(true);
			const double x2 = ivec[jj].date.getJulian(true);
			if (x1==x2) { //skip duplicate timestamps if any
				jj++;
				continue;
			}
			
			if ((x2-x1) > 2.*sampling_rate) { //only interpolate between nearby points, otherwise keep nodata
				vecResults.push_back( ProcessingBlock::offset_spec(dt, IOUtils::nodata) );
				continue;
			}
			
			const double y1 = ivec[jj-1](param);
			const double y2 = ivec[jj](param);
			if (y1!=IOUtils::nodata && y2!=IOUtils::nodata) {
				const double a = (y2 - y1) / (x2 - x1);
				const double b = y2 - a*x2;
				const double x = dt.getJulian(true);
				const double y = (a*x + b);
				
				vecResults.push_back( ProcessingBlock::offset_spec(dt, y) );
			} else {
				vecResults.push_back( ProcessingBlock::offset_spec(dt, IOUtils::nodata) );
			}
		} else { //we are at the start of the vector, there is no jj-1
			vecResults.push_back( ProcessingBlock::offset_spec(dt, IOUtils::nodata) );
		}
	}
	
	return vecResults;
}

/**
 * @brief Pearson's correlation coefficient between two datasets when one will be time shifted
 * It is assumed that the two datasets have the exact same timestamps. The sums for X are 
 * recomputed to make sure that if Y[ii] is nodata, no value is taken for the matching X[ii]
 * See https://en.wikipedia.org/wiki/Pearson_correlation_coefficient for more
 * @param vecX Reference values
 * @param vecY dataset to compare with
 * @param curr_idx position that should be attributed the Pearson's coefficient
 * @param offset index offset to apply to the vecY vector
 * @return Pearson's correlation coefficient
 */
double ProcShift::getPearson(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx, const int& offset) const
{
	//compute the required data window, 
	//accounting for offsets in vecY that could bring us outside vecY
	size_t startIdx = std::max(curr_idx-width_idx/2, static_cast<size_t>(0));
	if (-offset>(signed)startIdx) startIdx = static_cast<size_t>(-offset);
	size_t endIdx = std::min(curr_idx+width_idx/2, vecX.size());
	if (endIdx+static_cast<size_t>(offset)>vecX.size()) endIdx = vecX.size() - static_cast<size_t>(offset);
	
	size_t count=0;
	double sumX=0., sumX2=0., sumY=0., sumY2=0., sumXY=0.;
	for (size_t ii=startIdx; ii<endIdx; ii++) {
		const double valueX = vecX[ii].value;
		const double valueY = vecY[static_cast<size_t>((signed)ii+offset)].value;
		if (valueX!=IOUtils::nodata && valueY!=IOUtils::nodata) {
			sumX += valueX;
			sumX2 += valueX * valueX;
			
			sumY += valueY;
			sumY2 += valueY * valueY;
			
			sumXY += valueX * valueY;
			count++;
		}
	}
	
	const double Xcontribution = static_cast<double>(count)*sumX2 - sumX*sumX;
	const double Ycontribution = static_cast<double>(count)*sumY2 - sumY*sumY;
	if (Xcontribution<=0. || Ycontribution<=0.) return IOUtils::nodata; //either count==0, count==1 or all values are identical
	
	const double pearson = (static_cast<double>(count)*sumXY - sumX*sumY) / (sqrt(Xcontribution) * sqrt(Ycontribution));
	return pearson;
}

/**
 * @brief Find optimal time shift between two datasets to have the highest possible correlation coefficient
 * All possible time shift (expressed as index shift between the two vectors) are evaluated in order to find
 * the absolute maximum.
 * @param vecX Reference values
 * @param vecY dataset to compare with
 * @param curr_idx position that should be attributed the optimal time shift
 * @param range_min minimum time shift (expressed as index shift) to start the scan
 * @param range_max maximum time shift (expressed as index shift) to end the scan
 * @return index shift that leads to the highest correlation coefficient for the given position
 */
int ProcShift::getOffsetFullScan(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx, const int& range_min, const int& range_max) const
{
	static const unsigned int minPts = 10; //minimum number of valid points to consider
	int offset_at_max = IOUtils::inodata;
	double pearson_max = IOUtils::nodata;
	unsigned int count = 0;
	
	for (int offset=range_min; offset<=range_max; offset++) {
		const double pearson = getPearson(vecX, vecY, curr_idx, offset);
		if (pearson==IOUtils::nodata) continue;
		
		count++;
		if (pearson>pearson_max) {
			offset_at_max = offset;
			pearson_max = pearson;
		}
	}
	
	//we want to avoid looking for a maximum among too few points
	if (count<minPts || pearson_max==IOUtils::nodata)
		return IOUtils::inodata;
	
	return offset_at_max;
}

/**
 * @brief Find optimal time shift between two datasets to have the highest possible correlation coefficient
 * A golden section search is performed in order to find the optimal time shift (offset) that leads to the maximum
 * correlation coefficient between the reference data and the time-shifted one. In some cases, all possible
 * offsets are evaluated (for example if some nodata values prevent the evaluation of the correlation coefficient
 * as required by the golden section search algorithm).
 * See https://en.wikipedia.org/wiki/Golden-section_search for more
 * @param vecX Reference values
 * @param vecY dataset to compare with
 * @param curr_idx position that should be attributed the optimal time shift
 * @return index shift that leads to the highest correlation coefficient for the given position
 */
int ProcShift::getOffset(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx) const
{
	static const size_t max_count = 100; //after more than max_count iterations, we consider there is no convergence
	static const double r = Cst::phi - 1.;
	static const double c = 1. - r;
	
	int range_min = static_cast<int>(-0.5*offset_range/sampling_rate);
	int range_max = static_cast<int>(0.5*offset_range/sampling_rate);
	
	int offset_c = IOUtils::inodata, offset_r = IOUtils::inodata;
	double pearson_c = IOUtils::nodata, pearson_r = IOUtils::nodata;
	size_t count = 0;
	
	while (true) {
		if (offset_c==IOUtils::inodata) {
			offset_c = static_cast<int>( round(((double)range_max - (double)range_min)*c) ) + range_min;
			pearson_c = getPearson(vecX, vecY, curr_idx, offset_c);
		}
		
		if (offset_r==IOUtils::inodata) {
			offset_r = static_cast<int>( round(((double)range_max - (double)range_min)*r) ) + range_min;
			pearson_r = getPearson(vecX, vecY, curr_idx, offset_r);
		}
		
		if (pearson_c==IOUtils::nodata || pearson_r==IOUtils::nodata) { //we must do a full scan
			offset_c = getOffsetFullScan(vecX, vecY, curr_idx, range_min, range_max);
			break;
		}
		
		if (pearson_c > pearson_r) {
			range_max = offset_r;
			offset_r = offset_c;
			pearson_r = pearson_c;
			offset_c = IOUtils::inodata;
		} else {
			range_min = offset_c;
			offset_c = offset_r;
			pearson_c = pearson_r;
			offset_r = IOUtils::inodata;
		}
		
		//check convergence
		if ((range_max-range_min) < 2) 
			break; //convergence criteria: within 1 cell of optimum
		if (count>max_count) {
			std::cout << "No convergence at " << vecX[curr_idx].date.toString(Date::ISO) << " offset_c=" << offset_c*24.*60. << " offset_r=" << offset_r*24.*60. << "\n";
			break;
		}
		
		count++;
	}
	
	return offset_c;
}

void ProcShift::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	double TZ=0.;
	bool has_type=false, has_cst=false, has_TZ=false, has_sampling_rate=false;
	bool has_width=false, has_offset_range=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="EXTRACT_OFFSETS") {
			IOUtils::parseArg(vecArgs[ii], where, extract_offsets);
		} else if (vecArgs[ii].first=="TYPE") {
			has_type=true;
			const std::string type_str( IOUtils::strToUpper(vecArgs[ii].second) );
			if (type_str=="CST") {
				offsets_interp = cst;
			} /*else if (type_str=="LINEAR") { //not implemented for now
				offsets_interp = linear;
			} */else if (type_str=="STEPWISE") {
				offsets_interp = stepwise;
			} else
				throw InvalidArgumentException("Invalid type \""+type_str+"\" specified for "+where, AT);
		} else if (vecArgs[ii].first=="TZ") {
			has_TZ=true;
			IOUtils::parseArg(vecArgs[ii], where, TZ);
		} else if (vecArgs[ii].first=="CST") {
			has_cst=true;
			IOUtils::parseArg(vecArgs[ii], where, cst_offset);
			cst_offset /= 3600.*24.; //convert to days
		} else if (vecArgs[ii].first=="OFFSETS_FILE") { //both for input & output
			offsets_file = vecArgs[ii].second;
		} else if (vecArgs[ii].first=="REF_PARAM") {
			IOUtils::parseArg(vecArgs[ii], where, ref_param);
		} else if (vecArgs[ii].first=="SAMPLING_RATE") { //assumed to be in seconds
			has_sampling_rate=true;
			IOUtils::parseArg(vecArgs[ii], where, sampling_rate);
			sampling_rate /= 3600.*24.; //convert to days
		} else if (vecArgs[ii].first=="WIDTH") { //assumed to be in seconds
			has_width=true;
			IOUtils::parseArg(vecArgs[ii], where, width_d);
			width_d /= 3600.*24.; //convert to days
		} else if (vecArgs[ii].first=="OFFSET_RANGE") { //assumed to be in seconds
			has_offset_range=true;
			IOUtils::parseArg(vecArgs[ii], where, offset_range);
			offset_range /= 3600.*24.; //convert to days
		}
	}
	
	//check consistency of provided configuration options
	if (extract_offsets==false) { //correction mode
		if (!ref_param.empty() || has_sampling_rate || has_width || has_offset_range)
			throw InvalidArgumentException("Please only provide the interpolation type and the correction constant or file for "+where, AT);
		if (has_cst && !offsets_file.empty())
			throw InvalidArgumentException("It is not possible to provide both a correction constant and a correction file for "+where, AT);
		if (!has_cst && (offsets_file.empty() || !has_TZ))
			throw InvalidArgumentException("Please either provide a correction constant or a correction file with its default time zone TZ for "+where+" or set the EXTRACT_OFFSETS option to TRUE", AT);
		
		//now read the corrections file if necessary
		if (!offsets_file.empty()) {
			//if this is a relative path, prefix the path with the current path
			const std::string prefix = ( FileUtils::isAbsolutePath(offsets_file) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+offsets_file, true) );  //clean & resolve path
			offsets_file = path + "/" + FileUtils::getFilename(offsets_file);
			
			corrections = ProcessingBlock::readCorrections(block_name, offsets_file, TZ);
			std::sort(corrections.begin(), corrections.end());
		}
		
	} else { //extracting the offsets
		if (offsets_file.empty())
			throw InvalidArgumentException("Please set OFFSETS_FILE where to write the extracted offsets for "+where, AT);
		if (ref_param.empty())
			throw InvalidArgumentException("Please provide the parameter name REF_PARAM to use as reference to extract the time offsets for "+where, AT);
		if (has_type || has_cst)
			throw InvalidArgumentException("It is not possible to provide the interpolation type or the correction constant when extracting offsets for "+where, AT);
		if (sampling_rate!=IOUtils::nodata && 0.5*offset_range<=sampling_rate)
			throw InvalidArgumentException("It is not possible to provide the interpolation type or the correction constant when extracting offsets for "+where, AT);
		
		if (FileUtils::fileExists(offsets_file)) {
			if (remove( offsets_file.c_str() ) != 0)
				throw AccessException("File \""+offsets_file+"\" already exists and can not be deleted before running "+where, AT);
		}
	}
}

} //end namespace
