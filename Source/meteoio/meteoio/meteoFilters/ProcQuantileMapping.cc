/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcQuantileMapping.h>
#include <meteoio/FileUtils.h>
#include <meteoio/meteoStats/libinterpol1D.h>

#include <ctime>
#include <cstdlib>

using namespace std;

namespace mio {

ProcQuantileMapping::ProcQuantileMapping(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& i_root_path)
        : ProcessingBlock(vecArgs, name), quantiles(), corrections(), root_path(i_root_path), period_duration(IOUtils::nodata), type('N')
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcQuantileMapping::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	
	if (period_duration==IOUtils::nodata) { //we work on the full buffer at once
		correctPeriod(param, 0, ovec.size() - 1, ivec, ovec);
	} else { //we work on each period one by one
		const std::vector< std::pair<size_t, size_t> > starts( getStarts(ivec) );
		
		for (size_t ii=0; ii<starts.size(); ii++) 
			correctPeriod(param, starts[ii].first, starts[ii].second, ivec, ovec);
	}
}

void ProcQuantileMapping::correctPeriod(const unsigned int& param, const size_t& idx_start, const size_t& idx_end, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) const
{
	std::vector<double> vecX(idx_end - idx_start + 1);
	for (size_t ii=idx_start; ii<=idx_end; ii++) {
		const double val = ivec[ii]( param );
		if (val!=IOUtils::nodata) vecX[ii-idx_start] = val;
	}

	const std::vector<double> thresholds( Interpol1D::quantiles_core(vecX, quantiles) ); //compute the thresholds for each bin
	
	//correct each data point according to its rank (ie which threshold it belongs to)
	for (size_t ii=idx_start; ii<=idx_end; ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values
		const double corr = getCorrection(thresholds, tmp);

		if (type=='a') tmp += corr;
		else if (type=='m') tmp *= corr;
	}
}

//find which quantile a given value belongs to. Here we assume that there are few quantiles so there is no need for a binary search...
double ProcQuantileMapping::getCorrection(const std::vector<double>& thresholds, const double& value) const
{
	const size_t nrThresh = thresholds.size();
	
	for (size_t ii=0; ii<nrThresh; ii++) {
		if (value<=thresholds[ii]) {
			if (thresholds[ii]<0.5) return corrections[ii];
			else if (ii>0) return corrections[ii-1];
			else return 1.;
		}
	}
	
	return 1.;
}

std::vector< std::pair<size_t, size_t> > ProcQuantileMapping::getStarts(const std::vector<MeteoData>& ivec) const
{
	std::vector< std::pair<size_t, size_t> > results;

	if ((ivec.back().date.getJulian(true) - ivec.front().date.getJulian(true)) < .95*period_duration ) {
		std::cerr << "Not enough data in buffer to use Filters::"+block_name+"\n";
		return results;
	}

	//N periods will be found with a remainder until the end of the vector.
	//This remainder will be associated a period that starts one full period before the end of the vector.
	const double end_target = ivec.back().date.getJulian(true) - period_duration;
	size_t start_idx = 0, end_idx = IOUtils::npos; //end_idx is the index for the start of the end period
	bool found;
	do {
		found = false;
		const double target_jul = ivec[ start_idx ].date.getJulian(true) + period_duration;
		for (size_t jj = start_idx+1; jj<ivec.size(); jj++) {
			if (ivec[jj].date.getJulian(true) >= end_target && end_idx==IOUtils::npos)
				end_idx = jj;

			if (ivec[jj].date.getJulian(true) >= target_jul) {
				results.push_back( std::make_pair(start_idx, jj) );
				start_idx = jj;
				found = true;
				break;
			}
		}
	} while (found);

	if (end_idx!=IOUtils::npos) results.push_back( std::make_pair(end_idx, ivec.size()-1) );

	return results;
}

void ProcQuantileMapping::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	bool has_corr_type = false;
	std::string filename;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="CORRECTIONS") {
			//if this is a relative path, prefix the path with the current path
			const std::string in_filename( vecArgs[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			filename = path + "/" + FileUtils::getFilename(in_filename);
		} else if (vecArgs[ii].first=="PERIOD") {
			const std::string period_str( IOUtils::strToUpper(vecArgs[ii].second) );
			if (period_str=="YEARLY") {
				period_duration = 365.25;
			} else if (period_str=="MONTHLY") {
				period_duration = 365.25/12.;
			} else if (period_str=="DAILY") {
				period_duration = 1.;
			} else
				throw InvalidArgumentException("Invalid period \""+period_str+"\" specified for "+where, AT);
		} else if (vecArgs[ii].first=="TYPE") {
			const std::string type_str( IOUtils::strToUpper(vecArgs[ii].second) );
			if (type_str=="ADD") {
				type = 'a';
			} else if (type_str=="MULT") {
				type = 'm';
			} else
				throw InvalidArgumentException("Invalid type \""+type_str+"\" specified for "+where, AT);
			has_corr_type = true;
		}
	}

	if (!has_corr_type) throw InvalidArgumentException("Please provide a type for "+where, AT);
	if (filename.empty()) throw InvalidArgumentException("Please provide a correction file for "+where, AT);
	ProcessingBlock::readCorrections(getName(), filename, quantiles, corrections); //they must be in INCREASING order
	if (quantiles.empty()) throw InvalidArgumentException("Please provide a VALID correction file for "+where, AT);
	
	//check that the quantiles contain both 0 and 1 and that they are in increasing order
	bool status = true;
	if (quantiles.front()!=0.) status=false;
	if (quantiles.back()!=1.) status=false;
	double prevQ = 0.;
	for (size_t ii=1; ii<quantiles.size(); ii++) {
		if (quantiles[ii]<=prevQ) status = false;
		prevQ = quantiles[ii];
	}
	if (!status) throw InvalidArgumentException("Please provide a VALID correction file (quantiles from 0 to 1 in increasing order) for "+where, AT);
	
}

} //end namespace
