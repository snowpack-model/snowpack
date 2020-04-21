/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <algorithm>

#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

FilterSuppr::FilterSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ)
          : ProcessingBlock(vecArgs, name), suppr_dates(), range(IOUtils::nodata), type(NONE)
{
	const std::string where( "Filters::"+block_name );
	properties.stage = ProcessingProperties::first; //for the rest: default values
	bool has_type=false, has_range=false, has_file=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string type_str( IOUtils::strToUpper( vecArgs[ii].second ) );
			if (type_str=="FRAC") type=FRAC_SUPPR;
			else if (type_str=="FILE") type=FILE_SUPPR;
			else if (type_str=="ALL") type=ALL_SUPPR;
			else
				throw InvalidArgumentException("Unknown type '"+type_str+"' for " + where, AT);

			has_type = true;
		}
		if (vecArgs[ii].first=="FRAC") {
			if (!IOUtils::convertString(range, vecArgs[ii].second))
				throw InvalidArgumentException("Invalid range \""+vecArgs[ii].second+"\" specified for "+where, AT);
			if (range<0. || range>1.)
				throw InvalidArgumentException("Wrong range for " + where + ", it should be between 0 and 1", AT);
			has_range = true;
		}
		if (vecArgs[ii].first=="FILE") {
			const std::string in_filename( vecArgs[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );

			suppr_dates = ProcessingBlock::readDates(block_name, filename, TZ);
			has_file = true;
		}
	}
	
	if (!has_type) throw InvalidArgumentException("Please provide a TYPE for "+where, AT);
	if (type==FILE_SUPPR && (has_range || !has_file)) throw InvalidArgumentException("For "+where+" and type=FILE, please provide the FILE argument and no other", AT);
	if (type==FRAC_SUPPR && (!has_range || has_file)) throw InvalidArgumentException("For "+where+" and type=FRAC, please provide the RANGE argument and no other", AT);
}

void FilterSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ovec.empty()) return;
	
	if (type==FILE_SUPPR) { //remove periods provided in a file
		if (!suppr_dates.empty()) supprByDates(param, ovec);
	} else if (type==ALL_SUPPR) { //remove all
		for (size_t ii=0; ii<ovec.size(); ii++)
			ovec[ii](param) = IOUtils::nodata;
	} else { //only remove a given fraction
		supprFrac(param, ivec, ovec);
	}
}

//this assumes that the DATES_RANGEs in suppr_dates have been sorted by increasing starting dates
void FilterSuppr::supprByDates(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	const std::string station_ID( ovec[0].meta.stationID ); //we know it is not empty
	const std::map< std::string, std::vector<dates_range> >::const_iterator station_it( suppr_dates.find( station_ID ) );
	if (station_it==suppr_dates.end()) return;

	const std::vector<dates_range> &suppr_specs = station_it->second;
	const size_t Nset = suppr_specs.size();
	size_t curr_idx = 0; //we know there is at least one
	for (size_t ii=0; ii<ovec.size(); ii++) {
		if (ovec[ii].date<suppr_specs[curr_idx].start) continue;

		if (ovec[ii].date<=suppr_specs[curr_idx].end) { //suppress the interval
			ovec[ii](param) = IOUtils::nodata;
			continue;
		} else { //look for the next interval
			curr_idx++;
			if (curr_idx>=Nset) break; //all the suppression points have been processed
		}
	}
}

void FilterSuppr::supprFrac(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) const
{
	const size_t set_size = ovec.size();
	const size_t nrRemove = static_cast<size_t>( round( (double)set_size*range ) );

	srand( static_cast<unsigned int>(time(NULL)) );
	size_t ii=1;
	while (ii<nrRemove) {
		const size_t idx = (unsigned)rand() % set_size;
		if (ivec[idx](param)!=IOUtils::nodata && ovec[idx](param)==IOUtils::nodata) continue; //the point was already removed

		ovec[idx](param) = IOUtils::nodata;
		ii++;
	}
}

} //end namespace
