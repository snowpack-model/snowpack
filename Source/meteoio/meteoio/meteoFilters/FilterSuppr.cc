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
#include <errno.h>
#include <algorithm>

#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/FileUtils.h>

#include <fstream>

using namespace std;

namespace mio {

FilterSuppr::FilterSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ)
          : FilterBlock(vecArgs, name), suppr_dates(), range(IOUtils::nodata)
{
	const std::string where( "Filters::"+block_name );
	properties.stage = ProcessingProperties::first; //for the rest: default values
	const size_t nrArgs = vecArgs.size();

	if (nrArgs>1)
		throw InvalidArgumentException("Wrong number of arguments for " + where, AT);

	if (nrArgs==1) {
		if (vecArgs[0].first=="FRAC") {
			if (!IOUtils::convertString(range, vecArgs[0].second))
				throw InvalidArgumentException("Invalid range \""+vecArgs[0].second+"\" specified for "+where, AT);
			if (range<0. || range>1.)
				throw InvalidArgumentException("Wrong range for " + where + ", it should be between 0 and 1", AT);
		} else if (vecArgs[0].first=="SUPPR") {
			const std::string in_filename( vecArgs[0].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );

			fillSuppr_dates(filename, TZ);
		} else
			throw UnknownValueException("Unknown option '"+vecArgs[0].first+"' for "+where, AT);
	}
}

void FilterSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ovec.empty()) return;
	
	if (!suppr_dates.empty()) {
		supprByDates(param, ovec);
	} else if (range==IOUtils::nodata) { //remove all
		for (size_t ii=0; ii<ovec.size(); ii++)
			ovec[ii](param) = IOUtils::nodata;
	} else { //only remove a given fraction
		supprFrac(param, ivec, ovec);
	}
}

//this assumes that the SUPPR_SETs in suppr_dates have been sorted by increasing starting dates
void FilterSuppr::supprByDates(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	const std::string station_ID( ovec[0].meta.stationID ); //we know it is not empty
	const std::map< std::string, std::vector<suppr_spec> >::const_iterator station_it( suppr_dates.find( station_ID ) );
	if (station_it==suppr_dates.end()) return;

	const std::vector<suppr_spec> &suppr_specs = station_it->second;
	const size_t Nset = suppr_specs.size();
	size_t curr_idx=0; //we know there is at least one
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
		const size_t idx = rand() % set_size;
		if (ivec[idx](param)!=IOUtils::nodata && ovec[idx](param)==IOUtils::nodata) continue; //the point was already removed

		ovec[idx](param) = IOUtils::nodata;
		ii++;
	}
}

void FilterSuppr::fillSuppr_dates(const std::string& filename, const double& TZ)
{
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename, AT);
	if (!FileUtils::fileExists(filename)) throw NotFoundException(filename, AT);
	
	std::ifstream fin(filename.c_str());
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << block_name << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	
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
				throw InvalidFormatException("Invalid syntax for filter " + block_name + " in file \"" + filename + "\": expecting at least 2 arguments", AT);

			const std::string station_ID( vecString[0] );
			if (nrElems==2) {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				const suppr_spec specs(d1, d1);
				suppr_dates[ station_ID ].push_back( specs );
			} else if (nrElems==3) {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[2], TZ))
					throw InvalidFormatException("Could not process date "+vecString[2]+" in file \""+filename+"\"", AT);
				const suppr_spec specs(d1, d2);
				suppr_dates[ station_ID ].push_back( specs );
			} else if (nrElems==4 && vecString[2]=="-") {
				if (!IOUtils::convertString(d1, vecString[1], TZ))
					throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
				if (!IOUtils::convertString(d2, vecString[3], TZ))
					throw InvalidFormatException("Could not process date "+vecString[3]+" in file \""+filename+"\"", AT);
				const suppr_spec specs(d1, d2);
				suppr_dates[ station_ID ].push_back( specs );
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
	std::map< std::string, std::vector<suppr_spec> >::iterator station_it( suppr_dates.begin() );
	for (; station_it!=suppr_dates.end(); ++station_it) {
		std::sort(station_it->second.begin(), station_it->second.end());
	}
}

} //end namespace
