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

#include <meteoio/IOManager.h>

using namespace std;

namespace mio {

IOManager::IOManager(const std::string& filename_in) : cfg(filename_in), iohandler(cfg),
                                                       tsmanager(iohandler, cfg), gridsmanager(iohandler, cfg), interpolator(cfg, tsmanager, gridsmanager),
                                                       resampling(false)
{
	initIOManager();
}

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), iohandler(cfg),
                                            tsmanager(iohandler, cfg), gridsmanager(iohandler, cfg), interpolator(cfg, tsmanager, gridsmanager),
                                            resampling(false)
{
	initIOManager();
}

void IOManager::initIOManager()
{
	cfg.getValue("Resampling", "Input", resampling, IOUtils::nothrow);
	if (resampling) { //in this case, we do not want to re-apply the filters
		tsmanager.setProcessingLevel(IOUtils::resampled | IOUtils::generated);
		gridsmanager.setProcessingLevel(IOUtils::resampled | IOUtils::generated);
	}
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	if (!resampling) {
		tsmanager.setProcessingLevel(i_level);
		gridsmanager.setProcessingLevel(i_level);
	}
}

size_t IOManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (resampling) {
		return interpolator.getVirtualStationsMeta(date, vecStation);
	} else { //usual case
		return tsmanager.getStationData(date, vecStation);
	}
}

void IOManager::clear_cache()
{
	tsmanager.clear_cache();
	gridsmanager.clear_cache();
}

//data can be raw or processed (filtered, resampled)
size_t IOManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();
	
	if (resampling) 
		interpolator.pushVirtualMeteoData(i_date, tsmanager);

	tsmanager.getMeteoData(i_date, vecMeteo);
	return vecMeteo.size();
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result)
{
	std::string info_string;
	const bool status = getMeteoData(date, dem, meteoparam, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
	return status;
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, result, info_string);
	return (!result.empty());
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result)
{
	std::string info_string;
	interpolate(date, dem, meteoparam, in_coords, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, in_coords, result, info_string);
}

const std::string IOManager::toString() const {
	ostringstream os;
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &cfg << dec << "\n";
	os << iohandler.toString();
	os << tsmanager.toString();
	os << gridsmanager.toString();
	os << interpolator.toString();
	os << "Resampling = " << resampling << "\n";
	os << "</IOManager>\n";
	return os.str();
}

} //namespace
