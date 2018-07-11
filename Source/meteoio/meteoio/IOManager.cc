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
                                                       tsm1(iohandler, cfg), gdm1(iohandler, cfg), interpolator(cfg, tsm1, gdm1),
                                                       resampling(false)
{
	initIOManager();
}

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), iohandler(cfg),
                                            tsm1(iohandler, cfg), gdm1(iohandler, cfg), interpolator(cfg, tsm1, gdm1),
                                            resampling(false)
{
	initIOManager();
}

void IOManager::initIOManager()
{
	cfg.getValue("Resampling", "Input", resampling, IOUtils::nothrow);
	if (resampling) { //in this case, we do not want to re-apply the filters
		tsm1.setProcessingLevel(IOUtils::resampled | IOUtils::generated);
		gdm1.setProcessingLevel(IOUtils::resampled | IOUtils::generated);
	}
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	if (!resampling) {
		tsm1.setProcessingLevel(i_level);
		gdm1.setProcessingLevel(i_level);
	}
}

size_t IOManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (resampling) {
		return interpolator.getVirtualStationsMeta(date, vecStation);
	} else { //usual case
		return tsm1.getStationData(date, vecStation);
	}
}

size_t IOManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo) 
{
	return tsm1.getMeteoData(dateStart, dateEnd, vecVecMeteo);
}

void IOManager::clear_cache()
{
	tsm1.clear_cache();
	gdm1.clear_cache();
}

//data can be raw or processed (filtered, resampled)
size_t IOManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();
	
	if (resampling) 
		interpolator.pushVirtualMeteoData(i_date, tsm1);

	tsm1.getMeteoData(i_date, vecMeteo);
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

double IOManager::getAvgSamplingRate() const 
{
	return tsm1.getAvgSamplingRate();
}

void IOManager::add_to_points_cache(const Date& i_date, const METEO_SET& vecMeteo) 
{
	tsm1.add_to_points_cache(i_date, vecMeteo);
}

const std::string IOManager::toString() const {
	ostringstream os;
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &cfg << dec << "\n";
	os << iohandler.toString();
	os << tsm1.toString();
	os << gdm1.toString();
	os << interpolator.toString();
	os << "Resampling = " << resampling << "\n";
	os << "</IOManager>\n";
	return os.str();
}

} //namespace
