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
#include <meteoio/plugins/CsvIO.h>

using namespace std;

namespace mio {
/**
 * @page csvio CsvIO
 * @section csvio_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section csvio_units Units
 *
 *
 * @section csvio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double CsvIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

CsvIO::CsvIO(const std::string& configfile) : cfg(configfile)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	/* Example: how to read keys from the Config object*/
	/*const double factor = cfg.get("PLUGIN_FACTOR", "Input"); //if the key PLUGIN_FACTOR is not found in the [Input] section, an exception will be thrown
	 * 
	 * bool enable_feature = false;
	 * cfg.getValue("ENABLE_FEATURE", "Input", enable_feature, IOUtils::nothrow); //if the key is not found, it simply keeps its previous value
	 * 
	 * int parameter = 0;
	 * cfg.getValue("PLUGIN_NR_PARAMS", "Output", parameter); //if the key is not found, an exception will be thrown
	 * 
	 * //it is also possible to get all the keys starting with a given pattern at once and then loop through them:
	 * std::vector<std::string> vecFilenames;
	* cfg.getValues("STATION", "INPUT", vecFilenames);
	 */
}

CsvIO::CsvIO(const Config& cfgreader) : cfg(cfgreader)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

CsvIO::~CsvIO() throw()
{
	//if there is no need to cleanup some pointers before exiting, do not even declare a destructor!
}

void CsvIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
	
	/* Example: how to read coordinate in the input coordinate system and make them availabe to any other coordinate system*/
	/* //we assume we already have eatsing, northing and altitude in "double" variables
	* Coords point1(coordin, coordinparam);
	* point1.setXY(easting , northing., altitude);
	* 
	* //in order to get lat/lon:
	* const double lat = point1.getLat();
	* 
	* //if we have lat/lon instead:
	* point1.setLatLon(lat , lon., altitude);
	* 
	* //to set a coodrinate system by EPSG code:
	* point1.setEPSG(epsg);
	* 
	* //for more information: see the "Coords" class and the "coordinates.cc" example
	*/
	
	
	/* Example: how to copy the meteo data has read by the plugin into the vector of MeteoData*/
	/* 
	 * //for each station, repeat the following:
	 *  vector<MeteoData> timeseries;				//create a time vector of MeteoData for the station
	 * 
	 * //we read the station coordinates and put them in a Coords object (see above).
	 * //we also read a station id and station name (they might be identical if only one is available)
	 * StationData sd(point1, id, name);
	 *
	 * //then, we loop over all the available timestamps and fill the MeteoData
	 * for (size_t step=0; step<nr_steps; step++) {
	* 	//we read the date and put it in a Date object:
	* 	Date date(julian, TZ_in); 						//from julian date and input timezone
	* 
	* 	Date date;									//other possibility:
	* 	IOUtils::convertString(date, timestamp, TZ_in); 	//from a text timestamp and input timezone
	* 
	* 	MeteoData md(date, sd);						//create an empty MeteoData object initialized at a given date and station
	* 	md(MeteoData::TA) = my_ta;					//add each field
	* 	md(MeteoData::RH) = my_rh;					//of course, this could be done without relying on these "my_rh" intermediate variables!
	* 	//etc
	* 
	* 	timeseries.push_back( md );					//add the MeteoData to the time vector
	 * }
	 */
}

void CsvIO::cleanup() throw()
{
	//if there is nothing to cleanup, remove this method
}

} //namespace
