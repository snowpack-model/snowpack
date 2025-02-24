// SPDX-License-Identifier: LGPL-3.0-or-later
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

#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/IOHandler.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/DEMAlgorithms.h>
#include <algorithm>

namespace mio {

TauCLDGenerator::TauCLDGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ, const Config &i_cfg)
                              : GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), last_cloudiness(), masks(), horizons_outfile(), cfg(i_cfg), dem(), sun(), cloudiness_model(DEFAULT), use_rswr(false), use_rad_threshold(false), write_mask_out(), use_horizons(false), from_dem(false)
{
	const std::string where( section+"::"+algo );
	bool has_infile=false;
	
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="CLD_TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );
			
			if (user_algo=="LHOMME") cloudiness_model = CLF_LHOMME;
			else if (user_algo=="KASTEN") cloudiness_model = KASTEN;
			else if (user_algo=="CRAWFORD") cloudiness_model = CLF_CRAWFORD;
			else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for "+where, AT);
		}
		if (vecArgs[ii].first=="USE_RSWR") {
			IOUtils::parseArg(vecArgs[ii], where, use_rswr);
		}
		if (vecArgs[ii].first=="USE_RAD_THRESHOLD") {
			IOUtils::parseArg(vecArgs[ii], where, use_rad_threshold);
		}
		if (vecArgs[ii].first=="SHADE_FROM_DEM") {
			IOUtils::parseArg(vecArgs[ii], where, from_dem);
		}
		if (vecArgs[ii].first=="INFILE") {
			const std::string root_path( cfg.getConfigRootDir() );
			//if this is a relative path, prefix the path with the current path
			const std::string in_filename( vecArgs[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );
			masks = DEMAlgorithms::readHorizonScan(where, filename); //this mask is valid for ALL stations
			has_infile = true;
		} else if (vecArgs[ii].first=="OUTFILE") {
			IOUtils::parseArg(vecArgs[ii], where, horizons_outfile);
			write_mask_out = true;
		}
	}
	
	use_horizons = (from_dem || has_infile);
	if (write_mask_out && !use_horizons)
		throw InvalidArgumentException("In "+where+", please provide an INFILE and/or a DEM to compute the shading and write it out!", AT);
}

TauCLDGenerator::~TauCLDGenerator()
{
	if (from_dem && !masks.empty() && !horizons_outfile.empty())
		DEMAlgorithms::writeHorizons(masks, horizons_outfile);
}

/**
 * @brief Compute the clearness index from an atmospheric cloudiness value
 * @details
 * This is a convenience method that helps process the same way various types of inputs: if a cloudiness
 * is provided (which is quite rare), it can be converted to a clearness index (ie the ratio of the incoming 
 * short wave radiation over the ground potential radiation, projected on the horizontal) and then processed 
 * the same way as more traditional measurements (ie only ISWR provided) where it will be re-converted
 * to a cloudiness (thus falling abck to the same cloudiness as originally provided).

 * @param[in] cloudiness cloudiness (between 0 and 1)
 * @return clearness index (between 0 and 1)
 */
double TauCLDGenerator::getClearness(const double& cloudiness) const
{
	if (cloudiness_model==CLF_LHOMME) {
		return Atmosphere::Lhomme_cloudiness( cloudiness/8. );
	} else if (cloudiness_model==KASTEN || cloudiness_model==DEFAULT) {
		return Atmosphere::Kasten_cloudiness( cloudiness/8. );
	} else if (cloudiness_model==CLF_CRAWFORD) {
		return Atmosphere::Lhomme_cloudiness( cloudiness/8. );
	} else
		return IOUtils::nodata; //this should never happen
}

std::vector< std::pair<double,double> > TauCLDGenerator::computeMask(const DEMObject& i_dem, const StationData& sd)
{
	//compute horizon by "angularResolution" increments
	const double cellsize = i_dem.cellsize;
	const double angularResolution = (cellsize<=20.)? 5. : (cellsize<=100.)? 10. : 20.; 
	std::vector< std::pair<double,double> > o_mask( DEMAlgorithms::getHorizonScan(i_dem, sd.position, angularResolution) );
	if (o_mask.empty()) throw InvalidArgumentException( "In Generator, could not compute mask from DEM '"+i_dem.llcorner.toString(Coords::LATLON)+"'", AT);

	return o_mask;
}

//check if the station already has an associated mask, first by stationID then as wildcard, otherwise compute it from dem
double TauCLDGenerator::getHorizon(const MeteoData& md, const double& sun_azi)
{
	const std::string stationID( md.getStationID() );
	std::map< std::string , std::vector< std::pair<double,double> > >::const_iterator mask = masks.find( stationID );
	if (mask!=masks.end()) return DEMAlgorithms::getHorizon(mask->second, sun_azi);
	
	//now look for a wildcard fallback
	mask = masks.find( "*" );
	if (mask!=masks.end()) return DEMAlgorithms::getHorizon(mask->second, sun_azi);
	
	if (!from_dem)
		throw InvalidArgumentException("No horizon could be found for station "+stationID+", please either provide it in the horizon file or provide a DEM to compute the horizon", AT);

	//get the horizon from the DEM, save the computed mask
	if (dem.empty()) {
		IOHandler io(cfg);
		dem.setUpdatePpt( DEMObject::NO_UPDATE ); //we only need the elevations
		io.readDEM(dem);
	}

	masks[ stationID ] = computeMask(dem, md.meta);
	mask = masks.find( stationID );
	return DEMAlgorithms::getHorizon(mask->second, sun_azi);
}

double TauCLDGenerator::interpolateCloudiness(const std::string& station_hash, const double& julian_gmt) const
{
	const auto& cloudiness_point = last_cloudiness.find(station_hash); //we get a cloudCache object
	if (cloudiness_point==last_cloudiness.end()) {
		return IOUtils::nodata;
	}

	//TODO for now, exact same behavior as before. The goal is to implement better interpolations!
	const double last_cloudiness_julian = cloudiness_point->second.last_valid.first;
	const double last_cloudiness_value = cloudiness_point->second.last_valid.second;
	double cloudiness = IOUtils::nodata;
	if ((julian_gmt - last_cloudiness_julian) < 1.) cloudiness = last_cloudiness_value;

	return cloudiness;
}

/**
 * @brief Add a cloudiness value to the cache.
 * @param[in] julian_gmt Timestamp of the cloudiness value (in GMT)
 * @param[in] cloudiness cloudiness value (between 0 and 1)
 */
void TauCLDGenerator::cloudCache::addCloudiness(const double& julian_gmt, const double& cloudiness)
{
	last_valid = std::make_pair(julian_gmt, cloudiness);
}

/**
 * @brief Return the atmospheric cloudiness
 * @details
 * Based on the available data, it might be a direct copy of a cloudiness value available in the Meteodata, a parametrization
 * based on clearness index (ie the ratio of the incoming short wave radiation over the ground potential radiation, projected on the horizontal) or
 * a temporal interpolations (over the night or while the AWS is in the shade of the terrain). 
 * @param[in] md MeteoData
 * @return cloudiness (between 0 and 1)
 */
double TauCLDGenerator::getCloudiness(const MeteoData& md)
{
	const double TAU_CLD=md(MeteoData::TAU_CLD);
	const double CLD = (md.param_exists("CLD"))? md("CLD") : IOUtils::nodata;
	double cloudiness = (TAU_CLD!=IOUtils::nodata)? Atmosphere::Kasten_cloudiness( TAU_CLD ) : IOUtils::nodata;

	if (CLD!=IOUtils::nodata) {
		//Synop sky obstructed from view -> fully cloudy
		if (CLD>9. || CLD<0.) throw InvalidArgumentException("Cloud cover CLD should be between 0 and 8!", AT);
		cloudiness = std::max(std::min(CLD/8., 1.), 0.1);
	}

	const std::string station_hash( md.meta.stationID + ":" + md.meta.stationName );
	const double julian_gmt = md.date.getJulian(true);

	//try to get a cloudiness value
	if (cloudiness==IOUtils::nodata) {
		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return IOUtils::nodata;

		bool is_night;
		sun.setLatLon(lat, lon, alt); //there is some caching in Sun if lat/lon/alt are unchanged
		sun.setDate(julian_gmt, 0.);
		cloudiness = computeCloudiness(md, is_night);
		if (cloudiness==IOUtils::nodata && !is_night) return IOUtils::nodata;

		if (is_night) { //interpolate the cloudiness over the night but don't cache it!
			return interpolateCloudiness(station_hash, julian_gmt); //either interpolated cloudiness or nodata
		}
	}

	//save the last valid cloudiness
	last_cloudiness[station_hash] = cloudCache( julian_gmt, cloudiness );

	return cloudiness;
}


/**
 * @brief Compute the atmospheric cloudiness from the available measurements
 * @details
 * The clearness index (ie the ratio of the incoming short wave radiation over the ground potential radiation, projected on the horizontal) 
 * is computed and used to evaluate the cloudiness, based on the chosen parametrization.
 * @param[in] md MeteoData
 * be a soil ro a snow albedo
 * @param[out] is_night set to TRUE if it is night time
 * @return cloudiness (between 0 and 1) or IOUtils::nodata
 */
double TauCLDGenerator::computeCloudiness(const MeteoData& md, bool &is_night)
{
	//we know that TA and RH are available, otherwise we would not get called
	const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
	if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;		//NOTE: keep it for safety?
	const double HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR);
	double ISWR=md(MeteoData::ISWR);
	
	double sun_azi, sun_elev;
	sun.position.getHorizontalCoordinates(sun_azi, sun_elev);
	const double mask_elev = (use_horizons)? getHorizon(md, sun_azi) : 5.;

	//at sunrise or sunset, we might get very wrong results -> return nodata in order to use interpolation instead
	//obviously, when it is really night neither can we compute anything here...
	//we use 5 degrees to represent very low elevation rays of light that would not work well in the horizontal sensors
	is_night = (sun_elev <= mask_elev || sun_elev<5.);
	if (is_night) return IOUtils::nodata;
	
	double albedo = .5;
	if (RSWR!=IOUtils::nodata && ISWR!=IOUtils::nodata && RSWR>Atmosphere::day_iswr_thresh && ISWR>Atmosphere::day_iswr_thresh) {
		albedo = std::min( 0.99, std::max(0.01, RSWR / ISWR) );
	} else { //so some measurements are missing
		if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=Cst::snow_nosnow_thresh)? Cst::albedo_fresh_snow : Cst::albedo_short_grass;

		if (ISWR==IOUtils::nodata) { //ISWR is missing, trying to compute it
			if (!use_rswr) return IOUtils::nodata;
			if (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)
				ISWR = RSWR / albedo;
			else
				return IOUtils::nodata; //no way to get ISWR, aborting
		}
	}

	sun.calculateRadiation(TA, RH, albedo);
	double toa, direct, diffuse;
	sun.getHorizontalRadiation(toa, direct, diffuse);
	const double iswr_clear_sky = direct+diffuse;
	
	if (use_rad_threshold) {
		if (ISWR<Atmosphere::day_iswr_thresh || iswr_clear_sky<Atmosphere::day_iswr_thresh || iswr_clear_sky<ISWR) {
			is_night = true;
			return IOUtils::nodata;
		}
	}

	const double clearness_index = std::min(ISWR/iswr_clear_sky, 1.);
	if (cloudiness_model==CLF_LHOMME) {
		const double clf = Atmosphere::Lhomme_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else if (cloudiness_model==KASTEN || cloudiness_model==DEFAULT) {
		const double clf = Atmosphere::Kasten_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else if (cloudiness_model==CLF_CRAWFORD) {
		const double clf = Atmosphere::Lhomme_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else
		return IOUtils::nodata; //this should never happen
}

bool TauCLDGenerator::generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& /*vecMeteo*/)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		const double cloudiness = getCloudiness(md);
		value = 1.- cloudiness;
	}

	return true; //all missing values could be filled
}

bool TauCLDGenerator::create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool status = true;
	for (size_t ii=ii_min; ii<ii_max; ii++) {
		if (!generate(param, vecMeteo[ii], vecMeteo))
			status = false;
	}

	return status;
}

} //namespace
