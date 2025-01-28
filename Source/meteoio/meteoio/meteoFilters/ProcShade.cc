// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/meteoFilters/ProcShade.h>
#include <meteoio/meteoLaws/Sun.h>
#include <meteoio/IOHandler.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/DEMAlgorithms.h>

using namespace std;

namespace mio {

const double ProcShade::diffuse_thresh = 15.; //below this threshold, not correction is performed since it will only be diffuse

ProcShade::ProcShade(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& i_cfg)
        : ProcessingBlock(vecArgs, name, i_cfg), cfg(i_cfg), dem(), masks(), horizons_outfile(), has_dem(false)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

ProcShade::~ProcShade()
{
	if (has_dem && !masks.empty() && !horizons_outfile.empty())
		DEMAlgorithms::writeHorizons(masks, horizons_outfile);
}

void ProcShade::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ovec.empty()) return;
	
	SunObject Sun;
	
	//check if the station already has an associated mask, first by stationID then as wildcard
	const std::string stationID( ovec[0].getStationID() );
	std::map< std::string , std::vector< std::pair<double,double> > >::const_iterator mask = masks.find( stationID );
	if (mask==masks.end()) {
		//now look for a wildcard fallback
		mask = masks.find( "*" );
		if (mask==masks.end()) {
			masks[ stationID ] = computeMask(dem, ovec[0].meta);
			mask = masks.find( stationID);
		}
	}
	
	double Md_prev = IOUtils::nodata;
	double julian_prev = 0.;
	for (size_t ii=0; ii<ovec.size(); ii++) { //now correct all timesteps
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values
		if (tmp<diffuse_thresh) continue; //only diffuse radiation, there is nothing to correct

		const Coords position( ovec[ii].meta.position );
		Sun.setLatLon(position.getLat(), position.getLon(), position.getAltitude()); //if they are constant, nothing will be recomputed
		Sun.setDate(ovec[ii].date.getJulian(true), 0.); //quicker: we stick to gmt
		double sun_azi, sun_elev;
		Sun.position.getHorizontalCoordinates(sun_azi, sun_elev);
		
		const double mask_elev = DEMAlgorithms::getHorizon(mask->second, sun_azi);
		if (mask_elev>0 && mask_elev>sun_elev) { //the point is in the shade
			const double TA=ovec[ii](MeteoData::TA), RH=ovec[ii](MeteoData::RH), HS=ovec[ii](MeteoData::HS), RSWR=ovec[ii](MeteoData::RSWR);
			double ISWR=ovec[ii](MeteoData::ISWR);

			double albedo = .5;
			if (RSWR==IOUtils::nodata || ISWR==IOUtils::nodata || RSWR<=0 || ISWR<=0) {
				if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
					albedo = (HS>=Cst::snow_nosnow_thresh)? Cst::albedo_fresh_snow : Cst::albedo_short_grass;

				if (ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata))
					ISWR = RSWR / albedo;
			} else {
				albedo = RSWR / ISWR;
				if (albedo>=1.) albedo=0.99;
				if (albedo<=0.) albedo=0.01;
			}

			const bool has_potRad = (ISWR!=IOUtils::nodata && TA!=IOUtils::nodata && RH!=IOUtils::nodata);
			if (has_potRad) 
				Sun.calculateRadiation(TA, RH, albedo);
			else 
				if (ovec[ii].date.getJulian(true) - julian_prev > 1.) continue; //no way to get ISWR and/or potential radiation, previous Md is too old
			
			const double Md = (has_potRad)? Sun.getSplitting(ISWR) : Md_prev; //fallback: use previous valid value
			tmp *= Md; //only keep the diffuse radiation, either on RSWR or ISWR
			if (has_potRad) {
				Md_prev = Md;
				julian_prev = ovec[ii].date.getJulian(true);
			}
		}
	}
}

std::vector< std::pair<double,double> > ProcShade::computeMask(const DEMObject& i_dem, const StationData& sd)
{
	//compute horizon by "angularResolution" increments
	const double cellsize = i_dem.cellsize;
	const double angularResolution = (cellsize<=20.)? 5. : (cellsize<=100.)? 10. : 20.;
	std::vector< std::pair<double,double> > o_mask( DEMAlgorithms::getHorizonScan(i_dem, sd.position, angularResolution) );
	if (o_mask.empty()) throw InvalidArgumentException( "In filter 'SHADE', could not compute mask from DEM '"+i_dem.llcorner.toString(Coords::LATLON)+"'", AT);

	return o_mask;
}

void ProcShade::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	bool from_dem=true;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="INFILE") {
			const std::string root_path( cfg.getConfigRootDir() );
			//if this is a relative path, prefix the path with the current path
			const std::string in_filename( vecArgs[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			const std::string filename( path + "/" + FileUtils::getFilename(in_filename) );
			masks = DEMAlgorithms::readHorizonScan(getName(), filename); //this mask is valid for ALL stations
			from_dem = false;
		} else if (vecArgs[ii].first=="OUTFILE") {
			IOUtils::parseArg(vecArgs[ii], where, horizons_outfile);
		}
	}

	if (from_dem) {
		IOHandler io(cfg);
		dem.setUpdatePpt( DEMObject::NO_UPDATE ); //we only need the elevations
		io.readDEM(dem);
		has_dem=true;
	}
}

} //end namespace
