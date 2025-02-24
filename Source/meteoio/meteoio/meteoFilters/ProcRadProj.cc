// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2024 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcRadProj.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Sun.h>

using namespace std;

namespace mio {
	
RadProj::RadProj(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
          : ProcessingBlock(vecArgs, name, cfg), src_azi(0.), dest_azi(0.), src_angle(0.), dest_angle(0.) 
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void RadProj::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	//TODO: optimization: around solar noon, if elev is large enough, do not recompute more than twice an hour
	ovec = ivec;
	if (ivec.empty()) return;
	
	//using the "beta" naming in order to be the same as in Suntrajectory
	const double cos_beta_src  = cos( src_angle*Cst::to_rad );
	const double sin_beta_src  = sin( src_angle*Cst::to_rad );
	const double cos_beta_dest = cos( dest_angle*Cst::to_rad );
	const double sin_beta_dest = sin( dest_angle*Cst::to_rad );
	SunObject sun;
	
	for (auto& meteo : ovec) {
		double& iswr = meteo(param);
		if (iswr == IOUtils::nodata) continue; //preserve nodata values
		
		//get the Sun coordinates and cos/sin of angles
		sun.setLatLon( meteo.meta.position.getLat(), meteo.meta.position.getLon(), meteo.meta.position.getAltitude() );
		sun.setDate( meteo.date.getJulian(true), 0.);
		const double sun_elev = sun.position.getSolarElevation();
		const double sun_azi_rad = sun.position.getSolarAzimuth()*Cst::to_rad;
		const double cos_Z = cos( (90.-sun_elev)*Cst::to_rad );
		const double sin_Z = sin( (90.-sun_elev)*Cst::to_rad );
		
		//get the necessary forcings and compute the toa and clear sky radiation
		const double HS=meteo(MeteoData::HS), TA=meteo(MeteoData::TA), RH=meteo(MeteoData::RH), P=meteo(MeteoData::P);
		const double albedo = (HS>=Cst::snow_nosnow_thresh)? Cst::albedo_fresh_snow : Cst::albedo_short_grass;
		sun.calculateRadiation(TA, RH, P, albedo);
		double R_toa, R_direct, R_diffuse;
		sun.getBeamRadiation(R_toa, R_direct, R_diffuse);
		if (R_toa==IOUtils::nodata) { //could not compute the radiation, because of missing meteo forcings
			iswr = IOUtils::nodata;
			continue;
		}
		//the Sun is below the horizon, any radiation is diffuse and thus not reprojected
		if (R_toa==0.) continue;
		
		//get the projection coefficients from src slope to dest slope, knowing Sun zenith angle Z and azimuth
		const double cos_theta_src  = (cos_beta_src*cos_Z + sin_beta_src*sin_Z*cos((sun_azi_rad-src_azi)));
		const double cos_theta_dest = (cos_beta_dest*cos_Z + sin_beta_dest*sin_Z*cos((sun_azi_rad-dest_azi)));
		
		//split direct / diffuse and only reproject direct
		const double R_direct_slope = R_direct*cos_theta_src;
		if (R_direct_slope<=0.) { //we are in self shading, we can not reproject anything
			iswr = IOUtils::nodata;
			continue;
		}
		//R_direct and R_diffuse are at ground level, toa is not
		const double atm_losses = R_toa / (R_direct + R_diffuse);
		const double toa_slope  = (R_direct_slope + R_diffuse) * atm_losses;
		//splitting is normally performed on the horizontal. Hoping (!!) that by using toa on the 
		//slope and iswr on the slope, it would still work fine...
		const double Md = sun.getSplittingBoland(toa_slope, iswr, sun.position.getSolarTimeOfDay()*24.);	
		
		iswr = ((1.-Md) * iswr) * cos_theta_dest / cos_theta_src + Md * iswr; //only the direct component is reprojected
	}
}


void RadProj::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	bool has_src_azi=false, has_src_angle=false;
	//const bool has_src_azi=true, has_src_angle=true;
	bool has_dest_azi=false, has_dest_angle=false;

	for (const auto& arg : vecArgs) {
		if (arg.first=="SRC_AZI") {
			IOUtils::parseArg(arg, where, src_azi);
			has_src_azi = true;
		} else if (arg.first=="SRC_ANGLE") {
			IOUtils::parseArg(arg, where, src_angle);
			has_src_angle = true;
		} else if (arg.first=="DEST_AZI") {
			IOUtils::parseArg(arg, where, dest_azi);
			has_dest_azi = true;
		} else if (arg.first=="DEST_ANGLE") {
			IOUtils::parseArg(arg, where, dest_angle);
			has_dest_angle = true;
		}
	}

	if ((has_src_azi && !has_src_angle) || (!has_src_azi && has_src_angle)) throw InvalidArgumentException("Please provide both SRC_AZI and SRC_SLOPE for "+where, AT);
	if ((has_dest_azi && !has_dest_angle) || (!has_dest_azi && has_dest_angle)) throw InvalidArgumentException("Please provide both DEST_AZI and DEST_SLOPE for "+where, AT);
	if (!has_src_azi && !has_dest_azi) throw InvalidArgumentException("Please provide either the source slope/azimuth and/or the destination slope/azimuth for "+where, AT);
}

} //end namespace
