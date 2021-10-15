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
#include <fstream>
#include <cstring>
#include <meteoio/meteoFilters/ProcTransformWindVector.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/FileUtils.h>
#ifdef PROJ
	#include <proj_api.h>
#endif
#include <stdio.h>

namespace mio {

#ifndef PROJ
ProcTransformWindVector::ProcTransformWindVector(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
          : ProcessingBlock(vecArgs, name, cfg), s_coordparam(), t_coordparam(), RACMO2(false)
{
	throw IOException("ProcTransformWindVector requires PROJ library. Please compile MeteoIO with PROJ support.", AT);
}

ProcTransformWindVector::~ProcTransformWindVector() {}

void ProcTransformWindVector::process(const unsigned int&, const std::vector<MeteoData>&,
                        std::vector<MeteoData>&)
{
	throw IOException("ProcTransformWindVector requires PROJ library. Please compile MeteoIO with PROJ support.", AT);
}

#else

inline bool isEPSG(const std::string &c) {
	// Check if EPSG code, in which case the string should only contain numbers
	return c.find_first_not_of("0123456789") == std::string::npos;
}

ProcTransformWindVector::ProcTransformWindVector(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
          : ProcessingBlock(vecArgs, name, cfg), pj_src(nullptr), pj_dest(nullptr), vecArgs_i(vecArgs), name_i(name), cfg_i(cfg), s_coordparam(), t_coordparam(), RACMO2(false)
{
	parse_args(vecArgs, cfg);
	//the filters can be called at two points: before the temporal resampling (first stage, ProcessingProperties::first)
	//or after the temporal resampling (second stage, ProcessingProperties::second) or both (ProcessingProperties::both)
	//filters that do not depend on past data can safely use "both" (such as min/max filters) while
	//corrections should only use "first" (such as undercatch correction) in order to avoid double-correcting the data
	//filters relying on past data (through variance or other statistical parameters) usually should use "first"
	properties.stage = ProcessingProperties::first;
}

ProcTransformWindVector::~ProcTransformWindVector() {
	if (pj_src!=nullptr) pj_free(pj_src);
	if (pj_dest!=nullptr) pj_free(pj_dest);
}

ProcTransformWindVector::ProcTransformWindVector(const ProcTransformWindVector& c) :
	ProcessingBlock(c.vecArgs_i, c.name_i, c.cfg_i), pj_src(nullptr), pj_dest(nullptr), vecArgs_i(c.vecArgs_i), name_i(c.name_i), cfg_i(c.cfg_i), s_coordparam(c.s_coordparam), t_coordparam(c.t_coordparam), RACMO2(c.RACMO2)
{
	initPROJ();
}

ProcTransformWindVector& ProcTransformWindVector::operator=(const ProcTransformWindVector& source) {
	if (this != &source) {
		pj_src = nullptr;
		pj_dest = nullptr;
		vecArgs_i = source.vecArgs_i;
		name_i = source.name_i;
		cfg_i = source.cfg_i;
		t_coordparam = source.t_coordparam;
		initPROJ();
	}
	return *this;
}

void ProcTransformWindVector::initPROJ()
{
	std::string src_param;
	if (s_coordparam == "") {
		src_param = std::string("+proj=latlong +datum=WGS84 +ellps=WGS84");
	} else if (isEPSG(s_coordparam)) {
		src_param = std::string("+init=epsg:" + s_coordparam);
	} else {
		// Assume it is a proj string
		src_param = s_coordparam;
	}

	std::string dest_param;
	if (t_coordparam == "") {
		dest_param = std::string("+proj=latlong +datum=WGS84 +ellps=WGS84");
	} else if (isEPSG(t_coordparam)) {
		dest_param = std::string("+init=epsg:" + t_coordparam);
	} else {
		// Assume it is a proj string
		dest_param = t_coordparam;
	}

	if ( !(pj_src = pj_init_plus(src_param.c_str())) ) {
		throw InvalidArgumentException("Failed to initalize Proj with given arguments: "+src_param, AT);
	}
	if ( !(pj_dest = pj_init_plus(dest_param.c_str())) ) {
		throw InvalidArgumentException("Failed to initalize Proj with given arguments: "+dest_param, AT);
	}
}

void ProcTransformWindVector::TransformCoord(const double& X_in, const double& Y_in, double& X_out, double& Y_out)
{
	double x = (RACMO2) ? (Y_in) : (X_in);
	double y = (RACMO2) ? (X_in) : (Y_in);
	if (pj_is_latlong(pj_src)) {
		x *= Cst::to_rad;
		y *= Cst::to_rad;
	}
	const int p = pj_transform(pj_src, pj_dest, 1, 1, &x, &y, NULL );
	if (p!=0) throw ConversionFailedException("PROJ conversion failed: "+p, AT);
	X_out = x;
	Y_out = y;
	if (pj_is_latlong(pj_dest)) {
		X_out *= Cst::to_deg;
		Y_out *= Cst::to_deg;
	}
}

std::string ProcTransformWindVector::findUComponent(const MeteoData& md)
{
	if (md.param_exists("U")) return "U";
	else if (md.param_exists("VW_U")) return "VW_U";
	else if (md.param_exists("WIND_U")) return "WIND_U";
	
	return std::string();
}

std::string ProcTransformWindVector::findVComponent(const MeteoData& md)
{
	if (md.param_exists("V")) return "V";
	else if (md.param_exists("VW_V")) return "VW_V";
	else if (md.param_exists("WIND_V")) return "WIND_V";
	
	return std::string();
}

void ProcTransformWindVector::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	// Sanity check to see if the variable can be transformed
	if (ivec[0].getNameForParameter(param) != "DW"
	 && ivec[0].getNameForParameter(param) != "U"
	 && ivec[0].getNameForParameter(param) != "V"
	 && ivec[0].getNameForParameter(param) != "VW_U"
	 && ivec[0].getNameForParameter(param) != "VW_V"
	 && ivec[0].getNameForParameter(param) != "WIND_U"
	 && ivec[0].getNameForParameter(param) != "WIND_V") {
		throw InvalidArgumentException("Trying to use "+getName()+" filter on " + ivec[0].getNameForParameter(param) + " but it can only be applied to DW, or U and V, or VW_U and VW_V, or WIND_U and WIND_V!!", AT);
	}

	// Accuracy (length scale of projected vectors to determine wind direction)
	static const double eps = 1E-6;

	ovec = ivec;
	for (size_t ii=0; ii<ivec.size(); ii++){
		double u = IOUtils::nodata;
		double v = IOUtils::nodata;
		double uc = IOUtils::nodata;
		double vc = IOUtils::nodata;

		// Get coordinates
		const double lon = ivec[ii].meta.getPosition().getLon();
		const double lat = ivec[ii].meta.getPosition().getLat();

		// Check for wind speed components
		const std::string U_param( findUComponent(ivec[ii]) );
		const std::string V_param( findVComponent(ivec[ii]) );

		const bool has_components = !(U_param.empty() || V_param.empty());
		if (has_components) {
			if (ivec[ii](U_param) != IOUtils::nodata && ivec[ii](V_param) != IOUtils::nodata) {
				uc = ivec[ii](U_param);
				vc = ivec[ii](V_param);
			}
		}

		bool VWisZero = false;	// We use this flag to make sure that when DW is present and not nodata, but the wind speed is zero,
					// we can still transform, but not propagate an erroneous wind speed.
		const bool NotAtPoles = (lat > -90.+eps && lat < 90.-eps);
		if (param == MeteoData::DW && ivec[ii](MeteoData::DW) != IOUtils::nodata && NotAtPoles) {
			// If filter is applied on DW and DW is not nodata
			u = IOUtils::VWDW_TO_U(1., ivec[ii](MeteoData::DW));	// The filter may want to transform DW for VW == 0, so we calculate (u,v) assuming unity wind speed.
			v = IOUtils::VWDW_TO_V(1., ivec[ii](MeteoData::DW));
			if (ivec[ii](MeteoData::VW) != IOUtils::nodata) {
				if (ivec[ii](MeteoData::VW) != 0.) {
					u *= ivec[ii](MeteoData::VW);
					v *= ivec[ii](MeteoData::VW);
				} else {
					VWisZero = true;
				}
			}
		} else if (!has_components) {	// Otherwise, we try to see if both U and V exist
			if (param != MeteoData::DW) {
				// If both U and V do not exist, and the filter was applied on something else than DW, throw error
				throw InvalidArgumentException("Trying to use "+getName()+" filter on " + ivec[ii].getNameForParameter(param) + ", but not both components exist!!", AT);
			} else {
				if (!(NotAtPoles)) {
					throw ConversionFailedException("Trying to use "+getName()+" filter on " + ivec[ii].getNameForParameter(param) + " at latitude = " + IOUtils::toString(lat) + ", which leads to undefined results!!", AT);
				}
				// Otherwise we cannot do anything else
				continue;
			}
		} else {
			// We use the wind speed components when the filter is applied on them, or when MeteoData::DW is nodata *AND* wind speed components are available.
			u = uc;
			v = vc;
			if (u == 0. && v == 0.) continue;
		}
		if (u == IOUtils::nodata || v == IOUtils::nodata) continue;

		const double vw_old = u*u + v*v;	// For efficiency, we drop the sqrt.

		// Get easting and northing of point in target coordinate system (given by ivec[ii].meta.getPosition())
		double e0=0., n0=0.;
		// Note that we do not use the easting and northing from getPosition, since those may be a different coordinate system.
		TransformCoord(lon, lat, e0, n0);

		// Find ratio between m per degree latitude over m per degree longitude
		double et1=0., nt1=0., et2=0., nt2=0.;	// temporary variables to determine the ratio
		if (lat>0.) {
			TransformCoord(lon, lat-eps, et1, nt1);
		} else {
			TransformCoord(lon, lat+eps, et1, nt1);
		}
		if (lon>0.) {
			TransformCoord(lon-eps, lat, et2, nt2);
		} else {
			TransformCoord(lon+eps, lat, et2, nt2);
		}
		const double ratio = (et2!=e0 && nt2!=n0) ? (sqrt(  ((et1-e0)*(et1-e0) + (nt1-n0)*(nt1-n0)) / ((et2-e0)*(et2-e0) + (nt2-n0)*(nt2-n0))  )) : (1.);

		// Transform wind speed vector
		double e1, n1;		// end points of vector (start point is given by ivec[ii].meta.getPosition())
		double u_new, v_new;	// transformed wind speed components
		if (lat-(v*eps) >= -90. && lat-(v*eps) <= 90.) {
			if (lon-(u*eps*ratio)>=-360. && lon-(u*eps*ratio)<=360.) {
				TransformCoord(lon-(u*eps*ratio), lat-(v*eps), e1, n1);
				v_new = n0 - n1;
			} else {
				TransformCoord(lon+(u*eps*ratio), lat-(v*eps), e1, n1);
				v_new = n1 - n0;
			}
			u_new = e0 - e1;
		} else {
			if (lon-(u*eps*ratio)>=-360. && lon-(u*eps*ratio)<=360.) {
				TransformCoord(lon-(u*eps*ratio), lat+(v*eps), e1, n1);
				v_new = n0 - n1;
			} else {
				TransformCoord(lon+(u*eps*ratio), lat+(v*eps), e1, n1);
				v_new = n1 - n0;
			}
			u_new = e1 - e0;
		}

		// Calculate new wind direction
		const double dw_new = IOUtils::UV_TO_DW(u_new, v_new);

		// Assign transformed wind direction and wind speed components
		ovec[ii](MeteoData::DW) = dw_new;
		if (has_components) {
			const double scale = (VWisZero) ? (0.) : (sqrt(vw_old / (u_new * u_new + v_new * v_new)));
			u_new *= scale;
			v_new *= scale;
			ovec[ii](U_param) = u_new;
			ovec[ii](V_param) = v_new;
		}
	}
}

inline std::string readProjectionfromFile(const std::string in_filename, std::string root_path) {
	//if this is a relative path, prefix the path with the current path
	const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
	const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
	const std::string filename = path + "/" + FileUtils::getFilename(in_filename);
	const std::string str_out("");
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter TRANSFORMWINDVECTOR: ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	try {
		std::string line;
		getline(fin, line, eoln); //read complete line
		IOUtils::stripComments(line);
		IOUtils::trim(line);
		return line;
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) fin.close();
		throw;
	}
	return std::string("");
}

void ProcTransformWindVector::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg)
{
	const std::string where( "Filters::"+block_name );
	//parse the arguments (the keys are all upper case)
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="COORDPARAM") {
			IOUtils::parseArg(vecArgs[ii], where, t_coordparam);
			if (!isEPSG(t_coordparam)) {
				t_coordparam = readProjectionfromFile(t_coordparam, cfg.getConfigRootDir());
			}
		} else if (vecArgs[ii].first=="COORDPARAM_SRC") {
			IOUtils::parseArg(vecArgs[ii], where, s_coordparam);
			if (!isEPSG(s_coordparam)) {
				s_coordparam = readProjectionfromFile(s_coordparam, cfg.getConfigRootDir());
			}
		} else if (vecArgs[ii].first=="RACMO2") {
			IOUtils::parseArg(vecArgs[ii], where, RACMO2);
		}
	}

	initPROJ();
}
#endif
} //end namespace
