// SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/StationData.h>
#include <meteoio/IOUtils.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm> //for set_difference

#include <regex>

using namespace std;
namespace mio {

/************************************************************
 * static section                                           *
 ************************************************************/
const size_t MeteoGrids::nrOfParameters =  MeteoGrids::lastparam - MeteoGrids::firstparam + 1;
std::vector<std::string> MeteoGrids::paramname, MeteoGrids::description, MeteoGrids::units;
const bool MeteoGrids::__init = MeteoGrids::initStaticData();

bool MeteoGrids::initStaticData()
{
	//the order must be the same as in the enum
	paramname.push_back("TA");			description.push_back("Air Temperature");						units.push_back("K");
	paramname.push_back("RH");			description.push_back("Relative Humidity");					units.push_back("-");
	paramname.push_back("QI");			description.push_back("Specific Humidity");					units.push_back("g/g");
	paramname.push_back("TD");			description.push_back("Dew Point Temperature");			units.push_back("K");
	paramname.push_back("VW");			description.push_back("Wind Velocity");							units.push_back("m/s");
	paramname.push_back("DW");			description.push_back("Wind Direction");						units.push_back("°");
	paramname.push_back("VW_MAX");		description.push_back("Gust Wind Velocity");					units.push_back("m/s");
	paramname.push_back("ISWR");		description.push_back("Incoming Short Wave Radiation");				units.push_back("W/m2");
	paramname.push_back("RSWR");		description.push_back("Reflected Short Wave Radiation");				units.push_back("W/m2");
	paramname.push_back("ISWR_DIFF");	description.push_back("Incoming Diffuse Short Wave Radiation");	units.push_back("W/m2");
	paramname.push_back("ISWR_DIR");	description.push_back("Incoming Direct Short Wave Radiation");	units.push_back("W/m2");
	paramname.push_back("ILWR");		description.push_back("Incoming Long Wave Radiation");				units.push_back("W/m2");
	paramname.push_back("OLWR");		description.push_back("Outgoing Long Wave Radiation");				units.push_back("W/m2");
	paramname.push_back("LWR_NET");                 description.push_back("Net Long Wave Radiation");                               units.push_back("W/m2");
	paramname.push_back("TAU_CLD");		description.push_back("Atmospheric Transmissivity");					units.push_back("-");
	paramname.push_back("CLD");			description.push_back("Total cloud cover");					units.push_back("okta");
	paramname.push_back("HS");			description.push_back("Snow Height");							units.push_back("m");
	paramname.push_back("PSUM");		description.push_back("Precipitation Sum");					units.push_back("kg/m2");
	paramname.push_back("PSUM_PH");		description.push_back("Precipitation Phase");					units.push_back("-");
	paramname.push_back("PSUM_L");		description.push_back("Precipitation Sum of the Liquid Phase");	units.push_back("kg/m2");
	paramname.push_back("PSUM_LC");                  description.push_back("Precipitation Sum of the Liquid Convection"); units.push_back("kg/m2");
	paramname.push_back("PSUM_S");		description.push_back("Precipitation Sum of the Solid Phase");		units.push_back("kg/m2");
	paramname.push_back("TSG");			description.push_back("Ground Surface Temperature");					units.push_back("K");
	paramname.push_back("TSS");			description.push_back("Surface Temperature");				units.push_back("K");
	paramname.push_back("TSOIL");		description.push_back("Soil Temperature");					units.push_back("K");
	paramname.push_back("TSNOW");		description.push_back("Snow Temperature");					units.push_back("K");
	paramname.push_back("P");			description.push_back("Air Pressure");							units.push_back("Pa");
	paramname.push_back("P_SEA");		description.push_back("Air Pressure at Sea Level");		units.push_back("Pa");
	paramname.push_back("U");			description.push_back("Wind Velocity East Component");				units.push_back("m/s");
	paramname.push_back("V");			description.push_back("Wind Velocity North Component");			units.push_back("m/s");
	paramname.push_back("W");			description.push_back("Wind Velocity Vertical Component");			units.push_back("m/s");
	paramname.push_back("SWE");			description.push_back("Snow Water Equivalent");			units.push_back("kg/m2");
	paramname.push_back("RSNO");		description.push_back("Snow Mean Density");				units.push_back("kg/m3");
	paramname.push_back("ROT");			description.push_back("Total generated runoff");			units.push_back("gk/m2");
	paramname.push_back("ALB");			description.push_back("Albedo");									units.push_back("-");
	paramname.push_back("DEM");			description.push_back("Altitude above Sea Level");		units.push_back("m");
	paramname.push_back("SHADE");		description.push_back("Hillshade");								units.push_back("-");
	paramname.push_back("SLOPE");		description.push_back("Slope Angle");							units.push_back("degree");
	paramname.push_back("AZI");			description.push_back("Slope Aspect");							units.push_back("degree");

	return true;
}

const std::string MeteoGrids::getParameterName(const size_t& parindex)
{
	if (parindex >= MeteoGrids::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	return paramname[parindex];
}

const std::string MeteoGrids::getParameterDescription(const size_t& parindex, const bool& allow_ws)
{
	if (parindex >= MeteoGrids::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to get description for parameter that does not exist", AT);

	if (allow_ws) {
		return description[parindex];
	} else {
		std::string tmp( description[parindex] );
		IOUtils::replaceWhitespaces(tmp, '_');
		return tmp;
	}
}

const std::string MeteoGrids::getParameterUnits(const size_t& parindex)
{
	if (parindex >= MeteoGrids::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to get units for parameter that does not exist", AT);

	return units[parindex];
}

size_t MeteoGrids::getParameterIndex(const std::string& parname)
{
	for (size_t ii=0; ii<MeteoGrids::nrOfParameters; ii++) {
		if (paramname[ii] == parname)
			return ii;
	}

	return IOUtils::npos; //parameter not a part of MeteoGrids
}

/************************************************************
 * static section                                           *
 ************************************************************/
const double MeteoData::epsilon = 1e-5;
const size_t MeteoData::nrOfParameters =  MeteoData::lastparam - MeteoData::firstparam + 1;
std::vector<std::string> MeteoData::s_default_paramname;
const bool MeteoData::__init = MeteoData::initStaticData();
MeteoData::flag_field MeteoData::zero_flag; //for easy initialization

// Parameter helper
MeteoData::Parameters MeteoData::toParameter(const std::string& paramStr) {
	if (paramStr == "P") return MeteoData::Parameters::P;
	else if (paramStr == "TA") return MeteoData::Parameters::TA;
	else if (paramStr == "RH") return MeteoData::Parameters::RH;
	else if (paramStr == "QI") return MeteoData::Parameters::QI;
	else if (paramStr == "TSG") return MeteoData::Parameters::TSG;
	else if (paramStr == "TSOIL") return MeteoData::Parameters::TSOIL;
	else if (paramStr == "TSS") return MeteoData::Parameters::TSS;
	else if (paramStr == "TSNOW") return MeteoData::Parameters::TSNOW;
	else if (paramStr == "HS") return MeteoData::Parameters::HS;
	else if (paramStr == "VW") return MeteoData::Parameters::VW;
	else if (paramStr == "DW") return MeteoData::Parameters::DW;
	else if (paramStr == "VW_MAX") return MeteoData::Parameters::VW_MAX;
	else if (paramStr == "U") return MeteoData::Parameters::U;
	else if (paramStr == "V") return MeteoData::Parameters::V;
	else if (paramStr == "RSWR") return MeteoData::Parameters::RSWR;
	else if (paramStr == "ISWR") return MeteoData::Parameters::ISWR;
	else if (paramStr == "ILWR") return MeteoData::Parameters::ILWR;
	else if (paramStr == "OLWR") return MeteoData::Parameters::OLWR;
	else if (paramStr == "TAU_CLD") return MeteoData::Parameters::TAU_CLD;
	else if (paramStr == "PSUM") return MeteoData::Parameters::PSUM;
	else if (paramStr == "PSUM_PH") return MeteoData::Parameters::PSUM_PH;
	else if (paramStr == "PSUM_L") return MeteoData::Parameters::PSUM_L;
	else if (paramStr == "PSUM_LC") return MeteoData::Parameters::PSUM_LC;
	else if (paramStr == "PSUM_S") return MeteoData::Parameters::PSUM_S;
	else {
		// Handle the error case, could also throw an exception
		throw IOException("Unknown parameter: " + paramStr, AT);
	}
};

// Parameter helper
std::string MeteoData::parToString(const MeteoData::Parameters& param) {
	if (param == MeteoData::Parameters::P) return "P";
	else if (param == MeteoData::Parameters::TA) return "TA";
	else if (param == MeteoData::Parameters::RH) return "RH";
	else if (param == MeteoData::Parameters::QI) return "QI";
	else if (param == MeteoData::Parameters::TSG) return "TSG";
	else if (param == MeteoData::Parameters::TSOIL) return "TSOIL";
	else if (param == MeteoData::Parameters::TSS) return "TSS";
	else if (param == MeteoData::Parameters::TSNOW) return "TSNOW";
	else if (param == MeteoData::Parameters::HS) return "HS";
	else if (param == MeteoData::Parameters::VW) return "VW";
	else if (param == MeteoData::Parameters::DW) return "DW";
	else if (param == MeteoData::Parameters::VW_MAX) return "VW_MAX";
	else if (param == MeteoData::Parameters::U) return "U";
	else if (param == MeteoData::Parameters::V) return "V";
	else if (param == MeteoData::Parameters::RSWR) return "RSWR";
	else if (param == MeteoData::Parameters::ISWR) return "ISWR";
	else if (param == MeteoData::Parameters::ILWR) return "ILWR";
	else if (param == MeteoData::Parameters::OLWR) return "OLWR";
	else if (param == MeteoData::Parameters::TAU_CLD) return "TAU_CLD";
	else if (param == MeteoData::Parameters::PSUM) return "PSUM";
	else if (param == MeteoData::Parameters::PSUM_PH) return "PSUM_PH";
	else if (param == MeteoData::Parameters::PSUM_L) return "PSUM_L";
	else if (param == MeteoData::Parameters::PSUM_LC) return "PSUM_LC";
	else if (param == MeteoData::Parameters::PSUM_S) return "PSUM_S";
	else throw IOException("Unknown parameter: " + std::to_string(param), AT);
};

bool MeteoData::initStaticData()
{
	//Since the parameters enum starts at 0, this is enough to associate an index with its name
	s_default_paramname.push_back("P");
	s_default_paramname.push_back("TA");
	s_default_paramname.push_back("RH");
	s_default_paramname.push_back("QI");
	s_default_paramname.push_back("TSG");
	s_default_paramname.push_back("TSS");
	s_default_paramname.push_back("TSOIL");
	s_default_paramname.push_back("TSNOW");
	s_default_paramname.push_back("HS");
	s_default_paramname.push_back("VW");
	s_default_paramname.push_back("DW");
	s_default_paramname.push_back("VW_MAX");
	s_default_paramname.push_back("U");
	s_default_paramname.push_back("V");
	s_default_paramname.push_back("RSWR");
	s_default_paramname.push_back("ISWR");
	s_default_paramname.push_back("ILWR");
	s_default_paramname.push_back("OLWR");
	s_default_paramname.push_back("LWR_NET");
	s_default_paramname.push_back("TAU_CLD");
	s_default_paramname.push_back("PSUM");
	s_default_paramname.push_back("PSUM_PH");
	s_default_paramname.push_back("PSUM_L");
	s_default_paramname.push_back("PSUM_LC");
	s_default_paramname.push_back("PSUM_S");

	zero_flag.resampled = false; //init data qa flag that is initially set
	zero_flag.generated = false;
	zero_flag.filtered = false;
	zero_flag.created = false;
	zero_flag.offset = (float)IOUtils::nodata;

	return true;
}

const std::string& MeteoData::getParameterName(const size_t& parindex)
{
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);

	return MeteoData::s_default_paramname[parindex];
}

// check that the height is given in a proper format and return it as a number
static double parseInputHeight(const std::string& height_str)
{
	std::vector<std::string> num_parts(IOUtils::split(height_str, '.'));
	if (num_parts.size() > 2) throw InvalidArgumentException("Cannot parse decimal number: " + height_str);
	if (num_parts.size() == 2) {
		if (num_parts[1].size() != 2) throw InvalidArgumentException("Only exactly 2 decimal points is allowed, got: " + std::to_string(num_parts[1].size()) + ", in height: " + height_str);
		// Because we need to do back and forth conversion between string and number, and we cannot really differentiate between @1 and @1.00 then
		if (num_parts[1] == "00") throw InvalidArgumentException("Do not use .00 as decimals, as this will be ignored, and lead to potential bad behaviour. In " + height_str);
	}
	if (height_str.back() == '.' ) throw InvalidArgumentException("A height ending with \".\" is not allowed: " + height_str);
	return std::stod(height_str);
}

// parser for parameters of type TA@15 to be able to have multiple air temperatures...
// returns true if it is a known parameter that could be parsed
bool MeteoData::getTypeAndNo(const std::string& i_parname, std::string& parameter, double& number, const std::set<std::string>& additional_parameters)
{
	if (getStaticParameterIndex(i_parname) != IOUtils::npos) {
		parameter = i_parname;
		number = IOUtils::nodata;
		return true;
	} else if (additional_parameters.find(i_parname) != additional_parameters.end()) {
		parameter = i_parname;
		number = IOUtils::nodata;
		return true;
	} else if (i_parname.find("@") == std::string::npos) {
		// is not a known parameter without height, cannot be parsed
		parameter = i_parname;
		number = IOUtils::nodata;
		return false;
	}
	bool is_known_parameter = false;

	static const regex param_reg("^([A-Z]+)@(-?[0-9]{1,4}(\\.[0-9]{0,2})?)$", std::regex::icase | std::regex::optimize);
	std::smatch matches;

	parameter = "";
	number = IOUtils::nodata;

	if (std::regex_search(i_parname, matches, param_reg)) {
		const std::string param( matches[1].str() );
		if (getStaticParameterIndex(param) != IOUtils::npos) {
			is_known_parameter = true;
		} else if (additional_parameters.find(param) != additional_parameters.end()) {
			is_known_parameter = true;
		}
		parameter = param;
		number = parseInputHeight(matches[2].str());

		return is_known_parameter;
	}
	return is_known_parameter;
}

bool MeteoData::sameParameterType(const std::string& par1, const std::string& par2, const std::set<std::string>& additional_parameters)
{
	std::string parameter1;
	double number1;
	std::string parameter2;
	double number2;
	if (getTypeAndNo(par1, parameter1, number1, additional_parameters) && getTypeAndNo(par2, parameter2, number2, additional_parameters)) {
		return parameter1 == parameter2;
	}
	return false;
}

// casts to int or two decimal points when there is a decimal number, otherwise not allowed
std::string MeteoData::convertHeightToString(const double& height)
{
	std::string height_string;
	if (height == IOUtils::nodata) {
		std::cout << "Trying to convert NoData height to string, undefined behaviour!";
		return "";
	}
	if (height == static_cast<int>(height)) {
		// Height has no decimal part, cast to int
		height_string = std::to_string(static_cast<int>(height));
	} else {
		// Height has a decimal part, round to two decimal places
		std::stringstream stream;
		stream << std::fixed << std::setprecision(2) << height;
		height_string =  stream.str();
	}
	return height_string;
}

static bool comparePairOfHeightParam(const std::pair<double,std::string>& a, const std::pair<double,std::string>& b)
{
	return a.first < b.first;
}

std::vector<std::string> MeteoData::sortListByParams(const std::vector<std::string>& param_list, const std::set<std::string>& additional_parameters)
{
	std::map<std::string,std::vector<std::pair<double,std::string>>> params_with_heights;
	std::vector<std::string> sorted_list;
	for (const auto& par : param_list) {
		std::string parameter;
		double number;
		if (getTypeAndNo(par, parameter, number, additional_parameters)) {
			params_with_heights[parameter].push_back(std::make_pair(number, par));
		}
	}
	for (auto& entry : params_with_heights) {
		std::vector<std::pair<double,std::string>> par = entry.second;
		std::sort(par.begin(), par.end(), comparePairOfHeightParam);
		for (const auto& p : par) {
			sorted_list.push_back(p.second);
		}
	}
	return sorted_list;
}

std::vector<std::string> MeteoData::retrieveAllHeightsForParam(const std::set<std::string>& available_parameters, const std::string& param, const std::set<std::string>& additional_parameters)
{
	std::vector<std::string> param_with_heights;
	for (const auto& par : available_parameters) {
		if (sameParameterType(par, param, additional_parameters)) {
			param_with_heights.push_back(par);
		}
	}
	return sortListByParams(param_with_heights);
}

std::set<double> MeteoData::retrieveAllHeights(const std::set<std::string>& available_parameters, const std::set<std::string>& additional_parameters)
{
	std::set<double> heights;
	for (const auto& par : available_parameters) {
		std::string parameter;
		double number;
		if (getTypeAndNo(par, parameter, number, additional_parameters)) {
			heights.insert(number);
		}
	}
	return heights;
}

// no height: IOUtils::nodata
std::vector<std::string> MeteoData::retrieveAllParametersAtHeight(const std::set<std::string>& available_parameters, const double& height, const std::set<std::string>& additional_parameters)
{
	std::vector<std::string> param_at_height;
	for (const auto& par : available_parameters) {
		std::string parameter;
		double number;
		if (getTypeAndNo(par, parameter, number, additional_parameters)) {
			if (number == height) {
				param_at_height.push_back(par);
			}
		}
	}
	return param_at_height;
}

std::set<std::string> MeteoData::retriveUniqueParameters(const std::set<std::string>& available_parameters, const std::set<std::string>& additional_parameters)
{
	std::set<std::string> unique_parameters;
	for (const auto& par : available_parameters) {
		std::string parameter;
		double number;
		if (getTypeAndNo(par, parameter, number, additional_parameters)) {
			unique_parameters.insert(parameter);
		}
	}
	return unique_parameters;
}

size_t MeteoData::getStaticParameterIndex(const std::string& parname)
{
	for (size_t ii = 0; ii<MeteoData::nrOfParameters; ii++) {
		if (s_default_paramname[ii] == parname)
			return ii;
	}

	return IOUtils::npos; //parameter not a part of MeteoData
}

MeteoGrids::Parameters MeteoData::findGridParam(const MeteoData::Parameters& mpar)
{
	const std::string parname( MeteoData::getParameterName(mpar) );
	const size_t idx = MeteoGrids::getParameterIndex(parname);
	return static_cast<MeteoGrids::Parameters>(idx);
}

MeteoData::Parameters MeteoData::findMeteoParam(const MeteoGrids::Parameters& gpar)
{
	const std::string parname( MeteoGrids::getParameterName(gpar) );
	const size_t idx = MeteoData().getParameterIndex(parname);
	return static_cast<MeteoData::Parameters>(idx);
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

const std::string& MeteoData::getNameForParameter(const size_t& parindex) const
{
	if (parindex >= MeteoData::nrOfAllParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	if (parindex<MeteoData::nrOfParameters) return MeteoData::s_default_paramname[parindex];
	return extra_param_name[parindex-MeteoData::nrOfParameters];
}

bool MeteoData::param_exists(const std::string& i_paramname) const
{
	for (size_t ii = 0; ii<MeteoData::nrOfParameters; ii++) {
		if (s_default_paramname[ii] == i_paramname)
			return true;
	}

	for (size_t ii = 0; ii<extra_param_name.size(); ii++) {
		if (extra_param_name[ii] == i_paramname)
			return true;
	}

	return false;
}

size_t MeteoData::addParameter(const std::string& i_paramname)
{
	//check if name is already taken
	const size_t current_index = getParameterIndex(i_paramname);
	if (current_index != IOUtils::npos)
		return current_index; //do nothing, because parameter is already present

	//add parameter
	extra_param_name.push_back(i_paramname);
	data.push_back(IOUtils::nodata);
	flags.push_back(zero_flag);

	//Increase nrOfAllParameters
	nrOfAllParameters++;

	return (nrOfAllParameters - 1);
}

MeteoData::MeteoData()
         : date(0.0, 0.), meta(), extra_param_name(), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters),
           resampled(false), flags(MeteoData::nrOfParameters, zero_flag)
{ }

MeteoData::MeteoData(const Date& date_in)
         : date(date_in), meta(), extra_param_name(), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters),
           resampled(false), flags(MeteoData::nrOfParameters, zero_flag)
{ }

MeteoData::MeteoData(const StationData& meta_in)
         : date(0.0, 0.), meta(meta_in), extra_param_name(), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters),
           resampled(false), flags(MeteoData::nrOfParameters, zero_flag)
{ }

MeteoData::MeteoData(const Date& date_in, const StationData& meta_in)
         : date(date_in), meta(meta_in), extra_param_name(), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters),
           resampled(false), flags(MeteoData::nrOfParameters, zero_flag)
{ }

void MeteoData::reset()
{
	std::fill(data.begin(), data.end(), IOUtils::nodata);
}

/**
* @brief Standardize nodata values
* The plugin-specific nodata values are replaced by MeteoIO's internal nodata value
*/
void MeteoData::standardizeNodata(const double& plugin_nodata)
{
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		//loop through all meteo params and check whether they're nodata values
		if (data[ii] <= plugin_nodata)
			data[ii] = IOUtils::nodata;
	}
}

bool MeteoData::operator==(const MeteoData& in) const
{
	//An object is equal if the date is equal and all meteo parameters are equal
	if (date != in.date) return false;

	if (nrOfAllParameters != in.nrOfAllParameters) //the number of meteo parameters has to be consistent
		return false;

	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		if ( !IOUtils::checkEpsilonEquality(data[ii], in.data[ii], epsilon) ) return false;
	}

	return true;
}

double& MeteoData::operator()(const size_t& parindex)
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

const double& MeteoData::operator()(const size_t& parindex) const
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

double& MeteoData::operator()(const std::string& parname)
{
	const size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

const double& MeteoData::operator()(const std::string& parname) const
{
	const size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

size_t MeteoData::getParameterIndex(const std::string& parname) const
{
	for (size_t ii = 0; ii<MeteoData::nrOfParameters; ii++) {
		if (s_default_paramname[ii] == parname)
			return ii;
	}

	for (size_t ii=0; ii<extra_param_name.size(); ii++) {
		if (extra_param_name[ii] == parname)
			return ii+MeteoData::nrOfParameters;
	}

	return IOUtils::npos; //parameter not a part of MeteoData
}

// get all heights for parameter, default parameter is included by default, but given with nodata
std::vector<double> MeteoData::getHeightsForParameter(const std::string& in_parname, bool include_default, const std::set<std::string>& additional_parameters) const
{
	std::string parameter;
	double number;
	std::vector<double> heights;
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		getTypeAndNo(getNameForParameter(ii), parameter, number, additional_parameters);
		if (parameter == in_parname) {
			if (!include_default && number == IOUtils::nodata) {
				continue; // skip default parameters,
			}
			heights.push_back(number);
		}
	}
	return heights;
}



// counts the number of occurences of a reappearing parameter (default is excluded)
size_t MeteoData::getOccurencesOfParameter(const Parameters& par, const std::set<std::string>& additional_parameters) const
{
	const std::string in_parname( getParameterName(par) );
	std::string parameter;
	double number;
	size_t occurences = 0;
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		getTypeAndNo(getNameForParameter(ii), parameter, number, additional_parameters);
		if (parameter == in_parname) {
			if (number == IOUtils::nodata) {
				continue; // skip default parameters,
			}
			occurences++;
		}
	}
	return occurences;
}

bool MeteoData::isNodata() const
{
	for (size_t ii=0; ii<data.size(); ii++) {
		if (data[ii]!=IOUtils::nodata) return false;
	}

	return true;
}

void MeteoData::setFiltered(const size_t& param, const bool& in_filtered /* = true */)
{
	flags[param].filtered = in_filtered;
}

bool MeteoData::isFiltered(const size_t& param) const
{
	return flags[param].filtered;
}

void MeteoData::setGenerated(const size_t& param, const bool& in_generated /* = true */)
{
	flags[param].generated = in_generated;
}

bool MeteoData::isGenerated(const size_t& param) const
{
	return flags[param].generated;
}

void MeteoData::setResampledParam(const size_t& param, const bool& in_resampled /* = true */)
{
	flags[param].resampled = in_resampled;
}

bool MeteoData::isResampledParam(const size_t& param) const
{
	return flags[param].resampled;
}

const std::string MeteoData::toString(const FORMATS format) const {
	std::ostringstream os;

	if (format==DFLT) {
		os << "<meteo>\n";
		os << meta.toString();
		os << date.toString(Date::FULL) << "\n";
		os << setw(8) << nrOfAllParameters << " parameters\n";

		for (size_t ii=0; ii<nrOfAllParameters; ii++) {
			const double& value = operator()(ii);
			if (value != IOUtils::nodata)
				os << setw(8) << getNameForParameter(ii) << ":" << setw(15) << value << endl;
		}
		os << "</meteo>\n";
	} else if (format==FULL) {
		os << "<meteo>\n";
		os << meta.toString();
		os << date.toString(Date::FULL) << "\n";
		os << setw(8) << nrOfAllParameters << " parameters\n";

		for (size_t ii=0; ii<nrOfAllParameters; ii++) {
			os << setw(8) << getNameForParameter(ii) << ":" << setw(15) << operator()(ii) << endl;
		}
		os << "</meteo>\n";
	} else if (format==COMPACT) {
		os << "<meteo>\t";
		os << meta.stationID << " @ " << date.toString(Date::ISO) << " -> ";
		for (size_t ii=0; ii<nrOfAllParameters; ii++) {
			const double& value = operator()(ii);
			if (value != IOUtils::nodata)
				os <<  getNameForParameter(ii) << ":" << value << " ";
		}
		os << "</meteo>";
	}

	return os.str();
}

std::ostream& operator<<(std::ostream& os, const MeteoData& data) {
	os << data.date;
	os << data.meta;
	const size_t s_vector = data.extra_param_name.size();
	os.write(reinterpret_cast<const char*>(&s_vector), sizeof(size_t));
	for (size_t ii=0; ii<s_vector; ii++) {
		const size_t s_string = data.extra_param_name[ii].size();
		os.write(reinterpret_cast<const char*>(&s_string), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&data.extra_param_name[ii][0]), s_string*sizeof(data.extra_param_name[ii][0]));
	}

	const size_t s_data = data.data.size();
	os.write(reinterpret_cast<const char*>(&s_data), sizeof(s_data));
	os.write(reinterpret_cast<const char*>(&data.data[0]), s_data*sizeof(data.data[0]));

	os.write(reinterpret_cast<const char*>(&data.nrOfAllParameters), sizeof(data.nrOfAllParameters));
	os.write(reinterpret_cast<const char*>(&data.resampled), sizeof(data.resampled));
	return os;
}

std::istream& operator>>(std::istream& is, MeteoData& data) {
	is >> data.date;
	is >> data.meta;
	size_t s_vector;
	is.read(reinterpret_cast<char*>(&s_vector), sizeof(size_t));
	data.extra_param_name.resize(s_vector);
	for (size_t ii=0; ii<s_vector; ii++) {
		size_t s_string;
		is.read(reinterpret_cast<char*>(&s_string), sizeof(size_t));
		data.extra_param_name[ii].resize(s_string);
		is.read(reinterpret_cast<char*>(&data.extra_param_name[ii][0]), s_string*sizeof(data.extra_param_name[ii][0]));
	}

	size_t s_data;
	is.read(reinterpret_cast<char*>(&s_data), sizeof(size_t));
	data.data.resize(s_data);
	is.read(reinterpret_cast<char*>(&data.data[0]), s_data*sizeof(data.data[0]));

	is.read(reinterpret_cast<char*>(&data.nrOfAllParameters), sizeof(data.nrOfAllParameters));
	is.read(reinterpret_cast<char*>(&data.resampled), sizeof(data.resampled));
	return is;
}

MeteoData::Merge_Type MeteoData::getMergeType(std::string merge_type)
{
	IOUtils::toUpper( merge_type );
	if (merge_type=="STRICT_MERGE") return STRICT_MERGE;
	else if (merge_type=="EXPAND_MERGE") return EXPAND_MERGE;
	else if (merge_type=="FULL_MERGE") return FULL_MERGE;
	else if (merge_type=="WINDOW_MERGE") return WINDOW_MERGE;
	else
		throw UnknownValueException("Unknown merge type '"+merge_type+"'", AT);
}

MeteoData::Merge_Conflicts MeteoData::getMergeConflicts(std::string merge_conflicts)
{
	IOUtils::toUpper( merge_conflicts );
	if (merge_conflicts=="CONFLICTS_PRIORITY_FIRST") return CONFLICTS_PRIORITY_FIRST;
	else if (merge_conflicts=="CONFLICTS_PRIORITY_LAST") return CONFLICTS_PRIORITY_LAST;
	else if (merge_conflicts=="CONFLICTS_AVERAGE") return CONFLICTS_AVERAGE;
	else
		throw UnknownValueException("Unknown merge conflicts type '"+merge_conflicts+"'", AT);
}

/*
 * In the cases != STRICT_MERGE, it matters if vec2 is bigger than vec1. So we define the following indices
 * in order to store the information about the insertion positions:
 *
 * ----------------[-----------------]----------------------------	vec1
 *                 ↓                 ↓
 *             vec1_start        vec1_end
 *                 ↓                 ↓
 * ------[---------|-----------------|--------------]-------------	vec2
 */
size_t MeteoData::mergeTimeSeries(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const Merge_Type& strategy, const Merge_Conflicts& conflicts_strategy)
{
	if (vec2.empty()) return 0; //nothing to merge
	if ((strategy==STRICT_MERGE || strategy==WINDOW_MERGE) && vec1.empty()) return 0; //optimization for STRICT_MERGE
	if (vec1.empty()) {
		vec1 = vec2;
		return 0;
	}

	//adding the necessary extra parameters to vec1 elements, no matter which merge strategy
	if (!vec1.empty()) {
		const size_t nrExtra2 = vec2.back().nrOfAllParameters - nrOfParameters;
		for (size_t pp=0; pp<nrExtra2; pp++) {
			const std::string extra_name( vec2.back().extra_param_name[pp] );
			if (vec1.back().getParameterIndex(extra_name)==IOUtils::npos) {
				for (size_t ii=0; ii<vec1.size(); ii++) vec1[ii].addParameter( extra_name );
			}
		}
	}

	if (strategy==STRICT_MERGE || strategy==WINDOW_MERGE) { //optimization for STRICT_MERGE
		if (vec1.back().date<vec2.front().date) return 0; //vec1 is before vec2
		if (vec1.front().date>vec2.back().date) return 0; //vec1 is after vec2
	}

	size_t nr_conflicts = 0; //track the number of merge conflicts
	size_t vec1_start = 0; //the index in vec2 that matches the original start of vec1
	size_t vec1_end = 0; //the index in vec2 that matches the original end of vec1

	//filling data before vec1
	if (strategy!=STRICT_MERGE && strategy!=WINDOW_MERGE && vec1.front().date>vec2.front().date) {
		const Date start_date( vec1.front().date );
		vec1_start = vec2.size(); //if no overlap is found, take all vec2
		for(size_t ii=0; ii<vec2.size(); ii++) { //find the range of elements to add
			if (vec2[ii].date>=start_date) {
				vec1_start = ii;
				break;
			}
		}

		MeteoData md_pattern( vec1.front() ); //This assumes that station1 is not moving!
		md_pattern.reset(); //keep metadata and extra params
		vec1.insert(vec1.begin(), vec1_start, md_pattern);
		for (size_t ii=0; ii<vec1_start; ii++) {
			vec1[ii].date = vec2[ii].date;
			if (!vec1[ii].merge( vec2[ii], conflicts_strategy )) nr_conflicts++;
		}
	}

	//general case: merge one timestamp at a time
	if (strategy==FULL_MERGE || strategy==WINDOW_MERGE) {
		std::vector<MeteoData> tmp;
		tmp.reserve( vec1.size() + (vec2.size() - vec1_start)); //"worst case" scenario: all elements will be added
		MeteoData md_pattern( vec1.front() ); //This assumes that station1 is not moving!
		md_pattern.reset(); //keep metadata and extra params

		size_t idx2 = vec1_start; //all previous elements were handled before
		size_t last_v1 = vec1_start; //last element from vec1 that will have to be invalidated
		for(size_t ii=vec1_start; ii<vec1.size(); ii++) {
			const Date curr_date( vec1[ii].date );
			while ((idx2<vec2.size()) && (curr_date>vec2[idx2].date)) {
				tmp.push_back( md_pattern );
				tmp.back().date = vec2[idx2].date;
				if (!tmp.back().merge( vec2[idx2], conflicts_strategy )) nr_conflicts++; //so the extra params are properly handled
				idx2++;
			}
			if (idx2==vec2.size())  break; //nothing left to merge

			if (curr_date==vec2[idx2].date) {
				if (!vec1[ii].merge( vec2[idx2], conflicts_strategy )) nr_conflicts++;
				idx2++;
			}
			tmp.push_back( vec1[ii] );
			last_v1 = ii;
		}

		const size_t new_count = last_v1 - vec1_start + 1;
		if (new_count<tmp.size())
			vec1.insert( vec1.begin() + vec1_start, tmp.size()-new_count, tmp.front()); //so room for the extra params is allocated

		for(size_t ii=0; ii<tmp.size(); ii++)
			vec1[vec1_start+ii] = tmp[ii];

		vec1_end = idx2;
	} else {
		size_t idx2 = vec1_start;
		for (size_t ii=vec1_start; ii<vec1.size(); ii++) { //loop over the timestamps. If some elements were inserted, vec1 now starts at vec1_start. If not, vec1_start==0
			const Date curr_date( vec1[ii].date );
			while ((idx2<vec2.size()) && (curr_date>vec2[idx2].date)) idx2++;

			if (idx2==vec2.size()) return nr_conflicts; //nothing left to merge
			if (curr_date==vec2[idx2].date) { //merging
				if (!vec1[ii].merge( vec2[idx2], conflicts_strategy )) nr_conflicts++;
			}
		}
		vec1_end = idx2 + 1; //element at idx2 has already been merged
	}

	//filling data after vec1
	if (strategy!=STRICT_MERGE && strategy!=WINDOW_MERGE && vec1.back().date<vec2.back().date) {
		if (vec1_end<vec2.size()) {
			MeteoData md_pattern( vec1.back() ); //This assumes that station1 is not moving!
			md_pattern.reset(); //keep metadata and extra params
			for (size_t ii=vec1_end; ii<vec2.size(); ii++) {
				vec1.push_back( md_pattern );
				vec1.back().date = vec2[ii].date;
				if (!vec1.back().merge( vec2[ii], conflicts_strategy )) nr_conflicts++;
			}
		}
	}

	return nr_conflicts;
}

void MeteoData::merge(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const bool& simple_merge, const Merge_Conflicts& conflicts_strategy)
{
	if (vec2.empty()) return;

	if (simple_merge || vec1.empty()) {
		vec1.reserve( vec1.size()+vec2.size() );
		for (size_t ii=0; ii<vec2.size(); ii++) vec1.push_back( vec2[ii] );
	} else {
		for (size_t ii=0; ii<vec2.size(); ii++) merge(vec1, vec2[ii], conflicts_strategy);
	}
}

void MeteoData::merge(std::vector<MeteoData>& vec, const MeteoData& meteo2, const bool& simple_merge, const Merge_Conflicts& conflicts_strategy)
{
	if (!simple_merge) {
		for (size_t ii=0; ii<vec.size(); ii++) {
			//two stations are considered the same if they point to the same 3D position
			if (vec[ii].meta.position==meteo2.meta.position) {
				vec[ii].merge(meteo2, conflicts_strategy);
				return;
			}
		}
	}

	//the station was not found in the vector -> adding it
	vec.push_back( meteo2 );
}

void MeteoData::merge(std::vector<MeteoData>& vec, const Merge_Conflicts& conflicts_strategy)
{
	const size_t nElems = vec.size();
	if (nElems<2) return;

	std::vector<MeteoData> vecResult;
	std::vector<size_t> mergeIdx(nElems, 0);

	for (size_t ii=0; ii<nElems; ii++) {
		if (mergeIdx[ii]==IOUtils::npos) continue; //this element has already been merged, skip
		for (size_t jj=ii+1; jj<nElems; jj++) {
			if (vec[ii].meta.position==vec[jj].meta.position) {
				vec[ii].merge( vec[jj], conflicts_strategy );
				mergeIdx[jj]=IOUtils::npos; //this element will be skipped in the next loops
			}
		}
		vecResult.push_back( vec[ii] );
	}

	vec.swap( vecResult );
}

MeteoData MeteoData::merge(MeteoData meteo1, const MeteoData& meteo2, const Merge_Conflicts& conflicts_strategy)
{
	meteo1.merge(meteo2, conflicts_strategy);
	return meteo1;
}

bool MeteoData::merge(const MeteoData& meteo2, const Merge_Conflicts& conflicts_strategy)
{
	if (!date.isUndef() && !meteo2.date.isUndef() && date!=meteo2.date) {
		//the data must be time synchronized!
		std::ostringstream ss;
		ss << "Trying to merge MeteoData at " << date.toString(Date::ISO);
		ss << " with MeteoData at " << meteo2.date.toString(Date::ISO);
		throw InvalidArgumentException(ss.str(), AT);
	}

	meta.merge(meteo2.meta); //no brainer merging of metadata
	if (date.isUndef()) date=meteo2.date; //we don't accept different dates, see above
	if (meteo2.resampled==true ) resampled=true;

	if (conflicts_strategy==CONFLICTS_PRIORITY_FIRST) {
		//merge standard parameters
		for (size_t ii=0; ii<nrOfParameters; ii++) {
			if (data[ii]==IOUtils::nodata) {
				data[ii] = meteo2.data[ii];
				flags[ii] = meteo2.flags[ii];
			}
		}

		//for each meteo2 extra parameter, check if a matching parameter exist
		const size_t nrExtra2 = meteo2.nrOfAllParameters - nrOfParameters;
		for (size_t ii=0; ii<nrExtra2; ii++) {
			const std::string extra_name( meteo2.extra_param_name[ii] );
			const size_t extra_param_idx = getParameterIndex(extra_name);
			if (extra_param_idx==IOUtils::npos) { //no such parameter in current object
				const size_t new_idx = addParameter( extra_name );
				data[new_idx] = meteo2.data[nrOfParameters+ii];
				flags[new_idx] = meteo2.flags[nrOfParameters+ii];
			} else if (data[extra_param_idx]==IOUtils::nodata) {
				data[extra_param_idx] = meteo2.data[nrOfParameters+ii];
				flags[extra_param_idx] = meteo2.flags[nrOfParameters+ii];
			}
		}
		return true;
	} else if (conflicts_strategy==CONFLICTS_PRIORITY_LAST) {
		//merge standard parameters
		for (size_t ii=0; ii<nrOfParameters; ii++) {
			if (meteo2.data[ii]!=IOUtils::nodata) {
				data[ii] = meteo2.data[ii];
				flags[ii] = meteo2.flags[ii];
			}
		}

		//for each meteo2 extra parameter, check if a matching parameter exist
		const size_t nrExtra2 = meteo2.nrOfAllParameters - nrOfParameters;
		for (size_t ii=0; ii<nrExtra2; ii++) {
			const std::string extra_name( meteo2.extra_param_name[ii] );
			const size_t extra_param_idx = getParameterIndex(extra_name);
			if (extra_param_idx==IOUtils::npos) { //no such parameter in current object
				const size_t new_idx = addParameter( extra_name );
				data[new_idx] = meteo2.data[nrOfParameters+ii];
				flags[new_idx] = meteo2.flags[nrOfParameters+ii];
			} else if (meteo2.data[nrOfParameters+ii]!=IOUtils::nodata) {
				data[extra_param_idx] = meteo2.data[nrOfParameters+ii];
				flags[extra_param_idx] = meteo2.flags[nrOfParameters+ii];
			}
		}
		return true;
	} else if (conflicts_strategy==CONFLICTS_AVERAGE) {
		bool has_no_conflicts = true;
		//merge standard parameters
		for (size_t ii=0; ii<nrOfParameters; ii++) {
			if (data[ii]==IOUtils::nodata) {
				data[ii] = meteo2.data[ii];
				flags[ii] = meteo2.flags[ii];
			} else if (meteo2.data[ii]!=IOUtils::nodata && data[ii]!=meteo2.data[ii]) {
				data[ii] = .5 * (data[ii] + meteo2.data[ii]);
				flags[ii].resampled = true;
				has_no_conflicts = false;
			}
		}

		//for each meteo2 extra parameter, check if a matching parameter exist
		const size_t nrExtra2 = meteo2.nrOfAllParameters - nrOfParameters;
		for (size_t ii=0; ii<nrExtra2; ii++) {
			const std::string extra_name( meteo2.extra_param_name[ii] );
			const size_t extra_param_idx = getParameterIndex(extra_name);
			if (extra_param_idx==IOUtils::npos) { //no such parameter in current object
				const size_t new_idx = addParameter( extra_name );
				data[new_idx] = meteo2.data[nrOfParameters+ii];
				flags[new_idx] = meteo2.flags[nrOfParameters+ii];
			} else {
				if (data[extra_param_idx]==IOUtils::nodata) {
					data[extra_param_idx] = meteo2.data[nrOfParameters+ii];
					flags[extra_param_idx] = meteo2.flags[nrOfParameters+ii];
				} else if (meteo2.data[nrOfParameters+ii]!=IOUtils::nodata  && data[extra_param_idx]!=meteo2.data[nrOfParameters+ii]) {
					data[extra_param_idx] = .5 * (data[extra_param_idx] + meteo2.data[nrOfParameters+ii]);
					flags[nrOfParameters+ii].resampled = true;
					has_no_conflicts = false;
				}
			}
		}

		return has_no_conflicts;
	}

	return true;
}

bool MeteoData::hasConflicts(const MeteoData& meteo2) const
{
	if (!date.isUndef() && !meteo2.date.isUndef() && date!=meteo2.date) return true;

	//check for conflicts in standard parameters
	for (size_t ii=0; ii<nrOfParameters; ii++) {
		if (data[ii]!=IOUtils::nodata && meteo2.data[ii]!=IOUtils::nodata && data[ii]!=meteo2.data[ii])
			return true;
	}

	//check for conflicts in extra parameters
	const size_t nrExtra2 = meteo2.nrOfAllParameters - nrOfParameters;
	for (size_t ii=0; ii<nrExtra2; ii++) {
		const std::string extra_name = meteo2.extra_param_name[ii];
		const size_t extra_param_idx = getParameterIndex(extra_name);
		if (extra_param_idx==IOUtils::npos) continue;
		if (data[extra_param_idx]!=IOUtils::nodata && meteo2.data[ii]!=IOUtils::nodata && data[extra_param_idx]!=meteo2.data[ii])
			return true;
	}

	return false;
}

std::set<std::string> MeteoData::listAvailableParameters(const std::vector<MeteoData>& vecMeteo)
{
	std::set<std::string> results;

	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		for (const std::string& parname : MeteoData::s_default_paramname)
			if (vecMeteo[ii](parname) != IOUtils::nodata) results.insert( parname );

		for (const std::string& parname : vecMeteo[ii].extra_param_name)
			if (vecMeteo[ii](parname) != IOUtils::nodata) results.insert( parname );
	}

	return results;
}

size_t MeteoData::listUnknownParameters(const std::set<std::string>& additional_params) const {
	size_t result = 0;
	std::string paramter;
	double number;
	for (const std::string& parname : extra_param_name) {
		if (!getTypeAndNo(parname, paramter, number, additional_params))
			result++;
	}
	return result;
}

void MeteoData::unifyMeteoData(METEO_SET &vecMeteo)
{
	const size_t nElems = vecMeteo.size();
	if (nElems<2) return;

	std::vector<std::string> extra_params_ref( vecMeteo.front().extra_param_name );

	for (size_t ii=0; ii<nElems; ii++) {
		//easy case: exactly the same vectors
		if (vecMeteo[ii].extra_param_name == extra_params_ref) continue;

		//maybe the same parameters are present, but in a different order?
		std::set<std::string> ref_params(extra_params_ref.begin(), extra_params_ref.end());
		std::set<std::string> new_params(vecMeteo[ii].extra_param_name.begin(), vecMeteo[ii].extra_param_name.end());
		if (ref_params == new_params) { //yes, it is only in a different order
			//most probably the new order will remain from now on
			extra_params_ref = vecMeteo[ii].extra_param_name;
			continue;
		}

		//now comes the hard work: we need set the whole vector to the same extra_param_name
		//first, we add all new elements to the begining of the vector until now
		std::vector<std::string> new_elems;
		std::set_difference(new_params.begin(), new_params.end(), ref_params.begin(), ref_params.end(), std::inserter(new_elems, new_elems.end()));
		for (size_t jj=0; jj<ii; jj++) {
			for(const auto &new_param : new_elems) vecMeteo[jj].addParameter( new_param );
		}

		//then, we add all past elements to the rest of the vector starting now and until the end
		std::vector<std::string> past_elems;
		std::set_difference(ref_params.begin(), ref_params.end(), new_params.begin(), new_params.end(), std::inserter(past_elems, past_elems.end()));
		for (size_t jj=ii; jj<nElems; jj++) {
			for(const auto &new_param : past_elems) vecMeteo[jj].addParameter( new_param );
		}

		//reset the reference parameters list to contain all potentially new elements
		extra_params_ref = vecMeteo[ii].extra_param_name;
	}
}


} //namespace
