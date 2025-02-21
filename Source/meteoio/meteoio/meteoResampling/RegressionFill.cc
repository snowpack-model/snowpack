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

#include <meteoio/meteoResampling/RegressionFill.h>
#include <meteoio/meteoStats/libfit1D.h>

#include <sstream>

namespace mio {

RegressionFill::RegressionFill(const std::string& i_algoname, const std::string& i_parname, const double& dflt_max_gap_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
             : ResamplingAlgorithms(i_algoname, i_parname, dflt_max_gap_size, vecArgs), regression_coefficients(), verbose(false), reg_type(LINEAR)
{
	const std::string where( "Interpolations1D::"+i_parname+"::"+i_algoname );
	for (const auto& arg : vecArgs){
		if (arg.first == "VERBOSE") {
			IOUtils::parseArg(arg, where, verbose);
		} else if (arg.first == "TYPE") {
			std::string type;
			IOUtils::parseArg(arg, where, type);
			if (type == "LINEAR") {
				reg_type = LINEAR;
			} else {
				throw IOException("Regression type not implemented", AT);
			}
		} else if (arg.first.find("ADDITIONAL_STATIONS") != std::string::npos) continue;
		else {
			throw IOException("Unknown argument for RegressionFill", AT);
		}
	}
}

// ------------------------- HELPER -----------------------
std::string RegressionFill::toString() const
{
	//this should help when debugging, so output relevant parameters for your algorithm
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	ss << std::endl;
	ss << std::right << std::setw(10) << " " << std::left << std::setw(15) << "VERBOSE" << verbose << std::endl;
	ss << std::right << std::setw(10) << " " << std::left << std::setw(15) << "TYPE" << (reg_type == LINEAR ? "LINEAR" : "UNKNOWN") << std::endl;
	ss << std::right << std::setw(10) << " " << std::left << std::setw(15) << "Computed Coefficients" << regression_coefficients.size() << std::endl;
	return ss.str();
}

bool RegressionFill::findValueAt(const std::vector<MeteoData>& support_station, const Date& date, const size_t& paramindex, double& value)
{
	for (const MeteoData& md : support_station) {
		if ((md.date-date).getJulian(true)>max_gap_size || (md.date-date).getJulian(true)<(-1*max_gap_size)) continue;

		if (md.date == date) {
			value = md(paramindex);
			return true;
		}
	}
	return false;
}

void RegressionFill::getRegressionData(const size_t index, const size_t paramindex, const std::vector<MeteoData>& vecM, const std::vector<METEO_SET>& additional_stations,
                       std::vector<double>& x, std::vector<double>& y, std::vector<Date>& dates) {
    for (size_t i = 0; i < index; ++i) {
        if (vecM[i](paramindex) != IOUtils::nodata) {
            Date date = vecM[i].date;
            double x_val;
            if (findValueAt(additional_stations[0], date, paramindex, x_val)) {
                y.push_back(vecM[i](paramindex));
                x.push_back(x_val);
                dates.push_back(date);
            }
        }
    }
}

// ------------------------- REGRESSION -----------------------
double RegressionFill::linear(double julian_date, const std::vector<double>& coefficients)
{
	return coefficients[1] + coefficients[0] * julian_date;
}


// ------------------------- MAIN FUNCTION -----------------------
void RegressionFill::resample(const std::string& /* stationHash */, const size_t& /* index */, const ResamplingPosition& /* position */, const size_t& /* paramindex */,
                            const std::vector<MeteoData>& /* vecM */, MeteoData& /* md */) {
								return;
							}



// TODO: Possibility to cache the models directly, and easily support different kinds of regression
// TODO: Support giving multiple stations as support, and then use the one with the most data
void RegressionFill::resample(const std::string& /*stationHash*/, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md, const std::vector<METEO_SET>& additional_stations)
{
	if (index >= vecM.size()) throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value;
			return;
		}
	}

	if (additional_stations.empty()) throw IOException("The Regression Fill needs additional stations to work properly. Make sure to use the correct Station IDs", AT);
	if (additional_stations.size() != 1) throw IOException("The Regression Fill needs exactly one additional station to work properly", AT);

	if (verbose && !printed_info) {
		std::cout << "RegressionFill: Using station " << additional_stations[0].front().meta.getStationID() << " as support station for station " << vecM[index].meta.getStationID() << std::endl;
		printed_info = true;
	}


	double new_x_val;
	if (!findValueAt(additional_stations[0], md.date, paramindex, new_x_val)) {
		if (verbose) std::cout << "RegressionFill: could not find value at required timestamp in support station" << std::endl;
		return;
	}

	// if already fitted for parameter, use the cached coefficients
	if (regression_coefficients.find(paramindex) != regression_coefficients.end()) {
		md(paramindex) = linear(new_x_val, regression_coefficients[paramindex]); // TODO: do i need to convert to gmt?
		return;
	}

	// get the regression data before the missing value
	std::vector<double> x, y;
	std::vector<Date> dates;
	getRegressionData(index, paramindex, vecM, additional_stations, x, y, dates);

	if (x.size() < 2) {
		if (verbose) std::cout << "RegressionFill: Did not find enough data to perform regression" << std::endl;
		return;
	};

	// do the regression
	LinearLS model = LinearLS();
	model.setData(x, y);
	model.fit();

	// cache and fill
	regression_coefficients[paramindex] = model.getParams();

	md(paramindex) = model.f(new_x_val);
	return;
}

} //namespace
