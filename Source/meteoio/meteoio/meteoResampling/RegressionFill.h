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
#ifndef RESAMPLINGREGRESSIONFILL_H
#define RESAMPLINGREGRESSIONFILL_H

#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <unordered_map>

namespace mio {

/**
 * @brief Fill missing values from a station, using another station as reference.
 * 
 * @details
 * For this algorithm, a supportive station is needed for each station, which is used to fill missing data. 
 * Using the data from the support station, a linear regression is performed, and with this linear relation between the two stations, 
 * the missing values are filled using the data from the second station.
 * 
 * @note As the resampling algorithm is specified for each parameter, but not station, when the parameter of a station for which no support station 
 * is specified is tried to be resampled, it will not do anything, so it is advised to set linear as a backup algorithm in the resampling stack.
 * 
 * @section Parameters
 * - VERBOSE: If set to true, the algorithm will print out information about the state of the resampling.
 * - ADDITIONAL_STATIONS#: Specifiy which station should be used to support which station. In the form of StationID1--StationID2, where 2 will be used to fill in 1. To specify multiple stations, use ADDITIONAL_STATIONS1, ADDITIONAL_STATIONS2, etc.
 *                          It is also possible to use "CLOSEST" as the second station, which will use the closest station to the first one.
 * - TYPE: The type of regression to be used. Currently not implemented, i.e. only LINEAR is implemented.
 * 
 * @code
 * [Interpolations1D]
 * TA::resample         = REGFILL
 * TA::REGFILL::VERBOSE = true
 * TA::ADDITIONAL_STATIONS1 = StationID1--StationID2
 * TA::ADDITIONAL_STATIONS2 = StationID2--WFJ2
 * @endcode
 *
 * @author Patrick Leibersperger
 * @date 2024-03-20
 * 
 * @todo Implement other regression types.
 * @todo Make it possible to use multiple stations to support a single station.
 * 
 */
class RegressionFill : public ResamplingAlgorithms {
	public:
        enum RegressionType {
            LINEAR,
            QUADRATIC,
            CUBIC
        };

		RegressionFill(const std::string& i_algoname, const std::string& i_parname, const double& dflt_max_gap_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const std::string& stationHash, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md) override;
        void resample(const std::string& stationHash, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                const std::vector<MeteoData>& vecM, MeteoData& md, const std::vector<METEO_SET>& additional_stations);
		std::string toString() const override;

        bool findValueAt(const std::vector<MeteoData>& support_station, const Date& date, const size_t& paramindex, double& value);
        void getRegressionData(const size_t index, const size_t paramindex, const std::vector<MeteoData>& vecM, const std::vector<METEO_SET>& additional_stations, std::vector<double>& x, std::vector<double>& y, std::vector<Date>& dates);

        double linear(double julian_date, const std::vector<double>& coefficients);

    private: 
        std::unordered_map<size_t, std::vector<double>> regression_coefficients;
        bool verbose;
        RegressionType reg_type;

        // flag
        bool printed_info = false;
};

} //end namespace mio

#endif
