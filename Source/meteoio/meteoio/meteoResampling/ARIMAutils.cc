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

#include "ARIMAutils.h"
#include <cmath>
#include <unordered_map>
#include <numeric>
#include <algorithm>

static const double NODATA = 1e-20; // needs to be 0 so arima doesnt care

namespace mio {

namespace ARIMAutils {
const std::map<ObjectiveFunction,std::string> ObjectiveFunctionMap = {
            {CSS_MLE, "CSS_MLE"},
            {MLE, "MLE"},
            {CSS, "CSS"},
        };
const std::map<OptimizationMethod, std::string> OptimizationMethodMap = {
            {Nelder_Mead, "Nelder_Mead"},
            {Newton_Line_Search, "Newton_Line_Search"},
            {Newton_Trust_Region_Hook_Step, "Newton_Trust_Region_Hook_Step"},
            {Newton_Trust_Region_Double_Dog_Leg, "Newton_Trust_Region_Double_Dog_Leg"},
            {Conjugate_Gradient, "Conjugate_Gradient"},
            {BFGS, "BFGS"},
            {LBFGS, "LBFGS"},
            {BFGS_MTM, "BFGS_MTM"},
        };

Normalization::Normalization(): mean(0), std(0), min(0), max(0) {}

Normalization::Normalization(std::vector<double>& data): mean(calcVecMean(data)), std(stdDev(data)), min(*std::min_element(data.begin(), data.end())), max(*std::max_element(data.begin(), data.end())) {};

Normalization::Normalization(std::vector<double>& data, Mode new_mode): mean(calcVecMean(data)), std(stdDev(data)), min(*std::min_element(data.begin(),data.end())), max(*std::max_element(data.begin(), data.end())), mode(new_mode) {}

std::vector<double> Normalization::normalize(const std::vector<double>& data) {
	std::vector<double> normalizedData = data; // Create a copy of data
	if (mode == Mode::Nothing) {
		for (size_t i = 0; i < normalizedData.size(); i++) {
			if (normalizedData[i] == IOUtils::nodata)
				normalizedData[i] = 0.0;
		}
		return normalizedData;
	}
	if (mode == Mode::ZScore) {
		for (size_t i = 0; i < normalizedData.size(); i++) {
			if (normalizedData[i] != IOUtils::nodata)
				normalizedData[i] = (normalizedData[i] - mean) / std;
			else	
				normalizedData[i] = NODATA;
		}
	} else if (mode == Mode::MinMax) {
		for (size_t i = 0; i < normalizedData.size(); i++) {
			if (normalizedData[i] != IOUtils::nodata)
				normalizedData[i] = (normalizedData[i] - min) / (max - min);
			else
				normalizedData[i] = NODATA;
		}
	}
	return normalizedData; // Return the modified data
}

std::vector<double> Normalization::denormalize(const std::vector<double>& data) {
	if (std::all_of(data.begin(), data.end(), [](double val) { return val == 0.0; })) {
		return data; // we do not want to accidentally set a vector of zeros to the mean (we do not use random walks)
	}

	std::vector<double> denormalizedData = data; // Create a copy of data
	if (mode == Mode::ZScore) {
		for (size_t i = 0; i < denormalizedData.size(); i++) {
			if (denormalizedData[i] != NODATA)
				denormalizedData[i] = denormalizedData[i] * std + mean;
			else
				denormalizedData[i] = IOUtils::nodata;
		}
	} else if (mode == Mode::MinMax) {
		for (size_t i = 0; i < denormalizedData.size(); i++) {
			if (denormalizedData[i] != NODATA)
				denormalizedData[i] = denormalizedData[i] * (max - min) + min;
			else
				denormalizedData[i] = IOUtils::nodata;
		}
	}
	return denormalizedData; // Return the modified data
}


// slice a vector from start to start+N
std::vector<double> slice(const std::vector<double>& vec, size_t start, size_t N) { 
	assert(start + N < vec.size()); // Ensure the range is valid
	std::vector<double> vec_sliced;
	vec_sliced.assign(vec.begin() + start, vec.begin() + start + N);
	return vec_sliced;
}

// slice a vector from start to end
std::vector<double> slice(const std::vector<double>& vec, size_t start) { 
	assert(start < vec.size()); // Ensure the range is valid
	std::vector<double> vec_sliced;
	vec_sliced.assign(vec.begin() + start, vec.end());
	return vec_sliced;
}

// np.arange for c++
std::vector<double> arange(size_t start, size_t N) {
    std::vector<double> vec(N);
    for (size_t i = 0; i < N; i++) {
        vec[i] = static_cast<double>(start + i);
    }
    return vec;
}

//calculate the mean of a vector
double calcVecMean(const std::vector<double>& vec) {
	double sum = 0.0;
	int count = 0;
	for (const auto& val : vec) {
		if (val != IOUtils::nodata) {
			sum += val;
			++count;
		}
	}
	return count > 0 ? sum / static_cast<double>(count) : 0;
}

//calculate the standard deviation of a vector
double stdDev(const std::vector<double>& vec) {
	double mean_vec = calcVecMean(vec);
	double sum = 0;
	int count = 0;
	for (const auto& val : vec) {
		if (val != IOUtils::nodata) {
			sum += (val - mean_vec)*(val - mean_vec);
			++count;
		}
	}
	return count > 0 ? std::sqrt(sum/static_cast<double>(count)) : 0;
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(const std::vector<MeteoData>& vecM, const std::string &paramname) {
    size_t paramindex = vecM[0].getParameterIndex(paramname);
    std::vector<double> vec(vecM.size());
    std::transform(vecM.begin(), vecM.end(), vec.begin(), [paramindex](const MeteoData& data) {
        return data(paramindex);
    });
    return vec;
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(const std::vector<MeteoData>& vecM, const size_t &paramindex) {
    std::vector<double> vec(vecM.size());
    std::transform(vecM.begin(), vecM.end(), vec.begin(), [paramindex](const MeteoData& data) {
        return data(paramindex);
    });
    return vec;
}

// helper to parse direction argument for interpolarima
std::vector<double> decideDirection(const std::vector<double>& data, const std::string& direction, bool forward, size_t gap_loc, size_t length) {
	if (!forward) {
		return std::vector<double>(5, 0.0);
	}

	if (direction == "forward") {
		return slice(data, 0, gap_loc);
	} 

	if (direction == "backward") {
		auto reversedData = reverseVectorReturn(data);
		return slice(reversedData, 0, reversedData.size() - length);
	} 

	throw mio::IOException("Direction " + direction + " not recognized");
}

//return true if a valid point could be found backward from pos
size_t searchBackward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size)
{
	const Date windowStart( resampling_date - i_window_size );
	const bool knownGap = (!last_gap.startDate.isUndef() && !last_gap.endDate.isUndef());
	const bool currentInGap = (knownGap)? (resampling_date>=last_gap.startDate && resampling_date<=last_gap.endDate) : false;
	const bool windowStartInGap = (knownGap)? (windowStart>=last_gap.startDate && windowStart<=last_gap.endDate) : false;
	
	//the current point and window start are in a known gap, there is no hope
	if (currentInGap && windowStartInGap) return IOUtils::npos;
	
	//the current point is NOT in a known gap
	if (!currentInGap) { //or !knownGap
		const Date dateStart = (windowStartInGap)? last_gap.endDate : windowStart;
		size_t ii = pos; //because idx will get decremented right away
		for (; ii-- >0; ) {
			if (vecM[ii].date < dateStart) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) {
				last_gap.setStart(ii, vecM);
				last_gap.setEnd(pos, vecM);
				return ii;
			}
		}
		
		//no valid point found
		if (windowStartInGap) { //we can extend the current known gap
			last_gap.extend(pos, vecM);
		} else { //this is a new gap
			last_gap.setStart(ii, vecM);
			last_gap.setEnd(pos, vecM);
		}
		return IOUtils::npos;
	} else { //what's left: the current point is in a known gap, but there might be some data before
		size_t ii = last_gap.start; //start from the begining of the last known gap (and idx will be decremented right away)
		for (; ii-- >0; ) {
			if (vecM[ii].date < windowStart) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) { //there is some data before the gap!
				last_gap.extend(ii, vecM);
				return ii;
			}
		}
		
		last_gap.extend(ii, vecM);
		return IOUtils::npos;
	}
}

//return true if a valid point could be found forward from pos
size_t searchForward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, const size_t& indexP1)
{
	const Date windowEnd = (indexP1 != IOUtils::npos)? vecM[indexP1].date+i_window_size : resampling_date+i_window_size;
	const bool knownGap = (!last_gap.startDate.isUndef() && !last_gap.endDate.isUndef());
	const bool currentInGap = (knownGap)? (resampling_date>=last_gap.startDate && resampling_date<=last_gap.endDate) : false;
	const bool windowEndInGap = (knownGap)? (windowEnd>=last_gap.startDate && windowEnd<=last_gap.endDate) : false;
	
	//the current point and window start are in a known gap, there is no hope
	if (currentInGap && windowEndInGap) return IOUtils::npos;
	
	//the current point is NOT in a known gap
	if (!currentInGap) { //or !knownGap
		const Date dateEnd = (windowEndInGap)? last_gap.startDate : windowEnd;
		size_t ii = pos;
		for (; ii<vecM.size(); ++ii) {
			if (vecM[ii].date > dateEnd) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) {
				last_gap.setStart(pos, vecM);
				last_gap.setEnd(ii-1, vecM);
				return ii;
			}
		}

		if (ii == pos) {
			if (vecM[ii-1].date <= dateEnd && vecM[ii].date >dateEnd) {
				if (windowEndInGap) { //we can extend the current known gap
					last_gap.extend(pos-1, vecM);
				} else { //this is a new gap
					last_gap.setStart(pos-1, vecM);
					last_gap.setEnd(ii, vecM);
				}				
			}
		}
		
		if (windowEndInGap) { //we can extend the current known gap
			last_gap.extend(pos, vecM);
		} else { //this is a new gap
			last_gap.setStart(pos, vecM);
			last_gap.setEnd(ii-1, vecM);
		}
		return IOUtils::npos;
	} else { //what's left: the current point is in a known gap, but there might be some data after
		size_t ii = last_gap.end;
		for (; ii<vecM.size(); ++ii) { //start from the end of the last known gap
			if (vecM[ii].date > windowEnd) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) { //there is some data after the gap!
				last_gap.extend(ii-1, vecM);
				return ii;
			}
		}
		
		last_gap.extend(ii-1, vecM);
		return IOUtils::npos;
	}
}

bool requal(const Date &date1, const Date &date2) {
    double tolerance = DATE_TOLERANCE; // Define your tolerance level
    bool is_equal = std::abs((date1 - date2).getJulian(true)) <= tolerance;
	return is_equal;
}

static void adjustDataStartDate(ARIMA_GAP &last_gap, const std::vector<MeteoData>& vecM, const Date& resampling_date, Date& data_start_date, Date& data_end_date) {
    if (data_start_date < vecM[0].date && last_gap.startDate == resampling_date) {
        data_start_date = adjustStartDate(vecM, last_gap, resampling_date, data_end_date);
        data_start_date -= 1/last_gap.sampling_rate;
        last_gap.startDate = data_start_date;
    } else {
        data_start_date = adjustStartDate(vecM, last_gap, data_start_date, data_end_date);
    }
}

static void checkWindowSize(const ARIMA_GAP &last_gap, const std::vector<MeteoData>& vecM, const Date& resampling_date, Date& data_start_date, Date& data_end_date, const double& window_size) {
	if (data_start_date < vecM[0].date && last_gap.startDate != resampling_date) {
        data_start_date = vecM[0].date;
    } 
    if (data_end_date > vecM[vecM.size()-1].date) {
        data_end_date = vecM[vecM.size()-1].date;
    }
    if (data_end_date.getJulian(true) - data_start_date.getJulian(true) > window_size) {
        throw IOException("The data window needed to interpolate the gap " + last_gap.toString() + " is larger than the resampling window size");
    }
}


void computeARIMAGap(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date, size_t& indexP1, size_t& indexP2, double& before_window, double& after_window, double& window_size, Date& data_start_date, Date& data_end_date)
{
	indexP1 = searchBackward(last_gap, pos, paramindex, vecM, resampling_date, window_size);
	indexP2 = searchForward(last_gap, pos, paramindex, vecM, resampling_date, window_size, indexP1);

	if (indexP1 == IOUtils::npos  && indexP2 == 0) {
		last_gap.startDate = resampling_date;
		last_gap.start = 0;
		last_gap.end = 0;
		indexP1 = 0;
	} else if (indexP1 == vecM.size() - 2 && indexP2 == vecM.size() - 1) {
		last_gap.startDate = vecM[vecM.size() - 1].date;
		last_gap.start = vecM.size() - 1;
		last_gap.endDate = last_gap.startDate + window_size;
		last_gap.end = vecM.size() - 1;
		indexP1 = vecM.size() - 1;
		indexP2 = vecM.size() - 1;
	} else if (indexP1 == IOUtils::npos || indexP2 == IOUtils::npos) {
        std::cout << "Could not pinpoint the gap "<< last_gap.toString() << std::endl;
        std::cout << "Gap start: " << last_gap.startDate.toString(Date::ISO) << "("<< indexP1 <<")" << std::endl;
        std::cout << "Gap end: " << last_gap.endDate.toString(Date::ISO) << "("<< indexP2 <<")" << std::endl;
		std::cout << "Consider increasing the Window Size" << std::endl;
    }
	data_start_date = last_gap.startDate - before_window;
	data_end_date = last_gap.endDate + after_window;	


	checkWindowSize(last_gap, vecM, resampling_date, data_start_date, data_end_date, window_size);
	last_gap.sampling_rate = computeSamplingRate(data_start_date, data_end_date, vecM);
	adjustDataStartDate(last_gap, vecM, resampling_date, data_start_date, data_end_date);
	if ((indexP1 == vecM.size()-1) && (indexP2 == vecM.size()-1)){
		last_gap.endDate = last_gap.startDate + MAX_ARIMA_EXTRAPOLATION / last_gap.sampling_rate;
		data_end_date = last_gap.endDate;
	}
}

// returns the most often accuring value in a vector
double mostLikelyValue(const std::vector<double>& vec) {
    if (vec.empty()) {
        throw mio::IOException("Vector of sampling rates is empty");
    }
    std::unordered_map<double, int> counts;
    for (double num : vec) {
        counts[num]++;
    }

    return std::max_element(counts.begin(), counts.end(), [](const std::pair<double,int>& pair1, const std::pair<double,int>& pair2) {
        return pair1.second < pair2.second;
    })->first;
}

// compute the most often occuring sampling rate rounded to 1e-6
double computeSamplingRate(Date data_start_date, Date data_end_date, std::vector<MeteoData> vecM) {
    std::vector<double> time_diffs;
    for (size_t i = 0; i < vecM.size()-1; i++) {
        if (vecM[i].date >= data_start_date && vecM[i].date <= data_end_date) {
            double value = 1/std::abs((vecM[i+1].date - vecM[i].date).getJulian(true));
            time_diffs.push_back(std::round(value*100000/100000));
        }
    }
	double val = mostLikelyValue(time_diffs);
	if (val == 0) {
		throw IOException("Could not compute sampling rate");
	}
    return val;
}


static Date findFirstDateWithSamplingRate(const std::vector<MeteoData>& vecM, const double sampling_rate, const Date& data_start_date, const Date& data_end_date) {
    Date closestDate = data_start_date;
    double minDiff = std::numeric_limits<double>::max();

    for (const auto& data : vecM) {
        if (data.date < data_start_date || data.date > data_end_date) {
			continue;
        }
        double diff = std::abs((data.date - data_start_date).getJulian(true));
        if (diff <= 1.5 / sampling_rate) {
			
			minDiff = std::min(minDiff, diff);
            closestDate = data.date;
        }
    }
    return closestDate;
}

Date adjustStartDate(const std::vector<MeteoData>& vecM, const ARIMA_GAP& last_gap, Date data_start_date, const Date& data_end_date) {
    const Date neededDate = findFirstDateWithSamplingRate(vecM, last_gap.sampling_rate, data_start_date, data_end_date);
    if (data_start_date == neededDate) {
        return data_start_date;
    }
    double diff = std::abs((neededDate - data_start_date).getJulian(true));
    if (diff < last_gap.sampling_rate) {
        data_start_date = neededDate;
    } else {
        // closest date that can be reached with the sampling rate
        int num_steps = std::max(1, static_cast<int>(std::floor(diff*last_gap.sampling_rate)));
        data_start_date = neededDate - num_steps/last_gap.sampling_rate;
    }
    return data_start_date;
}

} // namespace ARIMAutils
} // namespace mio
