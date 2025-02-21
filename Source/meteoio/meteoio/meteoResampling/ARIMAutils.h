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
#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <meteoio/MeteoIO.h>
#include <vector>

static const double DATE_TOLERANCE = 1e-6;
static const int MIN_ARIMA_DATA_POINTS = 8;
static const int MAX_ARIMA_EXTRAPOLATION = 48; // is two days with hourly data

namespace mio {

    namespace ARIMAutils {

        enum ObjectiveFunction {
            CSS_MLE,
            MLE,
            CSS,
        };
        extern const std::map<ObjectiveFunction, std::string> ObjectiveFunctionMap;
        
        enum OptimizationMethod {
            Nelder_Mead,
            Newton_Line_Search,
            Newton_Trust_Region_Hook_Step,
            Newton_Trust_Region_Double_Dog_Leg,
            Conjugate_Gradient,
            BFGS,
            LBFGS,
            BFGS_MTM,
        };
        extern const std::map<OptimizationMethod, std::string> OptimizationMethodMap;

        // a struct to hold the coefficients for normalization and denormalization of a time series, cannot be used 
        // for multiple time series
        class Normalization {
            public:
                enum Mode {
                    MinMax,
                    ZScore,
                    Nothing,
                };
                Normalization(); 
                Normalization(std::vector<double>& data, Mode new_mode); 
                Normalization(std::vector<double>& data);
                
                void setMode(Mode new_mode) {mode = new_mode;}
                Mode getMode() {return mode;}
                std::vector<double> normalize(const std::vector<double>& data); 
                std::vector<double> denormalize(const std::vector<double>& data);

            private:
                double mean;
                double std;
                double min;
                double max;
                Mode mode = Nothing;
            };

        // slice a vector from start to start+N
        std::vector<double> slice(const std::vector<double> &vec, size_t start, size_t N);

        // slice a vector from start to end
        std::vector<double> slice(const std::vector<double> &vec, size_t start);

        // np.arange for c++
        std::vector<double> arange(size_t start, size_t N);

        template <typename T> T findMinMax(const std::vector<T> &vec, bool findMin) {
            assert(!vec.empty()); // Ensure the vector is not empty

            T extremeValue = vec[0];
            for (const auto &value : vec) {
                if (findMin ? value < extremeValue : value > extremeValue) {
                    extremeValue = value;
                }
            }
            return extremeValue;
        }

        // calculate the of a vector
        double calcVecMean(const std::vector<double> &vec);

        // calculate the standard deviation of a vector
        double stdDev(const std::vector<double> &vec);

        // reverse a vector in place
        template <typename T> void reverseVector(std::vector<T> &vec) {
            if (vec.empty()) {
                throw std::invalid_argument("Cannot reverse an empty vector");
            }
            int start = 0;
            int end = int(vec.size()) - 1;

            while (start < end) {
                std::swap(vec[start], vec[end]);
                start++;
                end--;
            }
        }

        // reverse a vector and return it
        template <typename T> std::vector<T> reverseVectorReturn(const std::vector<T> &vec) {
            if (vec.empty()) {
                throw std::invalid_argument("Cannot reverse an empty vector");
            }
            std::vector<T> reversed_vec = vec;
            int start = 0;
            int end = int(reversed_vec.size()) - 1;

            while (start < end) {
                std::swap(reversed_vec[start], reversed_vec[end]);
                start++;
                end--;
            }

            return reversed_vec;
        }

        // converts a vector of MeteoData to a vector of doubles
        std::vector<double> toVector(const std::vector<MeteoData> &vecM, const std::string &paramname);

        // converts a vector of MeteoData to a vector of doubles
        std::vector<double> toVector(const std::vector<MeteoData> &vecM, const size_t &paramindex);

        // helper to parse direction argument for interpolarima
        std::vector<double> decideDirection(const std::vector<double> &data, const std::string &direction, bool forward, size_t gap_loc,
                                            size_t length);

        // a struct to cache information about a gap
        struct ARIMA_GAP {
            ARIMA_GAP() : start(), end(), startDate(), endDate(), sampling_rate() {}
            void extend(const size_t &idx, const std::vector<MeteoData> &vecM) {
                if (idx < start)
                    setStart(idx, vecM);
                if (idx > end)
                    setEnd(idx, vecM);
            }
            void setStart(const size_t &idx, const std::vector<MeteoData> &vecM) {
                if (idx >= vecM.size())
                    return;
                start = idx;
                startDate = vecM[idx].date;
            }
            void setEnd(const size_t &idx, const std::vector<MeteoData> &vecM) {
                if (idx >= vecM.size())
                    return;
                end = idx;
                endDate = vecM[idx].date;
            }
            void reset() {
                start = IOUtils::npos;
                end = IOUtils::npos;
                startDate = Date();
                endDate = Date();
            }
            size_t start, end;
            Date startDate, endDate;
            double sampling_rate;
            bool isGap() {
                return (endDate - startDate).getJulian(true) * sampling_rate >= 2;
            } // TODO: should i always do arima prediction?
            std::string toString() const {
                std::ostringstream os;
                os << "ARIMA_GAP: {\n"
                   << "\tStart Date: " << startDate.toString(Date::ISO) << ",\n"
                   << "\tEnd Date: " << endDate.toString(Date::ISO) << ",\n"
                   << "\tSampling Rate: " << sampling_rate << ",\n"
                   << "}";
                return os.str();
            }
        };

        // return true if a valid point could be found backward from pos
        size_t searchBackward(ARIMA_GAP &last_gap, const size_t &pos, const size_t &paramindex, const std::vector<MeteoData> &vecM,
                              const Date &resampling_date, const double &i_window_size);

        // return true if a valid point could be found forward from pos
        size_t searchForward(ARIMA_GAP &last_gap, const size_t &pos, const size_t &paramindex, const std::vector<MeteoData> &vecM,
                             const Date &resampling_date, const double &i_window_size, const size_t &indexP1);

        void computeARIMAGap(ARIMA_GAP &last_gap, const size_t &pos, const size_t &paramindex, const std::vector<MeteoData> &vecM,
                             const Date &resampling_date, size_t &indexP1, size_t &indexP2, double &before_window, double &after_window,
                             double &window_size, Date &data_start_date, Date &data_end_date);

        // roughly equal between two dates, given a tolerance level
        bool requal(const Date &date1, const Date &date2);

        // returns the most often accuring value in a vector
        double mostLikelyValue(const std::vector<double> &vec);

        // compute the most often occuring sampling rate rounded to 1e-6
        double computeSamplingRate(Date data_start_date, Date data_end_date, std::vector<MeteoData> vecM);

        Date findFirstDateWithSamplingRate(const std::vector<MeteoData> &vecM, const double sampling_rate, const Date &data_start_date,
                                           Date &data_end_date);
        Date adjustStartDate(const std::vector<MeteoData> &vecM, const ARIMA_GAP &last_gap, Date data_start_date,
                             const Date &data_end_date);

        template <typename T> std::string convertVectorsToString(const std::vector<std::vector<T>> &vecs) {
            std::ostringstream oss;
            size_t maxSize = 0;
            for (const auto &vec : vecs) {
                maxSize = std::max(maxSize, vec.size());
            }

            // Print headers
            for (size_t i = 0; i < vecs.size(); i++) {
                oss << std::left << std::setw(10) << "Vector" + std::to_string(i + 1);
            }
            oss << std::endl;
            oss << std::string(vecs.size() * 10, '-') << std::endl;

            for (size_t i = 0; i < maxSize; i++) {
                for (const auto &vec : vecs) {
                    // Print elements from vec or "NaN" if out of range
                    if (i < vec.size()) {
                        oss << std::left << std::setw(10) << vec[i];
                    } else {
                        oss << std::left << std::setw(10) << "NaN";
                    }
                }
                oss << std::endl;
            }
            return oss.str();
        }

        template <typename T> void printVectors(const std::vector<std::vector<T>> &vecs) {
            size_t maxSize = 0;
            for (const auto &vec : vecs) {
                maxSize = std::max(maxSize, vec.size());
            }

            // Print headers
            for (size_t i = 0; i < vecs.size(); i++) {
                std::cout << std::left << std::setw(10) << "Vector" + std::to_string(i + 1);
            }
            std::cout << std::endl;
            std::cout << std::string(vecs.size() * 10, '-') << std::endl;

            for (size_t i = 0; i < maxSize; i++) {
                for (const auto &vec : vecs) {
                    // Print elements from vec or "NaN" if out of range
                    if (i < vec.size()) {
                        std::cout << std::left << std::setw(10) << vec[i];
                    } else {
                        std::cout << std::left << std::setw(10) << "NaN";
                    }
                }
                std::cout << std::endl;
            }
        }

        template <typename T> void printVectors(const std::vector<Date> &vec1, const std::vector<T> &vec2) {
            size_t maxSize = std::max(vec1.size(), vec2.size());

            // Print headers
            std::cout << std::left << std::setw(30) << "Date1"
                      << "| Date2" << std::endl;
            std::cout << "--------------------------------------------------" << std::endl;

            for (size_t i = 0; i < maxSize; i++) {
                // Print date from vec1 or "NaN" if out of range
                if (i < vec1.size()) {
                    std::cout << std::left << std::setw(30) << vec1[i].toString(Date::ISO) << "| ";
                } else {
                    std::cout << std::left << std::setw(30) << "NaN"
                              << "| ";
                }

                // Print date from vec2 or "NaN" if out of range
                if (i < vec2.size()) {
                    std::cout << vec2[i] << std::endl;
                } else {
                    std::cout << "NaN" << std::endl;
                }
            }
        }

    } // namespace ARIMAutils
} // namespace mio
#endif // UTILS_H
