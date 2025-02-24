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

#ifndef INTERPOLARIMA_H
#define INTERPOLARIMA_H

/* For peace of mind, disable strict ISO c++ compliance as it would emit warnings
 * because of ctsa's variable lenght arrays*/
#ifdef __GNUC__
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #include <meteoio/thirdParty/ctsa.h>
    #pragma GCC diagnostic pop
#elif defined __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wpedantic"
	#include <meteoio/thirdParty/ctsa.h>
	#pragma clang diagnostic pop
#else
	#include <meteoio/thirdParty/ctsa.h>
#endif

#include <iostream>
#include <map>
#include <meteoio/meteoResampling/ARIMAutils.h>
#include <string>
#include <vector>

namespace mio {

using namespace ARIMAutils;
/**
* @class InterpolARIMA
*
* @brief This class is used for interpolating or predicting missing data in a time series using the Auto ARIMA algorithm.
*
* @details Depending on the constructor that is used, the data and auto ARIMA models are set up for either interpolation or prediction.
* The interpolation is interpolate(). The prediction methods is predict().
*
* Interpolate will fill a gap in the data, whose start is specified by gap_loc
* and whose length is specified by N_gap. Data is assumed to be of equal time steps, and is split into two parts, data_before and
* data_after. So in the end data should be of size data_before + data_after + N_gap. The interpolation is done by fitting one ARIMA
* model to data_before and one to data_after. The ARIMA models are fitted using the auto.arima algorithm from the <a
* href="https://github.com/rafat/ctsa/tree/master">ctsa</a> (BSD-3 Clause, see below). The ARIMA models are then used to predict the
* missing data forward and backward in time. The final prediction is a weighted average of the two, where the weighting is done so more
* information comes from the closer data.
*
* Predict will predict the next n_steps values in the time series. It can either be forward in time (direction = "forward") or backward
* in time (direction = "backward"). For forward prediction data[0:gap_loc] is used to fit the ARIMA model, and for backward prediction
* data[gap_loc + N_gap:] is used.
*
* For more Information concerning ARIMA see, <a
* href="https://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average">Wikipedia</a>, and <a
* href="https://otexts.com/fpp2/arima.html">Forecast documentation</a>, and the final interpolation algorithm: <a
* href="https://www.tandfonline.com/doi/abs/10.1080/02664769624332?casa_token=fEVPFRYrr7sAAAAA:ozZFAcUWX4mKaUI8tvOn6R-3giOHefH0p8vaRDFCN1ORGy0d9evP7Hn9aLbMWsUQsIKrKEKxP-M">Time
* weighted average</a>
*
*
* @note Interpolate is meant to only be used, when there is actually backward data available. If there is no backward data, then
* predict should be used instead. Where predict is meant to be used in conjunction with the according constructor. Currently prediction
* forward or backward, when providing both is not implemented, but can be easily added on demand.
*
* @author Patrick Leibersperger
* @date 2024-01-25
*
*
* Copyright (c) 2014, Rafat Hussain
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions
* are met:
*
*  1.  Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
*
*  2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
*  	in the documentation and/or other materials provided with the distribution.
*
*  3.  The name of the author may not be used to endorse or promote products derived from this software without specific prior written
* permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
class InterpolARIMA {
	public:
		InterpolARIMA();
		InterpolARIMA(const std::vector<double>& data_in, const size_t& gap_loc, const size_t& N_gap, const int& s = 0);
		InterpolARIMA(const std::vector<double>& data_in, const size_t& gap_loc, const size_t& N_gap, const std::vector<double>& xreg_vec, const int& s = 0);
		InterpolARIMA(const std::vector<double>& data_in, const size_t& gap_loc, const size_t& n_predictions, const std::string& direction = "forward", const int& s = 0);

		// Setters
		void setAutoArimaMetaData(int max_p_param = 8, int max_d_param = 3, int max_q = 8, int start_p = 2, int start_q = 2, int max_P = 2,
									int max_D = 1, int max_Q = 2, int start_P = 1, int start_Q = 1, bool seasonal = true,
									bool stationary = false);
		void setOptMetaData(ObjectiveFunction method = CSS_MLE, OptimizationMethod opt_method = BFGS, bool stepwise = true,
							bool approximation = false, int num_models = 94);
		void setVerbose(bool verbose = false);
		void setNormalizationMode(Normalization::Mode mode);
		void setManualARIMA(int p, int d, int q, int P, int D, int Q, bool fill_backward);

		// Interpolation methods
		void fillGap();
		void fillGapManual();

		void interpolate();
		std::vector<double> predict(size_t n_steps = 0);
		// only do denormalization ARIMA predict! otherwise it will denorm twice
		std::vector<double> ARIMApredict(size_t n_steps);

		// Getters
		std::vector<double> getData() { return norm.denormalize(data); }
		std::vector<double> getForwardData() { return norm.denormalize(data_forward); }
		std::vector<double> getBackwardData() { return  norm.denormalize(data_backward); }
		std::vector<double> getInterpolatedData();

		// Copy constructor
		InterpolARIMA(const InterpolARIMA &other)
			: norm(other.norm), data(other.data), gap_loc(other.gap_loc), N_gap(other.N_gap), time(other.time), pred_forward(other.pred_forward),
			pred_backward(other.pred_backward), xreg_vec_f(other.xreg_vec_f), xreg_vec_b(other.xreg_vec_b),
			data_forward(other.data_forward), data_backward(other.data_backward), new_xreg_vec_f(other.new_xreg_vec_f),
			new_xreg_vec_b(other.new_xreg_vec_b), xreg_f((xreg_vec_f.empty()) ? nullptr : &xreg_vec_f[0]),
			xreg_b((xreg_vec_b.empty()) ? nullptr : &xreg_vec_b[0]), new_xreg_f((new_xreg_vec_f.empty()) ? nullptr : &new_xreg_vec_f[0]),
			new_xreg_b((new_xreg_vec_b.empty()) ? nullptr : &new_xreg_vec_b[0]), amse_forward(other.amse_forward),
			amse_backward(other.amse_backward), N_data_forward(other.N_data_forward), N_data_backward(other.N_data_backward),
			max_p(other.max_p), max_d(other.max_d), max_q(other.max_q), start_p(other.start_p), start_q(other.start_q),
			max_P(other.max_P), max_D(other.max_D), max_Q(other.max_Q), start_P(other.start_P), start_Q(other.start_Q), r(other.r),
			s(other.s), method(other.method), opt_method(other.opt_method), stepwise(other.stepwise), approximation(other.approximation),
			num_models(other.num_models), seasonal(other.seasonal), stationary(other.stationary),
			auto_arima_forward(auto_arima_copy(other.auto_arima_forward)), auto_arima_backward(auto_arima_copy(other.auto_arima_backward)), sarima_forward(other.sarima_forward) {
		}

    // Copy assignment operator
    InterpolARIMA &
    operator=(const InterpolARIMA &other) {
		// protect against invalid self-assignment
		if (this != &other) {
			auto_arima_forward = auto_arima_copy(other.auto_arima_forward);
			auto_arima_backward = auto_arima_copy(other.auto_arima_backward);

			// 3: copy all the other fields from the other object
			gap_loc = other.gap_loc;
			N_gap = other.N_gap;
			time = other.time;
			pred_forward = other.pred_forward;
			pred_backward = other.pred_backward;
			norm = other.norm;
			data = other.data;
			xreg_vec_f = other.xreg_vec_f;
			xreg_vec_b = other.xreg_vec_b;
			data_forward = other.data_forward;
			data_backward = other.data_backward;
			new_xreg_vec_f = other.new_xreg_vec_f;
			new_xreg_vec_b = other.new_xreg_vec_b;
			N_data_forward = other.N_data_forward;
			N_data_backward = other.N_data_backward;
			max_p = other.max_p;
			max_d = other.max_d;
			max_q = other.max_q;
			start_p = other.start_p;
			start_q = other.start_q;
			max_P = other.max_P;
			max_D = other.max_D;
			max_Q = other.max_Q;
			start_P = other.start_P;
			start_Q = other.start_Q;
			r = other.r;
			s = other.s;
			method = other.method;
			opt_method = other.opt_method;
			stepwise = other.stepwise;
			approximation = other.approximation;
			num_models = other.num_models;
			seasonal = other.seasonal;
			stationary = other.stationary;

			// 4: handle the pointers to the vectors
			xreg_f = (xreg_vec_f.empty()) ? nullptr : &xreg_vec_f[0];
			xreg_b = (xreg_vec_b.empty()) ? nullptr : &xreg_vec_b[0];
			new_xreg_f = (new_xreg_vec_f.empty()) ? nullptr : &new_xreg_vec_f[0];
			new_xreg_b = (new_xreg_vec_b.empty()) ? nullptr : &new_xreg_vec_b[0];
		}
		// by convention, always return *this
		return *this;
	}

	~InterpolARIMA() {
		auto_arima_free(auto_arima_forward);
		auto_arima_free(auto_arima_backward);
		delete xreg_f;
		delete xreg_b;
		delete new_xreg_f;
		delete new_xreg_b;
	}

	// info
	std::string toString();
	std::string autoArimaInfo(const auto_arima_object& obj);


private:
	// Interpolation variables
	Normalization norm;
	std::vector<double> data;
	size_t gap_loc;
	size_t N_gap;
	std::vector<double> time;
	std::vector<double> pred_forward, pred_backward;

	// Auto Arima variables
	// const doesnt work with c
	std::vector<double> xreg_vec_f, xreg_vec_b, data_forward, data_backward, new_xreg_vec_f, new_xreg_vec_b;
	double *xreg_f;
	double *xreg_b;
	double *new_xreg_f;
	double *new_xreg_b;
	std::vector<double> amse_forward, amse_backward;
	size_t N_data_forward, N_data_backward;
	int max_p = 8, max_d = 3, max_q = 8;
	int start_p = 2, start_q = 2;
	int max_P = 2, max_D = 1, max_Q = 2;
	int start_P = 1, start_Q = 1;
	int r = 0, s = 0;
	ObjectiveFunction method = CSS_MLE; OptimizationMethod opt_method = BFGS;
	bool stepwise = true, approximation = true;
	int num_models = 94;
	bool seasonal = true, stationary = false;

	bool consistencyCheck();
	auto_arima_object initAutoArima(const size_t& N_data);

	// last to be initialized
public:
	auto_arima_object auto_arima_forward;
	auto_arima_object auto_arima_backward;

private:
	// (S)ARIMA variables
	bool set_manual = false;
	bool fill_backward_manual = false;

public:
	sarima_object sarima_forward;
};

} // namespace mio

#endif // INTERPOLARIMA_H
