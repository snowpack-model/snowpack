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
#include <cmath>
#include <cstdlib> // for std::rand and std::srand
#include <cstring>
#include <iomanip>
#include <iostream>
#include <meteoio/MeteoIO.h>
#include <meteoio/meteoResampling/InterpolARIMA.h>
#include <sstream>
#include <unistd.h>

namespace mio {

// ------------------- Constructor ------------------- //

// Default constructor
InterpolARIMA::InterpolARIMA()
    : norm(data), data(5, 0.0), gap_loc(0), N_gap(5), time(), pred_forward(), pred_backward(), xreg_vec_f(), xreg_vec_b(), data_forward(), data_backward(), new_xreg_vec_f(), new_xreg_vec_b(),
        xreg_f(nullptr), xreg_b(nullptr), new_xreg_f(nullptr), new_xreg_b(nullptr), amse_forward(), amse_backward(), N_data_forward(5), N_data_backward(5),
        auto_arima_forward(initAutoArima(N_data_forward)), auto_arima_backward(initAutoArima(N_data_backward)), sarima_forward() 
{}

/**
    * @brief Main Constructor for an InterpolARIMA object. Used to fill 1 gap in the data.
    *
    * @param data_in A vector of double values representing the input data.
    * @param gap_location The location of the gap in the data.
    * @param gap_length The length of the gap in the data.
    * @param period (Optional) The period of the ARIMA model. Defaults to 0. Only needed when the period is known.
    */
InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, size_t gap_length, int period)
    : norm(data_in), data(norm.normalize(data_in)), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap), xreg_vec_f(0), xreg_vec_b(0),
        data_forward(slice(data, 0, gap_loc)), data_backward(slice(data, gap_loc + N_gap)), new_xreg_vec_f(0), new_xreg_vec_b(0), xreg_f(nullptr), xreg_b(nullptr), new_xreg_f(nullptr),
        new_xreg_b(nullptr), amse_forward(N_gap), amse_backward(N_gap), N_data_forward(data_forward.size()), N_data_backward(data_backward.size()), s(period),
        auto_arima_forward(initAutoArima(N_data_forward)), auto_arima_backward(initAutoArima(N_data_backward)), sarima_forward() 
{
    // reverse the backward data
    reverseVector(data_backward); // TODO: can be done with std::reverse in C++17
}

/**
    * @brief This constructor is used to initialize an InterpolARIMA object for filling a gap in the data with exogenous variables.
    *
    * @param data_in The input data for making predictions.
    * @param gap_location The starting location of the data gap.
    * @param gap_length The length of the data gap.
    * @param xreg_vec_in The exogenous inputs for the ARIMA model.
    * @param period The period for the ARIMA model.
    *
    *
    * @note This constructor is part of the [`InterpolARIMA`](meteoio/meteoResampling/InterpolARIMA.cc) class.
    */
InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, size_t gap_length, std::vector<double> xreg_vec_in, int period)
    : norm(data_in), data(data_in), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap), xreg_vec_f(slice(xreg_vec_in, 0, gap_loc)),
        xreg_vec_b(reverseVectorReturn(slice(xreg_vec_in, gap_loc + N_gap))), data_forward(slice(data, 0, gap_loc)), data_backward(slice(data, gap_loc + N_gap)),
        new_xreg_vec_f(xreg_vec_f.size() == 0 ? 0 : N_gap), new_xreg_vec_b(xreg_vec_b.size() == 0 ? 0 : N_gap), xreg_f(xreg_vec_f.size() == 0 ? nullptr : &xreg_vec_f[0]),
        xreg_b(xreg_vec_b.size() == 0 ? nullptr : &xreg_vec_b[0]), new_xreg_f(xreg_vec_f.size() == 0 ? nullptr : &new_xreg_vec_f[0]),
        new_xreg_b(xreg_vec_b.size() == 0 ? nullptr : &new_xreg_vec_b[0]), amse_forward(N_gap), amse_backward(N_gap), N_data_forward(data_forward.size()), N_data_backward(data_backward.size()),
        r(xreg_vec_in.size() == 0 ? 0 : static_cast<int>(xreg_vec_f.size() / N_data_forward)), s(period), auto_arima_forward(initAutoArima(N_data_forward)),
        auto_arima_backward(initAutoArima(N_data_backward)), sarima_forward() 
{
    // reverse the backward data
    reverseVector(data_backward); // TODO: Can be done with std::reverse in C++17
}

/**
    * @brief This constructor is used to initialize an InterpolARIMA object for making predictions ahead or backward in time.
    *
    * @param data_in The input data for making predictions.
    * @param data_end The end location of the data gap.
    * @param n_predictions The number of predictions to be made.
    * @param direction The direction of the prediction. Can be either "forward" or "past".
    * @param period The period for the ARIMA model.
    *
    * It also decides the direction of the data based on the `direction` parameter and the `gap_loc`:
    * - If the direction is "forward", the prediction is made based on the data from
    *   the beginning of the data set up to `gap_loc`.
    * - If the direction is "past", the prediction is made based on the data from `gap_loc`
    *   to the end of the data set.
    *
    * @note This constructor is part of the [`InterpolARIMA`](meteoio/meteoResampling/InterpolARIMA.cc) class in [InterpolARIMA.cc](meteoio/meteoResampling/InterpolARIMA.cc).
    */
InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t data_end, size_t n_predictions, std::string direction, int period)
    : norm(data_in), data(norm.normalize(data_in)), gap_loc(data_end), N_gap(n_predictions), time(arange(0, data.size())), pred_forward(n_predictions), pred_backward(n_predictions), xreg_vec_f(0),
        xreg_vec_b(0), data_forward(decideDirection(data, direction, true, gap_loc, n_predictions)), data_backward(decideDirection(data, direction, false, gap_loc, n_predictions)),
        new_xreg_vec_f(0), new_xreg_vec_b(0), xreg_f(nullptr), xreg_b(nullptr), new_xreg_f(nullptr), new_xreg_b(nullptr), amse_forward(N_gap), amse_backward(N_gap),
        N_data_forward(data_forward.size()), N_data_backward(data_backward.size()), s(period), auto_arima_forward(initAutoArima(N_data_forward)),
        auto_arima_backward(initAutoArima(N_data_backward)), sarima_forward() 
{}

// ------------------- Helper methods ------------------- //
auto_arima_object InterpolARIMA::initAutoArima(size_t N_data) {
    std::vector<int> pqdmax = {max_p, max_d, max_q};
    std::vector<int> PQDmax = {max_P, max_D, max_Q};
    return auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, static_cast<int>(N_data));
}

std::string InterpolARIMA::toString() {
    std::stringstream ss;

    ss << "\n---------------- Auto ARIMA Model Information ----------------\n";

#ifdef DEBUG
    // Base Arima variables
    ss << "\nBase Arima Variables:\n";
    ss << std::left << std::setw(10) << "Max p:" << max_p << "\n";
    ss << std::left << std::setw(10) << "Max d:" << max_d << "\n";
    ss << std::left << std::setw(10) << "Max q:" << max_q << "\n";
    ss << std::left << std::setw(10) << "Start p:" << start_p << "\n";
    ss << std::left << std::setw(10) << "Start q:" << start_q << "\n";

    ss << "\n-------------------------------------------------------------\n";

    // Seasonal Arima variables
    ss << "\nSeasonal Arima Variables:\n";
    ss << std::left << std::setw(10) << "Max P:" << max_P << "\n";
    ss << std::left << std::setw(10) << "Max D:" << max_D << "\n";
    ss << std::left << std::setw(10) << "Max Q:" << max_Q << "\n";
    ss << std::left << std::setw(10) << "Start P:" << start_P << "\n";
    ss << std::left << std::setw(10) << "Start Q:" << start_Q << "\n";
    ss << std::left << std::setw(10) << "r:" << r << "\n";
    ss << std::left << std::setw(10) << "s:" << s << "\n";

    ss << "\n-------------------------------------------------------------\n";
#endif

    // Data information
    ss << "\nData Information:\n";
    ss << std::left << std::setw(20) << "Data size:" << data.size() << "\n";
    ss << std::left << std::setw(20) << "Gap location:" << gap_loc << "\n";
    ss << std::left << std::setw(20) << "N_gap:" << N_gap << "\n";
    ss << std::left << std::setw(20) << "Data forward size:" << data_forward.size() << "\n";
    ss << std::left << std::setw(20) << "Data backward size:" << data_backward.size() << "\n";

    ss << "\n-------------------------------------------------------------\n";
    ss << "Forward Model:\n";
    ss << autoArimaInfo(auto_arima_forward);
    ss << "\n-------------------------------------------------------------\n";
    ss << "Backward Model:\n";
    ss << autoArimaInfo(auto_arima_backward);

    ss << "\n-------------------------------------------------------------\n";

    return ss.str();
}

std::string InterpolARIMA::autoArimaInfo(auto_arima_object obj) {
    int i, pq, t, ncxreg, mean;
    pq = obj->p + obj->q + obj->P + obj->Q + obj->M;
    mean = obj->M - obj->r;

    std::stringstream result;

    if (obj->method == 0 || obj->method == 1) {
        result << "\n\n Exit Status \n";
        result << "Return Code : " << obj->retval << "\n";
        result << "Exit Message : ";
        if (obj->retval == 0) {
            result << "Input Error";
        } else if (obj->retval == 1) {
            result << "Probable Success";
        } else if (obj->retval == 4) {
            result << "Optimization Routine didn't converge";
        } else if (obj->retval == 7) {
            result << "Exogenous Variables are collinear";
        } else if (obj->retval == 10) {
            result << "Nonstationary AR part";
        } else if (obj->retval == 12) {
            result << "Nonstationary Seasonal AR part";
        } else if (obj->retval == 15) {
            result << "Optimization Routine Encountered Inf/Nan Values";
        }
        result << "\n\n"; // Add a newline after the exit message
    }

    result << "  ARIMA Seasonal Order : ( " << obj->p << ", " << obj->d << ", " << obj->q << ") * (" << obj->P << ", " << obj->D << ", " << obj->Q << ")\n\n";

    result << std::setw(20) << "Coefficients" << std::setw(20) << "Value" << std::setw(20) << "Standard Error"
            << "\n\n";
    for (i = 0; i < obj->p; ++i) {
        result << "AR" << std::setw(15) << i + 1 << std::setw(20) << obj->phi[i] << std::setw(20) << sqrt(obj->vcov[i + pq * i]) << "\n";
    }
    for (i = 0; i < obj->q; ++i) {
        t = obj->p + i;
        result << "MA" << std::setw(15) << i + 1 << std::setw(20) << obj->theta[i] << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
    }
    for (i = 0; i < obj->P; ++i) {
        t = obj->p + obj->q + i;
        result << "SAR" << std::setw(14) << i + 1 << std::setw(20) << obj->PHI[i] << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
    }
    for (i = 0; i < obj->Q; ++i) {
        t = obj->p + obj->q + obj->P + i;
        result << "SMA" << std::setw(14) << i + 1 << std::setw(20) << obj->THETA[i] << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
    }
    result << "\n";
    t = obj->p + obj->q + obj->P + obj->Q;
    if (mean > 0) {
        result << std::setw(17) << "MEAN" << std::setw(20) << obj->mean << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
        t++;
    } else {
        result << std::setw(17) << "MEAN" << std::setw(20) << obj->mean << "\n";
    }
    ncxreg = 0;
    if (obj->idrift == 1) {
        result << std::setw(17) << "TREND" << std::setw(20) << obj->exog[0] << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
        t++;
        ncxreg++;
    } else {
        result << std::setw(17) << "TREND" << std::setw(20) << 0.0 << "\n";
    }
    for (i = ncxreg; i < obj->r; ++i) {
        result << std::setw(17) << "EXOG" << std::setw(20) << obj->exog[i] << std::setw(20) << sqrt(obj->vcov[t + pq * t]) << "\n";
        t++;
    }

    result << "\n" << std::setw(17) << "SIGMA^2" << std::setw(20) << obj->sigma2 << "\n\n";

    result << "ESTIMATION METHOD : ";
    if (obj->method == 0) {
        result << "CSS-MLE";
    } else if (obj->method == 1) {
        result << "MLE";
    } else if (obj->method == 2) {
        result << "CSS";
    }
    result << "\n\n";

    result << "OPTIMIZATION METHOD : ";
    if (obj->optmethod == 0) {
        result << "Nelder-Mead";
    } else if (obj->optmethod == 1) {
        result << "Newton Line Search";
    } else if (obj->optmethod == 2) {
        result << "Newton Trust Region - Hook Step";
    } else if (obj->optmethod == 3) {
        result << "Newton Trust Region - Double Dog-Leg";
    } else if (obj->optmethod == 4) {
        result << "Conjugate Gradient";
    } else if (obj->optmethod == 5) {
        result << "BFGS";
    } else if (obj->optmethod == 6) {
        result << "L-BFGS";
    } else if (obj->optmethod == 7) {
        result << "BFGS More-Thuente Line Search";
    }

    result << "\n\n";

    result << "AIC criterion : " << obj->aic << "\n\n";

    result << "BIC criterion : " << obj->bic << "\n\n";

    result << "AICC criterion : " << obj->aicc << "\n\n";

    if (obj->method == 0 || obj->method == 1 || obj->method == 2) {
        result << "Log Likelihood : " << obj->loglik << "\n\n";
    }

    result << "Auto ARIMA Parameters \n\n";

    result << "Approximation: " << (obj->approximation == 1 ? "TRUE" : "FALSE") << "\n";

    result << "Stepwise: " << (obj->stepwise == 1 ? "TRUE" : "FALSE");

    return result.str(); // Return the resulting string
}

// ------------------- Setters ------------------- //}
// Change the Normalization of the data
void InterpolARIMA::setNormalizationMode(Normalization::Mode mode) {
    if (mode == norm.getMode())
        return;
    data = norm.denormalize(data);
    data_backward = norm.denormalize(data_backward);
    data_forward = norm.denormalize(data_forward);
    norm.setMode(mode);
    data = norm.normalize(data);
    data_backward = norm.normalize(data_backward);
    data_forward = norm.normalize(data_forward);
    return;
}

// Set the metadata for the auto arima objects
void InterpolARIMA::setAutoArimaMetaData(int max_p_param, int max_d_param, int max_q_param, int start_p_param, int start_q_param, int max_P_param, int max_D_param, int max_Q_param,
                                            int start_P_param, int start_Q_param, bool seasonal_param, bool stationary_param) {
    this->max_p = max_p_param;
    this->max_d = max_d_param;
    this->max_q = max_q_param;
    this->start_p = start_p_param;
    this->start_q = start_q_param;
    this->max_P = max_P_param;
    this->max_D = max_D_param;
    this->max_Q = max_Q_param;
    this->start_P = start_P_param;
    this->start_Q = start_Q_param;
    this->seasonal = seasonal_param;
    this->stationary = stationary_param;
    auto_arima_backward->pmax = max_p;
    auto_arima_forward->pmax = max_p;
    auto_arima_backward->dmax = max_d;
    auto_arima_forward->dmax = max_d;
    auto_arima_backward->qmax = max_q;
    auto_arima_forward->qmax = max_q;
    auto_arima_backward->Pmax = max_P;
    auto_arima_forward->Pmax = max_P;
    auto_arima_backward->Dmax = max_D;
    auto_arima_forward->Dmax = max_D;
    auto_arima_backward->Qmax = max_Q;
    auto_arima_forward->Qmax = max_Q;
    auto_arima_backward->p_start = start_p;
    auto_arima_forward->p_start = start_p;
    auto_arima_backward->q_start = start_q;
    auto_arima_forward->q_start = start_q;
    auto_arima_backward->P_start = start_P;
    auto_arima_forward->P_start = start_P;
    auto_arima_backward->Q_start = start_Q;
    auto_arima_forward->Q_start = start_Q;
    auto_arima_backward->seasonal = seasonal;
    auto_arima_forward->seasonal = seasonal;
    auto_arima_backward->stationary = stationary;
    auto_arima_forward->stationary = stationary;
}

// Set the metadata for the auto arima objects optimization
// options for method: "css-mle", "ml", "css"
// options for opt_method: "Nelder-Mead", "Newton Line Search", "Newton Trust Region - Hook Step", "Newton Trust Region - Double
// Dog-Leg", "Conjugate Gradient", "BFGS", "Limited Memory BFGS", "BFGS Using More Thuente Method"
void InterpolARIMA::setOptMetaData(ObjectiveFunction method_param, OptimizationMethod opt_method_param, bool stepwise_param, bool approximation_param, int num_models_param) {
    this->method = method_param;
    this->opt_method = opt_method_param;
    this->stepwise = stepwise_param;
    this->approximation = approximation_param;
    this->num_models = num_models_param;
    auto_arima_backward->method = static_cast<int>(this->method);
    auto_arima_forward->method = static_cast<int>(this->method);
    auto_arima_backward->optmethod = static_cast<int>(this->opt_method);
    auto_arima_forward->optmethod = static_cast<int>(this->opt_method);
    auto_arima_backward->stepwise = stepwise;
    auto_arima_forward->stepwise = stepwise;
}

void InterpolARIMA::setVerbose(bool verbose) {
    auto_arima_backward->verbose = verbose;
    auto_arima_forward->verbose = verbose;
}

void InterpolARIMA::setManualARIMA(int p, int d, int q, int P, int D, int Q, bool fill_backward) {
    set_manual = true;
    fill_backward_manual = fill_backward;
    if (N_data_forward < 5) {
        throw NoDataException("Not enough data to set the ARIMA model manually");
    }
    sarima_forward = sarima_init(p, d, q, s, P, D, Q, static_cast<int>(N_data_forward));
}

// ------------------- Getters ------------------- //
// Get the interpolated data
std::vector<double> InterpolARIMA::getInterpolatedData() {
    std::vector<double> interpolated_data(data.begin() + gap_loc, data.begin() + gap_loc + N_gap);
    return norm.denormalize(interpolated_data);
}

// ------------------- Interpolation methods ------------------- //
// Simulate n_steps into the future
std::vector<double> InterpolARIMA::simulate(int n_steps, int seed) {
    std::vector<double> sim(n_steps);
    seed++;
    // use the equations to simulate with random errors
    std::cerr << "not implemented, and not needed for now\n";
    return sim;
}

std::vector<double> InterpolARIMA::ARIMApredict(size_t n_steps) {
    if (n_steps == 0) {
        n_steps = N_gap;
    } else {
        pred_forward.resize(n_steps);
        amse_forward.resize(n_steps);
    }

    if (set_manual) {
        sarima_exec(sarima_forward, data_forward.data());
        sarima_predict(sarima_forward, data_forward.data(), static_cast<int>(n_steps), pred_forward.data(), amse_forward.data());
    } else {
        auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
        // check if the models are valid (p and q should not be zero at the same time)
        if ((auto_arima_forward->p == 0 && auto_arima_forward->q == 0) && (auto_arima_forward->P == 0 && auto_arima_forward->Q == 0)) {
            bool current_stepwise = auto_arima_forward->stepwise;
            auto_arima_setStepwise(auto_arima_forward, !current_stepwise);
            auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
        }
        auto_arima_predict(auto_arima_forward, data_forward.data(), xreg_f, static_cast<int>(n_steps), new_xreg_f, pred_forward.data(), amse_forward.data());
    }
    return norm.denormalize(pred_forward);
}

// Check if the model is a random walk
static bool isRandomWalk(auto_arima_object model) 
{ 
    return (model->p == 0 && model->q == 0) && (model->P == 0 && model->Q == 0); 
}


// Fill the gap using the auto arima objects
void InterpolARIMA::fillGap() {
    bool isRandom_f = false;
    bool isRandom_b = false;
    for (int meth_Id = 0; meth_Id < 3; meth_Id++) {
        // fit the models
        auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
        auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);

        isRandom_b = isRandomWalk(auto_arima_backward);
        isRandom_f = isRandomWalk(auto_arima_forward);

        // refit but using full search
        if (isRandom_b && isRandom_f && meth_Id == 0) {
            isRandom_b = false;
            isRandom_f = false;
            const bool current_stepwise = auto_arima_forward->stepwise;
            auto_arima_setStepwise(auto_arima_forward, !current_stepwise);
            auto_arima_setStepwise(auto_arima_backward, !current_stepwise);
            auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
            auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);
            if ((auto_arima_forward->p == 0 && auto_arima_forward->q == 0) && (auto_arima_forward->P == 0 && auto_arima_forward->Q == 0)) {
                isRandom_f = true;
            }
            if ((auto_arima_backward->p == 0 && auto_arima_backward->q == 0) && (auto_arima_backward->P == 0 && auto_arima_backward->Q == 0)) {
                isRandom_b = true;
            }
        }
        if (isRandom_b && isRandom_f)
            break;
        // predict the gap
        // forward
        auto_arima_predict(auto_arima_forward, data_forward.data(), xreg_f, static_cast<int>(N_gap), new_xreg_f, pred_forward.data(), amse_forward.data());

        // backward
        auto_arima_predict(auto_arima_backward, data_backward.data(), xreg_b, static_cast<int>(N_gap), new_xreg_b, pred_backward.data(), amse_backward.data());
        // interpolate with the weighting according to
        assert(pred_forward.size() == pred_backward.size());
        assert(pred_forward.size() == N_gap);

        reverseVector(pred_backward); // TODO: can be done with std::reverse in C++17

        bool consistency = consistencyCheck();
        if (consistency)
            break;
    }

    // W1 = sqrt(T-t/T)
    // W2 = sqrt(t/T)
    for (size_t id = 0; id < N_gap; id++) {
        double weight_f, weight_b;

        if (isRandom_f) {
            weight_f = 0;
            weight_b = 1;
        } else if (isRandom_b) {
            weight_f = 1;
            weight_b = 0;
        } else {
            double N = static_cast<double>(N_gap);
            weight_f = std::sqrt((N - time[id]) / N);
            weight_b = std::sqrt(time[id] / N);
        }
        data[gap_loc + id] = weight_f * pred_forward[id] + weight_b * pred_backward[id];
    }
}

    // Fill the gap using the manually set arima model
void InterpolARIMA::fillGapManual() {
    if (!set_manual) {
        throw TimeOutException("No manual ARIMA model set");
    }

    sarima_exec(sarima_forward, data_forward.data());
    sarima_predict(sarima_forward, data_forward.data(), static_cast<int>(N_gap), pred_forward.data(), amse_forward.data());

    assert(pred_forward.size() == N_gap);

    if (fill_backward_manual) {
        bool isRandom_b = false;
        for (int meth_Id = 0; meth_Id < 3; meth_Id++) {
            // fit the models
            auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);

            isRandom_b = isRandomWalk(auto_arima_backward);

            // refit but using full search
            if (isRandom_b  && meth_Id == 0) {
                isRandom_b = false;
                const bool current_stepwise = auto_arima_forward->stepwise;
                auto_arima_setStepwise(auto_arima_backward, !current_stepwise);
                auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);
                if ((auto_arima_backward->p == 0 && auto_arima_backward->q == 0) && (auto_arima_backward->P == 0 && auto_arima_backward->Q == 0)) {
                    isRandom_b = true;
                }
            }
            if (isRandom_b)
                break;

            // backward
            auto_arima_predict(auto_arima_backward, data_backward.data(), xreg_b, static_cast<int>(N_gap), new_xreg_b, pred_backward.data(), amse_backward.data());
            // interpolate with the weighting according to
            assert(pred_forward.size() == pred_backward.size());

            reverseVector(pred_backward); // TODO: can be done with std::reverse in C++17

            const bool consistency = consistencyCheck();
            if (consistency)
                break;
        }
        
        // W1 = sqrt(T-t/T)
        // W2 = sqrt(t/T)
        for (size_t id = 0; id < N_gap; id++) {
            double weight_f, weight_b;

            if (isRandom_b) {
                weight_f = 1;
                weight_b = 0;
            } else {
                const double N = static_cast<double>(N_gap);
                weight_f = std::sqrt((N - time[id]) / N);
                weight_b = std::sqrt(time[id] / N);
            }
            data[gap_loc + id] = weight_f * pred_forward[id] + weight_b * pred_backward[id];
        }
    } else {
        for (size_t id = 0; id < N_gap; id++) {
            data[gap_loc + id] = pred_forward[id];
        }
    }
}

bool InterpolARIMA::consistencyCheck() {
    const double mean_before = calcVecMean(data_forward);
    const double std_before = stdDev(data_forward);
    const double mean_after = calcVecMean(data_backward);
    const double std_after = stdDev(data_backward);
    const double max_interpolated = findMinMax(slice(data, gap_loc, N_gap), false);

    // needs more checks
    if (max_interpolated > mean_before + 4 * std_before && max_interpolated > mean_after + 4 * std_after)
        return false;
    else if (max_interpolated < mean_before - 4 * std_before && max_interpolated < mean_after - 4 * std_after)
        return false;
    return true;
}

// ------------------- Wrappers ------------------- //
// wrapper for filling the gap
void InterpolARIMA::interpolate() {
    bool fit = true;
    if (N_data_backward == 0 || N_data_forward == 0) {
        throw NoDataException("No data to interpolate: forward datapoints " + std::to_string(N_data_forward) + ", backward datapoints " + std::to_string(N_data_backward) + "\n");
    }
    while (fit) {
        if (set_manual)
            fillGapManual();
        else
            fillGap();
            
        const int retval_f = set_manual ? sarima_forward->retval : auto_arima_forward->retval;
        const int retval_b = (set_manual && !fill_backward_manual) ? 1 : auto_arima_backward->retval;
        if (retval_f == 0 || retval_b == 0) {
            const std::string where = (retval_f == 0) ? "forward data" : "backward data";
            throw AccessException("Interpolation Input data is erroneous in " + where);
        } else if (retval_f == 15 || retval_b == 15) {
            const std::string where = (retval_f == 15) ? "forward data" : "backward data";
            throw InvalidFormatException("Interpolation Input data has Inf/Nan values in " + where);
        } else if (retval_f == 4 || retval_b == 4) {
            const std::string where = (retval_f == 4) ? "forward data" : "backward data";
            if (method != CSS_MLE && opt_method != BFGS) {
                throw IOException("Optimization of ARIMA did not converge in " + where + ".\n Please try another method and optimization method");
            } else {
                OptimizationMethod new_opt_method = Nelder_Mead;
                setOptMetaData(method, new_opt_method);
            }
        } else {
            fit = false;
        }
    }
    return;
}

// wrapper for predicting the gap
std::vector<double> InterpolARIMA::predict(size_t n_steps) {
    bool fit = true;
    std::vector<double> pred;
    if (N_data_backward == 0 || N_data_forward == 0) {
        throw NoDataException("No data to interpolate: forward datapoints " + std::to_string(N_data_forward) + ", backward datapoints " + std::to_string(N_data_backward) + "\n");
    }
    while (fit) {
        pred = ARIMApredict(n_steps);
        const int retval_f = set_manual ? sarima_forward->retval : auto_arima_forward->retval;

        if (retval_f == 0) {
            throw AccessException("Interpolation Input data is erroneous when trying to predict");
        } else if (retval_f == 15) {
            throw InvalidFormatException("Interpolation Input data has Inf/Nan values for prediction");
        } else if (retval_f == 4) {
            if (method != CSS_MLE && opt_method != BFGS) {
                throw IOException("Optimization of ARIMA did not converge for prediction.\n Please try another method and optimization method");
            } else {
                OptimizationMethod new_opt_method = Nelder_Mead;
                setOptMetaData(method, new_opt_method);
            }
        } else {
            fit = false;
        }
    }
    return pred;
}
} // end namespace mio
