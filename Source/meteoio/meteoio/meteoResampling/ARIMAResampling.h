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
#ifndef ARIMARESAMPLING_H
#define ARIMARESAMPLING_H

#include <meteoio/meteoResampling/ARIMAutils.h>
#include <meteoio/meteoResampling/InterpolARIMA.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <vector>

namespace mio {

    using namespace ARIMAutils;
    /**
     * @brief This class is designed to handle interpolation (resampling) of data using the ARIMA (AutoRegressive Integrated Moving Average)
     * model
     *
     * @details It uses the InterpolARIMA class to perform interpolation using ARIMA, with model selection and fitting done with the <a
     * href="https://github.com/rafat/ctsa/tree/master">ctsa</a> (BSD-3 Clause, see below) library. That implements the auto ARIMA algorithm
     * from <a href="https://www.jstatsoft.org/article/view/v027i03">Hyndman and Khandakar (2008)</a>.
     *
     * Gaps in the data are detected and interpolated using the ARIMA model. If available, data before and after the gap is used to fit an
     * ARIMA model each,
     * one model predicting the gap forward in time, the other backward in time. A weighted average of the two predictions is then used as
     * interpolated value.
     * If only data is available before or after the gap, only one model is fitted and used for a prediction, either forward or backward in
     * time.
     *
     * The ARIMA model needs constant sampling rates, therefore the most likely rate is calculated in the data used to fit, and the data is
     * resampled to that sampling rate. If a requested point falls in between available data, it will be linearly interpolated.
     *
     * Missing values in the data used to fit the ARIMA Model, will be linearly interpolated as well.
     *
     * A gap is defined as a period of missing data, that has at least 2 data points with the most likely sampling rate. Only 1 missing data
     * point is linearly interpolated (should maybe just return instead).
     *
     * As ARIMA is a stochastic process (white noise errors are an inherent part of the model), "prediction" needs to be specified. Here we
     * call it a simulation,
     * when we draw a random error for each timestep, and use the ARIMA equations (see below) to predict the proceeding values. A
     * "prediction" is
     * essentially the mean of many simulations, including a standard deviation (see Figure 1). As interpolated values, we only use the
     * mean, as
     * this is the most likely missing value, and not subject to a possible random divergence.
     * 
     * @note The automatic search for the correct ARIMA model can be varyingly successfull. If it does not work, try setting the parameters manually.
     * 
     * @section params Parameters
     * Mandatory parameters:
     * - `BEFORE_WINDOW` : The time before a gap that will be used to accumulate data to fit the ARIMA model.
     * - `AFTER_WINDOW` : The time after a gap that will be used to accumulate data to fit the ARIMA model.
     * (BEFORE_WINDOW + AFTER_WINDOW < window_size)
     *
     * Optional parameters:
     * - `MAX_P` : The maximum number of AR coefficients to use in the ARIMA model. Default: 8
     * - `MAX_D` : The maximum number of differences to use in the ARIMA model. Default: 3
     * - `MAX_Q` : The maximum number of MA coefficients to use in the ARIMA model. Default: 8
     * - `START_P` : The starting number of AR coefficients to use in the ARIMA model. Default: 2
     * - `START_Q` : The starting number of MA coefficients to use in the ARIMA model. Default: 2
     * - `MAX_P_SEASONAL` : The maximum number of seasonal AR coefficients to use in the ARIMA model. Default: 2
     * - `MAX_D_SEASONAL` : The maximum number of seasonal differences to use in the ARIMA model. Default: 1
     * - `MAX_Q_SEASONAL` : The maximum number of seasonal MA coefficients to use in the ARIMA model. Default: 2
     * - `START_P_SEASONAL` : The starting number of seasonal AR coefficients to use in the ARIMA model. Default: 1
     * - `START_Q_SEASONAL` : The starting number of seasonal MA coefficients to use in the ARIMA model. Default: 1
     * - `SEASONAL_PERIOD` : The period of the seasonal component. Default: 0s (no seasonal component)
     * - `LIKELIHOOD_METHOD` : The method used to fit the ARIMA model. Default: CSS-MLE\n
     *      Options are:
     *      - `CSS-MLE` : Conditional Sum of Squares - Maximum Likelihood Estimation
     *      - `ML` : Maximum Likelihood Estimation
     *      - `CSS` : Conditional Sum of Squares
     * - `OPTIMIZATION_METHOD` : The optimization method used to fit the ARIMA model. Default: BFGS\n
     *     Options are:
     *      - `BFGS` : Broyden–Fletcher–Goldfarb–Shanno algorithm
     *      - `Nelder-Mead` : Nelder-Mead method
     *      - `Newton_Line_Search` : Newton Line Search
     *      - `Newton_Trust_Region_Hook_Step` : Newton Trust Region Hook Step
     *      - `Newton_Trust_Region_Double_Dog_Leg` : Newton Trust Region Double Dog Leg
     *      - `Conjugate_Gradient` : Conjugate Gradient
     *      - `LBFGS` : Limited Memory BFGS
     *      - `BFGS_MTM` : BFGS Using More Thuente Method
     * - `STEPWISE` : Whether to use stepwise search of the best ARIMA model. Default: true (faster and more robust)
     * - `APPROXIMATION` : Whether to use approximation to determin the Information Criteria and the Likelihood. Default: true
     * - `NUM_MODELS` : The number of models to try when using stepwise search. Default: 94
     * - `SEASONAL` : Whether to use a seasonal component in the ARIMA model. Default: true
     * - `STATIONARY` : Whether to use a stationary ARIMA model. Default: false
     * - `NORMALIZATION` : The normalization method used to fit the ARIMA model. Default: Nothing\n
     *     Options are:
     *      - `MINMAX` : Min-Max Normalization
     *      - `ZSCORE` : Z-Score Normalization
     *      - `NOTHING` : No Normalization
     * - `VERBOSE` : Whether to print additional information. Default: false
     * 
     * It is also possible to set the parameters of the ARIMA model manually. However, only the forward part. It is optionally supported,
     * that a gap is filled in backwards with an auto arima model as well. Be careful when using this though, as it might lead to different results as you have 
     * been expecting. This can be done via these Keywords:
     * 
     * - `SET_MANUAL` : Whether to set the ARIMA parameters manually. Default: false
     * - `FILL_BACKWARD_MANUAL` : Whether to fill the gap backwards with an ARIMA model. Default: false
     * - `P` : The number of AR coefficients to use in the ARIMA model. Default: 0
     * - `D` : The number of differences to use in the ARIMA model. Default: 0
     * - `Q` : The number of MA coefficients to use in the ARIMA model. Default: 0
     * - `P_SEASONAL` : The number of seasonal AR coefficients to use in the ARIMA model. Default: 0
     * - `D_SEASONAL` : The number of seasonal differences to use in the ARIMA model. Default: 0
     * - `Q_SEASONAL` : The number of seasonal MA coefficients to use in the ARIMA model. Default: 0
     * 
     * @note When setting the parameters manually, please provide a seasonal period, as this will now not be estimated anymore, and per default assumed to be 0.
     * 
     * @note BEWARE: The algorithm will accumulate all available data in the before and after window, so if there is another gap in this window, which might also be outside of the 
     * specified resampling area, it will just linearly interpolate the data.
     * 
     * @note The accumulation window of the data, i.e. the BEFORE_WINDOW and AFTER_WINDOW, should be chosen carefully. It has a big impact on the performance, as this will be what the ARIMA model will see.
     *
     *
     * @code
     * [Interpolations1D]
     * TA::resample = ARIMA
     * TA::ARIMA::BEFORE_WINDOW = 86400
     * TA::ARIMA::AFTER_WINDOW = 86400
     *
     * @endcode
     *
     * @note In the case that only random/random walk arima models are found, the missing values will not be filled (It would be just the
     * mean otherwise)
     *
     * @section intro_arima Introduction to ARIMA
     *
     * Autoregressive Integrated Moving Average (ARIMA) is a method for forecasting on historic data. It assumes a stochastic process, with
     * errors that are uncorrelated and have a mean of zero. The ARIMA model is a generalization of an autoregressive moving average (ARMA)
     * which assumes
     * stationary (statistically stable over time) data. Autoregressive refers to a prediciton based on past values of the data. Integrated
     * refers to
     * the differencing of the data to make it stationary. Moving average refers to the use of past errors.
     *
     * To account for seasonality of the data this model can be extended to Seasonal ARIMA (SARIMA). See the <a
     * href="https://towardsdatascience.com/time-series-forecasting-with-arima-sarima-and-sarimax-ee61099e78f6#:~:text=SARIMA%20models%20allow%20for%20differencing,search%20frameworks%20such%20as%20pmdarina."
     * >following article</a> for more information.
     *
     * @subsection model Model Formulation
     *
     * @subsubsection AR Autoregressive (AR) Component:
     *
     *  The autoregressive component of order p (denoted as AR(p)) represents the correlation between the current observation and its
     *      \f$p\f$ past observations. The general formula for AR(\f$p\f$) is:
     *
     *      \f[
     *          Y_t = \phi_1 Y_{t-1} + \phi_2 Y_{t-2} + \ldots + \phi_p Y_{t-p} + a_t
     *      \f]
     *
     *  Here:
     *
     *  - \f$Y_t\f$ is the value at time t,
     *
     *  - \f$\phi_t\f$ are the autoregressive coefficients,
     *
     *  - \f$Y_{t-i}\f$ are the past observations,
     *
     *  - \f$a_t is\f$ a white noise term (random error) at time
     *
     * @subsubsection Int Integrated (I) Component:
     *
     *  The integrated component of order \f$d\f$ (denoted as I(\f$d\f$)) is responsible for differencing the time series data to achieve
     *      stationarity. The differenced series is denoted as \f$\mathbf{Y'}\f$ and is defined as:
     *
     *  \f[
     *      Y'_t = Y_t - Y_{t-d}
     *  \f]
     *
     *  Repeat differencing d times until stationarity is achieved.
     *
     * @subsubsection ma Moving Average (MA) Component:
     *
     *  The moving average component of order \f$q\f$ (denoted as MA(\f$q\f$)) represents the correlation between the current observation
     * and
     *      \f$q\f$ past white noise terms. The general formula for MA(\f$q\f$) is:
     *
     *  \f[
     *      Y_t = a_t - \theta_1 a_{t-1} - \theta_2 a_{t-2} - \ldots - \theta_q a_{t-q}
     *  \f]
     *
     *  Here:
     *
     *  - \f$\theta\f$ are the moving average coefficients,​
     *
     *  - \f$a_{t-q}\f$ are the past white noise terms.
     *
     *
     * @subsubsection arima ARIMA(p, d, q) Model:
     * Combining the AR, I, and MA components, an ARIMA(\f$p, d, q\f$) model is expressed as:
     *
     *  \f[
     *      Y'_t = \phi_1 Y'_{t-1} + \phi_2 Y'_{t-2} + \ldots + \phi_p Y'_{t-p} + a_t - \theta_1 a_{t-1} - \theta_2 a_{t-2} - \ldots -
     * \theta_q a_{t-q}
     *  \f]
     *
     * @subsection estimation Parameter Estimation
     *
     * The parameters of the model are estimated by either maximising the likelihood function or minimising the sum of squared errors.
     *  The likelihood function measures the probability of observing the given data given the model parameters. Optimization is done using
     *  optimization routines like <a href="https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm">BFGS</a>  or the <a
     *  href="https://en.wikipedia.org/wiki/Nelder–Mead_method">Nelder-Mead</a> method.
     *
     *  Maximum Likelihood estimation usually gives better results but is computationally more expensive. As default a mix between the two:
     *  CSS-MLE is used, with BFGS to optimize:
     *
     * @code
     * TA::ARIMA::LIK_METHOD = CSS-MLE
     * TA::ARIMA::OPT_METHOD = BFGS
     * @endcode
     *
     *
     * @subsection selection Model Selection
     * Usually, the parameter selection (p, d, q)(P,D,Q) is done by minimizing the Akaike Information Criterion <a
     * href="https://en.wikipedia.org/wiki/Akaike_information_criterion">AIC</a> or the Bayesian Information Criterion <a
     * href="https://en.wikipedia.org/wiki/Bayesian_information_criterion">BIC</a>.
     * Which parameterize the goodness of a fit. Either an extensive search is done, i.e. all possible combinations between p=0...max_P,
     * q=0...max_Q, ... are tried. Or a stepwise search is done, following the Implementation of <a
     * href="https://www.jstatsoft.org/article/view/v027i03">Hyndman and Khandakar (2008)</a>.
     *
     * As default a stepwise search is done. Additionally, an approximation of the ICs can be used to speed up the search, and the maxima of
     * the parameters can be chosen as well.
     *
     * @code
     * TA::ARIMA::STEPWISE = TRUE
     * TA::ARIMA::APPROXIMATION = TRUE
     * TA::ARIMA::MAX_P = 8
     * TA::ARIMA::MAX_D = 3
     * .
     * .
     * .
     * @endcode
     *
     *
     * In the case of a known seasonal period, e.g. yearly cycles, the seasonal component can be added to the model, by setting the seasonal
     * period (in seconds), if it is 0 it will be estimated automatically (not always reliable):
     * @code
     * TA::ARIMA::SEASONAL_PERIOD = 31536000
     * @endcode
     *
     * Not yet implemented:
     * If you are familiar with finding the best ARIMA parameters via the PACF and ACF plots, you can also set the parameters manually:
     * @code
     * TA::ARIMA::P = 2
     * TA::ARIMA::D = 1
     * TA::ARIMA::Q = 2
     * TA::ARIMA::P_SEASONAL = 1
     * TA::ARIMA::D_SEASONAL = 1
     * TA::ARIMA::Q_SEASONAL = 1
     * @endcode
     *
     * @subsection forecasting Forecasting and Interpolation
     *
     * Forecasting is then done by using the fitted model on the available data, and predicting the next time steps. As random errors are
     * needed, a prediction
     * will return the mean of predicted values, and not a single simulation (one draw of the random error).
     *
     * If data is available before and after a gap, two ARIMA models are fitted, one to predict forward the other to predict backward in
     * time. A weighted
     * average will then be computed and then used as actual prediction.
     *
     * As the ARIMA model only works with constant sampling rates, the data is resampled to the most likely sampling rate. If a requested
     * point falls in between available data, it will be linearly interpolated.
     *
     * \image html arima_simulation.png "Figure 1: Prediction and Simulation using an ARIMA model."
     * 
     * The cyan line is the prediction (mean)
     * with standard deviation (blue 1\f$\sigma\f$ and red 2\f$\sigma\f$). The Blue and Red lines show two simulations, each with a different
     * draw of random errors.
     * 
     * \image html arima_interpolation.png "Figure 2: Interpolation using an ARIMA model."
     * 
     * Green shows the data given to the Interpolation
     * Algorithm, in blue its output. And orange the original test data.
     *
     *
     * @note The performance of the ARIMA model can vary and also depends on the given data. In general the more data the better, and when
     * computing time is not
     * an issue an extensive search might be useful
     *
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
    class ARIMAResampling : public ResamplingAlgorithms {
    public:
        ARIMAResampling(const std::string &i_algoname, const std::string &i_parname, const double &dflt_window_size,
                        const std::vector<std::pair<std::string, std::string>> &vecArgs);

        // Performs ARIMA interpolation the first time it is called, as this is also the first time the data is available
        void resample(const std::string &stationHash, const size_t &index, const ResamplingPosition &position, const size_t &paramindex,
                      const std::vector<MeteoData> &vecM, MeteoData &md);
        std::string toString() const;

    private:
        bool verbose;
        // ARIMA related data
        std::vector<ARIMA_GAP> gap_data;
        std::vector<std::vector<double>> filled_data;
        std::vector<std::vector<Date>> all_dates;

        // Window parameters
        double before_window, after_window;

        // User defined Metadata
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        double period = 0;
        ObjectiveFunction method = CSS_MLE;
        OptimizationMethod opt_method = BFGS;
        bool stepwise = true, approximation = true;
        int num_models = 94;
        bool seasonal = true, stationary = false;
        Normalization::Mode normalize = Normalization::Mode::Nothing;

        bool set_arima_manual = false;
        bool fill_backward_manual = false;
        int p = 0, d = 0, q = 0;
        int P = 0, D = 0, Q = 0;

        // Flags
        bool is_zero_possible = false;
        bool checked_vecM = false;
        bool gave_warning_end = false;
        bool gave_warning_start = false;
        bool gave_warning_interpol = false;
        std::vector<bool> is_valid_gap_data;
        std::vector<bool> warned_about_gap;

        ARIMA_GAP newest_gap;
        // Private methods
        void setMetaData(InterpolARIMA &arima);
        std::vector<double> predictData(std::vector<double> &data, const std::string &direction, size_t startIdx_interpol,
                                        size_t length_gap_interpol, int sr_period);

        // Helper methods for resample
        void checkZeroPossibility(const std::vector<MeteoData> &vecM, size_t paramindex);
        bool processKnownGaps(const Date &resampling_date, const size_t paramindex,
                              const ResamplingAlgorithms::ResamplingPosition &position, const std::vector<MeteoData> &vecM, MeteoData &md);
        double interpolVecAt(const std::vector<MeteoData> &vecM, const size_t &idx, const Date &date, const size_t &paramindex);
        double interpolVecAt(const std::vector<double> &data, const std::vector<Date> &dates, const size_t &pos, const Date &date);
        double interpolVecAt(std::vector<MeteoData> &vecMet, const Date &date, const size_t &paramindex);
        void resampleInterpolationData(size_t &length_gap_interpol, size_t &endIdx_interpol, size_t &startIdx_interpol,
                                       const ARIMA_GAP &new_gap, const Date &data_start_date, const Date &data_end_date,
                                       std::vector<MeteoData> &data_vec_before, std::vector<MeteoData> &data_vec_after,
                                       bool has_data_before, bool has_data_after, size_t paramindex, std::vector<double> &data,
                                       std::vector<Date> &dates, size_t length);
        std::vector<double> getInterpolatedData(std::vector<double> &data, size_t size_before, size_t size_after, size_t startIdx_interpol,
                                                size_t length_gap_interpol, int period);
        void cacheGap(const std::vector<double> &interpolated_data, const std::vector<Date> &interpolated_dates, const ARIMA_GAP &new_gap);

        // info
        void infoARIMA(InterpolARIMA arima);
    };
} // end namespace mio

#endif
