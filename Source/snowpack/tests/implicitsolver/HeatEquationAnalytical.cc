/*
 * HeatEquationAnalytical.cpp
 *
 *  Created on: Jul 7, 2017
 *      Author: julien
 */

#include "HeatEquationAnalytical.h"

HeatEquationAnalytical::HeatEquationAnalytical(
    const size_t & nTestLayers,
    const double & totalThickness,
    const double &T0,
    const double &T1, const double &heatCapacity, const double &snowDensity,
    const double &thermalConductivity,
    const unsigned int &p)
    : _T0(T0),
      _T1(T1),
      _heatCapacity(heatCapacity),
      _snowDensity(snowDensity),
      _thermalConductivity(thermalConductivity),
      _p(p),
      _alpha(0.),
		  _totalThickness(totalThickness),
		  _nTestLayers(nTestLayers),
		  _layerTestThickness(totalThickness / (double) nTestLayers) {



  updateAlpha();




}

HeatEquationAnalytical::~HeatEquationAnalytical() {
  // TODO Auto-generated destructor stub
}

double HeatEquationAnalytical::getTemperatureAt(const size_t& n,
                                                const size_t& t) {

	const double x = (double) n * _layerTestThickness;

	if (t <= 0) {
    std::cerr << "t must be greater than 0" << std::endl;
    exit(1);
  }

  if (x < 0 || x > _totalThickness) {
    std::cerr
        << "x shall not me smaller than 0 or bigger than the total thickness "
        << _totalThickness << std::endl;
    exit(1);
  }

  double result = 0.0;

  for (unsigned int p = _p; p > 0; p--) {  // Going backward for numerical reasons...

    double tmp = 2 / M_PI * (_T1 - _T0) / p * (p % 2 ? -1.0 : 1.0);

    tmp *= exp(
		    -_alpha * M_PI * M_PI / _totalThickness / _totalThickness * p * p
		        * static_cast<double>(t));

    tmp *= sin(p * M_PI / _totalThickness * x);

    result += tmp;
  }

  result += (_T1 - _T0) * x / _totalThickness + _T0;

  return result;

}

inline void HeatEquationAnalytical::updateAlpha() {
  _alpha = _thermalConductivity / (_snowDensity * _heatCapacity);
}

bool HeatEquationAnalytical::computeErrors(
    const size_t &t,
    const size_t& n, const std::vector<NodeData>& nodesData,
    Error &errors) {

	errors.clear();

	std::vector<double> analytical = getTemperatureAtTestsPoints(t);
	std::vector<double> approximated = interpolate1D(n, nodesData);
	errors.temperature = approximated;

	double rmsError = 0.0;

	for (size_t N = 0; N < _nTestLayers + 1; N++) {

		double actualAbsError = fabs(analytical[N] - approximated[N]);


		errors.absError.push_back(actualAbsError);
		if (actualAbsError > errors.maxError) {
			errors.maxError = actualAbsError;
		}

		rmsError += pow(actualAbsError, 2.0);

	}

	rmsError /= static_cast<double>(_nTestLayers);

	errors.rmsError = sqrt(rmsError);

	return true;

}

std::vector<double> HeatEquationAnalytical::interpolate1D(
    const size_t& n_given, const std::vector<NodeData>& nodesData) {

	std::vector<double> yInterpolated;

	size_t N = 0;
	for (size_t n = 1; n <= n_given; n++) {

		const double xDiff = 1.0 / (double) n_given;
		const double xLower = xDiff * (double) (n - 1);

		const double yLower = nodesData[n - 1].T;

		const double slope = (nodesData[n].T - nodesData[n - 1].T) / xDiff;

		double xActual = (double) N / (double) _nTestLayers;

		while (xActual <= xLower + xDiff) {

			yInterpolated.push_back((xActual - xLower) * slope + yLower);
			xActual = (double) ++N / (double) _nTestLayers;
		}

	}

	return yInterpolated;

}

std::vector<double> HeatEquationAnalytical::getTemperatureAtTestsPoints(
    const size_t& t) {

	std::vector<double> benchmarkTemp;

	for (size_t n = 0; n <= _nTestLayers; n++) {

		benchmarkTemp.push_back(
		    this->getTemperatureAt(n, t));

	}

	return benchmarkTemp;

}
