/*
 * HeatEquationAnalytical.h
 *
 *  Created on: Jul 7, 2017
 *      Author: julien
 */

#ifndef TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_
#define TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include "snowpack/DataClasses.h"

struct Error {

		std::vector<double> absError, temperature;
		double rmsError, maxError;

		Error()
				: absError(),
				  temperature(),
				  rmsError(0.0),
				  maxError(0.0) {

		}

		void clear() {

			absError.clear();
			temperature.clear();
			rmsError = 0.0;
			maxError = 0.0;

		}

};



class HeatEquationAnalytical {
 public:
		explicit HeatEquationAnalytical(const size_t & nTestLayers,
		                                const double & totalThickness,
		                                const double &T0,
                                    const double &T1,
		                                const double &heatCapacity,
		                                const double &snowDensity = 400.0,
		                                const double &thermalConductivity = 1.0,
		                                const unsigned int &p = 1000);

		virtual ~HeatEquationAnalytical();

		bool computeErrors(const size_t &t, const size_t &n,
		                   const std::vector<NodeData> &nodesData,
		                   Error &errors);

 private:

  inline void updateAlpha();

  double _T0, _T1;

  double _heatCapacity;
  double _snowDensity;
  double _thermalConductivity;



  unsigned int _p;

		double _alpha;
		double _totalThickness;
		size_t _nTestLayers;
		double _layerTestThickness;

		double getTemperatureAt(const size_t& n, const size_t &t);

		std::vector<double> getTemperatureAtTestsPoints(const size_t &t);

		std::vector<double> interpolate1D(const size_t &n,
		                                  const std::vector<NodeData> &nodesData);


};

#endif /* TESTS_IMPLICITSOLVER_HEATEQUATIONANALYTICAL_H_ */
