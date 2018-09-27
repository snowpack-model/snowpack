/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef FilterDespikingPS_H
#define FilterDespikingPS_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <meteoio/meteoStats/libfit1D.h>

#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterDespikingPS
 * @brief Despiking filter using phase space.
 * @details
 * This filter iteratively finds and replaces spikes (outliers) in a data sequence until no more spikes are found. The filter was developed to despike Acoustic
 * Doppler Velocimeter (ADV) data, but might also be useful to remove outliers from other signals.
 *
 * \image html DespikingFilter_typicalFilterResult.png "Example of filter performance for one hour of ADV-data. Outliers are successfully removed by the filter."
 *
 * Two versions of this filter are implemented:
 *   - The "Goring"-version is implemented according to following paper:
 *      Goring, D.G. and Nikora, V.I. (2002). <i>"Despiking Acoustic Doppler Velocimeter Data"</i> J. Hydraul. Eng., 128(1): 117-126
 *   - And the "Mori"-version is implemented according to this paper:
 *      Mori, N., Suzuki, T. and Kakuno, S. (2007). <i>"Noise of Acoustic Doppler Velocimeter Data in Bubbly Flows"</i> J. Eng. Mech., 133(1):122-125
 *
 * Method:
 * - Detection of the spikes:
 *  - calculate first and second derivatives of the whole signal
 *  - the data points plotted in phase space are enclosed by an ellipsoid defined
 *     by the standard deviations and the universal threshold (=sqrt(2*ln(number of data points))) (see figure below)
 *  - points outside this ellipsoid are designated as spikes (the Goring-implementation uses three 2D-projections of the ellipsoid to identify
 *                                                            the outliers, while the Mori-implementation finds the outliers in 3D)
 *
 * - Replacement of the spikes:
 *  - find a cubic fit for 24 data points around the spike
 *  - replace the spike with a fitted value
 *
 * Parameters:
 * - Adjustable parameters:
 *  - The sensitivity parameter was added to be able to control the sensitivity of the filter. The universal threshold is divided by
 *      this value. Thus a sensitivity of 1 is the default value. The larger the sensitivity the smaller the threshold (=the smaller
 *      the ellipsoid) and thus the more spikes will be detected.
 *  - The method parameter decides which implementation to use: "Mori" or "Goring". The differences are small between both implementations.
 *    According to Mori the "Mori"-version performs slightly better. However, we suggest to use the Goring-method, because it is better tested.

 * - Hard-coded parameters:
 *    - the maximum number of iterations for the spike detection. This is set to 50.
 *    - the number of points used for the fitting (replacement) of the spikes. This is set to 24 (just as proposed in the paper by Goring).
 *    - the degree of the fit for the replacement. This is set to 3 (cubic fit) (as proposed in the paper by Goring).
 *
 * \image html DespikingFilter_phaceSpacePlots.png "2D projections of the phase space plots. Points outside the ellipses are spikes."
 *
 * @ingroup processing
 * @author Thiemo Theile
 * @date   2018-08-15
 * @code
 * VW::filter1	= despiking
 * VW::arg1::sensitivity = 1
 * VW::arg1::method = GORING
 * @endcode
 */


class FilterDespikingPS : public FilterBlock {
	public:
		FilterDespikingPS(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum IMPLEMENTATION_TYPE {
			GORING,
			MORI
		} implementation_type;

		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		std::vector<int> findSpikesGoring(const std::vector<double>& timeVec, const std::vector<double>& ivec, unsigned int& nNewSpikes);
		std::vector<int> findSpikesMori(const std::vector<double>& timeVec, const std::vector<double>& ivec, unsigned int& nNewSpikes);
		void replaceSpikes(const std::vector<double>& timeVec, std::vector<double>& uVec, std::vector<int>& spikesVec);

		std::vector<double> calculateDerivatives(const std::vector<double>& ivec, const std::vector<double>& timeVec);
		double calculateCrossCorrelation(const std::vector<double>& aVec,const std::vector<double>& bVec);
		unsigned int nNodataElements(const std::vector<double>& iVec);
		void findPointsOutsideEllipse(const std::vector<double>& xVec,const std::vector<double>& yVec,
                                      const double a,const double b,const double theta, std::vector<int>& outsideVec);
		void findPointsOutsideEllipsoid(const std::vector<double>& xVec,const std::vector<double>& yVec, const std::vector<double>& zVec,
                                                const double a,const double b,const double c, std::vector<int>& outsideVec);
		void getWindowForInterpolation(const size_t index,const std::vector<double>& timeVec, const std::vector<double>& uVec,
                                       const std::vector<int>& spikesVec, const unsigned int& windowWidth, std::vector<double>& xVec,
                                       std::vector<double>& yVec);
		bool checkIfWindowForInterpolationIsSufficient(const std::vector<double>& xVec,const double time,const unsigned int minPoints,
                                                       const bool avoidExtrapolation);
		//helper functions:
		void solve2X2LinearEquations(const double* a, const double* b, const double* c, double* x);
		const std::vector<double> helperGetDoubleVectorOutOfMeteoDataVector(const std::vector<MeteoData>& ivec, const unsigned int& param);
		const std::vector<double> helperGetTimeVectorOutOfMeteoDataVector(const std::vector<MeteoData>& ivec);
		void helperWriteDebugFile1DerivativesAndFittedEllipses (const std::vector<double>& uVec, const std::vector<double>& duVec,
                                                                const std::vector<double>& du2Vec,
                                                                double a1,double b1,double a2,double b2,double a3,double b3,double theta);
		void helperWriteDebugFile2Interpolation(const std::vector<double>& uVec, const std::vector<int>& spikesVec,
                                                const std::vector<double>& x, const std::vector<double>& y,
                                                const Fit1D& quadraticFit,unsigned int iiSpike);
		void helperWriteDebugFile3WindowForInterpolation(size_t ii,std::vector<double>& x,std::vector<double>& y);
		void helperWriteDebugFile4OriginalAndFinalSignal(std::vector<double>& ivec, std::vector<double>& ovec,
                                                         std::vector<int>& allSpikesVec);

		double sensitivityParam; //this parameter controls the sensitivity of the filter. standard value: 1
		implementation_type methodParam; //this parameter controls which implementation of the filter is used: according to Mori or Goring
		int nIterations;   //this counts the iterations
		int maxIterations; //this is a hard-coded parameter to stop the iteration
};

} //end namespace

#endif
