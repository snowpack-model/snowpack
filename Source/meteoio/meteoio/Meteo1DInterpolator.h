// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef METEO1DINTERPOLATOR_H
#define METEO1DINTERPOLATOR_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/Config.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <meteoio/meteoFilters/ProcessingBlock.h> //for ProcessingProperties

#include <string>
#include <vector>
#include <map>
#include <memory>

namespace mio {
/**
 * @page dev_1Dinterpol How to write a resampling algorithm
 * Timeseries can be resampled to any arbitrary timestep through the use of resampling algorithms. The user
 * defines for each parameter which algorithm should be used and whenever a requested timestep can not be found
 * in the original data, the user-defined resampling algorithms will be called for every parameter available in MeteoData.
 *
 * Please keep in mind that the input sampling rate can be anything (including variable) as well as the output
 * sampling rate!
 *
 * @section structure_1Dinterpol Structure
 * When a specific timestep is requested and can not be found in the original data (but lies within the time range of the original data), the
 * TimeSeriesManager calls the MeteoProcessor and requests a temporal interpolation. The request is forwarded to a Meteo1DInterpolator
 * that computes what should be the index of the element just after the one that needs to be computed as well as some flags to
 * know if is the first or the last element in the timeseries. Then the Meteo1DInterpolator calls for each meteorological parameters the
 * algorithm that is responsible for temporally interpolating this parameter with the above computed information. Depending on the
 * user configuration for the algorithm (that is parsed and setup in the constructor of the algorithm), the requested point is computed
 * (or the algorithm simply returns if it can not compute a value). Some helper functions are defined in the ResamplingAlgorithms
 * class that is inherited by every algorithm.
 *
 * @section implementation_1Dinterpol Implementation
 * Using the template.cc and template.h files in meteoio/meteoResampling, build your own algorithm:
 *  - rename template.cc and template.h into a proper name for your algorithm as well as all mentions of "TEMPLATE" in the files;
 *  - declare your cc file in meteoio/meteoResampling/CMakeLists.txt;
 *  - declare your class in meteoio/meteoResampling/ResamplingAlgorithms.cc (both as \#include and in the object factory
 * ResamplingAlgorithmsFactory::getAlgorithm);
 *  - implement your arguments parsing in your constructor;
 *  - implement the temporal interpolation in the resample() method. You receive the full vector of data, the index at or right after where
 * you should compute the interpolation, a hint about some specifics of this location and a single MeteoData object where you should fill
 * the value for the parameter given by paramindex;
 *  - implement an informative toString() method that could help for debugging.
 *
 * In order to debug and test your implementation, it is recommended to use a very small buffer (keys BUFFER_SIZE and BUFF_BEFORE)
 * so it is easy to follow what is going on. You can also always switch back and forth between <i>Enable_Resampling = true</i> and
 * <i>false</i> in the <i>[Interpolations1D]</i> section to see what your algorithm is doing. Please keep in mind that you can not make
 * any assumption regarding the sampling rate, so each neighbouring point has to be checked for validity.
 *
 * @section doc_1Dinterpol Documentation
 * The newly added resampling algorithm must be added to the list of available algorithms in
 * ResamplingAlgorithms.cc with a proper description and a link to its documentation.
 * Please feel free to add necessary bibliographic references to the bibliographic section!
 *
*/

class ResamplingStack {
    public:
		ResamplingStack();

        void addAlgorithm(std::shared_ptr<ResamplingAlgorithms> algo, const double& max_gap_size);
        std::vector<std::shared_ptr<ResamplingAlgorithms>> buildStack(const ResamplingAlgorithms::gap_info& gap) const;

		void resetResampling();
		void resample(const std::string &stationHash, const size_t &index, const ResamplingAlgorithms::ResamplingPosition elementpos, const size_t &par_idx, const std::vector<MeteoData> &vecM, MeteoData &md, const double& max_gap_size) const;
		std::string getStackStr() const;
		bool empty() const;

    private:
        std::vector<double> max_gap_sizes;
        std::vector<std::shared_ptr<ResamplingAlgorithms>> stack;
};


/**
 * @class Meteo1DInterpolator
 * @brief A class that can resample MeteoData objects
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-06-24
*/
class Meteo1DInterpolator {
	public:

		/**
		* @brief The default constructor
		* Set up the interpolation algorithm for each parameter
		* Init tasklist: a vector that holds one std::string for each parameter,
		*                representing the interpolation algorithm that will be executed
		*                for the respective parameter
		*                e.g. tasklist for TA: linear
		* taskargs:      a vector that holds the respective arguments for the algorithms
		*                as a std::vector<std::string>, so there can be multiple arguments
		*
		* @param[in] in_cfg Config object that holds the MeteoFilter configuration in the [Filters] section
		* @param[in] rank in case of multiple TimeSeriesManager, rank in the stack? (default: 1)
		* @param[in] mode spatial resampling operation mode (see IOUtils::OperationMode), default IOUtils::STD
		*/
		Meteo1DInterpolator(const Config& in_cfg, const char& rank=1, const IOUtils::OperationMode &mode=IOUtils::STD);
		Meteo1DInterpolator(const Meteo1DInterpolator& org) = default;

		/**
		 * @brief A function that executes all the resampling algorithms that have been setup in the constructor
		 * @param[in] date The requested date for a MeteoData object (to be resampled if not present)
		 * @param[in] stationHash A unique identifier for each timeseries (that could be used as an index for caching)
		 * @param[in] vecM A vector of MeteoData where the new object will be inserted if not present
		 * @param[in] md new MeteoData element, filled with the resampled values
		 * @return true if successfull, false if no resampling was possible (no element created)
		 */
		bool resampleData(const Date& date, const std::string& stationHash, const std::vector<MeteoData>& vecM, MeteoData& md);

		/**
		 * @brief Call each ResamplingAlgorithms to reset its cached data (as might be needed after a rebuffer)
		 */
		void resetResampling();

		void getWindowSize(ProcessingProperties& o_properties) const;

		Meteo1DInterpolator& operator=(const Meteo1DInterpolator&); ///<Assignement operator
		const std::string toString() const;

 	private:
		std::vector< std::pair<std::string, std::string> > getArgumentsForAlgorithm(const std::string& parname, const std::string& algorithm) const;
		std::string getAlgorithmsForParameter(const std::string& parname) const;

		void processAlgorithms(const std::string& parname, const std::vector<std::pair<int, std::string>>& vecAlgos, std::string base_parname="", const IOUtils::OperationMode& mode=IOUtils::STD, const char& rank=1);
		void createResamplingStacks(const IOUtils::OperationMode& mode, const char& rank);
		// resampling stack helpers
	    void addAlgorithmToStack(const std::string& parname,const std::string& algo_name ,const std::vector<std::pair<std::string, std::string>>& vecArgs, const double& i_max_gap_size);
    	void createDefaultAlgorithm(const std::string &parname);

		std::map< std::string, ResamplingStack > mapAlgorithms; //per parameter interpolation algorithms
		const Config& cfg;
		double max_gap_size; ///< In seconds
		bool enable_resampling, data_qa_logs; ///< easy way to turn resampling off
		std::string gap_size_key; // To support window size for now; TODO: remove at some point

	public:
		static const std::string interpol_section;
		static const std::string interpol_pattern;
		static const std::string arg_pattern;
};


} //end namespace

#endif
