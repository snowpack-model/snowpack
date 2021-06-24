// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCSHIFT_H
#define PROCSHIFT_H

#include <meteoio/meteoFilters/ProcessingBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcShift
 * @brief Time shifting filter for the selected meteo parameter
 * @ingroup processing
 * @details
 * This filter can correct the time base of a single parameter within a station, either with a constant offset or
 * by reading time varying offsets from a provided file. Compared to the TimeShift time filter, this filter operates
 * on only one parameter so it can shift parameters compared to each other within a given station dataset.
 * 
 * It can also compute the offset between two meteo parameter 
 * (one being used as reference parameter) by computing at each timestamp the offset that leads to the highest
 * Pearson's correlation coefficient and writing this timeseries of offsets to a file. This correlation 
 * coefficient is evaluated over a given time window and sampling rate (the data is first temporarily forced-resampled
 * at this constant sampling rate before computing the correlation coefficient). A large data window is more accurate
 * but will fail to properly capture sudden changes while a short data window might generate some spurious offsets.
 * In any case, it is advised to plot and tweak the generated timeseries of offsets first before using it to 
 * correct the meteo parameter (keep in mind that the offsets are given as the correction that should be applied to
 * temporally re-synchronize the meteo parameter of choice, so for example +600 means that the signal would need 
 * to be moved forward by 600 seconds in order to be synchronized with the reference parameter).
 * 
 * This filter takes the following arguments:
 *  - EXTRACT_OFFSETS: if set to TRUE, no correction is applied but the offsets between two meteo parameters
 *                     are computed (defaut: false);
 *      - OFFSETS_FILE: the file where the computed time offsets will be written (mandatory);
 *      - REF_PARAM: the name of the reference meteo parameter (mandatory);
 *      - SAMPLING_RATE: the sampling rate to use when computing the time offsets (default: automatically extracted from the data)
 *      - WIDTH: data window width (in seconds) over which to compute each correlation coefficient (default: 2 days);
 *      - OFFSET_RANGE: range of allowed variation (in seconds) for the offset when searching for the maximum correlation (default: 1 day, so looking for 
 *                      a maximum between 1/2 a day before and after the current point);
 *  - if EXTRACT_OFFSETS is set to false, the meteo parameter will be corrected.
 *      - TYPE: the type of correction to apply, either CST (constant over the whole dataset) or STEPWISE 
 *              (keeping the last correction until finding a new correction). Default is CST;
 *      - CST: when using the CST type, the offset value (in seconds, mandatory in this case);
 *      - OFFSETS_FILE: for other correction types, a timeseries of offsets (in seconds, mandatory in this case);
 *      - TZ: the time zone for the timestamps given in the OFFSETS_FILE (mandatory).
 * 
 * A typical workflow for this filter would be:
 *    1. run on the data with EXTRACT_OFFSETS = TRUE in order to compute the temporal offsets between the
 *       meteo parameter of interest and a reference one;
 *    2. plot the resulting offsets timeseries and manually pickup the points to keep in the correction file
 *       (the computed offsets might be noisy, specially around periods of nodata). In doubt, look back at the
 *       original data to decide if you want to keep or reject a given offset.
 *    3. run on the data with EXTRACT_OFFSETS = FALSE, providing your correction file.
 * 
 * Example of configuration to compute the time offset between TA_2 and TA_1 used as reference and write the
 * results in the TA2_extracted.dat file (located in the current working directory):
 * @code
 * TA_2::FILTER5 = SHIFT
 * TA_2::ARG5::EXTRACT_OFFSETS = TRUE
 * TA_2::ARG5::OFFSETS_FILE = TA2_extracted.dat
 * TA_2::ARG5::REF_PARAM = TA_1
 * @endcode
 * 
 * Applying the corrections provided in the TA2_corrections.dat file (located in the same directory as the ini file):
 * @code
 * TA_2::FILTER5 = SHIFT
 * TA_2::ARG5::TYPE = STEPWISE
 * TA_2::ARG5::OFFSETS_FILE = TA2_offsets.dat
 * TA_2::ARG5::TZ = 1
 * @endcode
 * The TA2_corrections.dat file gives a an hour forward correction from 2018-12-03T00:00, then 2 
 * hours backward from 2019-03-12T12:00 until 2019-10-01T15:30 when there is no correction any more:
 * @code
 * 2018-12-03T00:00 3600
 * 2019-03-12T12:00 -7200
 * 2019-10-01T15:30 0
 * @endcode
 */

class ProcShift : public ProcessingBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
	public:
		ProcShift(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum INTERPOL_TYPE {
			cst,
			stepwise,
			linear
		} interpol_type;
		
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void writeOffsets(const unsigned int& param, const std::vector<MeteoData>& ivec);
		void correctOffsets(const unsigned int& param, std::vector<MeteoData>& ovec);
		
		static bool isAllNodata(const std::vector<ProcessingBlock::offset_spec>& vecX, const size_t& startIdx, const size_t& endIdx);
		static double getMedianSampling(const size_t& param, const std::vector<MeteoData>& ivec);
		
		std::vector<offset_spec> resampleVector(const std::vector<MeteoData>& ivec, const size_t& param) const;
		double getPearson(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx, const int& offset) const;
		int getOffsetFullScan(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx, const int& range_min, const int& range_max) const;
		int getOffset(const std::vector<ProcessingBlock::offset_spec>& vecX, const std::vector<ProcessingBlock::offset_spec>& vecY, const size_t& curr_idx) const;
		
		std::vector<ProcessingBlock::offset_spec> corrections; ///< Corrections to apply to the data, read from the user provided file
		std::string ref_param; ///< The reference parameter to compare to
		std::string root_path;
		std::string offsets_file; ///< File name that contains the extracted offsets or the correction offsets
		double cst_offset; ///< Constant correction offset, to be provided by the user
		double sampling_rate; ///< dataset sampling rate to use, either automatically extracted or provided by the user
		double offset_range; ///< range of time offsets to consider, in days
		double width_d; ///< size of the data window in days
		size_t width_idx; ///< size of the data window over which to compute the correlation
		interpol_type offsets_interp; ///< type of interpolation to use to interpolate the provided offsets
		bool extract_offsets; ///< do not apply any correcxtion but extract the time-varying offsets from the data
};

} //end namespace

#endif
