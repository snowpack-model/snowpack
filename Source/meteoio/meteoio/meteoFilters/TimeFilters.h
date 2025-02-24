// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2017 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef TIMEFILTERS_H
#define TIMEFILTERS_H

#include <meteoio/IOUtils.h>
#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <vector>
#include <string>

namespace mio {
/**
 * @class  TimeProcStack
 * @brief Since the time filters are quite specific to TIME (and need to be applied before), they have their own
 * ProcessingStack.
 */
class TimeProcStack {
	public:
		TimeProcStack(const Config& i_cfg);
		TimeProcStack(const TimeProcStack& source);
		~TimeProcStack() {for (size_t ii=0; ii<filter_stack.size(); ii++) delete filter_stack[ii];}
		
		void process(std::vector< std::vector<MeteoData> >& ivec);
		void process(Date &dateStart, Date &dateEnd);

		TimeProcStack& operator=(const TimeProcStack& source);
		static void checkUniqueTimestamps(std::vector<METEO_SET> &vecVecMeteo);
		static const std::string timeParamName;
		
	private:
		std::vector<ProcessingBlock*> filter_stack; //for now: strictly linear chain of processing blocks
		bool enable_time_filtering;
};

/**
 * @class  TimeSuppr
 * @ingroup processing
 * @brief Timesteps suppression filter.
 * @details
 * This filter deletes some timesteps based on the provided arguments:
 *  - TYPE: defines the strategy to delete timesteps. It is one of:
 *       - CLEANUP: merge duplicated and out-of-order timestamps (a warning will be emitted anyway for the problematic timestamp). Merge are only performed when they don't conflict (two identical fields with different, non-nodata value is a conflict)
 *       - BYDATES: delete specific timesteps
 *       - FRAC: suppress a given fraction of the data at random
 *  - FILE: when type=BYDATE, a file that contains a list of station ID's and timesteps that should be suppressed;
 *  - FRAC: when type=FRAC, the fraction of data to suppress. For example, <i>0.5</i> would ensure that at least <i>50%</i> of the
 * data set's points are deleted.
 *  - WIDTH: when type=FRAC, the width of data gaps to create (in seconds). If not set, individual data points are deleted.
 *
 * @code
 * TIME::filter1     = suppr
 * TIME::arg1::type = cleanup
 * 
 * TIME::filter2     = suppr
 * TIME::arg2::type = bydates
 * TIME::arg2::file = ./input/meteo/suppr.dat
 * @endcode
 * 
 * The file <i>suppr.dat</i> would look like this (the time is given in the timezone declared in Input::TIME_ZONE):
 * @code
 * *WFJ 2015-10-01T12:00
 * *DAV 2015-10-02T15:00
 * *WFJ 2015-11-10T06:00
 * *WFJ 2015-12-25T01:00 2015-12-27T13:30
 * *WFJ 2015-09-01T07:15 - 2015-09-10T20:30
 * STB2 2015-10-01T21:30
 * @endcode
 * Time ranges are declared by providing two dates on the same line. For more visibility, the said two dates can be separated by " - " (which a white
 * space character on both sides, as shown in the example above).
 */
class TimeSuppr : public ProcessingBlock {
	public:
		TimeSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) override;

	private:
		//possible modes of operation
		typedef enum SUPPR_MODE {
		            NONE,
		            BYDATES,
		            FRAC,
		            CLEANUP
		} suppr_mode;
		
		void supprByDates(std::vector<MeteoData>& ovec) const;
		void supprFrac(std::vector<MeteoData>& ovec) const;
		void supprInvalid(std::vector<MeteoData>& ovec) const;
		
		std::map< std::string, std::vector<DateRange> > suppr_dates;
		double range, width;
		suppr_mode op_mode;
};

/**
 * @class  TimeShift
 * @ingroup processing
 * @brief Time corrections.
 * @details
 * This filter add offsets to the time at specified timestamps. This is for example usefull to remove the 
 * <A HREF="https://en.wikipedia.org/wiki/Daylight_saving_time">Daylight Saving Time</A> in order to bring 
 * the timestamps back to Winter time only (or "Standard Time", as it should always be!) or to correct 
 * a clock that has drifted in a data logger. 
 * 
 * This takes either a fixed offset or a correction file that contains on each line an ISO formatted 
 * timestamp as well as an offset (in seconds) to apply to the timestamps starting at the provided time
 * (always assumed to be in the input timezone + a potential previous offset). So the supported arguments 
 * are:
 *    + CORRECTIONS: the correction file (optional);
 *    + OFFSET: the fixe offset (in seconds) to add to the time (optional).
 *
 * @code
 * TIME::filter1     = SHIFT
 * TIME::arg1::CORRECTIONS = ./input/meteo/dst.dat
 * @endcode
 * 
 * The file <i>dst.dat</i> would look like this (the time is given in the timezone declared in 
 * Input::TIME_ZONE, first without offset and with -3600s offset for the second one):
 * @code
 * 2016-03-27T02:00 -3600
 * 2016-10-30T03:00 0
 * @endcode
 * 
 * @note The current implementation does not handle data conflicts if a shifted timestep conflicts with another timestep. In
 * such a case, an error message "[E] timestamp error for station..." will be returned and you will have to manually decide 
 * in your input data which timestep to keep.
 * 
 * @note The correction file supports comments (everything behind ';' or '#' is a comment) and empty lines.
 */
class TimeShift : public ProcessingBlock {
	public:
		TimeShift(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) override;

	private:
		void process_corrections(std::vector<MeteoData>& ovec) const;
		void process_offset(std::vector<MeteoData>& ovec) const;
		std::vector<offset_spec> dst_changes;
		double offset;
		bool has_offset;
};

/**
 * @class  TimeSort
 * @ingroup processing
 * @brief Sort out of order timesteps.
 * @details This filter guarantees that the timestamps are sorted in increasing order. This is convenient to fix known problems, but
 * please do not use this filter blindly: it is most probably a better idea to know that some timestamps are not in increasing
 * order than to always resort them without knowing if there are underlying problems with the dataset...
 * 
 * @code
 * TIME::filter1 = sort
 * @endcode
 * 
 */
class TimeSort : public ProcessingBlock {
	public:
		TimeSort(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) override;
};

/**
 * @class  TimeLoop
 * @ingroup processing
 * @brief Loops over a specific time period.
 * @details
 * This filter repeats a given time period over an over (for example, for simulation spin-up). It takes the following arguments (all dates are 
 * ISO formatted date where the time can be ommitted to mean "00:00:00"):
 *    - REF_START: start of the repeated period, mandatory argument;
 *    - REF_END: end of the repeated period, mandatory argument;
 *    - MATCH_DATE: date in the output data that matches the start date of the repeated period, mandatory argument;
 *
 * The MATCH_DATE argument is used to define a "synchronization point": from this point on, the reference data will be copied over and over 
 * (please keep in mind that copying exactly a month of data will still lead to some drift over the years as not every month has the same
 * number of days).
 * 
 * @code
 * TIME::filter1     = TimeLoop
 * TIME::arg1::ref_start = 2018-01-01	;this assumes 00:00:00
 * TIME::arg1::ref_end = 2018-02-01
 * TIME::arg1::match_date = 2019-03-01	;this assumes that the data from ref_start will be copied as data for 2019-03-01
 * @endcode
 * 
 * The example above copies over and over the timestamps from 2018-01-01T00:00:00 to 2018-02-01T00:00:00 while making sure that the timesteps 
 * of 2018-01-01T00:00:00 are copied to 2019-03-01T00:00:00. To do a ten years simulation spin up, then run your simulation for example as
 * "my_model 2000-01-01 2010-01-01". 
 */
class TimeLoop : public ProcessingBlock {
	public:
		TimeLoop(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec) override;
		void process(Date &dateStart, Date &dateEnd) override;

	private:
		Date req_start, req_end;
		Date match_date, ref_start, ref_end;
};

} //end namespace

#endif
