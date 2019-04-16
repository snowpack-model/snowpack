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
		TimeProcStack(const Config& cfg);
		~TimeProcStack() {for (size_t ii=0; ii<filter_stack.size(); ii++) delete filter_stack[ii];}
		
		void process(std::vector< std::vector<MeteoData> >& ivec, const bool& second_pass);
		static void checkUniqueTimestamps(const std::vector<METEO_SET>& vecVecMeteo);
		static const std::string timeParamName;
		
	private:
		std::vector<ProcessingBlock*> filter_stack; //for now: strictly linear chain of processing blocks
};


/**
 * @class  TimeSuppr
 * @ingroup processing
 * @brief Timesteps suppression filter.
 * @details
 * This filter deletes some timesteps based on the provided arguments:
 *  - SUPPR: provide a file that contains a list of station ID's and timesteps that should be suppressed;
 *  - FRAC: suppress a given fraction of the data at random. For example, <i>0.5</i> would ensure that at least <i>50%</i> of the
 * data set's points are deleted.
 *
 * @code
 * TIME::filter1     = suppr
 * TIME::arg1::suppr = ./input/meteo/suppr.dat
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
		TimeSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec);

	private:
		void supprByDates(std::vector<MeteoData>& ovec) const;
		void supprFrac(std::vector<MeteoData>& ovec) const;
		
		std::map< std::string, std::vector<dates_range> > suppr_dates;
		double range;
};

/**
 * @class  TimeUnDST
 * @ingroup processing
 * @brief Timesteps Daylight Saving Time correction.
 * @details
 * This filter removes the <A HREF="https://en.wikipedia.org/wiki/Daylight_saving_time">Daylight Saving Time</A> in order to bring the timestamps back to Winter time only (or "Standard Time", as it should always be!). In order to
 * do so, a correction file has to be provided that contains on each line an ISO formatted timestamp as well as an offset (in seconds) to apply
 * to the timestamps starting at the provided time. 
 *
 * @code
 * TIME::filter1     = UnDST
 * TIME::arg1::CORRECTIONS = ./input/meteo/dst.dat
 * @endcode
 * 
 * The file <i>dst.dat</i> would look like this (the time is given in the timezone declared in Input::TIME_ZONE, with or without the DST):
 * @code
 * 2016-03-27T02:00 -3600
 * 2016-10-30T03:00 0
 * @endcode
 * Where the first change point is in Winter time while the second one is in Summer time.
 * 
 */
class TimeUnDST : public ProcessingBlock {
	public:
		TimeUnDST(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec);

	private:
		std::vector<offset_spec> dst_changes;
};

} //end namespace

#endif
