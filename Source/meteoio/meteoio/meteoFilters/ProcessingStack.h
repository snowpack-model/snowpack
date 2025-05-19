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
#ifndef PROCESSINGSTACK_H
#define PROCESSINGSTACK_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/Config.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcessingStack
 * @brief This builds and runs through a filter stack for filtering a given parameter.
 * @author Thomas Egger
 * @date   2011-01-11
 */
class ProcessingStack {
	public:
		/**
		 * @brief Constructor parses cfg and builds up a filter stack for param_name
		 */
		ProcessingStack(const Config& cfg, const std::string& param_name);
		virtual ~ProcessingStack() {for (size_t ii=0; ii<filter_stack.size(); ii++) delete filter_stack[ii];}

		void process(const std::vector< std::vector<MeteoData> >& ivec,
		             std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass=false);
		void getWindowSize(ProcessingProperties& o_properties) const;
		const std::string toString() const;
		
		static const std::string filter_section, filter_pattern, arg_pattern;
		
	private:
		virtual bool applyFilter(const size_t& param, const size_t& jj, const std::vector<MeteoData>& ivec, std::vector<MeteoData> &ovec);
		virtual bool filterStation(std::vector<MeteoData> ivec, std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass, const size_t& stat_idx);
		
		std::vector<ProcessingBlock*> filter_stack; //for now: strictly linear chain of processing blocks
		const std::string param_name;
		bool data_qa_logs;
};

} //end namespace

#endif
