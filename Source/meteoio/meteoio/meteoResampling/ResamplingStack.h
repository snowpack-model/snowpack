// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2026 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef RESAMPLINGSTACK_H
#define RESAMPLINGSTACK_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <meteoio/Config.h>

#include <string>
#include <vector>
#include <memory>

namespace mio {

class ResamplingStack {
    public:
		ResamplingStack();

		void addAlgorithm(std::shared_ptr<ResamplingAlgorithms> algo, const double& max_gap_size);
		std::vector<std::shared_ptr<ResamplingAlgorithms>> buildStack(const ResamplingAlgorithms::gap_info& gap) const;

		void resetResampling();
		void resample(const std::string &stationHash, const size_t &index, const ResamplingAlgorithms::ResamplingPosition elementpos, const size_t &par_idx, const std::vector<MeteoData> &vecM, MeteoData &md, const double& max_gap_size) const;
		std::string getStackStr() const;
		bool empty() const;

		// Stack configuration methods
		void addAlgorithmToStack(const std::string& parname, const std::string& algo_name, 
		                         const std::vector<std::pair<std::string, std::string>>& vecArgs, 
		                         const double& i_max_gap_size, const Config& cfg);
		void createDefaultAlgorithm(const std::string &parname, const double& i_max_gap_size, const Config& cfg);

    private:
        std::vector<double> max_gap_sizes;
        std::vector<std::shared_ptr<ResamplingAlgorithms>> stack;
};

} //end namespace

#endif
