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
#include <meteoio/meteoResampling/ResamplingStack.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>

using namespace std;

namespace mio {

ResamplingStack::ResamplingStack() : max_gap_sizes(), stack() {}

void ResamplingStack::addAlgorithm(std::shared_ptr<ResamplingAlgorithms> algo, const double &max_gap_size) 
{
	stack.push_back(algo);
	max_gap_sizes.push_back(max_gap_size);
}

std::vector<std::shared_ptr<ResamplingAlgorithms>> ResamplingStack::buildStack(const ResamplingAlgorithms::gap_info &gap) const
{
	std::vector<std::shared_ptr<ResamplingAlgorithms>> res;
	for (size_t ii = 0; ii < stack.size(); ii++) {
		if (max_gap_sizes[ii] == IOUtils::nodata)
			res.push_back(stack[ii]);
		else if (gap.size() <= max_gap_sizes[ii]) {
			res.push_back(stack[ii]);
		}
	}
	return res;
}

void ResamplingStack::resetResampling()
{
	for (size_t ii = 0; ii < stack.size(); ii++) {
		stack[ii]->resetResampling();
	}
}

//TODO The current implementation is not efficient as well as does not handle exact matches properly.
//We should build once and for all the resampling stack (for every parameter) in the constructor
//(consider having a '*' parameter as default for parameters that have not been explicitly configured)
//and then just loop over the resampling algos with the gap that we have, asking who can take care of it.
//Thus each specific algo could recognize an exact match and handle it...
void ResamplingStack::resample(const std::string &stationHash, const size_t &index, const ResamplingAlgorithms::ResamplingPosition elementpos, const size_t &par_idx, const std::vector<MeteoData> &vecM,
                                MeteoData &md, const double &i_max_gap_size) const
{
	const ResamplingAlgorithms::gap_info gap( ResamplingAlgorithms::findGap(index, par_idx, vecM, md.date, i_max_gap_size) );
	const std::vector<std::shared_ptr<ResamplingAlgorithms>> resampling_stack = buildStack(gap);

	for (size_t jj = 0; jj < resampling_stack.size(); jj++) {
		resampling_stack[jj]->resample(stationHash, index, elementpos, par_idx, vecM, md);
		if (jj > 0 && (index != IOUtils::npos) && vecM[index](par_idx) != md(par_idx)) {
			break;
		} else if (ResamplingAlgorithms::exact_match == elementpos && vecM[index](par_idx) == md(par_idx)) {
			break;
		}
	}
}

bool ResamplingStack::empty() const { return stack.empty(); }

std::string ResamplingStack::getStackStr() const
{
	ostringstream os;
	os << "[";
	for (size_t ii = 0; ii < stack.size(); ii++) {
		os << stack[ii]->getAlgo();
		if (ii < stack.size() - 1)
			os << ", ";
	}
	os << "]";
	return os.str();
}

// Methods moved from Meteo1DInterpolator to ResamplingStack
void ResamplingStack::addAlgorithmToStack(const std::string& parname, const std::string& algo_name, 
                                        const std::vector<std::pair<std::string, std::string>>& vecArgs, 
                                        const double& i_max_gap_size, const Config& cfg)
{
	const std::shared_ptr<ResamplingAlgorithms> algo_ptr( ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, i_max_gap_size, vecArgs, cfg) );
	addAlgorithm(algo_ptr, i_max_gap_size);
}

void ResamplingStack::createDefaultAlgorithm(const std::string &parname, const double& i_max_gap_size, const Config& cfg)
{
	std::vector<std::pair<std::string, std::string>> vecArgs; // is empty anyways as we already iterated through all specified algorithms
	const std::shared_ptr<ResamplingAlgorithms> algo_ptr( ResamplingAlgorithmsFactory::getAlgorithm("LINEAR", parname, i_max_gap_size, vecArgs, cfg) );
	addAlgorithm(algo_ptr, i_max_gap_size);
}

} // namespace
