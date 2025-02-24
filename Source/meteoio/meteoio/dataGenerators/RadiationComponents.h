// SPDX-License-Identifier: LGPL-3.0-or-later
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
#ifndef RADCOMPONENTS_H
#define RADCOMPONENTS_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class RadiationComponents
 * @brief Compute the global radiation (ISWR) from the direct and diffuse components.
 * @details The split radiation must be nammed ISWR_DIR and ISWR_DIFF to be recognized.
 * @code
 * [Generators]
 * ISWR::generator1 = RadComponents
 * @endcode
 */
class RadiationComponents : public GeneratorAlgorithm {
	public:
		RadiationComponents(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ)
			: GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), where( section+"::"+algo ) { parse_args(vecArgs); }
		bool generate(const size_t& param, MeteoData& md, const std::vector<MeteoData>& vecMeteo) override;
		bool create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo) override;
	private:
		const std::string where;
};

} //end namespace mio

#endif
