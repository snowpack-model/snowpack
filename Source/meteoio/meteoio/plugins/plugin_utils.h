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
#ifndef PLUGIN_UTILS_H
#define PLUGIN_UTILS_H

#include <meteoio/FileUtils.h>
#include <meteoio/Config.h>

#include <vector>

namespace mio
{
namespace PLUGIN
{

    std::vector<std::string> getFilesWithPaths(const std::vector<std::string> &vecFilenames, const std::string &inpath, const std::string& dflt_extension);
    void scanMeteoPath(const Config &cfg, const std::string &inpath, std::vector<std::string> &vecFilenames, const std::string& pattern);

    
} // namespace PLUGIN
} // namespace mio


#endif // PLUGIN_UTILS_H