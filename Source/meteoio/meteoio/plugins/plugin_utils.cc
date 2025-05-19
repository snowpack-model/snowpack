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

#include <meteoio/plugins/plugin_utils.h>
#include <meteoio/FileUtils.h>

namespace mio
{
    
namespace PLUGIN
{
    
/**
* @brief Get the file names together with paths and check for validity
*
* @param vecFilenames The vector of filenames to check.
* @param inpath The input path
* @param dflt_extension The extension to add to the filename if it is missing
* @return std::vector<std::string> The vector of filenames together with paths
* @throw InvalidNameException if a filename is invalid
*
*/
std::vector<std::string> getFilesWithPaths(const std::vector<std::string> &vecFilenames, const std::string &inpath, const std::string &dflt_extension) {
    std::vector<std::string> all_files_and_paths;
    
    for (const std::string& filename : vecFilenames) {
        const std::string extension( FileUtils::getExtension(filename) );
#if defined _WIN32 || defined __MINGW32__
        const std::string file_and_path = (!extension.empty()) ? inpath + "\\" + filename : inpath + "\\" + filename + dflt_extension;
#else
        const std::string file_and_path = (!extension.empty()) ? inpath + "/" + filename : inpath + "/" + filename + dflt_extension;
#endif

        if (!FileUtils::validFileAndPath( file_and_path )) // Check whether filename is valid
            throw InvalidNameException(file_and_path, AT);
        all_files_and_paths.push_back( file_and_path );
    }
    
    return all_files_and_paths;
}

void scanMeteoPath(const Config &cfg, const std::string &inpath, std::vector<std::string> &vecFilenames, const std::string &pattern) {
    bool is_recursive = false;
    cfg.getValue("METEOPATH_RECURSIVE", "Input", is_recursive, IOUtils::nothrow);
    
    std::list<std::string> dirlist( FileUtils::readDirectory(inpath, pattern, is_recursive) );
    dirlist.sort();
    vecFilenames.reserve( dirlist.size() );
    std::copy(dirlist.begin(), dirlist.end(), std::back_inserter(vecFilenames));
}
} // namespace PLUGIN
} // namespace mio
