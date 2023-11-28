// SPDX-License-Identifier: LGPL-3.0-or-later
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

#ifndef FSTREAM_H
#define FSTREAM_H

#include <sstream>
#include <fstream>
#include <meteoio/IOManager.h>

namespace mio {

/**
* @class ofilestream
* @brief A class that extends std::ofstream, adding some output functionality. Limiting the write access of the software, 
*        writing non-existing output directories, and adding a timestamp to output filenames.
*
* @author Patrick Leibersperger
* @date   2023-11-21
*/
class ofilestream : public std::ofstream
{
    public:
        ofilestream() {}
        ofilestream(const char* filename, std::ios_base::openmode mode = std::ios_base::out);
        ofilestream(const char* filename, const Config& cfgreader, std::ios_base::openmode mode = std::ios_base::out);
        ofilestream(const std::string filename, std::ios_base::openmode mode = std::ios_base::out);
        ofilestream(const std::string filename, const Config& cfgreader, std::ios_base::openmode mode = std::ios_base::out);


        void open(const char* filename, std::ios_base::openmode mode = std::ios_base::out);
        static void createDirectoriesOfFile(const char* filename);

        static bool getDefault();
        static std::string getLimitBaseDir();

    private:
        static std::string initializeFilesystem(const char* filename, const Config& cfgreader);
        static std::string initializeFilesystem(const char* filename, bool write_directories);
        static bool write_directories_default;
        static bool keep_old_files;
        friend void IOManager::setOfstreamDefault(const Config& i_cfg);

        static std::string cutPathToLimitDir(const std::string &path);
        static std::string limitAccess(std::string path, const bool& write_directories);

        static bool warn_abs_path;
};
}

#endif
