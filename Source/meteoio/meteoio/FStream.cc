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

#include <sys/stat.h>
#include <sstream>
#include <regex>
#include <fstream>
#include <iostream>

#include <meteoio/FileUtils.h>
#include <meteoio/FStream.h>

namespace mio {

std::string ofilestream::cutPathToCWD(const std::string &path)
{
	std::string outpath;
	if (FileUtils::isAbsolutePath(path)) {
		std::stringstream cwds( FileUtils::getCWD() );
		std::vector<std::string> cwd_parts;
		std::string item;
		while (std::getline(cwds, item, '/')) {
			cwd_parts.push_back(item);
		}

		std::stringstream ps(path);
		while (std::getline(ps,item,'/')) {
			if (std::find(cwd_parts.begin(), cwd_parts.end(), item)!=cwd_parts.end()) {
				continue;
			} else {
				outpath += item;
				outpath += "/";
			}
		}
		if (outpath.empty()) return "./";
		else if (outpath.back()=='/') outpath.pop_back();
	} else {
		std::regex e("\\.\\.\\/");
		outpath = std::regex_replace(path, e, "");
		std::regex e2("\\.\\/");
		outpath = std::regex_replace(outpath, e2,"");
	}
	outpath = "./"+ outpath;
	return outpath;
}

std::string ofilestream::limitAccess(std::string path, const bool& write_directories)
{
	if (write_directories) {
#ifdef LIMIT_WRITE_ACCESS
		if (FileUtils::isAbsolutePath(path) && warn_abs_path) {
			std::cerr << "Output path is absolute, i.e. trying to access home directory or similar, which is not allowed."<<std::endl;
			std::cerr << "Creating directory at " << cutPathToCWD(FileUtils::cleanPath(path,true)) << std::endl;
			warn_abs_path = false;
		}

		if (FileUtils::directoryExists(path)) {
			path = cutPathToCWD(FileUtils::cleanPath(path,true));
		} else {
			path = cutPathToCWD(path);
		}

#endif
		if (!FileUtils::directoryExists(path)) {
			FileUtils::createDirectories(path, false);
		}
	} else {
#ifdef LIMIT_WRITE_ACCESS
		const std::string cwd( FileUtils::getCWD() );
		const std::string clean_path( FileUtils::cleanPath(path,true) );
		if (clean_path.find(cwd)==std::string::npos || (cwd.substr(0,4)!=clean_path.substr(0,4))) {
			std::cerr <<"Write access was restricted, but directories are not supposed to be created. Making it impossible to use the specified directory" << std::endl;
			throw IOException("Unqualified directory path "+path+"\n Please make sure you are not trying to access outside of the directory, or set WRITE_DIRECTORIES to true in the configuration file");
		}
#endif
	}

    return path;
}


std::string ofilestream::initialize(const char* filename, const Config& cfgreader)
{
	const bool write_directories = cfgreader.get("WRITE_DIRECTORIES", "Output", true);
	return initialize(filename, write_directories);
}

std::string ofilestream::initialize(const char* filename, bool write_directories)
{
	const std::string path( FileUtils::getPath(filename) );
	std::string file( FileUtils::getFilename(filename) );
	if (keep_old_files) {
		const std::string extension( FileUtils::getExtension(file) );
		if (!extension.empty()) {
			file = FileUtils::removeExtension(file) + "_" + FileUtils::getDateTime() +extension;
		} else {
			file = file + "_" + FileUtils::getDateTime();
		}
	}
#if !defined _WIN32 && !defined __MINGW32__
	if (FileUtils::isWindowsPath(path)) throw IOException("Windows paths are not allowed for UNIX systems");
#endif
	const std::string FILE( limitAccess(path, write_directories) + "/" + file );
	warn_abs_path = false;
	return FILE;
}

void ofilestream::open(const char* filename, std::ios_base::openmode mode)
{
	std::ofstream::open(initialize(filename, write_directories_default).c_str(), mode);
}

ofilestream::ofilestream(const char* filename, std::ios_base::openmode mode) : std::ofstream(initialize(filename, write_directories_default).c_str(), mode), warn_abs_path(true)
{}

ofilestream::ofilestream(const std::string filename, std::ios_base::openmode mode) : std::ofstream(initialize(filename.c_str(), write_directories_default).c_str(), mode), warn_abs_path(true)
{}

ofilestream::ofilestream(const char* filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initialize(filename, cfgreader).c_str(),mode), warn_abs_path(true)
{}

ofilestream::ofilestream(const std::string filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initialize(filename.c_str(), cfgreader).c_str(), mode), warn_abs_path(true)
{}

bool ofilestream::write_directories_default = true;
bool ofilestream::keep_old_files = false;

bool ofilestream::getDefault()
{
	return ofilestream::write_directories_default;
}

}
