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
/**
 * @page ofstream_wrapper ofilestream
 * @section opt Output options
 * This wrapper provides some custom functionality for the writing of files. It is based on std::ofstream.
 * 
 * It is possible to limit the write access of MeteoIO. This is done by setting the cmake macro RESTRICT_WRITE_ACCESS to ON. (Default: off)
 * By default the access is limited to the current working directory This can be changed by setting the cmake macro RESTRICTED_BASE_DIR to the desired directory. (Default: ".")
 * 
 * The following configuration keys are available to control the output:
 * - WRITE_DIRECTORIES: If set to true, the software will create non-existing directories in the output path. If set to false, the software will only write to existing directories. Default: true
 * - KEEP_OLD_FILES: If set to true, the software will not overwrite existing files, but create a new file with a timestamp in the filename. Default: false
 *
 * @section Example
 * @code
 * [Output]
 * WRITE_DIRECTORIES = false
 * KEEP_OLD_FILES = true
 * @endcode
 */


/**
* @brief Adjust the path, so that it starts from the current working directory.
* @details
* The path is checked, whether it points to a directory outside of the current working directory. If so, the path is adjusted to point to the current working directory.
* If an absolute path is given, it will be prefixed by ./, possibly creating a lot of directories, if an absolute path to a restricted area was given.
* @param path The output path that needs to be checked.
* @return outpath The adjusted path.
*/
std::string ofilestream::cutPathToLimitDir(const std::string &path)
{
	static const std::regex regex_up_path("\\.\\.\\/", std::regex::optimize);
	static const std::regex regex_cwd_path("\\.\\/", std::regex::optimize);
	
	std::string outpath;
	if (FileUtils::isAbsolutePath(path)) { //processing absolute path
		std::stringstream lim_dirs( FileUtils::cleanPath(getRestrictedBaseDir(), true) );
		std::vector<std::string> lim_dir_parts;
		std::string item;
		while (std::getline(lim_dirs, item, '/')) {
			lim_dir_parts.push_back(item);
		}

		std::stringstream ps(path);
		size_t ii = 0;
		while (std::getline(ps,item,'/')) {
			if (ii >= lim_dir_parts.size() || lim_dir_parts[ii]!=item) {
				outpath += item;
				outpath += "/";
			}
			ii++;
		}
		if (outpath.empty()) return getRestrictedBaseDir();
		else if (outpath.back()=='/') outpath.pop_back();
	} else { //processing relative path
		outpath = path;
		std::replace(outpath.begin(), outpath.end(), '\\', '/');
		outpath = std::regex_replace(path, regex_up_path, "");
		outpath = std::regex_replace(outpath, regex_cwd_path,"");
	}
	outpath = getRestrictedBaseDir() +"/"+ outpath;
	return outpath;
}


/**
* @brief Handles the limited access depending on different options and input
* @details
* Throws an error if directories are not supposed to be created, but the given path is outside of the permitted area.
* @param path output filename and path
* @param write_directories whether directories are supposed to be created
* @return The adjusted path.
*/
std::string ofilestream::limitAccess(std::string path, const bool& write_directories)
{
	path = FileUtils::cleanPath(path, false);
	if (write_directories) {
#ifdef RESTRICT_WRITE_ACCESS
		if (!FileUtils::directoryExists(getRestrictedBaseDir())) {
			throw IOException("Write permissions are restricted to the '"+getRestrictedBaseDir()+"' base directory, but it does not exist");
		}
		if (FileUtils::isAbsolutePath(path)) {
			if (warn_abs_path) {
				std::cerr << "[W] This version of MeteoIO has been compiled with restricted write permissions, absolute output paths are not allowed. ";
				std::cerr << "[W] Instead creating directory at " << cutPathToLimitDir(FileUtils::cleanPath(path, true, true)) << std::endl;
			}
			path = cutPathToLimitDir(FileUtils::cleanPath(path, true, true));
		}
		else {
			path = cutPathToLimitDir(path);
		}

#endif
		if (!FileUtils::directoryExists(path)) {
			FileUtils::createDirectories(path);
		} 
	} else { //creating new directories is not allowed
#ifdef RESTRICT_WRITE_ACCESS
		if (!FileUtils::directoryExists(getRestrictedBaseDir())) {
			throw IOException("Write permissions are restricted to the '"+getRestrictedBaseDir()+"' directory, but it does not exist...");
		}

		if (!FileUtils::isSubDirOf(path, getRestrictedBaseDir())) {
			std::ostringstream ss;
			ss << "Could not write to '" << path << "'. This version of MeteoIO has been compiled with restricted write permissions, it is therefore not allowed to write outside of '" << getRestrictedBaseDir() << "'";
			throw IOException(ss.str(), AT);
		}

		if (!FileUtils::directoryExists(path)) {
			throw IOException("Attempting to create '" + path + "', but the creation of new directories is not allowed, set WRITE_DIRECTORIES to true in the configuration file to change this behavior", AT);
		}
#endif
	}

    return path;
}


std::string ofilestream::initializeFilesystem(const char* filename, const Config& cfgreader)
{
	const bool write_directories = cfgreader.get("WRITE_DIRECTORIES", "Output", true);
	return initializeFilesystem(filename, write_directories);
}

std::string ofilestream::initializeFilesystem(const char* filename, bool write_directories)
{
	const std::string path( FileUtils::getPath(filename) );
	std::string file( FileUtils::getFilename(filename) );
	if (keep_old_files) {
		const std::string extension( FileUtils::getExtension(file) );
		if (!extension.empty()) {
			file = FileUtils::removeExtension(file) + "_" + FileUtils::getDateTime()+"." +extension;
		} else {
			file = file + "_" + FileUtils::getDateTime();
		}
	}
#if !defined _WIN32 && !defined __MINGW32__
	if (FileUtils::isWindowsPath(path)) throw IOException("Windows paths are not allowed for UNIX systems");
#endif
	const std::string fileAndPath( limitAccess(path, write_directories) + "/" + file );
	warn_abs_path = false;
	return fileAndPath;
}


/**
* @brief The actual writing function
* @details
* works the same as std::ofstream::open, but with the additional functionality
* @param filename file to write
* @param mode mode to open the file in
*/
void ofilestream::open(const char* filename, std::ios_base::openmode mode)
{
	std::ofstream::open(initializeFilesystem(filename, write_directories_default).c_str(), mode);
}

/**
* @brief Constructor
* @param filename file to write
* @param mode mode to open the file in
*/
ofilestream::ofilestream(const char* filename, std::ios_base::openmode mode) : std::ofstream(initializeFilesystem(filename, write_directories_default).c_str(), mode)
{}

ofilestream::ofilestream(const std::string filename, std::ios_base::openmode mode) : std::ofstream(initializeFilesystem(filename.c_str(), write_directories_default).c_str(), mode)
{}

/**
* @brief Constructor
* @param filename file to write
* @param cfgreader instance of Config, to read the inishell config keywords
* @param mode mode to open the file in
*/
ofilestream::ofilestream(const char* filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initializeFilesystem(filename, cfgreader).c_str(),mode)
{}

ofilestream::ofilestream(const std::string filename, const Config& cfgreader, std::ios_base::openmode mode) : std::ofstream(initializeFilesystem(filename.c_str(), cfgreader).c_str(), mode)
{}

bool ofilestream::write_directories_default = true;
bool ofilestream::keep_old_files = false;
bool ofilestream::warn_abs_path = true;


bool ofilestream::getDefault()
{
	return ofilestream::write_directories_default;
}

std::string ofilestream::getRestrictedBaseDir()
{
#if defined  RESTRICTED_BASE_DIR
	return RESTRICTED_BASE_DIR;
#else
	return ".";
#endif
}

void ofilestream::createDirectoriesOfFile(const char* filename)
{
	initializeFilesystem(filename, write_directories_default);
}

}
