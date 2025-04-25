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
#include "IOUtils.h"
#include "meteoio/IOExceptions.h"
#include <cstddef>
#include <meteoio/Config.h>
#include <meteoio/FileUtils.h>
#include <meteoio/FStream.h>
#include <meteoio/thirdParty/tinyexpr.h>

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <regex>

using namespace std;

namespace mio {

static const std::string defaultSection = "GENERAL";
static const char NUM[] = "0123456789";

//these are not defined as part of any class in order to avoid dragging a <regex> include into the Config.h header
//regex to try sorting the keys by numerical order of their index, so STATION2 comes before STATION10
static const std::regex index_regex("([^0-9]+)([0-9]{1,9})$", std::regex::optimize); //limit the number of digits so it fits within an "int" for the call to atoi()
static const std::regex section_regex(R"(^\[((?:\w|\-)+)\].*$)", std::regex::optimize); //valid chars for section: letters, numbers, _ and -

//Constructors
Config::Config() : properties(), sections(), sourcename(), configRootDir() {}

Config::Config(const std::string& i_filename) : properties(), sections(), sourcename(i_filename), configRootDir(FileUtils::getPath(i_filename, true))
{
	addFile(i_filename);
}

const ConfigProxy Config::get(const std::string& key, const std::string& section) const
{
	return ConfigProxy(*this, key, section);
}

template <typename T> T Config::get(const std::string& key, const std::string& section, const T& dflt) const
{
	if (keyExists(key, section)) {
		T tmp;
		getValue(key, section, tmp);
		return tmp;
	} else {
		return dflt;
	}
}

std::string Config::get(const std::string& key, const std::string& section, const std::string& dflt) const
{
	if (keyExists(key, section)) {
		std::string tmp;
		getValue(key, section, tmp);
		return tmp;
	} else {
		return dflt;
	}
}

std::string Config::get(const std::string& key, const std::string& section, const char dflt[]) const
{
	if (keyExists(key, section)) {
		std::string tmp;
		getValue(key, section, tmp);
		return tmp;
	} else {
		return std::string(dflt);
	}
}

double Config::get(const std::string& key, const std::string& section, const double& dflt) const
{
	if (keyExists(key, section)) {
		double tmp;
		getValue(key, section, tmp);
		return tmp;
	} else {
		return dflt;
	}
}

bool Config::get(const std::string& key, const std::string& section, const bool& dflt) const
{
	if (keyExists(key, section)) {
		bool tmp;
		getValue(key, section, tmp);
		return tmp;
	} else {
		return dflt;
	}
}

//Populating the property map
void Config::addFile(const std::string& i_filename)
{
	if (configRootDir.empty()) configRootDir = FileUtils::getPath(i_filename, true);
	sourcename = i_filename;
	const ConfigParser parser( i_filename, properties, sections );
}

void Config::addKey(std::string key, std::string section, const std::string& value)
{
	IOUtils::toUpper(section);
	IOUtils::toUpper(key);
	properties[ section + "::" + key ] = value;
}

void Config::deleteKey(std::string key, std::string section)
{
	IOUtils::toUpper(section);
	IOUtils::toUpper(key);
	properties.erase( section + "::" + key );
}

void Config::deleteKeys(std::string keymatch, std::string section, const bool& anywhere)
{
	IOUtils::toUpper(section);
	IOUtils::toUpper(keymatch);
	const size_t section_len = section.length();

	//Loop through keys, look for match - delete matches
	if (anywhere) {
		std::map<std::string,std::string>::iterator it = properties.begin();
		while (it != properties.end()) {
			const size_t section_start = (it->first).find(section, 0);

			if ( section_start==0 && (it->first).find(keymatch, section_len)!=string::npos )
				properties.erase( it++ ); // advance before iterator become invalid
			else //wrong section or no match
				++it;
		}

	} else {
		keymatch = section + "::" + keymatch;
		std::map<std::string,std::string>::iterator it = properties.begin();
		while (it != properties.end()) {
			if ( (it->first).find(keymatch, 0)==0 ) //match found at start
				properties.erase( it++ ); // advance before iterator become invalid
			else
				++it;
		}
	}
}

bool Config::keyExists(std::string key, std::string section) const
{
	IOUtils::toUpper(section);
	IOUtils::toUpper(key);
	const std::string full_key( section + "::" + key );
	const std::map<string,string>::const_iterator it = properties.find(full_key);
	return (it!=properties.end());
}

// will do a regex search with the given key
bool Config::keyExistsRegex(std::string key_regex, std::string section) const
{
	IOUtils::toUpper(section);
	std::regex full_pattern( section + "::" + key_regex, std::regex::icase | std::regex::optimize);
	std::map<string,string>::const_iterator it = properties.begin();
	for (; it!=properties.end(); ++it) {
		if (std::regex_match(it->first, full_pattern)) return true;
	}
	return false;
}

bool Config::sectionExists(std::string section) const
{
	IOUtils::toUpper( section );
	std::set<std::string>::const_iterator it = sections.begin();

	for (; it!=sections.end(); ++it) {
		if (*it==section) return true;
	}

	return false;
}

void Config::moveSection(std::string org, std::string dest, const bool& overwrite)
{
	IOUtils::toUpper( org );
	IOUtils::toUpper( dest );
	std::map<string,string>::iterator it = properties.begin();

	//delete all current keys in "dest" if overwrite==true
	if (overwrite) {
		while (it != properties.end()) {
			const std::string::size_type pos = it->first.find("::");
			if (pos!=std::string::npos) {
				const std::string sectionname( IOUtils::strToUpper(it->first.substr(0, pos)) );
				if (sectionname==dest) properties.erase( it++ ); // advance before iterator become invalid
				else ++it;
			} else {
				++it;
			}
		}
	}

	//move the keys from org to dest
	it = properties.begin();
	while (it != properties.end()) {
		const std::string::size_type pos = it->first.find("::");
		if (pos!=std::string::npos) {
			const std::string sectionname( IOUtils::strToUpper(it->first.substr(0, pos)) );
			if (sectionname==org) {
				const std::string key( it->first.substr(pos) );
				properties[ dest+key ] = it->second;
				properties.erase( it++ ); // advance before iterator become invalid
			} else {
				++it;
			}
		} else {
			++it;
		}
	}
}

std::vector< std::pair<std::string, std::string> > Config::getValuesRegex(const std::string& regex_str, std::string section) const
{
	//regex for selecting the keys
	static const std::regex user_regex(regex_str, std::regex::icase | std::regex::optimize);
	std::smatch index_matches;
	
	IOUtils::toUpper(section);
	std::vector< std::pair<std::string, std::string> > vecResult;
	std::map< std::pair<int, std::string>, std::pair<std::string, std::string> > keyMap;
	bool indexed_keys = true;
	
	for (const auto& prop : properties) {
		const size_t section_start = (prop.first).find(section, 0);
		if (section_start==0) { //found the section!
			const size_t section_len = section.length();
			const std::string key_no_section( prop.first.substr(section_len+2) ); //account for "::" between section and key
			
			if (std::regex_match(key_no_section, user_regex)) {
				//we want to figure out of the keys are all indexed, ie like {some prefix}{some integral number}
				if (indexed_keys) {
					if (std::regex_match(key_no_section, index_matches, index_regex)) { //retrieve the key index
						const std::string key_root( index_matches.str(1) );
						const int index = atoi( index_matches.str(2).c_str() ); //we take the first capture group, guaranteed to fit in an int
						keyMap[ make_pair(index, key_root) ] = make_pair(key_no_section, prop.second);
					} else {
						indexed_keys = false;
					}
				}
				
				//keys are not indexed, store them directly in the results vector
				//if indexed_keys just got toggled above, we recover already processed keys
				if (!indexed_keys) {
					//the keys are not indexed, move the processed keys into the results vector
					for (const auto& key_record : keyMap) vecResult.push_back( key_record.second );
					keyMap.clear();
					vecResult.push_back( make_pair(key_no_section, prop.second) ); //push the current, unprocessed key
				}
			}
		}
	}
	
	if (indexed_keys && !keyMap.empty()) {
		for (const auto& key_record : keyMap) {
			vecResult.push_back( key_record.second );
		}
	}

	return vecResult;
}

std::vector< std::pair<std::string, std::string> > Config::getValues(std::string keymatch, std::string section, const bool& anywhere) const
{
	std::smatch index_matches;
	IOUtils::toUpper(section);
	IOUtils::toUpper(keymatch);
	
	std::vector< std::pair<std::string, std::string> > vecResult;
	//stores (index, key_root) as map index and (key, value) as map value
	//although the key could be rebuilt from (index, key_root), storing it makes conversion 
	//to the vector of pair to be returned easier and robust
	std::map< std::pair<int, std::string>, std::pair<std::string, std::string> > keyMap;
	bool indexed_keys = true;

	//Loop through keys, look for match - push it into vecResult
	if (anywhere) {
		for (const auto& prop : properties) {
			const size_t section_start = (prop.first).find(section, 0);
			if (section_start==0) { //found the section!
				const size_t section_len = section.length();
				const size_t found_pos = (prop.first).find(keymatch, section_len);
				if (found_pos!=string::npos) { //found it!
					const std::string key( (prop.first).substr(section_len + 2) ); //from pos to the end
					
					//we want to figure out of the keys are all indexed, ie like {some prefix}{some integral number}
					if (indexed_keys) {
						if (std::regex_match(key, index_matches, index_regex)) { //retrieve the key index
							const std::string key_root( index_matches.str(1) );
							const int index = atoi( index_matches.str(2).c_str() ); //we take the first capture group, guaranteed to fit in an int
							keyMap[ make_pair(index, key_root) ] = make_pair(key, prop.second);
						} else {
							indexed_keys = false;
						}
					}
					
					//keys are not indexed, store them directly in the results vector
					//if indexed_keys just got toggled above, we recover already processed keys
					if (!indexed_keys) {
						//the keys are not indexed, move the processed keys into the results vector
						for (const auto& key_record : keyMap) vecResult.push_back( key_record.second );
						keyMap.clear();
						vecResult.push_back( make_pair(key, prop.second) ); //push the current, unprocessed key
					}
				}
			}
		}
	} else {
		keymatch = section + "::" + keymatch;
		for (const auto& prop : properties) {
			const size_t found_pos = (prop.first).find(keymatch, 0);
			if (found_pos==0) { //found it starting at the begining
				const size_t section_len = section.length();
				const std::string key( (prop.first).substr(section_len + 2) ); //from pos to the end
				
				//we want to figure out of the keys are all indexed, ie like {some prefix}{some integral number}
				if (indexed_keys) {
					if (std::regex_match(key, index_matches, index_regex)) { //retrieve the key index
						const std::string key_root( index_matches.str(1) );
						const int index = atoi( index_matches.str(2).c_str() ); //we take the first capture group, guaranteed to fit in an int
						keyMap[ make_pair(index, key_root) ] = make_pair(key, prop.second);
					} else {
						indexed_keys = false;
					}
				}
				
				//keys are not indexed, store them directly in the results vector
				//if indexed_keys just got toggled above, we recover already processed keys
				if (!indexed_keys) {
					//the keys are not indexed, move the processed keys into the results vector
					for (const auto& key_record : keyMap) vecResult.push_back( key_record.second );
					keyMap.clear();
					vecResult.push_back( make_pair(key, prop.second) ); //push the current, unprocessed key
				}
			}
		}
	}
	
	if (indexed_keys && !keyMap.empty()) {
		for (const auto& key_record : keyMap) {
			vecResult.push_back( key_record.second );
		}
	}

	return vecResult;
}


std::vector<std::string> Config::getKeysRegex(const std::string& regex_str, std::string section) const
{
	const std::vector< std::pair<std::string, std::string> > vecKeys( getValuesRegex(regex_str, section) );
	std::vector<std::string> vecResult;
	
	for (const auto& key_record : vecKeys) {
		vecResult.push_back( key_record.first );
	}
	
	return vecResult;
}

std::vector<std::string> Config::getKeys(std::string keymatch,
                        std::string section, const bool& anywhere) const
{
	
	const std::vector< std::pair<std::string, std::string> > vecKeys( getValues(keymatch, section, anywhere) );
	std::vector<std::string> vecResult;
	
	for (const auto& key_record : vecKeys) {
		vecResult.push_back( key_record.first );
	}
	
	return vecResult;
}

void Config::write(const std::string& filename) const
{
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename,AT);
	ofilestream fout(filename.c_str(), ios::out);
	if (fout.fail()) throw AccessException(filename, AT);

	try {
		std::string current_section;
		unsigned int sectioncount = 0;
		for (const auto& prop : properties) {
			const std::string key_full( prop.first );
			const std::string section( ConfigParser::extract_section(key_full, true) );

			if (current_section != section) {
				current_section = section;
				if (sectioncount != 0)
					fout << endl;
				sectioncount++;
				fout << "[" << section << "]" << endl;
			}

			const size_t key_start = key_full.find_first_of(":");
			const std::string value( prop.second );
			if (value.empty()) continue;

			if (key_start!=string::npos) //start after the "::" marking the section prefix
				fout << key_full.substr(key_start+2) << " = " << value << endl;
			else //every key should have a section prefix, but just in case...
				fout << key_full << " = " << value << endl;
		}
	} catch(...) {
		if (fout.is_open()) //close fout if open
			fout.close();

		throw;
	}

	if (fout.is_open()) //close fout if open
		fout.close();
}

const std::string Config::toString() const {
	std::ostringstream os;
	os << "<Config>\n";
	os << "Source: " << sourcename << "\n";
	for (auto it = properties.begin(); it != properties.end(); ++it){
		os << (*it).first << " -> " << (*it).second << "\n";
	}
	os << "</Config>\n";
	return os.str();
}

std::ostream& operator<<(std::ostream& os, const Config& cfg) {
	const size_t s_source = cfg.sourcename.size();
	os.write(reinterpret_cast<const char*>(&s_source), sizeof(size_t));
	os.write(reinterpret_cast<const char*>(&cfg.sourcename[0]), s_source*sizeof(cfg.sourcename[0]));
	
	const size_t s_rootDir = cfg.configRootDir.size();
	os.write(reinterpret_cast<const char*>(&s_rootDir), sizeof(size_t));
	os.write(reinterpret_cast<const char*>(&cfg.configRootDir[0]), s_rootDir*sizeof(cfg.configRootDir[0]));

	const size_t s_map = cfg.properties.size();
	os.write(reinterpret_cast<const char*>(&s_map), sizeof(size_t));
	for (map<string,string>::const_iterator it = cfg.properties.begin(); it != cfg.properties.end(); ++it){
		const string& key = (*it).first;
		const size_t s_key = key.size();
		os.write(reinterpret_cast<const char*>(&s_key), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&key[0]), s_key*sizeof(key[0]));

		const string& value = (*it).second;
		const size_t s_value = value.size();
		os.write(reinterpret_cast<const char*>(&s_value), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&value[0]), s_value*sizeof(value[0]));
	}
	
	const size_t s_set = cfg.sections.size();
	os.write(reinterpret_cast<const char*>(&s_set), sizeof(size_t));
	for (set<string>::const_iterator it = cfg.sections.begin(); it != cfg.sections.end(); ++it){
		const string& value = *it;
		const size_t s_value = value.size();
		os.write(reinterpret_cast<const char*>(&s_value), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&value[0]), s_value*sizeof(value[0]));
	}

	return os;
}

std::istream& operator>>(std::istream& is, Config& cfg) {
	size_t s_source;
	is.read(reinterpret_cast<char*>(&s_source), sizeof(size_t));
	cfg.sourcename.resize(s_source);
	is.read(reinterpret_cast<char*>(&cfg.sourcename[0]), s_source*sizeof(cfg.sourcename[0]));
	
	size_t s_rootDir;
	is.read(reinterpret_cast<char*>(&s_rootDir), sizeof(size_t));
	cfg.configRootDir.resize(s_rootDir);
	is.read(reinterpret_cast<char*>(&cfg.configRootDir[0]), s_rootDir*sizeof(cfg.configRootDir[0]));

	cfg.properties.clear();
	size_t s_map;
	is.read(reinterpret_cast<char*>(&s_map), sizeof(size_t));
	for (size_t ii=0; ii<s_map; ii++) {
		size_t s_key, s_value;
		is.read(reinterpret_cast<char*>(&s_key), sizeof(size_t));
		string key;
		key.resize(s_key);
		is.read(reinterpret_cast<char*>(&key[0]), s_key*sizeof(key[0]));

		is.read(reinterpret_cast<char*>(&s_value), sizeof(size_t));
		string value;
		value.resize(s_value);
		is.read(reinterpret_cast<char*>(&value[0]), s_value*sizeof(value[0]));

		cfg.properties[key] = value;
	}
	
	cfg.sections.clear();
	size_t s_set;
	is.read(reinterpret_cast<char*>(&s_set), sizeof(size_t));
	for (size_t ii=0; ii<s_set; ii++) {
		size_t s_value;
		is.read(reinterpret_cast<char*>(&s_value), sizeof(size_t));
		string value;
		value.resize(s_value);
		is.read(reinterpret_cast<char*>(&value[0]), s_value*sizeof(value[0]));

		cfg.sections.insert( value );
	}
	return is;
}

unsigned int Config::getCommandNr(const std::string& section, const std::string& cmd_pattern, const std::string& cmd_key)
{
	//extract the cmd number and perform basic checks on the syntax
	const size_t end_cmd = cmd_key.find(cmd_pattern);
	if (end_cmd==std::string::npos) 
		throw InvalidArgumentException("Can not parse command number in "+section+"::"+cmd_key, AT);

	const size_t start_cmd_nr = cmd_key.find_first_of(NUM, end_cmd+cmd_pattern.length());
	const size_t end_cmd_nr = cmd_key.find_first_not_of(NUM, end_cmd+cmd_pattern.length());
	if (start_cmd_nr==std::string::npos || end_cmd_nr!=std::string::npos) 
		throw InvalidArgumentException("Can not parse command number in "+section+"::"+cmd_key, AT);

	unsigned int cmd_nr;
	const std::string cmd_nr_str( cmd_key.substr(start_cmd_nr) );
	if ( !IOUtils::convertString(cmd_nr, cmd_nr_str) ) InvalidArgumentException("Can not parse command number in "+cmd_key, AT);
	return cmd_nr;
}

std::vector< std::pair<std::string, std::string> > Config::parseArgs(const std::string& section, const std::string& cmd_id, const unsigned int& cmd_nr, const std::string& arg_pattern) const
{
	//read the arguments and clean them up (ie get all key/values matching {cmd_id}::{arg_pattern}#:: and remove this prefix)
	std::ostringstream arg_str;
	arg_str << cmd_id << arg_pattern << cmd_nr << "::";
	std::vector< std::pair<std::string, std::string> > vecArgs( getValues(arg_str.str(), section) );
	for (auto& arg : vecArgs) {
		const size_t beg_arg_name = arg.first.find_first_not_of(":", arg_str.str().length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+arg.first+"'", AT);
		arg.first = arg.first.substr(beg_arg_name);
	}
	
	return vecArgs;
}

std::vector< std::pair<std::string, std::string> > Config::getArgumentsForAlgorithm(const std::string& parname, const std::string& algorithm, const std::string& section) const
{
	const std::string key_prefix( parname+"::"+algorithm+"::" );
	std::vector< std::pair<std::string, std::string> > vecArgs( getValues(key_prefix, section) );

	//clean the arguments up (ie remove the {Param}::{algo}:: in front of the argument key itself)
	for (auto& arg : vecArgs) {
		const size_t beg_arg_name = arg.first.find_first_not_of(":", key_prefix.length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+arg.first+"'", AT);
		arg.first = arg.first.substr(beg_arg_name);
	}

	return vecArgs;
}

std::vector< std::pair<std::string, std::string> > Config::getArgumentsForAlgorithm(const std::string& parname, const std::string& algorithm, const size_t& algo_index, const std::string& section) const
{
	if (algo_index==IOUtils::npos) return getArgumentsForAlgorithm(parname, algorithm, section);
	
	return getArgumentsForAlgorithm(parname, algorithm+std::to_string(algo_index), section);
}


///////////////////////////////////////////////////// ConfigParser helper class //////////////////////////////////////////
//the local values must have priority -> we initialize from the given i_properties, overwrite from the local values 
//and swap the two before returning i_properties
ConfigParser::ConfigParser(const std::string& filename, std::map<std::string, std::string> &i_properties, std::set<std::string> &i_sections) : properties(i_properties), imported(), sections(), vars(), sourcename(filename)
{
	parseFile( fileProperties(filename, sourcename) );
	
	//expand all variables that might be used in key/values. 
	//Instead of relying on a dependency tree between keys, we just keep on running on all keys that require expansion until 
	//they are all done. If absolutely no replacement is done in a round, this means that there is a circular dependency...
	while (!vars.empty()) {
		bool hasSomeSuccesses = false; //it will be set to true once at least one variable could be expanded
		
		//loop over all key/values that contain variables to expand (and are therefore not yet part of the 'properties' map
		for (std::map<std::string, std::vector<variable> >::iterator it = vars.begin(); it!=vars.end(); ) {
			const std::string section( extract_section(it->first) );
			std::vector<variable>& vecVars( it->second ); //for clarity
			const bool status = expandVarsForKey(vecVars, section, hasSomeSuccesses);
			
			if (status) { //all expansions could be done, the root of the tree contains the final string, remove the key from 'vars'
				properties[it->first] = vecVars[0].value; //save the key/value pair
				it = vars.erase(it);
			} else {
				++it;
			}
		}
		
		if (!hasSomeSuccesses) { //not a single variable could be expanded
			std::string msg("In file "+filename+", the following keys could not be resolved (circular dependency? undefined variable name?):");
			for (const auto& var : vars) msg.append( " "+var.first );
			throw InvalidArgumentException(msg, AT);
		}
	}
	
	//swap the caller and local properties before returning, so the caller gets the new version
	std::swap(properties, i_properties);
	i_sections.insert(sections.begin(), sections.end());
}

ConfigParser::FILE_PPT::FILE_PPT(const std::string& filename, const std::string& sourcename)
             : original_name(filename), restrict_section(), clean_name()
{
	//extract a section name appended to the filename, if any
	const size_t section_delim = filename.find_last_of(':');
	if (section_delim!=std::string::npos) {
		original_name = filename.substr(0, section_delim);
		restrict_section = filename.substr(section_delim+1);
		IOUtils::toUpper( restrict_section );
	}

	//resolve symlinks, resolve relative path w/r to the path of the current ini file
	//if this is a relative path, prefix the import path with the current path
	const std::string prefix = ( FileUtils::isAbsolutePath(original_name) )? "" : FileUtils::getPath(sourcename, true)+"/";
	const std::string path( FileUtils::getPath(prefix+original_name, true) );  //clean & resolve path
	const std::string clean_filename( FileUtils::getFilename(original_name) );
	clean_name = path + "/" + clean_filename;
}


bool ConfigParser::onlyOneEqual(const std::string& str)
{
	const size_t firstEqualPos = str.find('=');
	if (firstEqualPos != std::string::npos) {
		const size_t secondEqualPos = str.find('=', firstEqualPos + 1);
		if (secondEqualPos != std::string::npos && str[secondEqualPos - 1] == '\\') {
			return onlyOneEqual(str.substr(secondEqualPos + 1));
		}
		return secondEqualPos == std::string::npos;
	}
	
	return true;
}

/**
* @brief Parse the whole file, line per line
* @param[in] filename file to parse
*/
void ConfigParser::parseFile(const fileProperties& iniFile)
{
	if (!FileUtils::validFileAndPath(iniFile.original_name)) throw InvalidNameException("Invalid configuration file name '"+iniFile.original_name+"'",AT);
	if (!FileUtils::fileExists(iniFile.original_name)) throw NotFoundException("Configuration file '"+iniFile.original_name+"' not found", AT);

	//Open file
	std::ifstream fin(iniFile.original_name.c_str(), ifstream::in);
	if (fin.fail()) throw AccessException(iniFile.original_name, AT);
	imported.insert( iniFile ); //keep track of this file being processed to prevent circular IMPORT directives
	
	std::string section( defaultSection );
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	unsigned int linenr = 1;
	std::vector< fileProperties > import_after; //files to import after the current one
	bool accept_import_before = true;

	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			parseLine(linenr++, iniFile.restrict_section, import_after, accept_import_before, line, section);
		} while (!fin.eof());
		fin.close();
	} catch(const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}

	std::reverse(import_after.begin(), import_after.end());
	while (!import_after.empty()) {
		parseFile( import_after.back());
		import_after.pop_back();
	}
	
	imported.erase( iniFile );
}

bool ConfigParser::processSectionHeader(const std::string& line, std::string &section, const unsigned int& linenr)
{
	std::smatch section_matches;
	
	if (std::regex_match(line, section_matches, section_regex)) {
		//valid chars for section: letters, numbers, _ and -. Any number of whitespaces after the section header, potentially comments too
		static const std::regex sectionValidation_regex(R"(^\[((?:\w|\-)+)\]\s*((;|#).*)*$)", std::regex::optimize); 
		std::smatch sectionValidation_matches;
		if (!std::regex_match(line, sectionValidation_matches, sectionValidation_regex))
			throw IOException("Section header corrupt at line " + IOUtils::toString(linenr), AT);
		
		section = section_matches.str(1);
		IOUtils::toUpper(section);
		sections.insert( section );
		return true;
	}

	return false;
}

bool ConfigParser::processImports(const std::string& key, const std::string& value, std::vector< fileProperties > &import_after, const bool &accept_import_before)
{
	if (key=="IMPORT_BEFORE") {
		const fileProperties imported_file( value, sourcename );
		if (!accept_import_before)
			throw IOException("Error in \""+sourcename+"\": IMPORT_BEFORE key MUST occur before any other key!", AT);
		if (imported.count( imported_file ) != 0)
			throw IOException("IMPORT Circular dependency with \"" + value + "\"", AT);
		parseFile( imported_file );
		return true;
	}
	if (key=="IMPORT_AFTER") {
		const fileProperties imported_file( value, sourcename );
		if (imported.count( imported_file ) != 0)
			throw IOException("IMPORT Circular dependency with \"" + value + "\"", AT);
		import_after.push_back(imported_file);
		return true;
	}

	return false;
}

void ConfigParser::handleNonKeyValue(const std::string& line_backup, const std::string& section, const unsigned int& linenr, bool &accept_import_before)
{
	std::string key, value;
	if (IOUtils::readKeyValuePair(line_backup, "=", key, value, true)) {
		if (value==";" || value=="#") { //so we can accept the comments char if given only by themselves (for example, to define a CSV delimiter)
			properties[section+"::"+key] = value; //save the key/value pair
			accept_import_before = false; //this is not an import, so no further import_before allowed
			return;
		}
	}

	const std::string key_msg = (key.empty())? "" : "key "+key+" ";
	const std::string key_value_link = (key.empty() && !value.empty())? "value " : "";
	const std::string value_msg = (value.empty())? "" : value+" " ;
	const std::string keyvalue_msg = (key.empty() && value.empty())? "key/value " : key_msg+key_value_link+value_msg;
	const std::string section_msg = (section.empty())? "" : "in section "+section+" ";
	const std::string source_msg = (sourcename.empty())? "" : "from \""+sourcename+"\" at line "+IOUtils::toString(linenr);

	throw InvalidFormatException("Error reading "+keyvalue_msg+section_msg+source_msg, AT);
}

/**
* @brief Parse the given INI line
* @param[in] linenr line number in the current file for better error messages
* @param[in] restrict_section if not an empty string, only the keys from this section will be processed (useful to only import a given section from another ini file)
* @param[out] import_after vector to keep track of all the files that should be imported after the parsing of the current file
* @param accept_import_before if set to false, any IMPORT directive will trigger an error (since imports must be done before any other key is declared)
* @param[in] line The line to parse
* @param section The current section (either returned when a new section header has been encountered or used to attribute the current key to the right section)
*/
void ConfigParser::parseLine(const unsigned int& linenr, const std::string& restrict_section, std::vector< fileProperties > &import_after, bool &accept_import_before, std::string line, std::string &section)
{
	const std::string line_backup( line ); //this might be needed in some rare cases
	//First thing cut away any possible comments (may start with "#" or ";")
	IOUtils::stripComments(line);
	IOUtils::trim(line);    //delete leading and trailing whitespace characters
	if (line.empty()) return;//ignore empty lines

	//if this is a section header, read it and return
	if (processSectionHeader(line, section, linenr)) return;

	//if we should only process a given section and this is not it, return
	if (!restrict_section.empty() && section!=restrict_section) return;

	//first, we check that we don't have two '=' chars in one line (this indicates a missing newline)
	if (!onlyOneEqual(line)) {
		const std::string source_msg = (sourcename.empty())? "" : " in \""+sourcename+"\"";
		throw InvalidFormatException("Error reading line "+IOUtils::toString(linenr)+source_msg, AT);
	}
	
	IOUtils::cleanEscapedCharacters(line, std::vector<char>({' ', '=', ';'}));
	
	//this can only be a key value pair...
	std::string key, value;
	if (IOUtils::readKeyValuePair(line, "=", key, value, true)) {
		//if this is an import, process it and return
		if (processImports(key, value, import_after, accept_import_before)) return;
		
		if (value.find("${")==std::string::npos) { //normal key/value pair
			properties[section+"::"+key] = value; //save the key/value pair
		} else { //variables are handled separately
			vars[ section+"::"+key ] = parseVariable(value); //it will be parsed and expanded later
		}
		
		accept_import_before = false; //this is not an import, so no further import_before allowed
	} else {
		handleNonKeyValue(line_backup, section, linenr, accept_import_before);
	}
}

//extract the section name from a section+"::"+key value
std::string ConfigParser::extract_section(const std::string& key, const bool& provide_default)
{
	const std::string::size_type pos = key.find("::");

	if (pos != string::npos){
		const std::string sectionname( key.substr(0, pos) );
		return sectionname;
	}
	if (!provide_default)
		return std::string();
	else
		return defaultSection;
}


bool ConfigParser::expandVar(const variable& var, const std::string& section, std::string &replacement) const
{
	if (var.type==ENV) { //environment variable
		char *tmp = getenv( var.value.c_str() );
		if (tmp==nullptr) 
			throw InvalidNameException("Environment variable '"+var.value+"' declared in ini file could not be resolved", AT);
		replacement = std::string(tmp);
		return true;
		
	} else if (var.type==EXPR) { //arithmetic expression
		int status;
		const double val = te_interp(var.value.c_str(), &status);
		if (status!=0)
			throw InvalidNameException("Arithmetic expression '"+var.value+"' declared in ini file could not be evaluated", AT);
		replacement = IOUtils::toString(val);
		return true;
	} else if (var.type==REF) { //reference to another variable
		//reference to another key/value
		std::string value_Up( IOUtils::strToUpper( var.value ) );
		std::string var_section( extract_section( value_Up ) ); //extract the section name out of the variable name. If none is found, it will be empty
		//since sections' headers are declared before the variables, every extracted section should be available in our sections
		//list, so we can validate that it is indeed a proper section. if not, prepend with the section given as argument
		const bool extracted_section_exists = sections.find(var_section) != sections.end();
		if (!extracted_section_exists && !section.empty()) value_Up = section+"::"+value_Up;
		if (properties.count( value_Up )!=0) {
			replacement = properties.find(value_Up)->second;
			return true;
		}
		if (vars.count( value_Up)==0) {
			throw InvalidNameException("Reference to key '"+value_Up+"' declared in ini file does not exist", AT);
		}
		//if the key/value referred to is not ready yet (ie. not expanded), return false
		return false;
	}
	
	return false; //ROOT type does not require replacements
}

/** 
 * @brief Tries to expand all the variables that are contained in the value string of a key/value
 * @details
 * The expansion of the variables is attempted in reverse order (starting from the end of the root string). When all the
 * children of a particular node have been expanded, the value of the node itself will also be expanded again in reverse
 * order. This guarantees that the positions in the string remain valid and don't have to be recomputed.
 * 
 * Two flags are returned: "hasSomeSuccesses" tells if at least one variable could be expanded (that had not been expanded before)
 * and the return value tells if all the variables in the root string could be replaced in the root string.
 * @param[in] vecVars all the variables (in a tree) that are present in the root string (including the root string at position zero)
 * @param[in] section section of this key/value pair (for error messages)
 * @param[in] hasSomeSuccesses set to true if at least one new variable vould be expanded
 * @return true if all variables could be expanded and replaced in the root string
 */
bool ConfigParser::expandVarsForKey(std::vector<variable>& vecVars, const std::string& section, bool& hasSomeSuccesses)
{
	bool allSuccess = true;
	
	//we go through the vector representation of the tree in reverse order so we start with the children...
	//for (size_t ii=(vecVars.size()-1); ii>0; ii--) { //at 0, there is the root expression that should not be called here
	for (size_t ii=vecVars.size(); ii-- > 0; ) {
		variable& var = vecVars[ii];
		if (var.isExpanded) continue; //this variable has already been processed before
		if (var.children_not_ready!=0) continue; //it still has at least one un-expanded child
				
		//for nodes having children that are all ready, perform the replacements, starting from the last so the positions remain correct
		if (!var.children_idx.empty()) {
			for (std::set<size_t>::const_reverse_iterator rit = var.children_idx.rbegin(); rit != var.children_idx.rend(); ++rit) {
				const size_t child_idx = *rit;
				size_t pos_offset = 0;
				if (var.type!=ROOT) {
					//parent position contains the variable markers which are missing in its "value" member
					const size_t marker_offset = (var.type==EXPR)? 2 : (var.type==ENV)? 5 : 1; 
					pos_offset = var.pos_start + marker_offset + 1;
				}
				const size_t relative_start_idx = vecVars[child_idx].pos_start - pos_offset;
				const size_t length = vecVars[child_idx].pos_end - vecVars[child_idx].pos_start + 1;
				var.value.replace(relative_start_idx, length, vecVars[child_idx].value);
			}
			var.children_idx.clear();
		}
		
		if (var.type==ROOT) return allSuccess; //in ROOT, we perform replacements but no expansions
		
		//attempt to do variable expansion
		std::string replacement;
		const bool status = expandVar(var, section, replacement);
		if (status) {
			hasSomeSuccesses = true; //at least one variable could be expanded
			var.value = replacement;
			vecVars[ var.parent_idx ].children_not_ready--;
			var.isExpanded = true;
		} else {
			allSuccess = false;
		}
	}
	
	return allSuccess;
}

/** 
 * @brief Parse a value string to extract all variables that will require expansion
 * @details
 * Variables can be either environment variables ${env:XXX}, airthmetic expressions ${{XXX}} or references to other
 * keys ${XXX}. Recursion is supported, such as ${{ ${SAMPLING_RATE}*60 + ${env:${MYOFFSET}_${MYVAR}} }}.
 * The parsed variables are returned in the binary tree representation of a m-ary tree: each node may 
 * have a parent, and multiple children. 
 * 
 * As a side note, having a different marker between env / vars and arithmetic expressions makes 
 * parsing much more unpleasant...
 * @param[in] value string to parse for variables
 * @return vector of variables, organized as binary tree
 */
std::vector<ConfigParser::variable> ConfigParser::parseVariable(const std::string &value)
{
	const size_t max_len = value.length();
	size_t start_pos = 0;
	std::stack<size_t> variables_stack; //keep the indices of parents during recursion
	std::vector<variable> variables_tree; //tree of the whole set of expressions in the value string
	
	//the root node contains the original value string associated with the key
	variable tmp_root(value, 0, value.size(), ROOT);
	variables_tree.push_back( tmp_root );
	variables_stack.push( 0 );
	
	do {
		const size_t pos_open = value.find("${", start_pos);
		size_t pos_close = value.find("}", start_pos);	//for arithmetic expressions, we will increment it to contain the last '}'
		const bool hasOpen = (pos_open!=std::string::npos);
		const bool hasClose = (pos_close!=std::string::npos);
		
		if (!hasClose) {
			if (hasOpen) { //we shall not open a new expression if there won't be any closing
				throw InvalidFormatException("Invalid configuration key '"+value+"': the expressions will not be properly closed", AT);
			} else { //no open, no close, there shall be no remaining expressions in the stack
				if (variables_stack.size()>1) throw InvalidFormatException("Invalid configuration key '"+value + "' starting at position "+IOUtils::toString(variables_tree[variables_stack.top()].pos_start)+": remaining expressions are not properly closed", AT);
				break; //all good: no open, no close, no stack
			}
		}
		
		//closing one remaining expression from the stack
		if (!hasOpen || pos_open>pos_close) {
			if (variables_stack.size()<2) throw InvalidFormatException("Invalid configuration key '"+value+"': closing an expression that has not been opened", AT);
			const size_t parent_idx = variables_stack.top();
			if (variables_tree[parent_idx].type==EXPR) {
				const bool closeArithmExpr = (pos_close+1 < max_len && value[pos_close+1]=='}')? true : false;
				if (!closeArithmExpr) throw InvalidFormatException("Invalid configuration key '"+value+"': not closing the current arithmetic expression", AT);
				pos_close++; //to capture the additional '}'
			}
			
			const size_t pos_start = variables_tree[parent_idx].pos_start;
			const size_t pos_end = (variables_tree[parent_idx].pos_end!=IOUtils::npos) ? variables_tree[parent_idx].pos_end : pos_close; //it contains the additional '}'
			const size_t content_begin = (variables_tree[parent_idx].type==EXPR ? pos_start+3 : (variables_tree[parent_idx].type==ENV ? pos_start+6 : pos_start+2));
			const size_t content_end = (variables_tree[parent_idx].type==EXPR ? pos_end-2 : pos_end-1);
			const std::string parsedVar( value.substr(content_begin, content_end-content_begin+1) );
			variables_tree[parent_idx].finalize(parsedVar, pos_end);
			start_pos = pos_end+1;
			
			variables_stack.pop();
			continue;
		}
	
		// from now one, pos_open<pos_close, starting a new expression
		// identify if this is an arithmetic expression or Not
		const bool isEnv = (pos_open+5 < max_len && value.substr(pos_open, 6)=="${env:")? true : false;
		const bool isArithmExpr = (pos_open+2 < max_len && value[pos_open+2]=='{')? true : false;
		const size_t pos_next_open = value.find("${", pos_open+2);
		
		// easy case: opening and closing, no recursion
		if (pos_next_open==std::string::npos || (pos_next_open>pos_close)) {
			// check that the closing matches the opening
			if (isArithmExpr) {
				const bool closingArithmExpr = (pos_close+1 < max_len && value[pos_close+1]=='}')? true : false;
				if (!closingArithmExpr) throw InvalidFormatException("Invalid configuration key '"+value+"': the arithmetic expression is not properly closed", AT);
				pos_close++; //to capture the additional '}'
			}
			const size_t content_begin = (isArithmExpr ? pos_open+3 : (isEnv ? pos_open+6 : pos_open+2));
			const size_t content_end = (isArithmExpr ? pos_close-2 : pos_close-1);
			const std::string parsedVar( value.substr(content_begin, content_end-content_begin+1) );
			const VarType type = (isEnv)? ENV : (isArithmExpr)? EXPR : REF;
			variable tmp_var(parsedVar, pos_open, pos_close, type);
			
			// set relationships
			if (variables_stack.size()>1) {
				const size_t parent_idx = variables_stack.top();
				tmp_var.setParent( parent_idx );
				variables_tree[parent_idx].setChild( variables_tree.size() ); //add this child
			} else {
				std::size_t parent_idx = 0;
				tmp_var.setParent( parent_idx );
				variables_tree[parent_idx].setChild( variables_tree.size() ); //this will be its own index as soon as the push is done
			}
			
			variables_tree.push_back( tmp_var );
			start_pos = pos_close+1;
		} else { //recursion...
			const VarType type = (isEnv)? ENV : (isArithmExpr)? EXPR : REF;
			variable tmp_var(pos_open, type);
			std::size_t parent_idx = variables_stack.top();
			tmp_var.setParent( parent_idx );
			variables_tree[parent_idx].setChild( variables_tree.size() ); //this will be its own index as soon as the push is done
			
			variables_tree.emplace_back( tmp_var );
			variables_stack.push( variables_tree.size()-1 ); //keep track of the parent of upcoming variables
			start_pos = (isArithmExpr)? pos_open+3 : pos_open+2;
		}
	} while (true);
	
	return variables_tree;
}


} //end namespace
