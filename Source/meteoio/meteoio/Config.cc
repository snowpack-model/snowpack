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
//arithmetic expression syntax, delimited by ${{ and }} and contains anything in between including sub-variables
static const std::regex expr_regex(R"(\$\{\{(.+?)\}\})", std::regex::optimize);
//environment variable syntax, delimited by ${env: and } and not sub-variables / expressions
static const std::regex envvar_regex(R"(\$\{env:([^\$\{\}]+?)\})", std::regex::optimize);
//variable refering to another key, delimited by ${ and } and can contain sub-variables / expressions
static const std::regex var_regex(R"(\$\{((?!\{).*?)\}(?!\}))", std::regex::optimize);
static const std::regex sub_expr_regex(R"(.*\$\{((?!\{).*?)\}(?!\}).*)", std::regex::optimize); //same as above but allowing match anywhere in the string
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
ConfigParser::ConfigParser(const std::string& filename, std::map<std::string, std::string> &i_properties, std::set<std::string> &i_sections) : properties(i_properties), imported(), sections(), deferred_vars(), deferred_exprs(), sourcename(filename)
{
	parseFile(filename);

	//process potential deferred vars and exprs, at least one must be successfuly replaced / evaluated at each iteration, 
	//otherwise we consider that there is a circular dependency
	size_t deferred_count = deferred_vars.size() + deferred_exprs.size();
	while (deferred_count>0) {
		//try to perform variable expansion
		for (std::set<std::string>::iterator it = deferred_vars.begin(); it!=deferred_vars.end(); ) {
			const std::string section( extract_section( *it, true ) );
			std::string value( properties[ *it ] );
			const bool status = processVars(value, section);
			
			properties[ *it ] = value; //save the key/value pair
			if (status) //the variable could be fully expanded
				deferred_vars.erase( it++ );
			else
				++it;
		}

		//try to arithmetic expression evaluation
		for (std::set<std::string>::iterator it = deferred_exprs.begin(); it!=deferred_exprs.end(); ) {
			std::string value( properties[ *it ] );
			const bool status = processExpr(value);

			properties[ *it ] = value; //save the key/value pair
			if (status) //the variable could be fully expanded
				deferred_exprs.erase( it++ );
			else
				++it;
		}

		const size_t new_deferred_count = deferred_vars.size() + deferred_exprs.size();
		if (new_deferred_count==deferred_count) {
			std::string msg("In file "+filename+", the following keys could not be resolved (circular dependency? invalid variable name? syntax error?):");
			for (const std::string& var : deferred_vars) 
				msg.append( " "+var );
			throw InvalidArgumentException(msg, AT);
		}
		deferred_count = new_deferred_count;
	}
	
	//swap the caller and local properties before returning, so the caller gets the new version
	std::swap(properties, i_properties);
	i_sections.insert(sections.begin(), sections.end());
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
void ConfigParser::parseFile(const std::string& filename)
{
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException("Invalid configuration file name '"+filename+"'",AT);
	if (!FileUtils::fileExists(filename)) throw NotFoundException("Configuration file '"+filename+"' not found", AT);

	//Open file
	std::ifstream fin(filename.c_str(), ifstream::in);
	if (fin.fail()) throw AccessException(filename, AT);
	imported.insert( FileUtils::cleanPath(filename, true) ); //keep track of this file being processed to prevent circular IMPORT directives
	
	std::string section( defaultSection );
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	unsigned int linenr = 1;
	std::vector<std::string> import_after; //files to import after the current one
	bool accept_import_before = true;

	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			parseLine(linenr++, import_after, accept_import_before, line, section);
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
	
	imported.erase( filename );
}

bool ConfigParser::processSectionHeader(const std::string& line, std::string &section, const unsigned int& linenr)
{
	std::smatch section_matches;
	
	if (std::regex_match(line, section_matches, section_regex)) {
		static const std::regex sectionValidation_regex(R"(^\[((?:\w|\-)+)\]\s*((;|#).*)*$)", std::regex::optimize); //valid chars for section: letters, numbers, _ and -. Any number of whitespaces after the section header, potentially comments too
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

/**
* @brief Process keys that refer to environment variables
* @details Environment variables are used as key = ${env:envvar}.
*/
void ConfigParser::processEnvVars(std::string& value)
{
	std::smatch matches;
	
	std::string::const_iterator searchStart( value.cbegin() );
    while( regex_search( searchStart, value.cend(), matches, envvar_regex ) ) { //there might be several replacements to perform in a single value string
		const std::string envVar( matches[1].str() );
		if (envVar.empty())
			throw InvalidFormatException("Wrong syntax for environment variable: '"+envVar+"'", AT);
		
		char *tmp = getenv( envVar.c_str() );
		if (tmp==nullptr) 
			throw InvalidNameException("Environment variable '"+envVar+"' declared in ini file could not be resolved", AT);
		
		const size_t match_pos = std::distance(static_cast<std::string::const_iterator>(value.begin()), searchStart) + matches.position();
		value.replace(match_pos, matches.length(), std::string(tmp));
		//we don't increment searchStart as after replacement, we have to search again from the same point
	}
	
}

/**
* @brief Process keys that are arithmetic expressions
* @details Variables are used as key = ${{expr}} and can refere to keys that will be defined later.
* @return false if the key has to be processed later (as in the case of part of the expression refering to a key
* that has not yet been read)
*/
bool ConfigParser::processExpr(std::string& value)
{
	//process env. variables
	static const std::regex var_ref_regex(R"(.*\$\{((?!\{).*?)\}(?!\}).*)", std::regex::optimize); //to be able to catch sub-expressions referring to variables
	std::smatch matches;
	
	//the replacements won't be done in place but by building a separate string for efficiency
	bool var_fully_parsed = true;
	
	std::string::const_iterator searchStart( value.cbegin() );
    while( regex_search( searchStart, value.cend(), matches, expr_regex ) ) { //there might be several replacements to perform in a single value string
		const std::string expression( matches[1].str() );
		if (expression.empty())
			throw InvalidFormatException("Wrong syntax for arithmetic expression: '"+expression+"'", AT);

		if (std::regex_match(expression, var_ref_regex, std::regex_constants::match_any)) {
				//the expression contains a ${VAR}, so it still requires further replacement before evaluation
				var_fully_parsed = false;
				searchStart = matches.suffix().first;
				continue;
		}
		int status_code;
		const double val = te_interp(expression.c_str(), &status_code);
		if (status_code!=0)
			throw InvalidNameException("Arithmetic expression '"+expression+"' declared in ini file could not be evaluated", AT);
		
		
		const size_t match_pos = std::distance(static_cast<std::string::const_iterator>(value.begin()), searchStart) + matches.position();
		value.replace(match_pos, matches.length(), IOUtils::toString(val));
		//we don't increment searchStart as after replacement, we have to search again from the same point
	}
	
	if (!var_fully_parsed) return false;
	return true;
}

/**
* @brief Process keys that are standard variables
* @details Variables are used as key = ${var} and can refere to keys that will be defined later.
* @return false if the key has to be processed later (as in the case of standard variables refering to a key
* that has not yet been read)
*/
bool ConfigParser::processVars(std::string& value, const std::string& section) const
{
	//process env. variables
	std::smatch matches;
	bool var_fully_parsed = true;

	std::string::const_iterator searchStart( value.cbegin() );
    while( std::regex_search( searchStart, value.cend(), matches, var_regex ) ) { //there might be several replacements to perform in a single value string
		std::string var( matches[1].str() );
		if (var.empty()) throw InvalidFormatException("Wrong syntax for variable: '"+var+"'", AT);
		IOUtils::toUpper( var );
	
		const std::string var_section( extract_section( var ) ); //extract the section name out of the variable name. If none is found, it will be empty
		//since sections' headers are declared before the variables, every extracted section should be available in our sections
		//list, so we can validate that it is indeed a proper section. if not, prepend with the section given as argument
		const bool extracted_section_exists = sections.find(var_section) != sections.end();
		if (!extracted_section_exists && !section.empty()) var = section+"::"+var;
		if (properties.count( var )!=0) {
			const std::string replacement( properties.find(var)->second );
			if (std::regex_match(replacement, sub_expr_regex, std::regex_constants::match_any)) {
				//the replacement contains a ${VAR}, so it still requires further replacement
				var_fully_parsed = false;
			}
			
			const size_t match_pos = std::distance(static_cast<std::string::const_iterator>(value.begin()), searchStart) + matches.position();
			value.replace(match_pos, matches.length(), replacement);
			//we don't increment searchStart as after replacement, we have to search again from the same point
		} else {
			var_fully_parsed = false;
			searchStart = matches.suffix().first; //so we can look for the next match
		}
	}
	
	if (!var_fully_parsed) return false;
	return true;
}

//resolve symlinks, resolve relative path w/r to the path of the current ini file
std::string ConfigParser::clean_import_path(const std::string& in_path) const
{
	//if this is a relative path, prefix the import path with the current path
	const std::string prefix = ( FileUtils::isAbsolutePath(in_path) )? "" : FileUtils::getPath(sourcename, true)+"/";
	const std::string path( FileUtils::getPath(prefix+in_path, true) );  //clean & resolve path
	const std::string filename( FileUtils::getFilename(in_path) );

	return path + "/" + filename;
}

bool ConfigParser::processImports(const std::string& key, const std::string& value, std::vector<std::string> &import_after, const bool &accept_import_before)
{
	if (key=="IMPORT_BEFORE") {
		const std::string file_and_path( clean_import_path(value) );
		if (!accept_import_before)
			throw IOException("Error in \""+sourcename+"\": IMPORT_BEFORE key MUST occur before any other key!", AT);
		if (imported.count( file_and_path ) != 0)
			throw IOException("IMPORT Circular dependency with \"" + value + "\"", AT);
		parseFile( file_and_path );
		return true;
	}
	if (key=="IMPORT_AFTER") {
		const std::string file_and_path( clean_import_path(value) );
		if (imported.count( file_and_path ) != 0)
			throw IOException("IMPORT Circular dependency with \"" + value + "\"", AT);
		import_after.push_back(file_and_path);
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

void ConfigParser::parseLine(const unsigned int& linenr, std::vector<std::string> &import_after, bool &accept_import_before, std::string &line, std::string &section)
{
	const std::string line_backup( line ); //this might be needed in some rare cases
	//First thing cut away any possible comments (may start with "#" or ";")
	IOUtils::stripComments(line);
	IOUtils::trim(line);    //delete leading and trailing whitespace characters
	if (line.empty()) return;//ignore empty lines

	//if this is a section header, read it and return
	if (processSectionHeader(line, section, linenr)) return;

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

		if (value.find("}}}")!=std::string::npos) throw InvalidFormatException("Invalid configuration key '"+section+"::"+key+"': please use a space to separate the variables from the arithmetic expression delimiters, like in ${{ 10*${SAMPLING_RATE_MIN} }} ", AT);
		
		processEnvVars( value );
		if (!processExpr( value )) deferred_exprs.insert( section+"::"+key );
		if (!processVars(value, section)) deferred_vars.insert( section+"::"+key );
		if (key.find(' ')!=std::string::npos) throw InvalidFormatException("Invalid configuration key '"+section+"::"+key+"': keys can not contain spaces", AT);
		properties[section+"::"+key] = value; //save the key/value pair
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
		//key.erase(key.begin(), key.begin() + pos + 2); //delete section name
		return sectionname;
	}
	if (!provide_default)
		return std::string();
	else
		return defaultSection;
}

} //end namespace
