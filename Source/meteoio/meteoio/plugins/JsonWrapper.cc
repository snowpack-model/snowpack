// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 SLF                                                                                                                                */
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
#include <meteoio/plugins/JsonWrapper.h>
#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <sstream>
#include <iostream>
#include <cstring>

using namespace std;

namespace mio {

const int JsonWrapper::http_timeout_dflt = 60; // seconds until connect time out for libcurl

////////////////////////////////////// Wrapper calls around libcurl
JsonWrapper::JsonWrapper()
            : json_tree(), curl(nullptr), proxy_url(), http_timeout(http_timeout_dflt), debug(false), proxy(false)
{
	curl_global_init(CURL_GLOBAL_ALL); //NOT thread-safe
}

JsonWrapper::JsonWrapper( const JsonWrapper& c)
            : json_tree(c.json_tree), curl(nullptr), proxy_url( c.proxy_url ), http_timeout( c.http_timeout ), debug( c.debug ), proxy( c.proxy ) {}

JsonWrapper& JsonWrapper::operator=(const mio::JsonWrapper& c)
{
	if (this != &c) {
		json_tree = c.json_tree;
		curl = nullptr;
		proxy_url = c.proxy_url;
		http_timeout = c.http_timeout;
		debug = c.debug;
		proxy = c.proxy;
	}
	
	return *this;
}

JsonWrapper::~JsonWrapper()
{
	if (curl != nullptr) curl_easy_cleanup(curl);
}

/**
* @brief Set the connection parameters for CURL
* @param[in] i_proxy_url optional SOCKS5 proxy URL
* @param[in] i_http_timeout timeout in seconds for the connections
* @param[in] i_debug enable debug outputs?
*/
void JsonWrapper::setConnectionParams(const std::string& i_proxy_url, const int& i_http_timeout, const bool& i_debug)
{
	proxy_url = i_proxy_url;		//an optional SOCKS5 proxy for the connection
	http_timeout = i_http_timeout;
	debug = i_debug;
	proxy = !proxy_url.empty();
}

////////////////////////////////////// Wrapper calls around picoJson
std::string JsonWrapper::indent(const unsigned int& depth)
{
	std::ostringstream ss;
	for (unsigned int jj=0; jj<depth; jj++) 
		ss << "\t";
	return ss.str();
}

void JsonWrapper::printJSON(const picojson::value& v, const unsigned int& depth)
{
	if (v.is<picojson::null>()) {
		std::cout << indent(depth) << "NULL\n";
		return;
	}

	if (v.is<picojson::object>()) {
		//NOTE: v.get() will only be called once, see
		//https://stackoverflow.com/questions/15766020/does-a-c11-range-based-for-loop-condition-get-evaluated-every-cycle
		for (const auto& it : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
			std::cout << indent(depth) << it.first << "\n";
			printJSON(it.second, depth+1);
		}
	} else if (v.is<std::string>()){
		std::cout << indent(depth) << v.get<std::string>() << "\n";
	} else if (v.is<double>()){
		std::cout << indent(depth) << v.get<double>() << "\n";
	} else if (v.is<bool>()){
		std::cout << indent(depth) << std::boolalpha << v.get<bool>() << "\n";
	} else if (v.is<picojson::array>()){ //ie vector<picojson::value>
		const auto& array = v.get<picojson::array>();
		std::cout << indent(depth) << "array " << array.size() << "\n";
		for (const auto& vec_elem : array)
			printJSON(vec_elem, depth+1);
	}
}

picojson::value JsonWrapper::goToJSONPath(const std::string& path) const
{
	return goToJSONPath(path, json_tree);
}

picojson::value JsonWrapper::goToJSONPath(const std::string& path, const picojson::value& v)
{
	if (v.is<picojson::null>()) return picojson::value();
	if (!v.is<picojson::object>()) return picojson::value();

	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	for (auto& keyvalue : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
		if (keyvalue.first==local_path) {
			if (!remaining_path.empty())
				goToJSONPath(remaining_path, keyvalue.second);
			else
				return keyvalue.second;
		}
	}

	return picojson::value();
}

std::vector<picojson::value> JsonWrapper::JSONQuery(const std::string& path) const
{
	std::vector<picojson::value> results;
	JSONQuery(path, json_tree, results);
	return results;
}

void JsonWrapper::JSONQuery(const std::string& path, const picojson::value& v, std::vector<picojson::value>& results)
{
	if (v.is<picojson::null>()) return;
	if (!v.is<picojson::object>()) return;

	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	for (const auto& keyvalue : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
		if (keyvalue.first==local_path) {
			if (!remaining_path.empty()) {
				if (keyvalue.second.is<picojson::array>()) { //ie vector<picojson::value>
					for (const auto& vec_elem : keyvalue.second.get<picojson::array>())
						JSONQuery(remaining_path, vec_elem, results);
				} else
					JSONQuery(remaining_path, keyvalue.second, results);
			} else {
				results.push_back( keyvalue.second );
			}
		}
	}
}

std::string JsonWrapper::getString(const std::string& path) const
{
	const std::vector<picojson::value> results( JSONQuery(path) );
	
	if (!results.empty()) {
		if (! results.front().is<picojson::null>() && results.front().is<std::string>()) 
			return results.front().get<std::string>();
	}

	return std::string();
}

std::vector<std::string> JsonWrapper::getStrings(const std::string& path) const
{
	const std::vector<picojson::value> results( JSONQuery(path) );

	std::vector<std::string> vecString;
	for (const auto& result : results) {
		 if (result.is<picojson::array>()) {
			for (const auto& vec_elem : result.get<picojson::array>()) { //loop over vector<picojson::value>
				if (!vec_elem.is<picojson::null>() && vec_elem.is<std::string>())
					vecString.push_back( vec_elem.get<std::string>() );
			}
		} else
			if (!result.is<picojson::null>() && result.is<std::string>())
				vecString.push_back( result.get<std::string>() );
	}

	return vecString;
}

double JsonWrapper::getDouble(const std::string& path) const
{
	const std::vector<picojson::value> results( JSONQuery(path) );
	
	if (!results.empty()) {
		if (!results.front().is<picojson::null>() && results.front().is<double>())
			return results.front().get<double>();
	}

	return IOUtils::nodata;
}

std::vector<double> JsonWrapper::getDoubles(const std::string& path) const
{
	const std::vector<picojson::value> results( JSONQuery(path) );

	std::vector<double> vecDouble;
	for (const auto& result : results) {
		 if (result.is<picojson::array>()) {
			for (const auto& vec_elem : result.get<picojson::array>()) {//loop over vector<picojson::value>
				if (!vec_elem.is<picojson::null>() && vec_elem.is<double>())
					vecDouble.push_back( vec_elem.get<double>() );
			}
		} else
			if (!result.is<picojson::null>() && result.is<double>())
				vecDouble.push_back( result.get<double>() );
	}

	return vecDouble;
}

void JsonWrapper::readAndParse(const std::string& request, const std::string& where)
{
	std::stringstream ss;
	const bool read_status = curl_read(request, ss);
	if (!read_status) {
		if (debug) std::cout << "****\nRequest: " << request << "\n****\n";
		throw IOException(where+", could not request data from server", AT);
	}
	
	if (ss.str().empty()) throw IOException(where+", no data returned from server", AT);
	
	const std::string err( picojson::parse(json_tree, ss) );
	if (!err.empty()) throw IOException(where+", error while parsing JSON: "+err, AT);
}

////////////////////////////////////// Wrapper calls around libcurl
size_t JsonWrapper::data_write(void* buf, const size_t size, const size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const std::streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool JsonWrapper::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	if (curl == nullptr) curl = curl_easy_init(); //calls curl_global_init() if not already done, so it is NOT thread-safe
	if (!curl) {
		if (debug)
			std::cout << "****\nRequest: " << url_query << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}
	
	if (proxy && !(
		CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_PROXY, proxy_url.c_str()))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_PROXYTYPE, CURLPROXY_SOCKS5_HOSTNAME)) ))
	{
		if (debug)
			std::cout << "****\nRequest: " << url_query << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	if (CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, http_timeout))
		&& CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url_query.c_str())) )
	{
		code = curl_easy_perform(curl);
	}

	if (code!=CURLE_OK) {
		if (debug)
			std::cout << "****\nRequest: " << url_query << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	return (code==CURLE_OK);
}

} //namespace
