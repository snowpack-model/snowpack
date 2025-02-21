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
#ifndef JSONWRAPPER_H
#define JSONWRAPPER_H

#include <meteoio/IOInterface.h>
#include <meteoio/thirdParty/picojson.h>

#include <curl/curl.h>
#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

/**
 * @class JsonWrapper
 * @brief This is a wrapper class around picoJson and Curl.
 *
 * This class relies on an embedded version of <A HREF="https://github.com/kazuho/picojson/">picojson</A> for reading and parsing
 * <A HREF="https://en.wikipedia.org/wiki/JSON">JSON</A> data. Picojson is released under a
 * <A HREF="https://opensource.org/licenses/BSD-2-Clause">2-Clause BSD License</A>. Please find here below
 * the full license agreement for picojson:
 *
 * @code
 * Copyright 2009-2010 Cybozu Labs, Inc.
 * Copyright 2011-2014 Kazuho Oku
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * @endcode
 * @ingroup plugins
 * @date   2021-11-23
 */

class JsonWrapper {
	public:
		JsonWrapper();
		JsonWrapper( const JsonWrapper& );
		~JsonWrapper();
		JsonWrapper& operator=(const mio::JsonWrapper&);

		void setConnectionParams(const std::string& i_proxy_url, const int& i_http_timeout, const bool& i_debug=false);
		void readAndParse(const std::string& request, const std::string& where);
		
		static void printJSON(const picojson::value& v, const unsigned int& depth);
		std::vector<picojson::value> JSONQuery(const std::string& path) const;
		picojson::value goToJSONPath(const std::string& path) const;
		
		std::string getString(const std::string& path) const;
		std::vector<std::string> getStrings(const std::string& path) const;
		double getDouble(const std::string& path) const;
		std::vector<double> getDoubles(const std::string& path) const;
		
	private:
		static std::string indent(const unsigned int& depth);
		static void JSONQuery(const std::string& path, const picojson::value& v, std::vector<picojson::value>& results);
		static picojson::value goToJSONPath(const std::string& path, const picojson::value& v);
		static size_t data_write(void* buf, const size_t size, const size_t nmemb, void* userp);
		bool curl_read(const std::string& url, std::ostream& os);
		
		picojson::value json_tree;
		CURL *curl;
		std::string proxy_url;
		int http_timeout; //time out for http connections
		bool debug, proxy;
		static const int http_timeout_dflt;
};

} //end namespace mio

#endif
