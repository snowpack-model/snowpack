// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef LIBMYSQLWRAPPER_H
#define LIBMYSQLWRAPPER_H

#ifdef _WIN32
	#include <winsock.h>
#endif // _WIN32

#include <mysql.h>
#include <cstring>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/IOUtils.h>

static const int MYSQL_STRING_SIZE = 50;

class SQL_FIELD {
	public:
		#if LIBMYSQL_VERSION_ID > 80001
			typedef bool BOOL_TYPE;
		#else
			typedef my_bool BOOL_TYPE;
		#endif
		
		enum unitsConversions {
			NONE=0,
			C_TO_K,
			CM_TO_M,
			NORMALIZE_PC,
			HPA_TO_PA
		};
	
		SQL_FIELD(const std::string& i_param, const enum_field_types &type, const unsigned int &i_processing=0);
		SQL_FIELD(const std::string& i_str, const std::string& i_param="", const unsigned int &i_processing=0);
		SQL_FIELD(const mio::Date& i_dt, const std::string& i_param="", const unsigned int &i_processing=0);
		
		void resetDate();
		void reset();
		void setFromDate(const mio::Date& i_dt, MYSQL_TIME &ts);
		void setString(const std::string& i_str);
		void setDate(const mio::Date& i_dt);
		void setDouble(const double& i_val);
		mio::Date getDate(const double& TZ) const;
		
		//several members could be const, but we need a working '=' operator in the plugin for easier code...
		std::string param;               ///< the parameter name for MeteoIO
		char str[MYSQL_STRING_SIZE];     ///< for MySQL to store string data
		MYSQL_TIME dt;                   ///< for MySQL to store datetime data
		unsigned long int str_len;       ///< for MySQL, length of a string
		unsigned long int buffer_len;    ///< for MySQL, allocated data buffer length
		double val;                      ///< for MySQL to store double data
		unsigned int processing;         ///< what kind of corrections to apply to the raw data (such as C to K conversion, etc
		BOOL_TYPE is_null, error;        ///< for MySQL ro report a NULL value or errors
		bool isDate;                     ///< for MeteoIO to quickly identify datetime fields
		enum_field_types MysqlType;      ///< for MySQL, data type of the field (see mysql/field_types.h)
};

namespace mysql_wrp {
	enum MysqlOptions {
		COMPRESSION	= 1,
		ENCRYPTION	= 1 << 1
	};
	
	MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb, const unsigned int& options=0);
	MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count);
	
	void bindParams(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &params_fields);
	void bindResults(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &result_fields);
	double retrieveData(const SQL_FIELD &field, const unsigned int& conversion=SQL_FIELD::NONE);
}

#endif
