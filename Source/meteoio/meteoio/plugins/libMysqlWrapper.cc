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
#include <meteoio/plugins/libMysqlWrapper.h>
#include <meteoio/IOExceptions.h>

#include <stdio.h>

using namespace std;
using namespace mio;

SQL_FIELD::SQL_FIELD(const std::string& i_param, const enum_field_types &type, const unsigned int &i_processing) 
          : param(i_param), str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), processing( i_processing ), is_null(false), error(false), isDate(param=="DATETIME"), MysqlType(type) 
{}

SQL_FIELD::SQL_FIELD(const std::string& i_str, const std::string& i_param, const unsigned int &i_processing) 
          : param(i_param), str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), processing( i_processing ), is_null(false), error(false), isDate(param=="DATETIME"), MysqlType(MYSQL_TYPE_STRING) 
{
	setString(i_str);
}

SQL_FIELD::SQL_FIELD(const mio::Date& i_dt, const std::string& i_param, const unsigned int &i_processing) 
          : param(i_param), str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), processing( i_processing ), is_null(false), error(false), isDate(param=="DATETIME"), MysqlType(MYSQL_TYPE_DATETIME) 
{
	setFromDate(i_dt, dt);
}

void SQL_FIELD::resetDate() 
{
	dt.year = 0; 
	dt.month = 0; 
	dt.day = 0; 
	dt.hour = 0; 
	dt.minute = 0; 
	dt.second = 0; 
	dt.second_part = 0;
}

void SQL_FIELD::reset() 
{
	str[0] = '\0'; 
	resetDate(); 
	str_len = 0; 
	val = mio::IOUtils::nodata; 
	MysqlType=MYSQL_TYPE_NULL;
}

void SQL_FIELD::setFromDate(const mio::Date& i_dt, MYSQL_TIME &ts) 
{ 
	int year, month, day, hour, minute; 
	double second;
	i_dt.getDate(year, month, day, hour, minute, second);
	
	ts.year = static_cast<unsigned int>( year );
	ts.month = static_cast<unsigned int>( month );
	ts.day = static_cast<unsigned int>( day );
	ts.hour = static_cast<unsigned int>( hour );
	ts.minute = static_cast<unsigned int>( minute );
	ts.second = static_cast<unsigned int>( floor(second) );
	ts.second_part = static_cast<unsigned long int>( floor((second - floor(second))*1e6) );
}

void SQL_FIELD::setString(const std::string& i_str) 
{
	reset(); 
	strncpy(str, i_str.c_str(), std::min(static_cast<int>(i_str.size()), MYSQL_STRING_SIZE)); 
	str_len = strlen(str); 
	MysqlType = MYSQL_TYPE_STRING;
}

void SQL_FIELD::setDate(const mio::Date& i_dt) 
{
	reset(); 
	setFromDate(i_dt, dt); 
	MysqlType = MYSQL_TYPE_DATETIME;
}

void SQL_FIELD::setDouble(const double& i_val) 
{
	reset(); 
	val = i_val; 
	MysqlType = MYSQL_TYPE_DOUBLE;
}

mio::Date SQL_FIELD::getDate(const double& TZ) const 
{
	const int year = static_cast<int>(dt.year);
	const int month = static_cast<int>(dt.month);
	const int day = static_cast<int>(dt.day);
	const int hour = static_cast<int>(dt.hour);
	const int minute = static_cast<int>(dt.minute);
	const double seconds = static_cast<double>(dt.second) + static_cast<double>(dt.second_part)*1e-6;
	
	mio::Date o_dt(year, month, day, hour, minute, seconds, TZ); 
	return o_dt;
}

namespace mysql_wrp {

MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb, const unsigned int& options)
{
	MYSQL *mysql = mysql_init(nullptr);
	
	//set some options
	const unsigned int timeout = 2; // in seconds
	mysql_options(mysql, MYSQL_OPT_CONNECT_TIMEOUT, &timeout);
	if ((COMPRESSION & options) == COMPRESSION)
		mysql_options(mysql, MYSQL_OPT_COMPRESS, 0);
	if ((ENCRYPTION & options) == ENCRYPTION) {
		const unsigned int enforce_ssl = SSL_MODE_REQUIRED;
		mysql_options(mysql, MYSQL_OPT_SSL_MODE, &enforce_ssl);
	}
	
	if (!mysql_real_connect(mysql, mysqlhost.c_str(), mysqluser.c_str(), mysqlpass.c_str(), mysqldb.c_str(), 0, NULL, 0))
		throw AccessException("Could not initiate connection to Mysql server "+mysqlhost+": "+std::string(mysql_error(mysql)), AT);

	return mysql;
}

MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count)
{
	MYSQL_STMT* stmt = mysql_stmt_init(*mysql);
	if (!stmt) throw IOException("Could not allocate memory for mysql statement", AT);

	if (mysql_stmt_prepare(stmt, query.c_str(), query.size())) {
		throw IOException("Error preparing mysql statement, please check for syntax and typos in field names", AT);
	} else {
		const long unsigned int param_count = mysql_stmt_param_count(stmt);
		if (param_count!=ref_param_count) throw InvalidArgumentException("Wrong number of parameters in mysql statement", AT);
	}
	
	return stmt;
}

void bindParams(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &params_fields)
{
	const size_t params_count = params_fields.size();
	MYSQL_BIND *stmtParams = (MYSQL_BIND*)calloc(params_count, sizeof(MYSQL_BIND));
	if (stmtParams==nullptr) throw IOException("Could not allocate memory for parameter binding to Mysql query", AT);
	
	for (size_t ii=0; ii<params_count; ++ii) {
		stmtParams[ii].buffer_type = params_fields[ii].MysqlType;
		stmtParams[ii].is_null = nullptr;
		
		if (params_fields[ii].MysqlType==MYSQL_TYPE_STRING) {	//strings
			stmtParams[ii].buffer = (char *)params_fields[ii].str;
			stmtParams[ii].buffer_length = MYSQL_STRING_SIZE;
			stmtParams[ii].length = &params_fields[ii].str_len;
		} else if(params_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {	//doubles
			stmtParams[ii].buffer = (char *)&params_fields[ii].val;
		} else if(params_fields[ii].MysqlType==MYSQL_TYPE_DATETIME) {	//dates
			stmtParams[ii].buffer = (char *)&params_fields[ii].dt;
		}
	}
	
	if (mysql_stmt_bind_param(*stmt, stmtParams)) {
		free( stmtParams );
		throw IOException("Error binding parameters", AT);
	}
	free( stmtParams );
}

void bindResults(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &result_fields)
{
	MYSQL_RES *prepare_meta_result = mysql_stmt_result_metadata(*stmt);
	if (!prepare_meta_result) throw IOException("Error executing meta statement", AT);
	const size_t column_count = static_cast<size_t>( mysql_num_fields(prepare_meta_result) );
	if (column_count!=result_fields.size()) throw InvalidArgumentException("Wrong number of columns returned", AT);
	mysql_free_result(prepare_meta_result);
	
	MYSQL_BIND *result = (MYSQL_BIND*)calloc(column_count, sizeof(MYSQL_BIND));
	if (result==nullptr) throw IOException("Could not allocate memory for results binding to Mysql query", AT);
	
	for (size_t ii=0; ii<column_count; ++ii) {
		result[ii].buffer_type = result_fields[ii].MysqlType;
		if (result_fields[ii].MysqlType==MYSQL_TYPE_STRING) {	//strings
			result[ii].buffer = (char *)result_fields[ii].str;
			result[ii].buffer_length = MYSQL_STRING_SIZE;
		} else if(result_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {	//doubles
			result[ii].buffer_type = MYSQL_TYPE_DOUBLE;
			result[ii].buffer = (char *)&result_fields[ii].val;
		} else if(result_fields[ii].MysqlType==MYSQL_TYPE_DATETIME) {	//dates
			result[ii].buffer_type = MYSQL_TYPE_DATETIME;
			result[ii].buffer = (char *)&result_fields[ii].dt;
		}
		
		result[ii].is_null = &result_fields[ii].is_null;
		result[ii].length = &result_fields[ii].buffer_len;
		result[ii].error = &result_fields[ii].error;
	}
	
	if (mysql_stmt_bind_result(*stmt, result)) throw IOException("Error binding results", AT);
	if (mysql_stmt_store_result(*stmt)) throw IOException("mysql_stmt_store_result failed", AT);
	free( result );
}

double retrieveData(const SQL_FIELD &field, const unsigned int& conversion)
{
	const double val = field.val;
	if (field.is_null==1) return IOUtils::nodata;
	
	if (conversion==SQL_FIELD::C_TO_K) return IOUtils::C_TO_K( val );
	if (conversion==SQL_FIELD::NORMALIZE_PC || conversion==SQL_FIELD::CM_TO_M) return val / 100.;
	if (conversion==SQL_FIELD::HPA_TO_PA ) return val * 100.;
	
	return val;
}

}
