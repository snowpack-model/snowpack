/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef LIBACDD_H
#define LIBACDD_H

#include <meteoio/IOUtils.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/Grid2DObject.h>

#include <string>
#include <vector>

namespace mio {
/**
 * @class ACDD
 * @brief This class contains and handles NetCDF Attribute Conventions Dataset Discovery attributes (see 
 * <A href="http://wiki.esipfed.org/index.php?title=Category:Attribute_Conventions_Dataset_Discovery">ACDD</A>).
 * @details A few attributes can get their default value automatically from the data. For the others, some "best efforts" are made in order to keep
 * the whole process as simple as possible. It is however possible to provide some of these attributes from the configuration file, using the
 * following keys:
 *  - ACDD_CREATOR: the name of the creator of the data set (default: login name);
 *  - ACDD_EMAIL: the email of the creator;
 *  - ACDD_KEYWORDS: a list of AGU Index Terms (default: hard-coded list);
 *  - ACDD_TITLE: a short title for the data set;
 *  - ACDD_INSTITUTION: the institution providing the data set (default: domain name);
 *  - ACDD_PROJECT: the scientific project that created the data;
 *  - ACDD_PROGRAM: The overarching program(s) of which the dataset is a part;
 *  - ACDD_ID: an identifier for the data set, provided by and unique within its naming authority. Example: DOI, URL, text string, but without white spaces
 *  - ACDD_NAMING_AUTHORITY: The organization that provides the initial id (see above) for the dataset
 *  - ACDD_PROCESSING_LEVEL: a textual description of the processing level
 *  - ACDD_SUMMARY: a paragraph describing the dataset;
 *  - ACDD_SUMMARY_FILE: a file containing a description of the dataset, it overwrites the value of ACDD_SUMMARY if present;
 *  - ACDD_COMMENT: miscellaneous informartion about the dataset;
 *  - ACDD_ACKNOWLEDGEMENT: acknowledgement for the various types of support for the project that produced this data;
 *  - ACDD_METADATA_LINK: A URL/DOI that gives more complete metadata;
 *  - ACDD_LICENSE: describes the license applicable to the dataset;
 *  - ACDD_PRODUCT_VERSION: Version identifier of the data file or product as assigned by the data creator (default: 1.0).
*/
class ACDD {
	public:
		enum Mode {MERGE, REPLACE, APPEND};
		
		/**
		* @brief Constructor, the argument allows the object to know if the acdd metadata should be written out or not
		* @param[in] set_enable enable ACDD support?
		*/
		ACDD(const bool& set_enable) : name(), cfg_key(), value(), enabled(set_enable) {defaultInit();}
		
		/**
		* @brief Set an internal boolean as a helper for the caller to know if ACDD support should be enabled or not
		* @param[in] i_enable enable ACDD support?
		*/
		void setEnabled(const bool& i_enable)  {enabled=i_enable;}
		void setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line=true);
		
		void addAttribute(const std::string& att_name, const std::string& att_value, const std::string& att_cfg_key="", Mode mode=MERGE);
		void addAttribute(const std::string& att_name, const double& att_value, const std::string& att_cfg_key="", const Mode& mode=MERGE);

		/**
		* @brief Get an internal boolean as a helper for the caller to know if ACDD support should be enabled or not
		* @return enable ACDD support from the caller side?
		*/
		bool isEnabled() const {return enabled;}
		void getAttribute(const size_t ii, std::string &att_name, std::string & att_value) const;
		size_t getNrAttributes() const {if(enabled) return name.size(); else return 0;}
		
		void setGeometry(const mio::Grid2DObject& grid, const bool& isLatLon);
		void setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon);
		void setGeometry(const mio::Coords& location, const bool& isLatLon);
		void setGeometry(const std::vector< mio::Coords >& vecLocation, const bool& isLatLon);
		
		void setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo);
		void setTimeCoverage(const std::vector<mio::MeteoData>& vecMeteo);
		void setTimeCoverage(const std::vector<std::string>& vec_timestamp);
		
	private:
		void defaultInit();
		size_t find(const std::string& search_name) const;
		
		std::vector<std::string> name, cfg_key, value;
		bool enabled; //helper boolean for callers to know if this object should be used or not
};

} //namespace
#endif
