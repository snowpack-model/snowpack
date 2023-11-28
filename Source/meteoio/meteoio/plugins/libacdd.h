// SPDX-License-Identifier: LGPL-3.0-or-later
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
 * the whole process as simple as possible. It is however possible to provide some of these attributes in the INI configuration file, using the
 * following keys (all in the [Output] section):
 *  - Overview of the dataset
 *     - ACDD_TITLE: a short title for the data set;
 *     - ACDD_SUMMARY: a paragraph describing the dataset;
 *     - ACDD_SUMMARY_FILE: a file containing a description of the dataset, it overwrites the value of ACDD_SUMMARY if present;
 *     - ACDD_COMMENT: miscellaneous informartion about the dataset;
 *     - ACDD_ACKNOWLEDGEMENT: acknowledgement for the various types of support for the project that produced this data;
 *     - ACDD_KEYWORDS: a list of AGU Index Terms (default: hard-coded list);
 *     - ACDD_KEYWORDS_VOCABULARY: the unique name or identifier of the vocabulary from which keywords are taken (default: AGU);
 *  - Linking the dataset to other resources and metadata
 *     - ACDD_ID: an identifier for the data set, provided by and unique within its naming authority. Example: DOI, URL, text string, but without white spaces
 *     - ACDD_NAMING_AUTHORITY: The organization that provides the initial id (see above) for the dataset
 *     - ACDD_METADATA_LINK: A URL/DOI that gives more complete metadata;
 *     - ACDD_REFERENCES: Published or web-based references that describe the data or methods used to produce it;
 *     - WIGOS_ID: although this is not an ACDD key, it can be very useful in linking datasets together through their <a href="https://oscar.wmo.int/surface/#/">WIGOS ID</a>.
 *  - Origin of the data
 *     - ACDD_PROJECT: the scientific project that created the data;
 *     - ACDD_PROGRAM: The overarching program(s) of which the dataset is a part;
 *     - ACDD_SOURCE: The method of production of the original data;
 *     - ACDD_ACTIVITY_TYPE: Activity types are used to identify the origin of the dataset. Pick one from this <a href="https://htmlpreview.github.io/?https://github.com/metno/mmd/blob/master/doc/mmd-specification.html#activity-type">controlled vocabulary</a>;
 *  - Contact information regarding the data set
 *     - ACDD_INSTITUTION: the institution providing the data set (default: domain name);
 *     - The following can be a list of comma-delimited values but must be kept in-sync (ie if providing two creators, then two emails must also be provided, etc):
 *        - ACDD_CREATOR: the name of the creator of the data set (default: login name);
 *        - ACDD_CREATOR_EMAIL: the email of the creator;
 *        - ACDD_CREATOR_INSTITUTION: the institution of the creator; should uniquely identify the creator's institution;
 *        - ACDD_CREATOR_URL: the URL of the creator principally responsible for creating this data;
 *        - ACDD_CREATOR_TYPE: either person, group, institution, or position (default: person);
 *     - The following can be a list of comma-delimited values but must be kept in-sync (ie if providing two publishers, then two emails must also be provided, etc):
 *        - ACDD_PUBLISHER: the name of the person / entity responsible for publishing the data file or product to users, with its current metadata and format;
 *        - ACDD_PUBLISHER_EMAIL: the email of the person / entity responsible for publishing the data file or product to users;
 *        - ACDD_PUBLISHER_URL: the url of the person / entity responsible for publishing the data file or product to users;
 *        - ACDD_PUBLISHER_TYPE: either person, group, institution, or position (default: person);
 *  - Miscellaneous
 *     - ACDD_PROCESSING_LEVEL: a textual description of the processing level
 *     - ACDD_LICENSE: describes the license applicable to the dataset;
 *     - ACDD_PRODUCT_VERSION: Version identifier of the data file or product as assigned by the data creator (default: 1.0);
 *     - ACDD_OPERATIONAL_STATUS: The current operational status of the product. Choose from the <a href="https://htmlpreview.github.io/?https://github.com/metno/mmd/blob/master/doc/mmd-specification.html#operational-status">controlled vocabulary</a>;
 *
 * 
 * This list contains all mandatory ACDD fields as listed at the <a href="https://adc.met.no/node/4">Arctic Data Centre</a> as well as most of the 
 * optional fields (please note that the geospatial and time coverage are automatically generated based on the data itself and the history is 
 * also automatically handled).
 * 
 * Example of ACDD configuration for a NetCDF file generated for the <a href="https://public.wmo.int/en">WMO</a>'s 
 * <a href="https://globalcryospherewatch.org/">Global Cryosphere Watch</a> (GCW) <a href="https://gcw.met.no/metsis/search">data portal</a>:
 * @code
 * [Output]
 * ACDD_WRITE = TRUE
 * ACDD_CREATOR = Mathias Bavay
 * ACDD_CREATOR_EMAIL = bavay@slf.ch
 * ACDD_CREATOR_INSTITUTION = SLF
 * ACDD_CREATOR_URL = https://slf.ch
 * ACDD_CREATOR_TYPE = person
 * ACDD_PUBLISHER = Mathias Bavay
 * ACDD_PUBLISHER_EMAIL = bavay@slf.ch
 * ACDD_PUBLISHER_URL = https://slf.ch
 * ACDD_PUBLISHER_TYPE = person
 * ACDD_INSTITUTION = SLF
 * ACDD_ID = 5LAR1_MeteoBase
 * ACDD_NAMING_AUTHORITY = SLF
 * ACDD_TITLE = Meteo station at Laret val/cal site
 * ACDD_SUMMARY = Meteo station at Laret val/cal site
 * ACDD_LICENSE = CC-BY-NC
 * ACDD_PROCESSING_LEVEL = raw
 * ACDD_ACTIVITY_TYPE = in situ land-based station
 * ACDD_OPERATIONAL_STATUS = scientific
 * @endcode
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
		void setTimeCoverage(const std::vector<std::string>& vec_timestamp, const double& TZ);
		
	private:
		void defaultInit();
		size_t find(const std::string& search_name) const;
		void checkMultiValueConsistency();
		
		std::vector<std::string> name, cfg_key, value;
		bool enabled; //helper boolean for callers to know if this object should be used or not
};

} //namespace
#endif
