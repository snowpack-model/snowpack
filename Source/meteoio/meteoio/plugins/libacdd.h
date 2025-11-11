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
 * @details The final value of any acdd field is provided by three sources, any attempts
 * stops as soon as some value has been found: 1. the INI configuration key; 2. reading the matching environment variable; 3. some hard-coded
 * defaults (whenever it makes sense and is possible). It is therefore possible to define relevant default values system-wide by setting an environment variable of the
 * same name as the INI configuration key and individual configuration files might still overwrite this default by explicitly setting the configuration key.
 *
 * The following keys are supported to define %ACDD attributes (all in the [Output] section):
 *  - Overview of the dataset
 *     - ACDD_TITLE: a short title for the data set;
 *     - ACDD_SUMMARY: a paragraph describing the dataset;
 *     - ACDD_SUMMARY_FILE: a file containing a description of the dataset, it overwrites the value of ACDD_SUMMARY if present;
 *     - ACDD_COMMENT: miscellaneous information about the dataset;
 *     - ACDD_ACKNOWLEDGEMENT: acknowledgment for the various types of support for the project that produced this data;
 *     - ACDD_KEYWORDS: a list of comma delimited keywords, for example from the <a href="https://gcmd.earthdata.nasa.gov/KeywordViewer/">GCMD Science Keywords</a> (GCMDSK) with their full path (default: hard-coded list);
 *     - ACDD_KEYWORDS_VOCABULARY: the unique name or identifier of the vocabulary from which keywords are taken (default: GCMDSK);
 *  - Linking the dataset to other resources and metadata
 *     - ACDD_ID: an identifier for the data set, provided by and unique within its naming authority. Example: DOI, URL, text string, but without white spaces
 *     - ACDD_NAMING_AUTHORITY: The organization that provides the initial dataset id (see above) for the dataset
 *     - ACDD_METADATA_LINK: A URL/DOI that gives more complete metadata (if a WIGOS_ID is provided, the matching <a href="https://oscar.wmo.int/surface/#/">Oscar surface</a> metadata page is used as default);
 *     - ACDD_REFERENCES: Published or web-based references that describe the data or methods used to produce it;
 *     - WIGOS_ID: although this is not an ACDD key, it can be very useful in linking datasets together through their <a href="https://oscar.wmo.int/surface">WIGOS ID</a>.
 *  - Origin of the data
 *     - ACDD_PROJECT: the scientific project that created the data;
 *     - ACDD_PROGRAM: The overarching program(s) of which the dataset is a part;
 *     - ACDD_SOURCE: The method of production of the original data;
 *     - ACDD_COVERAGE_CONTENT_TYPE: the type of the source data (physicalMeasurement, modelResult... as ISO 19115-1 code);
 *     - ACDD_ACTIVITY_TYPE: Activity types are used to identify the origin of the dataset. Pick one from this <a href="https://html-preview.github.io/?url=https://github.com/metno/mmd/blob/master/doc/mmd-specification.html#activity-type">controlled vocabulary</a>;
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
 *     - The following can be a list of comma-delimited values but must be kept in-sync (ie if providing two contributors, then two roles must also be provided, etc):
 *        - ACDD_CONTRIBUTOR: the name of the individuals, institution or projects that have contributed to the data set (default: login name);
 *        - ACDD_CONTRIBUTOR_ROLE: the role of the contributor;
 *  - Miscellaneous
 *     - ACDD_PROCESSING_LEVEL: a textual description of the processing level
 *     - ACDD_LICENSE: describes the license applicable to the dataset;
 *     - ACDD_PRODUCT_VERSION: Version identifier of the data file or product as assigned by the data creator (default: 1.0);
 *     - ACDD_OPERATIONAL_STATUS: The current operational status of the product. Choose from the <a href="https://html-preview.github.io/?url=https://github.com/metno/mmd/blob/master/doc/mmd-specification.html#operational-status">controlled vocabulary</a>;
 *     - ACDD_DATASET_PRODUCTION_STATUS: the production status of the product. Choose from the <a href="https://html-preview.github.io/?url=https://github.com/metno/mmd/blob/master/doc/mmd-specification.html#dataset-production-status-types">controlled vocabulary</a>;
 *
 * 
 * This list contains all mandatory ACDD fields as listed at the <a href="https://adc.met.no/node/4">Arctic Data Centre</a> as well as most of the 
 * optional fields (please note that the geospatial and time coverage are automatically generated based on the data itself and the history is 
 * also automatically handled). 
 * 
 * \note It is possible to write an <a href="https://docs.unidata.ucar.edu/netcdf-java/5.6/userguide/ncml_overview.html">NcML</a> 
 * file containing the ACDD metadata alongside each output data file. This is mostly aimed at data catalogs (such as for the meteoio
 * webservice) that can therefore easily retrieve metadata without even having to open the data file itself (which would also require
 * them to have a parser for each supported file format). This is controlled by the ACDD_WRITE_NCML boolean configuration key that can
 * also be forced to TRUE with the FORCE_NCML compilation flag (in order to enforce the creation of NcML file on a server for example).
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
 * ACDD_KEYWORDS = EARTH SCIENCE > CRYOSPHERE > SNOW/ICE, EARTH SCIENCE > ATMOSPHERE
 * ACDD_KEYWORDS_VOCABULARY = GCMDSK
 * ACDD_LICENSE = CC-BY-NC
 * ACDD_PROCESSING_LEVEL = raw
 * ACDD_ACTIVITY_TYPE = in situ land-based station
 * ACDD_OPERATIONAL_STATUS = scientific
 * @endcode
*/
class ACDD {
	public:
		enum Mode {MERGE, REPLACE, APPEND};
		struct ACDD_ATTR; //forward declaration
		struct VARS_ATTR; //forward declaration
		typedef struct ACDD_ATTR acdd_attrs;
		typedef struct VARS_ATTR vars_attr;
		
		/**
		* @brief Constructor, the argument allows the object to know if the acdd metadata should be written out or not
		* @param[in] set_enable enable ACDD support?
		*/
		ACDD(const bool& set_enable);
		
		//defining some iterators so the callers can loop over all available attributes
		using const_iterator = std::map<std::string, acdd_attrs>::const_iterator;
		const_iterator cbegin() const noexcept { return attributes.cbegin(); }
		const_iterator cend() const noexcept { return attributes.cend(); }
		
		/**
		* @brief Set an internal boolean as a helper for the caller to know if ACDD support should be enabled or not. Moreover, it
		* initializes the attributes map if not already done.
		* @param[in] i_enable enable ACDD support?
		*/
		void setEnabled(const bool& i_enable);
		
		void setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line=true);
		
		void addAttribute(const std::string& att_name, const std::string& att_value, const Mode& mode=MERGE);
		void addAttribute(const std::string& att_name, const double& att_value, const Mode& mode=MERGE);

		/**
		* @brief Delete an attribute (this is required in some rare cases if an attribute is handled separately, such as "history"
		* @param[in] att_name name of the attibute to delete
		*/
		void deleteAttribute(const std::string& att_name) {attributes.erase( att_name );}

		/**
		* @brief Get an internal boolean as a helper for the caller to know if ACDD support should be enabled or not
		* @return enable ACDD support from the caller side?
		*/
		bool isEnabled() const {return enabled;}

		std::string getAttribute(std::string &att_name) const;
		
		void setGeometry(const mio::Grid2DObject& grid, const bool& isLatLon);
		void setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon);
		void setGeometry(const mio::Coords& location, const bool& isLatLon);
		void setGeometry(const std::set< mio::Coords >& vecLocation, const bool& isLatLon);
		
		void setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo);
		void setTimeCoverage(const std::vector<mio::MeteoData>& vecMeteo);
		void setTimeCoverage(const std::vector<std::string>& vec_timestamp, const double& TZ);

		
		//support for NcML
		/**
		* @brief Should an <a href="https://docs.unidata.ucar.edu/netcdf-java/5.6/userguide/ncml_overview.html">NcML</a> file be written
		* alongside the data file?  
		* \note Please notice that this can be forced set to true by a compilation flag independently of the argument
		* provided to this call (as is necessary for the meteoio webservice).
		* @param[in] i_enable set to true if an NcML file should be produced (see caveat above).
		*/
		void setEnableNcML(const bool& i_enable);
		bool enableNcML() const {return enable_ncml;}
		void writeNcML(const std::string& data_filename) const;
		void addDimension( const std::string& var_name, const std::string& var_long_name, const size_t& length);
		void addVariable( const std::string& var_name, const std::string& var_long_name, const std::string& var_units);
		
		std::string toString() const;
		
	private:
		static std::map<std::string, acdd_attrs> initAttributes();
		static std::set< std::pair< std::string, std::set<std::string> > > initLinks();
		static size_t countCommas(const std::string& str);
		void checkLinkedAttributes();
		
		static bool isWigosID(const std::string& str);
		std::map<std::string, acdd_attrs> attributes; //all the ACDD attributes with their properties
		std::set< std::pair< std::string, std::set<std::string> > > linked_attributes; //attribute names that are linked together, ie must have the same number of sub-elements (comma delimited)
		std::set< vars_attr > dimensions, variables;
		bool enabled; //helper boolean for callers to know if this object should be used or not
		bool enable_ncml; //helper boolean for callers to trigger the writing of NcML files together with the data file
};

struct ACDD::VARS_ATTR {	//TODO at some point, this should be stored in MeteoData as timeseries metadata
	VARS_ATTR(const std::string& i_name, const std::string& i_standard_name, const std::string& i_units) : name(i_name), standard_name(i_standard_name), units(i_units), length(IOUtils::npos) {}
	VARS_ATTR(const std::string& i_name, const std::string& i_standard_name, const size_t& i_length) : name(i_name), standard_name(i_standard_name), units(), length(i_length) {}
	
	//"units" is ignored here, since with the same name and standard_name, they would have the same units
	bool operator<(const ACDD::VARS_ATTR& other) const {
        return (name < other.name) ||
               (name == other.name && standard_name < other.standard_name) ||
               (name == other.name && standard_name == other.standard_name && length < other.length);
    }
	
	std::string name, standard_name, units;
	size_t length;
};


/**
 * @struct ACDD::ACDD_ATTR
 * @brief This structure provides low level functions to handle and store individual ACDD fields.
 * @details
 * It is designed to retrieve ACDD values from configuration keys (including from a file pointed to by such
 * a key), from a default value provided to the constructor or as provided by the user. Several modes
 * of combining these different sources can be chosen from.
 */
struct ACDD::ACDD_ATTR {
	ACDD_ATTR() : name(), value(), cfg_key(), default_value(), Default(true) {}

	ACDD_ATTR(const std::string& i_name, const std::string& i_cfg_key, const std::string& i_default_value="") : name(i_name), value(), cfg_key(i_cfg_key), default_value(i_default_value), Default(i_name.empty()) {}

	ACDD_ATTR(const std::string& i_name, const std::string& i_value, const std::string& i_cfg_key, const std::string& i_default_value) : name(i_name), value(i_value), cfg_key(i_cfg_key), default_value(i_default_value), Default(false) {}

	std::string getValue() const {return value;}
	std::string getName() const {return name;}
	void setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line);
	void setValue(const std::string& i_value, const Mode& mode=MERGE);
	bool isDefault() const {return Default;}

private:
	static void readFromFile(std::string& value, const mio::Config& cfg, const std::string& cfg_key, const std::string& section, const bool& allow_multi_line);

	std::string name, value, cfg_key, default_value;
	bool Default;
};

} //namespace
#endif
