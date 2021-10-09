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
#ifndef METEODATA_H
#define METEODATA_H

#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/StationData.h>

#include <string>
#include <vector>
#include <set>

namespace mio {

class MeteoData; //forward declaration
typedef std::vector<MeteoData> METEO_SET;

/**
 * @class MeteoGrids
 * @brief A class to represent the meteorological parameters that could be contained in a grid.
 * This should be very close to MeteoData with a few additions (like the wind u,v,w)
 * @ingroup data_str
 * @author Mathias Bavay
 * @date   2011-12-22
 */

class MeteoGrids {
	public:
		/// \anchor meteogrids this enum provides names for possible meteogrids (from an ARPS file, etc)
		enum Parameters {firstparam=0,
				TA=firstparam, ///< Air temperature
				RH, ///< Relative humidity
				QI, ///< Specific humidity
				TD, ///< Dew Point temperature
				VW, ///< Wind velocity
				DW, ///< Wind direction
				VW_MAX, ///< Maximum wind velocity
				ISWR, ///< Incoming short wave radiation
				RSWR, ///< Reflected short wave radiation
				ISWR_DIFF, ///< Incoming short wave, diffuse
				ISWR_DIR, ///< Incoming short wave, direct
				ILWR, ///< Incoming long wave radiation
				OLWR, ///< Outgoing long wave radiation
				LWR_NET, ///< Net long wave radiation
				TAU_CLD, ///< Cloud transmissivity or ISWR/ISWR_clear_sky
				HS, ///< Height of snow
				PSUM, ///< Water equivalent of precipitations, either solid or liquid
				PSUM_PH, ///<  Precipitation phase, between 0 (fully solid) and 1 (fully liquid)
				PSUM_L, ///< Water equivalent of liquid precipitation
				PSUM_LC, ///< Water equivalent of liquid convective precipitation
				PSUM_S, ///< Water equivalent of solid precipitation
				TSG, ///< Temperature ground surface
				TSS, ///< Temperature snow surface
				TSOIL, ///< Temperature within the soil, at a given depth
				P, ///< Air pressure
				P_SEA, ///< Sea level air pressure
				U, ///< East component of wind
				V, ///< North component of wind
				W, ///< Vertical component of wind
				SWE, ///< Snow Water Equivalent
				RSNO, ///< Snow mean density
				ROT, ///< Total generated runoff
				ALB, ///< Albedo
				DEM, ///< Digital Elevation Model
				SHADE, ///< Hillshade
				SLOPE, ///< DEM slope angle
				AZI, ///< DEM slope azimuth
				lastparam=AZI};

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static size_t getParameterIndex(const std::string& parname);
		static const std::string getParameterName(const size_t& parindex);
		static const std::string getParameterDescription(const size_t& parindex, const bool& allow_ws=true);
		static const std::string getParameterUnits(const size_t& parindex);

	private:
		//static methods
		static std::vector<std::string> paramname, description, units;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static vector paramname
};

/**
 * @class MeteoData
 * @brief A class to represent a singular measurement received from one station at a certain time (represented by the Date object)
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-12-05
 */

class MeteoData {
	public:
		///Keywords for selecting the toString formats
		typedef enum {
			DFLT, ///< Shows detailed information, skipping nodata fields
			FULL, ///< Shows detailed information, including nodata fields
			COMPACT ///< output optimized to print the content of vector<MeteoData>
		} FORMATS;
		
		/** @brief Available %MeteoData merging strategies.
		* When the two stations both have data at a given time step, only the parameters that are *not* present
		* in station1 will be taken from station2 (ie. station1 has priority).
		*
		* \image html merging_strategies.png "Merging strategies for two stations with different sampling rates"
		* \image latex merging_strategies.eps "Merging strategies for two stations with different sampling rates" width=0.9\textwidth
		* @note Keep in mind that if a station is moving (ie. if its location might change in time) merge strategies other than STRICT_MERGE
		* will introduce potentially invalid metadata (since the new position can not be reconstructed).
		*/
		typedef enum MERGE_TYPE {
				STRICT_MERGE=0, ///< Station1 receives data from station2 only for common timestamps
				EXPAND_MERGE=1, ///< If station2 can provide some data before/after station1, this extra data is added to station1
				FULL_MERGE=2, ///< All timestamps from station2 are brought into station1 even if the timestamps don't match
				WINDOW_MERGE=3 ///< Similar to a full merge but within the time range of station1
		} Merge_Type;
		
		/** @brief Available %MeteoData conflict resolution strategies.
		 * When two stations are merged, there is a merge conflict when they have different values for the same parameter 
		 * at the same timestamp (please note that if one of them has nodata, it is not a conflict: the valid value will be 
		 * used to replace the nodata). How to handle conflicts is defined by the choice of conflicts resolution strategy. 
		 */
		typedef enum MERGE_CONFLICTS {
				CONFLICTS_PRIORITY_FIRST=0, ///< Station1 has priority over station 2
				CONFLICTS_PRIORITY_LAST=1, ///< Station2 has priority over station 1
				CONFLICTS_AVERAGE=2 ///< The merged value is the average of station1 and station2
		} Merge_Conflicts;

		/// \anchor meteoparam this enum provides indexed access to meteorological fields
		enum Parameters {firstparam=0,
		                 P=firstparam, ///< Air pressure
		                 TA, ///< Air temperature
		                 RH, ///< Relative humidity
		                 QI, ///< Specific humidity
		                 TSG, ///< Temperature of the ground surface
		                 TSS, ///< Temperature of the snow surface
		                 HS, ///< Height of snow
		                 VW, ///< Wind velocity
		                 DW, ///< Wind direction
		                 VW_MAX, ///< Maximum wind velocity
		                 U, ///< Wind velocity (west-east component)
		                 V, ///< Wind velocity (sout-north component)
		                 RSWR, ///< Reflected short wave radiation
		                 ISWR, ///< Incoming short wave radiation
		                 ILWR, ///< Incoming long wave radiation (downwelling)
		                 OLWR, ///< Outgoing long wave radiation
		                 LWR_NET, ///< Net long wave radiation
		                 TAU_CLD, ///< Cloud transmissivity or ISWR/ISWR_clear_sky
		                 PSUM, ///< Water equivalent of precipitations, either solid or liquid
		                 PSUM_PH, ///< Precipitation phase: between 0 (fully solid) and 1(fully liquid)
		                 PSUM_L, ///< Water equivalent of liquid precipitation
		                 PSUM_LC, ///< Water equivalent of liquid convective precipitation
		                 PSUM_S, ///< Water equivalent of solid precipitation
		                 lastparam=PSUM_S};

		static const std::string& getParameterName(const size_t& parindex);

		static size_t getStaticParameterIndex(const std::string& parname);

		/**
		 * @brief The default constructor initializing every double attribute to nodata and the Date to julian==0.0
		 */
		MeteoData(void);

		/**
		* @brief A constructor that sets the measurment time
		* @param in_date A Date object representing the time of the measurement
		*/
		MeteoData(const Date& in_date);

		/**
		* @brief A constructor that sets the measurment time and meta data
		* @param date_in A Date object representing the time of the measurement
		* @param meta_in A StationData object containing the meta data
		*/
		MeteoData(const Date& date_in, const StationData& meta_in);

		/**
		* @brief A setter function for the measurement date
		* @param in_date A Date object representing the time of the measurement
		*/
		void setDate(const Date& in_date) {date = in_date;}

		/**
		* @brief Add another variable to the MeteoData object,
		*        a double value will be added and the nrOfParameters increased
		* @param i_paramname A parameter name, e.g. "VSWR"
		* @return A size_t denoting the index of the the parameter added
		*/
		size_t addParameter(const std::string& i_paramname);

		/**
		* @brief Check whether a certain parameter is a part of this MeteoData instance
		* @param parname A string parameter, representing a meteo parameter, e.g. "VSWR"
		* @return A boolean indicating whether the parameter is a part of the object
		*/
		bool param_exists(const std::string& parname) const;

		/**
		 * @brief Resets all the meteo parameters to IOUtils::nodata
		 *        NOTE: member vars date and resampled are not affected
		 */
		void reset();
		
		/**
		 * @brief Are all the fields set to nodata?
		 * @return true if no meteo field has a value
		 */
		bool isNodata() const;

		bool isResampled() const {return resampled;}
		void setResampled(const bool& in_resampled) {resampled = in_resampled;}

		bool isFiltered(const size_t& param) const;
		void setFiltered(const size_t& param, const bool& in_filtered = true);
		bool isGenerated(const size_t& param) const;
		void setGenerated(const size_t& param, const bool& in_generated = true);
		bool isResampledParam(const size_t& param) const;
		void setResampledParam(const size_t& param, const bool& in_resampled = true);

		void standardizeNodata(const double& plugin_nodata);

		double& operator()(const size_t& parindex);
		const double& operator()(const size_t& parindex) const;
		double& operator()(const std::string& parname);
		const double& operator()(const std::string& parname) const;

		const std::string& getNameForParameter(const size_t& parindex) const;
		size_t getParameterIndex(const std::string& parname) const;
		size_t getNrOfParameters() const {return nrOfAllParameters;}

		static MeteoGrids::Parameters findGridParam(const Parameters& mpar);
		static Parameters findMeteoParam(const MeteoGrids::Parameters& gpar);

		/**
		 * @brief Simple merge strategy for two vectors containing meteodata time series for two stations.
		 * If some fields of the MeteoData objects given in the first vector are nodata, they will be
		 * filled by the matching field from the MeteoData objects given in the second vector (if the
		 * same timestamp exist).
		 * @note Only timestamps common to both data sets will be merged!
		 * @param vec1 reference vector, highest priority
		 * @param[in] vec2 extra vector to merge, lowest priority
		 * @param[in] strategy how should the merge be done? (default: STRICT_MERGE)
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 * @return number of merge conflicts
		 */
		static size_t mergeTimeSeries(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const Merge_Type& strategy=STRICT_MERGE, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Simple merge strategy for vectors containing meteodata for a given timestamp.
		 * If some fields of the MeteoData objects given in the first vector are nodata, they will be
		 * filled by the matching field from the MeteoData objects given in the second vector (if the
		 * same location exist). Stations only occuring in the second vector will be appended to the
		 * first vector.
		 * @note two stations are considered to be identical if they fit within a 5m 3D box
		 * @note the vectors are supposed to contain data at a given time stamp. If both vectors don't match a
		 * common time stamp, nothing is done
		 * @param vec1 reference vector, highest priority
		 * @param[in] vec2 extra vector to merge, lowest priority
		 * @param[in] simple_merge if set to true, assume all stations are unique (ie. simply append vec2 to vec1), default: false
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 */
		static void merge(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const bool& simple_merge=false, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Simple merge strategy for vectors containing meteodata for a given timestamp.
		 * If some fields of the MeteoData objects given in the first vector are nodata, they will be
		 * filled by the matching field from the MeteoData object given in the second argument (if the
		 * same location exist). If meteo2 does not describe a station already in vec, it will simply be appended.
		 * @note two stations are considered to be identical if they fit within a 5m 3D box
		 * @note the datasets are supposed to contain data at a given time stamp. If vec1 and meteo2 don't match a
		 * common time stamp, nothing is done
		 * @param vec reference vector, highest priority
		 * @param[in] meteo2 extra MeteoData object to merge, lowest priority
		 * @param[in] simple_merge if set to true, assume all stations are unique (ie.simply append meteo2 to vec), default: false
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 */
		static void merge(std::vector<MeteoData>& vec, const MeteoData& meteo2, const bool& simple_merge=false, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Simple merge strategy within a vector of MeteoData.
		 * All stations that can be considerd as identical (see note) will be merged using the defined conflict resolution
		 * policy (see Merge_Conflicts).
		 * The priority goes to the stations at the beginning of the vector. For example, if vec[0] has TA but no HS and
		 * vec[3] has TA and HS, then vec[0] will <i>keep</i> its TA and get HS from vec[3]. If vec[2] is further away than
		 * 5m from vec[0], then it can not contribute to vec[0].
		 * @note two stations are considered to be identical if they fit within a 5m 3D box
		 * @note the datasets are supposed to contain data at a given time stamp. If the stations don't match a
		 * common time stamp, nothing is done
		 * @param vec vector of stations
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 */
		static void merge(std::vector<MeteoData>& vec, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Simple merge strategy.
		 * If some fields of the object given as first argument are nodata, they will be filled by the matching field from the
		 * provided argument.
		 * @note no check on the location is performed, ie. it can merge data from stations kilometers away...
		 * @param meteo1 reference MeteoData, highest priority
		 * @param[in] meteo2 extra MeteoData to merge, lowest priority
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 * @return new MeteoData object resulting from the merge
		 */
		static MeteoData merge(MeteoData meteo1, const MeteoData& meteo2, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Simple merge strategy.
		 * If some fields of the current object are nodata, they will be filled by the matching field from the
		 * provided argument.
		 * @note no check on the location is performed, ie. it can merge data from stations kilometers away...
		 * @param[in] meteo2 extra MeteoData to merge, lowest priority
		 * @param[in] conflicts_strategy In case of conflicts, define how to resolve them (see Merge_Conflicts), default: CONFLICTS_PRIORITY_FIRST
		 * @return true if no conflicts were found
		 */
		bool merge(const MeteoData& meteo2, const Merge_Conflicts& conflicts_strategy=CONFLICTS_PRIORITY_FIRST);

		/**
		 * @brief Check for data conflicts between two MeteoData objects
		 * Conflicts are defined as two identical data fields having different 
		 * non-nodata values. The metadata are NOT checked for conflict.
		 * @param[in] meteo2 other MeteoData to compare to
		 * @return true if meteo2 has conflicts with the current MeteoData object
		 */
		bool hasConflicts(const MeteoData& meteo2) const;
		
		/**
		 * @brief Parse a string containing a merge type and return the proper enum member for it.
		 * @param[in] merge_type
		 * @return Merge_Type
		 */
		static MeteoData::Merge_Type getMergeType(std::string merge_type);
		
		/**
		 * @brief Parse a string containing a merge conflcits type and return the proper enum member for it.
		 * @param[in] merge_conflicts
		 * @return Merge_Conflicts
		 */
		static MeteoData::Merge_Conflicts getMergeConflicts(std::string merge_conflicts);

		/**
		 * @brief List the parameters that have a least one valid value in a vector of MeteoData.
		 * @param[in] vecMeteo vector to read the data from
		 * @return parameters that have at least one valid value
		 */
		static std::set<std::string> listAvailableParameters(const std::vector<MeteoData>& vecMeteo);

		/**
		 * @brief Print the content of the current object
		 * @param[in] format select the preferred output format
		 * @return string containing a human-readable representation of the content of the object
		 */
		const std::string toString(const FORMATS format=DFLT) const;
		friend std::ostream& operator<<(std::ostream& os, const MeteoData& data);
		friend std::istream& operator>>(std::istream& is, MeteoData& data);

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		inline bool operator!=(const MeteoData& in) const {return !(*this==in);} ///<Operator that tests for inequality
		inline bool operator<(const MeteoData& in) const {return (this->date < in.date);} ///<so vectors can be sorted by timestamps
		inline bool operator>(const MeteoData& in) const {return (this->date > in.date);} ///<so vectors can be sorted by timestamps

		//direct access allowed
		Date date; ///<Timestamp of the measurement
		StationData meta; ///<The meta data of the measurement

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData

		const std::string getStationID() const {return meta.stationID;}

	private:

		struct flag_field { //object to hold all data quality / datapoint meta flags
			bool filtered : 1;
			bool resampled : 1;
			bool generated : 1;
			bool created : 1; //TODO: not yet filled in the creators!
			bool : 1; //leave one empty for future flag
			unsigned int extra_flag : 3; //rest of byte space could be any 3-bit-value
			float offset; //an offset that has been applied to the data at this time step
		};

		//static methods
		static std::vector<std::string> s_default_paramname; ///<Associate a name with meteo parameters in Parameters
		static const double epsilon; ///<for comparing fields
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static vector s_default_paramname

		//private data members, please keep the order consistent with declaration lists and logic!
		std::vector<std::string> extra_param_name;
		std::vector<double> data;

		size_t nrOfAllParameters;

		//data qa containers:
		bool resampled; ///<set this to true if MeteoData is result of resampling
		std::vector<flag_field> flags; ///<Per-parameter data quality flags BUG these are currently not serialized!!
		static flag_field zero_flag; //TODO: is there a way to make this const?
};

} //end namespace

#endif
