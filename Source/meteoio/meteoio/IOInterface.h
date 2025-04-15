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
#ifndef IOINTERFACE_H
#define IOINTERFACE_H

#include <meteoio/Config.h> //so the plugins can get access to Config for their constructor
#include <meteoio/dataClasses/Array2D.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/dataClasses/Grid3DObject.h>
#include <meteoio/dataClasses/MeteoData.h>

#include <vector>
#include <string>

namespace mio {

/**
 * @class LinesRange
 * @brief A class to represent and handle ranges of lines. They can be sorted, 
 * checked for uniqueness and a line number can be compared to the range (is it 
 * before or after?).
 *
 * @author Mathias Bavay
 */
class LinesRange {
	public:
		LinesRange() : start(), end() {}
		LinesRange(const size_t& l1, const size_t& l2) : start(l1), end(l2) {}
		
		/**
		 * @brief Is the provided line number within the current range?
		 * @param[in] ll line number to check
		 * @return true if the line number is within the current range, false otherwise
		 */
		bool in(const size_t& ll) const {
			return (ll >= start && ll <= end);
		}
		
		/**
		 * @brief Is the provided line number before the *end* of the range?
		 * @param[in] ll line number to check
		 * @return true if the line number is less than the end of the current range, false otherwise
		 */
		bool operator<(const size_t& ll) const {
			return end < ll;
		}
		
		/**
		 * @brief Is the provided line number after the *start* of the range?
		 * @param[in] ll line number to check
		 * @return true if the line number is greater than the end of the current range, false otherwise
		 */
		bool operator>(const size_t& ll) const {
			return start > ll;
		}
		
		bool operator<(const LinesRange& ll) const { //needed for "sort"
			if (start==ll.start) return end < ll.end;
			return start < ll.start;
		}
		
		bool operator==(const LinesRange& ll) const { //needed to check for uniqueness
			return (start==ll.start) && (end==ll.end);
		}
		
		const std::string toString() const {std::ostringstream os; os << "[" << start << " - " << end << "]"; return os.str();}
		
		size_t start, end;
};

/**
 * @class IOInterface
 * @brief A class representing the IO Layer of the software Alpine3D. For each type of IO (File, DB, Webservice, etc)
 * a derived class is to be created that holds the specific implementation of the appropriate virtual methods.
 * The IOHandler class is a wrapper class that is able to deal with all above implementations of the IOInterface abstract base class.
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2009-01-08
 */
class IOInterface {
	public:
		virtual ~IOInterface() {}

		/**
		* @brief Return the list of grids within a given time period that could be read by the plugin, if requested.
		* @details This call should be implemented by all plugins reading grids, so the GridsManager can perform all kinds of
		* advanced features with grids (computing a parameter from other ones, temporally interpolating, etc)
		* @param[in] start the start of the time interval
		* @param[in] end the end of the interval
		* @param[out] list the list of grids. The first key is the date, then a set of parameters (from MeteoGrids::Parameters)
		* @return true if the list could be filled (even if empty), false if such as list can not be filled by the plugin (for example, it does not have
		* the data before actually reading it)
		*/
		virtual bool list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list);

		/**
		* @brief A generic function for parsing 2D grids into a Grid2DObject. The string parameter shall be used for addressing the
		* specific 2D grid to be parsed into the Grid2DObject, relative to GRID2DPATH for most plugins.
		* @param grid_out A Grid2DObject instance
		* @param parameter A std::string representing some information for the function on what grid to retrieve
		*/
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");

		/**
		* @brief Read the given meteo parameter into a Grid2DObject.
		* Each plugin has its own logic for finding the requested meteo parameter grid relative to GRID2DPATH for most plugins
		* @param grid_out A Grid2DObject instance
		* @param parameter The meteo parameter grid type to return (ie: air temperature, wind component, etc)
		* @param date date of the data to read
		*/
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		/**
		* @brief Read the given meteo parameter into a vector for a list of points.
		* Each plugin has its own logic for finding the requested meteo parameter grid relative to GRID2DPATH for most plugins
		* @param data A double vector to hold the data
		* @param parameter The meteo parameter grid type to return (ie: air temperature, wind component, etc)
		* @param date date of the data to read
		* @param Pts vector of points to read from the grid
		*/
		virtual void readPointsIn2DGrid(std::vector<double>& data, const MeteoGrids::Parameters& parameter, const Date& date, const std::vector< std::pair<size_t, size_t> >& Pts);

		/**
		* @brief A generic function for parsing 3D grids into a Grid3DObject. The string parameter shall be used for addressing the
		* specific 3D grid to be parsed into the Grid3DObject, relative to GRID3DPATH for most plugins.
		* @param grid_out A Grid3DObject instance
		* @param parameter A std::string representing some information for the function on what grid to retrieve
		*/
		virtual void read3DGrid(Grid3DObject& grid_out, const std::string& parameter="");

		/**
		* @brief Read the given meteo parameter into a Grid3DObject.
		* Each plugin has its own logic for finding the requested meteo parameter grid relative to GRID3DPATH for most plugins
		* @param grid_out A Grid3DObject instance
		* @param parameter The meteo parameter grid type to return (ie: air temperature, wind component, etc)
		* @param date date of the data to read
		*/
		virtual void read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		/**
		* @brief Parse the DEM (Digital Elevation Model) into the Grid2DObject
		*
		* Example Usage:
		* @code
		* Grid2DObject dem;
		* IOHandler io1("io.ini");
		* io1.readDEM(dem);
		* @endcode
		* @param dem_out A Grid2DObject that holds the DEM
		*/
		virtual void readDEM(DEMObject& dem_out);

		/**
		* @brief Parse the landuse model into the Grid2DObject
		*
		* Example Usage:
		* @code
		* Grid2DObject landuse;
		* IOHandler io1("io.ini");
		* io1.readLanduse(landuse);
		* @endcode
		* @param landuse_out A Grid2DObject that holds the landuse model
		*/
		virtual void readLanduse(Grid2DObject& landuse_out);

		/**
		* @brief Parse the input glacier grid into the Grid2DObject
		*
		* Example Usage:
		* @code
		* Grid2DObject glacier;
		* IOHandler io1("io.ini");
		* io1.readGlacier(glacier);
		* @endcode
		* @param glacier_out A Grid2DObject that holds the glacier height
		*/
		virtual void readGlacier(Grid2DObject& glacier_out);

		/**
		* @brief Fill vecStation with StationData objects for a certain date of interest
		*
		* Example Usage:
		* @code
		* vector<StationData> vecStation;  //empty vector
		* Date d1(2008,06,21,11,0, 1.);       //21.6.2008 11:00 UTC+1
		* IOHandler io1("io.ini");
		* io1.readStationData(d1, vecStation);
		* @endcode
		* @param date A Date object representing the date for which the meta data is to be fetched
		* @param vecStation  A vector of StationData objects to be filled with meta data
		*/
		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);

		/**
		* @brief Fill vecMeteo with a time series of objects
		* corresponding to the interval indicated by dateStart and dateEnd.
		*
		* Matching rules:
		* - if dateStart and dateEnd are the same: return exact match for date
		* - if dateStart > dateEnd: return first data set with date > dateStart
		* - read in all data starting with dateStart until dateEnd
		* - if there is no data at all then the vectors will be empty, no exception will be thrown
		*
		* Example Usage:
		* @code
		* vector< vector<MeteoData> > vecMeteo;      //empty vector
		* Date d1(2008,06,21,11,0, 1);       //21.6.2008 11:00 UTC+1
		* Date d2(2008,07,21,11,0, 1);       //21.7.2008 11:00 UTC+1
		* IOHandler io1("io.ini");
		* io1.readMeteoData(d1, d2, vecMeteo);
		* @endcode
		* @param dateStart   A Date object representing the beginning of an interval (inclusive)
		* @param dateEnd     A Date object representing the end of an interval (inclusive)
		* @param vecMeteo    A vector of vector<MeteoData> objects to be filled with data
		*/
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

		/**
		* @brief Write vecMeteo time series to a certain destination
		*
		* Example Usage:
		* Configure the io.ini to use a certain plugin for the output:
		* @code
		* METEODEST     = GEOTOP
		* METEODESTPATH = /tmp
		* METEODESTSEQ  = Iprec SWglobal
		* @endcode
		* An example implementation (reading and writing):
		* @code
		* vector< vector<MeteoData> > vecMeteo;      //empty vector
		* Date d1(2008,06,21,11,0, 1.);       //21.6.2008 11:00 UTC+1
		* Date d2(2008,07,21,11,0, 1.);       //21.7.2008 11:00 UTC+1
		* IOHandler io1("io.ini");
		* io1.readMeteoData(d1, d2, vecMeteo);
		* io1.writeMeteoData(vecMeteo)
		* @endcode
		* @param vecMeteo    A vector of vector<MeteoData> objects to be filled with data
		* @param name        (optional string) Identifier useful for the output plugin (it could become part
		*                    of a file name, a db table, etc)
		*/
		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		/**
		* @brief Parse the assimilation data into a Grid2DObject for a certain date represented by the Date object
		*
		* Example Usage:
		* @code
		* Grid2DObject adata;
		* Date d1(2008,06,21,11,0, 1.);       //21.6.2008 11:00 UTC+1
		* IOHandler io1("io.ini");
		* io1.readAssimilationData(d1, adata);
		* @endcode
		* @param date_in A Date object representing the date of the assimilation data
		* @param da_out  A Grid2DObject that holds the assimilation data for every grid point
		*/
		virtual void readAssimilationData(const Date& date_in, Grid2DObject& da_out);

		/**
		* @brief Read a list of points by their grid coordinates
		* This allows for example to get a list of points where to produce more detailed outputs.
		* @param pts (std::vector<Coords>) A vector of points coordinates
		*/
		virtual void readPOI(std::vector<Coords>& pts);

		/**
		* @brief Write a Grid2DObject
		* The filename is specified relative to GRID2DPATH for most plugins
		* @param grid_out (Grid2DObject) The grid to write
		* @param options (string) Identifier useful for the output plugin (it could become part of a file name, a db table, etc)
		*/
		virtual void write2DGrid(const Grid2DObject& grid_out, const std::string& options="");

		/**
		* @brief Write a Grid2DObject containing a known meteorological parameter
		* A filename is built relative to GRID2DPATH for most plugins.
		* @param grid_out (Grid2DObject) The grid to write
		* @param parameter The meteo parameter grid type of the provided grid object (ie: air temperature, wind component, etc)
		* @param date date of the data to write
		*/
		virtual void write2DGrid(const Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		/**
		* @brief Write a Grid3DObject
		* The filename is specified relative to GRID3DPATH for most plugins
		* @param grid_out (Grid3DObject) The grid to write
		* @param options (string) Identifier useful for the output plugin (it could become part of a file name, a db table, etc)
		*/
		virtual void write3DGrid(const Grid3DObject& grid_out, const std::string& options="");

		/**
		* @brief Write a Grid3DObject comtaining a known meteorological parameter
		* A filename is build relative to GRID3DPATH for most plugins.
		* @param grid_out (Grid3DObject) The grid to write
		* @param parameter The meteo parameter grid type of the provided grid object (ie: air temperature, wind component, etc)
		* @param date date of the data to write
		*/
		virtual void write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		static void set2DGridLatLon(Grid2DObject &grid, const double& i_ur_lat, const double& i_ur_lon);
		static double computeGridXYCellsize(const std::vector<double>& vecX, const std::vector<double>& vecY);
		
		/**
		 * @brief built the set of line ranges to read or skip.
		 * @details Then each plugin is responsible to call this method if necessary and implement the lines skipping if necessary.
		 * Obviously this can not be implemented by every plugin! The line ranges are given as a comma delimited list of 
		 * either single line numbers or ranges (line numbers delimited by a "-" character). Extra spaces can be given for more clarity
		 * in the input.
		 * @param[in] args the textual representation of the line ranges or lines to parse
		 * @param[in] where informative string to describe which component it is in case of error messages (ex. "CSV plugin")
		 * @param[in] negate take the negation of the provided ranges (converting a "ONLY" statement into an "EXCLUDE" statement)
		 * @return set of line ranges
		 */
		static std::vector< LinesRange > initLinesRestrictions(const std::string& args, const std::string& where, const bool& negate);
		
	protected:
		//this enum defines the different options to add some version information into a file name
		typedef enum VERSIONING_TYPE {
		            NO_VERSIONING, ///< no type selected
					STRING, ///< fixed string, user provided
		            NOW, ///< creation time
		            DATA_START, ///< date of the start of the data
		            DATA_END, ///< date of the end of the data
		            DATA_YEARS ///< start and end year of the data (if they are the same, it is not repeated)
		} VersioningType;

		/**
		 * @brief Build a version identification to add to the output file name
		 * @param[in] versioning versioning type
		 * @param[in] vecMeteo the whole metadata that will be written to an output file (to extract start/end dates)
		 * @param[in] tz timezone to use for the dates
		 * @param[in] versioning_str fixed string to use as version (to be provided by the user)
		 * @return string version to integrate into the output file name
		 */
		static std::string buildVersionString(const VersioningType& versioning, const std::vector< std::vector<MeteoData> >& vecMeteo, const double& tz, const std::string& versioning_str);

		static std::string buildVersionString(const VersioningType& versioning, const std::vector<MeteoData>& vecMeteo, const double& tz, const std::string& versioning_str);

		/**
		 * @brief Merge potentially overlaping line ranges
		 * @param[in] lines_specs sorted, unique and non-overlapping set of line ranges
		 */
		static void mergeLinesRanges(std::vector< LinesRange >& lines_specs);
};

} //end namespace

#endif
