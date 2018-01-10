/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef METEO2DINTERPOLATOR_H
#define METEO2DINTERPOLATOR_H

#include <meteoio/TimeSeriesManager.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/Buffer.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

#include <vector>
#include <map>

namespace mio {

class InterpolationAlgorithm;

/**
 * @page dev_2Dinterpol How to write a spatial interpolation algorithm
 * Point measurements can be spatially interpolated by MeteoIO, through the use of interpolation
 * algorithms. The user will then be able to choose for each meteorological parameter which
 * interpolations could be applicable and the system will choose (through a heuristic) which one
 * to apply at each time step (depending on the conditions of the moment, like the number of measurements).
 *
 * @section structure_2Dinterpol Structure
 * The selection of which interpolation algorithm to use at any given time step, for a given parameter is
 * performed by the Meteo2DInterpolator class. This class provides the interface for the spatial
 * interpolations. The interpolation algorithms themselves derive from the
 * InterpolationAlgorithm class that standardizes their public API (which would not be used by anything
 * but a Meteo2DInterpolator object). It contains a getQualityRating() method that must return a quality
 * index between 0 (algorithm not applicable) and 1 (perfect result if using this algorithm). This is currently
 * only based on extremely simple heuristics, using general knowledge about the applicability of the various
 * spatial interpolation methods depending on some obvious factors (number of measurement points, etc). The
 * Meteo2DInterpolator object will call this method from all the algorithms listed by the user (in his io.ini
 * configuration file) and keep the one that gets the highest score for interpolating the current parameter
 * at the current time step. The interpolation is then done calling the algorithm's calculate method.
 *
 * @section implementation_2Dinterpol Implementation
 * It is therefore necessary to create in InterpolationAlgorithms.cc (and declared in the .h) a new class,
 * nammed after the algorithm that will be implemented and inheriting InterpolationAlgorithm. Three methods need
 * to be implemented:
 * - a constructor that does the arguments parsing (if any);
 * - double getQualityRating()
 * - void calculate(Grid2DObject& grid)
 *
 * The getQualityRating() method takes the meteorological parameter that will be interpolated and set the param
 * private member to it. It then computes the private member nrOfMeasurments that contains the number of
 * stations that have this meteorological parameter available by either calling getData(param, vecData, vecMeta), which
 * also fills the vectors vecData and vecMeta with the available data (as double) and metadata (as StationData) or
 * directly filling veMeteo and vecMeta. Custom data preparation can obviously be done in this method.
 *
 * The calculate method must properly erase and reste the grid that it receives before filling it. If necessary,
 * (as is the case for precipitation, relative humidity and snow height, for example) the grid can be checked for min/max by
 * calling checkMinMax() at the end of Meteo2DInterpolator::interpolate.It can also add extra information about the
 * interpolation process (such as a regression coefficient or error estimate) to the InterpolationAlgorithm::info
 * stringstream (which will be made available to external programs, such as GUIs).
 *
 * The new class and its associated end user key must be used and its constructor called in AlgorithmFactory::getAlgorithm.
 * It is recommended that any generic statistical
 * spatial processing be implemented as a static class in libinterpol2D.cc so that it could be reused by other
 * algorithms (see for example Interpol2D::IDW and IDWCore). In any case, proper doxygen documentation
 * must be written alongside the implementation.
 *
 * @section doc_2Dinterpol Documentation
 * The newly added interpolation algorithm must be added to the list of available algorithms in
 * InterpolationAlgorithms.h with a proper description. An example can also be given in the example section
 * of the same file. Please feel free to add necessary bibliographic references to the bibliographic section!
 *
*/

/**
 * @class Meteo2DInterpolator
 * @brief A class to spatially interpolate meteo parameters. For more, see \ref interpol2d
 *
 * @ingroup stats
 */

class Meteo2DInterpolator {
	public:
		/**
		* @brief Constructor.
		*/
		Meteo2DInterpolator(const Config& i_cfg, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		Meteo2DInterpolator(const Meteo2DInterpolator&);
		Meteo2DInterpolator& operator=(const Meteo2DInterpolator&); ///<Assignement operator

		~Meteo2DInterpolator();

		///Keywords for virtual stations strategy
		enum RESAMPLING_STRATEGY {
			NONE, ///< default: no resampling
			VSTATIONS, ///< extract virtual stations as specified in the ini file
			GRID_EXTRACT, ///< extract data from grids at locations provided in the ini file
			GRID_ALL, ///< extract all grid points from a provided grid
			GRID_SMART ///< extract all relevant grid points from a provided grid
		};

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result);

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 * @param InfoString some information about the interpolation process (useful for GUIs)
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result, std::string& InfoString);

		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string);

		/**
		 * @brief Retrieve the arguments vector for a given interpolation algorithm
		 * @param[in] parname the meteorological parameter that is concerned
		 * @param[in] algorithm the desired algorithm
		 * @param[in] section the section into which to look for the arguments
		 * @return a vector containing the arguments
		 */
		std::vector< std::pair<std::string, std::string> > getArgumentsForAlgorithm(const std::string& parname,
		                                const std::string& algorithm, const std::string& section) const;

		/**
		 * @brief Returns the metadata associated with the configured virtual stations
		 * @param date when to extract the virtual stations' metadata
		 * @param vecStation a vector of stationdata for the configured virtual stations
		 */
		size_t getVirtualStationsMeta(const Date& date, STATIONS_SET& vecStation);

		/**
		 * @brief Compute point measurements from grids following a given computing strategy
		 * @param i_date when to compute the virtual stations
		 * @param vecMeteo a vector of meteodata for the configured virtual stations
		 */
		size_t getVirtualMeteoData(const Date& i_date, METEO_SET& vecMeteo);
		
		/**
		 * @brief Push virtual points resampling into the provided tsmanager, so it appears as if they were coming from real data
		 * @details This call garantees thatif no point exactly matches the requested date, the nearest points around will be
		 * provided so the TimeSeriesManager can perform a temporal interpolation.
		 * @param i_date when to compute the virtual stations
		 * @param user_tsmanager TimeSeriesManager where to push the data
		 */
		void pushVirtualMeteoData(const Date& i_date, TimeSeriesManager &user_tsmanager);

		const std::string toString() const;

	private:
		static Config stripVirtualConfig(const Config& cfg);
		static void checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj);
		static void check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo);
		static std::set<std::string> getParameters(const Config& cfg);
		static std::vector<std::string> getAlgorithmsForParameter(const Config& cfg, const std::string& parname);

		size_t getVirtualStationsData(const Date& i_date, METEO_SET& vecMeteo);
		size_t getVirtualStationsFromGrid(const Date& i_date, METEO_SET& vecMeteo);
		void setAlgorithms();
		void initVirtualStations(const bool& adjust_coordinates);
		void initVirtualStationsAtAllGridPoints();

		const Config& cfg; ///< Reference to Config object, initialized during construction
		TimeSeriesManager *tsmanager; ///< Reference to TimeSeriesManager object, used for callbacks, initialized during construction
		GridsManager *gridsmanager; ///< Reference to GridsManager object, used for callbacks, initialized during construction
		GridBuffer grid_buffer;
		DEMObject internal_dem; ///< With virtual stations & resampling, we must have a DEM, so we keep an internal copy

		std::map< std::string, std::vector<InterpolationAlgorithm*> > mapAlgorithms; //per parameter interpolation algorithms

		std::vector<size_t> v_params; ///< Parameters for virtual stations
		std::vector<Coords> v_coords; ///< Coordinates for virtual stations
		std::vector<StationData> v_stations; ///< metadata for virtual stations
		unsigned int vstations_refresh_rate, vstations_refresh_offset; ///< when using virtual stations, how often should the data be spatially re-interpolated?
		RESAMPLING_STRATEGY resampling_strategy; ///< Should we perform resampling and with which strategy?
		
		bool algorithms_ready; ///< Have the algorithms objects been constructed?
		bool use_full_dem; ///< use full dem for point-wise spatial interpolations
};

} //end namespace

#endif
