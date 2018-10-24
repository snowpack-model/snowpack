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
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>

#include <meteoio/meteoResampling/Accumulate.h>
#include <meteoio/meteoResampling/DailyAverageResampling.h>
#include <meteoio/meteoResampling/DailySolarResampling.h>
#include <meteoio/meteoResampling/LinearResampling.h>
#include <meteoio/meteoResampling/NearestNeighbour.h>
#include <meteoio/meteoResampling/NoResampling.h>
#include <meteoio/meteoResampling/SolarResampling.h>

namespace mio {
/**
 * @page resampling Resampling overview
 * The resampling infrastructure is described in ResamplingAlgorithms (for its API).
 * The goal of this page is to give an overview of the available resampling algorithms and their usage.
 *
 * @section resampling_section Resampling section
 * The resampling is specified for each parameter in the [Interpol1D] section. This section contains
 * a list of the various meteo parameters with their associated choice of resampling algorithm and
 * optional parameters. If a meteo parameter is not listed in this section, a linear resampling would be
 * assumed. An example of such section is given below:
 * @code
 * [Interpolations1D]
 * WINDOW_SIZE     = 86400
 * TA::resample    = linear
 *
 * RH::resample            = linear
 * RH::linear::window_size = 172800
 *
 * VW::resample             = nearest
 * VW::nearest::extrapolate = true
 *
 * PSUM::resample           = accumulate
 * PSUM::accumulate::period = 3600
 * @endcode
 *
 * Most of the resampling algorithms allow you to define per-meteo parameter and per-algorithm the WINDOW_SIZE. Otherwise, the section's WINDOW_SIZE is
 * used as default window size. This represents the biggest gap that can be interpolated (in seconds). Therefore if two valid points are less than
 * WINDOW_SIZE seconds apart, points in between will be interpolated. If they are further apart, all points in between will remain IOUtils::nodata.
 * If using the "extrapolate" optional argument, points at WINDOW_SIZE distance of only one valid point will be extrapolated, otherwise they will remain
 * IOUtils::nodata. Please keep in mind that allowing extrapolated values can lead to grossly out of range data: using the slope
 * between two hourly measurements to extrapolate a point 10 days ahead is obviously risky!
 *
 * By default, WINDOW_SIZE is set to 2 days. This key has a <b>potentially large impact on run time/performance</b>.
 *
 * @section algorithms_available Available Resampling Algorithms
 * Several algorithms for the resampling are implemented:
 * - none: do not perform resampling, see NoResampling
 * - nearest:  nearest neighbor data resampling, see NearestNeighbour
 * - linear: linear data resampling, see LinearResampling
 * - accumulate: data re-accumulation as suitable for precipitations, see Accumulate
 * - solar: resample solar radiation by interpolating an atmospheric loss factor, see Solar
 * - daily_solar: generate solar radiation (ISWR or RSWR) from daily sums, see Daily_solar
 * - daily_avg: generate a sinusoidal variation around the measurement taken as daily average and of a given amplitude, see DailyAverage
 *
 * By default a linear resampling will be performed. It is possible to turn off all resampling by setting the *Enable_Resampling* key
 * to *false* in the [Interpolations1D] section.
 */

ResamplingAlgorithms* ResamplingAlgorithmsFactory::getAlgorithm(const std::string& i_algoname, const std::string& parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string algoname( IOUtils::strToUpper(i_algoname) );

	if (algoname == "NONE" || algoname == "NO") {
		return new NoResampling(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "LINEAR") {
		return new LinearResampling(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "NEAREST") {
		return new NearestNeighbour(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "ACCUMULATE") {
		return new Accumulate(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "SOLAR") {
		return new Solar(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "DAILY_SOLAR") {
		return new Daily_solar(algoname, parname, dflt_window_size, vecArgs);
	} else if (algoname == "DAILY_AVG") {
		return new DailyAverage(algoname, parname, dflt_window_size, vecArgs);
	} else {
		throw IOException("The resampling algorithm '"+algoname+"' is not implemented" , AT);
	}
}

const double ResamplingAlgorithms::soil_albedo = .23; ///< grass albedo
const double ResamplingAlgorithms::snow_albedo = .85; ///< snow albedo
const double ResamplingAlgorithms::snow_thresh = .1; ///< if snow height greater than this threshold -> snow albedo

//compute the partial accumulation at the std::left of curr_date within a sampling interval
double ResamplingAlgorithms::partialAccumulateAtLeft(const std::vector<MeteoData>& vecM, const size_t& paramindex,
                                                     const size_t& pos, const Date& curr_date)
{
	const size_t end = pos+1;
	if (end>=vecM.size()) return IOUtils::nodata; //reaching the end of the input vector

	const double valend = vecM[end](paramindex);
	if (valend == IOUtils::nodata) return IOUtils::nodata;

	const double jul1 = vecM[pos].date.getJulian(true);
	const double jul2 = vecM[end].date.getJulian(true);
	const double jul_curr = curr_date.getJulian(true);

	if (jul_curr>jul2) return IOUtils::nodata; //jul1 and jul2 should  frame jul_curr, if not some data points are missing
	const double left_accumulation = linearInterpolation(jul1, 0., jul2, valend, jul_curr);
	return left_accumulation;
}

//compute the partial accumulation at the std::right of curr_date within a sampling interval
double ResamplingAlgorithms::partialAccumulateAtRight(const std::vector<MeteoData>& vecM, const size_t& paramindex,
                                                      const size_t& pos, const Date& curr_date)
{
	const size_t end = pos+1;
	if (end>=vecM.size()) return IOUtils::nodata; //reaching the end of the input vector

	const double valend = vecM[end](paramindex);
	if (valend == IOUtils::nodata) return IOUtils::nodata;

	const double jul1 = vecM[pos].date.getJulian(true);
	const double jul2 = vecM[end].date.getJulian(true);
	const double jul_curr = curr_date.getJulian(true);

	if (jul_curr<jul1) return IOUtils::nodata; //jul1 and jul2 should  frame jul_curr, if not some data points are missing
	const double left_accumulation = linearInterpolation(jul1, 0., jul2, valend, curr_date.getJulian(true));

	return valend - left_accumulation;
}

/**
 * @brief This function returns the last and next valid points around a given position
 * @param pos current position (index)
 * @param paramindex meteo parameter to use
 * @param vecM vector of MeteoData
 * @param resampling_date date to resample
 * @param i_window_size size of the search window
 * @param indexP1 index of point before the current position (IOUtils::npos if none could be found)
 * @param indexP2 index of point after the current position (IOUtils::npos if none could be found)
 */
void ResamplingAlgorithms::getNearestValidPts(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, size_t& indexP1, size_t& indexP2)
{
	indexP1=IOUtils::npos;
	indexP2=IOUtils::npos;

	const Date dateStart( resampling_date - i_window_size );
	for (size_t ii=pos; ii-- >0; ) { //because idx gets decremented right away
		if (vecM[ii].date < dateStart) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP1 = ii;
			break;
		}
	}

	//make sure the search window remains i_window_size
	const Date dateEnd = (indexP1 != IOUtils::npos)? vecM[indexP1].date+i_window_size : resampling_date+i_window_size;
	for (size_t ii=pos; ii<vecM.size(); ++ii) {
		if (vecM[ii].date > dateEnd) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP2 = ii;
			break;
		}
	}
}

/**
 * @brief This function solves the equation y = ax + b for two given points and returns y for a given x
 * @param x1 x-coordinate of first point
 * @param y1 y-coordinate of first point
 * @param x2 x-coordinate of second point
 * @param y2 y-coordinate of second point
 * @param x x-coordinate of desired point
 * @return y-coordinate of desired point
 */
double ResamplingAlgorithms::linearInterpolation(const double& x1, const double& y1,
                                       const double& x2, const double& y2, const double& x)
{
	if (x1 == x2) throw IOException("Attempted division by zero", AT);

	//Solving y = ax + b
	const double a = (y2 - y1) / (x2 - x1);
	const double b = y2 - a*x2;
	return (a*x + b);
}

/**
 * @brief For a given date, find the start of the day, considering that for midnight we return the day before!
 * (as is necessary for daily averages, sums, etc that can be provided at midnight for the day before)
 * @param resampling_date current date
 * @return start of the day or start of the day before in case of midnight
 */
Date ResamplingAlgorithms::getDailyStart(const Date& resampling_date)
{
	Date Start( resampling_date );
	Start.rnd(24*3600, Date::DOWN);
	if (Start==resampling_date) //if resampling_date=midnight GMT, the rounding lands on the exact same date
		Start -= 1.;
	
	return Start;
}

/**
 * @brief Find a unique value in a given time interval. 
 * This is useful for retrieving a unique daily average, daily sum, etc
 * @param vecM meteo data time series for the station
 * @param paramindex index of the meteo parameter to process
 * @param pos index of the current time step
 * @param intervalStart start of the interval within which we will look for a valid value (often, start of the day)
 * @param intervalEnd end of the interval within which we will look for a valid value (often, end of the day)
 * @return index of the parameter's valid value
 */
size_t ResamplingAlgorithms::getDailyValue(const std::vector<MeteoData>& vecM, const size_t& paramindex, size_t pos, const Date& intervalStart, const Date& intervalEnd)
{
	size_t indexP1=IOUtils::npos;
	size_t indexP2=IOUtils::npos;
	
	//if pos was not properly pre-positioned, do it
	if (vecM[pos].date<intervalStart) {
		for (; pos<vecM.size(); pos++)
			if (vecM[pos].date>=intervalStart) break;
	}
	if (vecM[pos].date>intervalEnd) {
		for (size_t ii = pos; ii-- >0; ) {
			pos=ii; //because ii gets corrupted at the final iteration if going all the way down
			if (vecM[pos].date<=intervalEnd) break;
		}
	}

	//look for daily sum before the current point
	for (size_t ii=pos; ii-- >0; ) {
		if (vecM[ii].date < intervalStart) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP1 = ii;
			break;
		}
	}

	//look for daily sum after the current point
	for (size_t ii=pos; ii<vecM.size(); ++ii) {
		if (vecM[ii].date > intervalEnd) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP2 = ii;
			break;
		}
	}
	
	if (indexP1!=IOUtils::npos && indexP2!=IOUtils::npos) {
		//if the data was provided at 00:00, a sum for each day has been found
		//we only keep the later one as being the sum of the past day
		int hour1, min1;
		vecM[indexP1].date.getTime(hour1, min1);
		if (hour1==0 && min1==0)
			indexP1=IOUtils::npos;
		else { //otherwise, this means multiple daily sums have been found for the same day
			const std::string msg = "More than one daily sum of solar radiation found between "+intervalStart.toString(Date::ISO)+" and "+intervalEnd.toString(Date::ISO);
			throw IOException(msg, AT);
		}
	}

	if (indexP1!=IOUtils::npos) return indexP1;
	if (indexP2!=IOUtils::npos) return indexP2;
	return IOUtils::npos;
}

} //namespace

