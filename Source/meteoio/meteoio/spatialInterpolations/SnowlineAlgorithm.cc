/***********************************************************************************/
/*  Copyright 2020 ALPsolut.eu                                                     */
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

#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/spatialInterpolations/SnowlineAlgorithm.h>   

#include <fstream>
#include <limits>

namespace mio {

SnowlineAlgorithm::SnowlineAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs,
    const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
    GridsManager& i_gdm, Meteo2DInterpolator& i_mi) : 
    InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), gdm_(i_gdm), mi_(i_mi), base_alg_("IDW_LAPSE"),
    input_args_(vecArgs), where_("Interpolations2D::" + algo),
    assimilateFunction_(&SnowlineAlgorithm::assimilateCutoff), smoothing_(false),
    snowlines_(), snowlines_file_(),
    enforce_positive_rate_(false), calc_base_rate_(false), fallback_rate_(IOUtils::nodata),
    cutoff_val_(0.), band_height_(10.), band_no_(10), formula_(std::string()),
    has_warned_deduced_trend_(false), has_warned_fixed_trend_(false), verbose_(true)
{
	std::string algo_info( "CUTOFF" );
	for (size_t ii = 0; ii < vecArgs.size(); ii++) {
		if (vecArgs[ii].first == "BASE") {
			base_alg_ = IOUtils::strToUpper(vecArgs[ii].second);
		} else if (vecArgs[ii].first == "SNOWLINE") {
			aspect sl_aspect;
			IOUtils::parseArg(vecArgs[ii], where_, sl_aspect.sl_elevation);
			snowlines_.push_back(sl_aspect);
		} else if (vecArgs[ii].first == "SNOWLINEFILE") {
			snowlines_file_ = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "METHOD") {
			const std::string mode( IOUtils::strToUpper(vecArgs[ii].second) );
			if (mode == "CUTOFF")
				assimilateFunction_ = &SnowlineAlgorithm::assimilateCutoff;
			else if (mode == "BANDS")
				assimilateFunction_ = &SnowlineAlgorithm::assimilateBands;
			else if (mode == "FORMULA")
				assimilateFunction_ = &SnowlineAlgorithm::assimilateFormula;
			else
				throw InvalidArgumentException("Snowline assimilation mode \"" + mode +
				    "\" supplied for " + where_ + " not known.", AT);
			algo_info = mode;
		} else if (vecArgs[ii].first == "SMOOTHING") {
			IOUtils::parseArg(vecArgs[ii], where_, smoothing_);
		} else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where_, verbose_);
		} else if (vecArgs[ii].first == "SET") {
			IOUtils::parseArg(vecArgs[ii], where_, cutoff_val_);
		} else if (vecArgs[ii].first == "ENFORCE_POSITIVE_RATE") {
			IOUtils::parseArg(vecArgs[ii], where_, enforce_positive_rate_);
		} else if (vecArgs[ii].first == "CALC_BASE_RATE") {
			IOUtils::parseArg(vecArgs[ii], where_, calc_base_rate_);
		} else if (vecArgs[ii].first == "FALLBACK_RATE") {
			IOUtils::parseArg(vecArgs[ii], where_, fallback_rate_);
		/* args of method BANDS */
		} else if (vecArgs[ii].first == "BAND_HEIGHT") {
			IOUtils::parseArg(vecArgs[ii], where_, band_height_);
		} else if (vecArgs[ii].first == "BAND_NO") {
			IOUtils::parseArg(vecArgs[ii], where_, band_no_);
		/* args of method FORMULA */
		} else if (vecArgs[ii].first == "FORMULA") {
			formula_ = vecArgs[ii].second;
		}
	}
	info << "method: " << algo_info << ", ";
}

double SnowlineAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date; //this is always called before calculate() and sets the date!
	nrOfMeasurments = getData(date, param, vecData, vecMeta); //more initialization
	if (nrOfMeasurments == 0)
		return 0.0;
	if (!snowlines_.empty()) //TODO: possible to propagate from base alg?
		return 0.8;
	else
		return 0.7;
}

void SnowlineAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	getSnowlines();
	if (snowlines_.empty()) { //we already gave notice for this
		baseInterpol(IOUtils::nodata, dem, grid);
		return; 
	}

	std::vector<Grid2DObject> azi_grids;
	for (size_t ii = 0; ii < snowlines_.size(); ++ii) {
		Grid2DObject grid_copy(grid, IOUtils::nodata); //copy geo metadata
		const double snowline = snowlines_[ii].sl_elevation;
		baseInterpol(snowline, dem, grid_copy);
		(this->*assimilateFunction_)(snowline, dem, grid_copy);
		azi_grids.push_back(grid_copy);
	}
	grid = mergeSlopes(dem, azi_grids);


	has_warned_deduced_trend_ = false; //reset warning flags in case Algorithm lives on
	has_warned_fixed_trend_ = false;
}

void SnowlineAlgorithm::baseInterpol(const double& snowline, const DEMObject& dem, Grid2DObject& grid)
{
	std::string base_alg_fallback( base_alg_ ); //fallback algorithm if only 1 station is used
	if (nrOfMeasurments == 1) {
		base_alg_fallback = "AVG";
		msg("[W] Falling back to \"AVG\" for " + where_ + " (insufficient number of stations on " +
		    date.toString(Date::ISO_DATE) + ").");
	}

	//read and adjust parameters for base algorithm:
	std::vector< std::pair<std::string, std::string> > vecArgs( prepareBaseArgs(snowline, base_alg_fallback) );
	
	InterpolationAlgorithm* algorithm(
	    AlgorithmFactory::getAlgorithm(base_alg_fallback, mi_, vecArgs, tsmanager, gdm_, param) );
	//apply base interpolation algorithm to whole grid:
	algorithm->getQualityRating(date); //set date
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
	delete algorithm;
}

void SnowlineAlgorithm::assimilateCutoff(const double& snowline, const DEMObject& dem, Grid2DObject& grid)
{ //set everything below snowline elevation to fixed value
	for (size_t ii = 0; ii < grid.getNx(); ++ii) {
		for (size_t jj = 0; jj < grid.getNy(); ++jj) {
			if (dem(ii, jj) == IOUtils::nodata)
				continue;
			if (dem(ii, jj) < snowline)
				grid(ii, jj) = cutoff_val_;
		}
	}
}

void SnowlineAlgorithm::assimilateBands(const double& snowline, const DEMObject& dem, Grid2DObject& grid)
{ //multiply elevation bands above snowline with factors from 0 to 1
	for (size_t ii = 0; ii < grid.getNx(); ++ii) {
		for (size_t jj = 0; jj < grid.getNy(); ++jj) {
			if (dem(ii, jj) == IOUtils::nodata) {
				continue;
			} else if (dem(ii, jj) > snowline + band_no_ * band_height_) {
				continue;
			} else if (dem(ii, jj) < snowline) {
				grid(ii, jj) = cutoff_val_;
				continue;
			}
			for (unsigned int bb = 0; bb < band_no_; ++bb) { //bin DEM into bands
				if ( (dem(ii, jj) >= snowline + bb * band_height_) && (dem(ii, jj) < snowline + (bb + 1.) * band_height_) )
					grid(ii, jj) = grid(ii, jj) * bb / band_no_;
			}
		} //endfor jj
	} //endfor ii
}

void SnowlineAlgorithm::assimilateFormula(const double& snowline, const DEMObject& dem, Grid2DObject& grid)
{ //set to result of formula evaluated at grid points
	std::vector< std::pair<std::string, double> > substitutions; //fixed memory for tinyexpr

	static const size_t nr_sub = 3;
	static const std::string sub[nr_sub] = {"snowline", "altitude", "param"};
	for (size_t ii = 0; ii < nr_sub; ++ii) {
		std::pair<std::string, double> item(sub[ii], IOUtils::nodata);
		substitutions.push_back(item);
	}
	substitutions[0].second = snowline; //this is the same for all points

	te_variable *te_vars = new te_variable[substitutions.size()];
	initExpressionVars(substitutions, te_vars); //build te_variables from substitution vector
	te_expr *expr_formula = compileExpression(formula_, te_vars, nr_sub);

	for (size_t ii = 0; ii < grid.getNx(); ++ii) {
		for (size_t jj = 0; jj < grid.getNy(); ++jj) {
			if (dem(ii, jj) == IOUtils::nodata)
				continue;
			if (dem(ii, jj) < snowline) {
				grid(ii, jj) = cutoff_val_;
				continue;
			}

			substitutions[1].second = dem(ii, jj); //point-dependent substitutions
			substitutions[2].second = grid(ii, jj);
			grid(ii, jj) = te_eval(expr_formula);
		}
	}

	te_free(expr_formula);
	delete[] te_vars;
}

Grid2DObject SnowlineAlgorithm::mergeSlopes(const DEMObject& dem, const std::vector<Grid2DObject>& azi_grids)
{ //merge separate slope aspects back together

	DEMObject dem_copy(dem);
	dem_copy.setUpdatePpt(DEMObject::SLOPE);
	dem_copy.update(DEMObject::HORN);
	dem_copy.sanitize(); //set DEM to nodata at points where we have no slope/curvature
	Grid2DObject outgrid(dem_copy, IOUtils::nodata); //copy geo metadata
	Grid2DObject azi_classes(dem_copy); //for external visualization we would need the cellsize
	azi_classes.grid2D = dem_copy.azi; //container for binned azimuth classification

	std::vector<double> azi_thresholds;
	std::vector<double> azi_ids;
	azi_ids.push_back(static_cast<double>(snowlines_.size()) - 1.); //wrap back around North (0)
	for (size_t ii = 0; ii < snowlines_.size(); ++ii) {
		azi_thresholds.push_back(snowlines_[ii].deg_beg);
		azi_ids.push_back(static_cast<double>(ii));
	}
	azi_classes.binning(azi_thresholds, azi_ids);

	/*
	 * Now we have a matrix containing indices which correspond to a single slope aspect.
	 * Thus, this is a lookup table for which slope aspect matrix to choose the points from.
	 */

	for (size_t ii = 0; ii < outgrid.getNx(); ++ii) {
		for (size_t jj = 0; jj < outgrid.getNy(); ++jj) {
			if (dem_copy(ii, jj) != IOUtils::nodata) {
				const size_t idx = static_cast<size_t>(azi_classes(ii, jj));
				outgrid(ii, jj) = azi_grids.at(idx)(ii, jj);
			}
		}
	} //endfor ii

	if (smoothing_)
		throw InvalidArgumentException("Border smoothing not implemented yet in " + where_, AT);

	return outgrid;
}

void SnowlineAlgorithm::initExpressionVars(const std::vector< std::pair<std::string, double> >& substitutions, te_variable* vars) const
{ //build a substitutions expression for tinyexpr
	size_t cc = 0;
	for (std::vector< std::pair<std::string, double> >::const_iterator it = substitutions.begin();
	    it != substitutions.end(); ++it) {
		vars[cc].name = it->first.c_str();
		vars[cc].address = &it->second;
		vars[cc].type = 0;
		vars[cc].context = 0;
		cc++;
	}
}

te_expr* SnowlineAlgorithm::compileExpression(const std::string& expression, const te_variable* te_vars, const size_t& sz) const
{ //ready the lazy expressions (with syntax check)
	int te_err;
	te_expr *expr = te_compile(expression.c_str(), te_vars, static_cast<int>(sz), &te_err);
	if (!expr)
		throw InvalidFormatException("Arithmetic expression \"" + expression +
		        "\" could not be evaluated for " + where_ + "; parse error at " + IOUtils::toString(te_err), AT);
	return expr;
}

/**
 * @brief Read in snowline elevation information from a textfile.
 * @details Example file:
 * #Snowline elevation estimated from satellite data split up into slope aspects.
 * #Format: aspect_#,beg_azimuth end_azimuth snowline
 * #Snowline in meters, cloud coverage in percent, azimuth in degrees.
 * image_dates,2020-10-08T10:20:31
 * cloud_coverage,46.3
 * aspect_1,45 135 2130
 * aspect_2,135 225 2230
 * aspect_3,225 315 2130
 * aspect_4,315 45 2130
 * @return Vector with all parsed aspects.
 */
std::vector<SnowlineAlgorithm::aspect> SnowlineAlgorithm::readSnowlineFile()
{
	std::ifstream fin( snowlines_file_.c_str() );
	if (fin.fail())
		return std::vector<SnowlineAlgorithm::aspect>();
	
	std::vector<SnowlineAlgorithm::aspect> snowline_elevations;
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	try {
		do {
			aspect sl_aspect;
			std::string line;
			getline(fin, line, eoln);
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty())
				continue; //skip comments
			std::istringstream iss( line );

			std::string field_name; //each line has a field name followed by a comma...
			std::string field_value; //... and then an arbitrary number of value fields.
			if (std::getline(iss, field_name, ','))
				std::getline(iss, field_value, ',');
			else
				throw InvalidFormatException(
				    "The snowline elevations file is ill-formatted (no field/value separator in line).", AT);

			if (field_name == "image_dates") {
				//nothing yet
			} else if (field_name == "cloud_coverage") {
				//nothing yet
			} else if (field_name.substr(0, 7) == "aspect_") { //aspect_#,azi1 azi2 snowline
				std::istringstream ias( field_value );
				ias.setf(std::ios::fixed);
				ias.precision(std::numeric_limits<double>::digits10);
				ias >> std::skipws >> sl_aspect.deg_beg;
				ias >> std::skipws >> sl_aspect.deg_end;
				ias >> std::skipws >> sl_aspect.sl_elevation;
				if (!ias)
					throw InvalidFormatException(
					    "The snowline elevations file is ill-formatted (could not parse all 3 aspect parameters).", AT);
				snowline_elevations.push_back(sl_aspect);
			}
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&) {
		if (fin.is_open())
			fin.close();
		throw;
	}
	return snowline_elevations;
}

void SnowlineAlgorithm::getSnowlines()
{ //read snowline information according to input method and validate it
	if (!snowlines_file_.empty()) {
		if (!snowlines_.empty())
			msg("[i] Ignoring additional SNOWLINE argument since SNOWFILE is given for " + where_ + ".");
		snowlines_ = readSnowlineFile();
		if (snowlines_.empty()) {
			msg("[W] No valid snowline elevation can be read from SNOWFILE \"" + snowlines_file_ + "\" for " +
			    where_ + ". Continuing without...");
			return;
		}
	} else {
		if (snowlines_.empty()) {
			msg("[W] No numeric value found for SNOWLINE, no SNOWFILE provided either for " +
			    where_ + ". Continuing without...");
			return;
		}
	}

	if (!isSorted(snowlines_))
		throw InvalidFormatException("The snowline elevations file is ill-formatted (azimuth values must be in ascending order).", AT);
	if (snowlines_.at(0).deg_beg != snowlines_.at(snowlines_.size() - 1).deg_end)
		throw InvalidFormatException("The snowline elevations file is ill-formatted (the start azimuth must be the end azimuth).", AT);
	for (std::vector<SnowlineAlgorithm::aspect>::iterator it = snowlines_.begin(); it != snowlines_.end(); ++it) {
		    if ((*it).deg_beg > 360. || (*it).deg_end > 360. || (*it).deg_beg < 0. || (*it).deg_end < 0.)
			throw InvalidArgumentException("The snowline elevations file is ill-formatted (degree value out of range).", AT);
	}
}

double SnowlineAlgorithm::probeTrend()
{
	//calculate model parameters - we don't want to change data yet:
	Trend trend(input_args_, algo, param);
	std::vector<double> vecData_copy( vecData );
	std::vector<StationData> vecMeta_copy( vecMeta );
	trend.detrend(vecMeta_copy, vecData_copy);
	
	/*
	 * Regression parameters are available only in low-level function calls when
	 * trending/detrending is being performed.
	 * But: they are propagated in the info string, so we parse this here.
	 */
	const std::string str_tag( "Model parameters:" );
	std::string fit_info( trend.toString() );
	const size_t beg = fit_info.find(str_tag) + str_tag.length();
	const size_t end = fit_info.find("\n", beg);
	fit_info = IOUtils::trim(fit_info.substr(beg, end - beg));
	const size_t space = fit_info.find(" ");
	double slope;
	const bool success = IOUtils::convertString(slope, fit_info.substr(space + 1));
	if (!success)
		throw ConversionFailedException(where_ + " could not extract internal fit parameters of " + fit_info, AT);
	return slope;
}

double SnowlineAlgorithm::calculateRate(const double& snowline)
{ //calculate lapse rate from snowline elevation and highest available station
	double max_alt = -999.;
	double HS = IOUtils::nodata;

	for (size_t ii = 0.; ii < vecData.size(); ++ii) {
		if (vecMeta[ii].position.getAltitude() > max_alt) {
			max_alt = vecMeta[ii].position.getAltitude();
			HS = vecData[ii];
		}
	}
	return HS / (max_alt - snowline);
}

std::vector< std::pair<std::string, std::string> > SnowlineAlgorithm::prepareBaseArgs(const double& snowline, const std::string& base_alg_fallback)
{ //inject calculated lapse rate into base algorithm arguments
	std::vector< std::pair<std::string, std::string> > vecArgs( mi_.getArgumentsForAlgorithm(param, base_alg_fallback, "Interpolations2D") );
	if (enforce_positive_rate_ || calc_base_rate_) {
		const double fit_slope = probeTrend();
		if (snowline != IOUtils::nodata && (fit_slope > 0. || calc_base_rate_)) {
			const std::string reason( (fit_slope > 0.? "Reverse trend" : "Trend") );
			if (!has_warned_deduced_trend_) {
				msg("[i] " + reason + " in data is being substituted with rate deduced from snowline elevation.");
				has_warned_deduced_trend_ = true;
			}
			if (snowline != IOUtils::nodata) {
				for(std::vector< std::pair<std::string, std::string> >::iterator it = vecArgs.begin(); it != vecArgs.end(); ++it) {
					if ((*it).first == "RATE" || (*it).first == "FRAC")
						vecArgs.erase(it);
				}
				const double snowline_rate = calculateRate(snowline);
				vecArgs.push_back(std::pair<std::string, std::string>( "RATE", IOUtils::toString(snowline_rate) ));
				vecArgs.push_back(std::pair<std::string, std::string>( "FRAC", "FALSE" ));
			}
		} else if (snowline == IOUtils::nodata && fit_slope > 0. && fallback_rate_ != IOUtils::nodata) {
			//reversed lapse rate but no snowline: switch to fixed fallback rate
			for(std::vector< std::pair<std::string, std::string> >::iterator it = vecArgs.begin(); it != vecArgs.end(); ++it) {
				if ((*it).first == "RATE" || (*it).first == "FRAC")
					vecArgs.erase(it);
			}
			vecArgs.push_back(std::pair<std::string, std::string>( "RATE", IOUtils::toString(fallback_rate_) ));
			vecArgs.push_back(std::pair<std::string, std::string>( "FRAC", "FALSE" ));
			if (!has_warned_fixed_trend_) {
				msg("[i] Reverse trend in data is being substituted with a fixed fallback rate.");
				has_warned_fixed_trend_ = true;
			}
		}
	}
	return(vecArgs);
}

void SnowlineAlgorithm::msg(const std::string& message)
{
	if (verbose_)
		std::cerr << message << std::endl;
}

} //namespace
