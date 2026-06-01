/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file Meteo.h
 * @author Charles Fierz
 * @version 9.x
 * @date -
 */

#ifndef METEO_H
#define METEO_H

#include <meteoio/MeteoIO.h>

#include <snowpack/SnowpackConfig.h>
#include <snowpack/snowpackCore/Canopy.h>
#include <snowpack/DataClasses.h>

class Meteo {
	public:
		//except Richardson and Neutral, all are standard MO iterations with the specified stability correction. In unstable conditions, they use Paulson and Stearns & Weidner, 1993, for scalars
		typedef enum {
			RICHARDSON,  ///< Simplified Richardson number stability correction
			NEUTRAL,  ///< Assume neutral stratification
			MO_LOG_LINEAR, ///< Simple log-linear model
			MO_HOLTSLAG, ///< Holtslag and DeBruin (1988) prepared from Ed Andreas
			MO_STEARNS, ///< Stearns C. and Weidner G., <i>"sensible and latent heat flux estimates in antarctica"</i>, Antarctic meteorology and climatology: studies based on automatic weather stations, Antarctic Research Series, <b>61</b>, pp 190--138, 1993
			MO_MICHLMAYR, ///< Stearns & Weidner, 1993 modified by Michlmayr, 2008
			MO_SCHLOEGL_UNI, ///< Schloegl univariate, see Schlögl et al. <i>"How do stability corrections perfom over snow in the stable boundary layer?"</i>, Boundary-Layer Meteorol., in review, 2016
			MO_SCHLOEGL_MULTI, ///< Schloegl multivariate without offset
			MO_SCHLOEGL_MULTI_OFFSET ///<Schloegl multivariate with offset
		} ATM_STABILITY;

		Meteo(const SnowpackConfig& i_cfg);

		void compMeteo(CurrentMeteo &Mdata, SnowStation &Xdata, const bool runCanopyModel);
		void compMeteo(CurrentMeteo &Mdata, SnowStation &Xdata, const bool runCanopyModel, const bool i_adjust_height_of_wind_value);
		static double windspeedProfile(const CurrentMeteo& Mdata, const double& target_z, const double& source_vw = -1.);
		void setStability(const ATM_STABILITY& i_stability);
		static ATM_STABILITY getStability(const std::string& stability_model);
		ATM_STABILITY getStability() const;

 	private:
		void MicroMet(const SnowStation& Xdata, CurrentMeteo& Mdata, const double& roughness_length, const bool& adjust_VW_height=true);
		static double getParameterAverage(mio::IOManager& io, const mio::MeteoData::Parameters& param,
		                                  const mio::Date& current_date, const int& time_span, const int& increment);
		static void RichardsonStability(const double& ta_v, const double& t_surf_v, const double& zref,
		                                const double& vw, const double& z_ratio, double &ustar, double &psi_s);
		static void MOStability(const ATM_STABILITY& use_stability, const double& ta_v, const double& t_surf_v, const double& t_surf,
		                                       const double& zref, const double& vw, const double& z_ratio, double &ustar, double &psi_s, double &psi_m);
		double compZ0(const std::string& model, const SnowStation& Xdata, CurrentMeteo& Mdata);
		double parametrizeSnowZ0(const SnowStation& Xdata, CurrentMeteo& Mdata) const;
		
		Canopy canopy;
		std::string roughness_length_parametrization;
		const double roughness_length_dflt; ///< Initial estimate of the snow roughness length for the site; will be adjusted iteratively, default value and operational mode: 0.002 m
		double height_of_wind_value; ///< Define the heights of the meteo measurements above ground (m). Required for surface energy exchange computation and for drifting and blowing snow.
		ATM_STABILITY stability; ///< Atmospheric stability correction model
		bool research_mode; ///< Either research mode or operational mode (this is some kind of "meta setting" that toggle some defaults values)
		bool useCanopyModel; ///< Defines whether the canopy model is used. OUT_CANOPY must also be set to dump canopy parameters to file; see Constants_local.h
		bool adjust_height_of_wind_value; ///< Adjust the wind speed measurement height above the surface as snow height increases
};

#endif //END of Meteo.h
