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
#include <snowpack/libsnowpack.h>
#include <meteoio/MeteoIO.h>

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

#ifdef _MSC_VER
/*
 This software contains code under BSD license (namely, getopt for Visual C++).
 Therefore, this product includes software developed by the University of
 California, Berkeley and its contributors when compiling with Visual C++.
 */
#include "getopt.h"
#else
//#include <unistd.h> //for getopt
#include <getopt.h> //for getopt_long
#endif

using namespace std;
using namespace mio;

#ifdef DEBUG_ARITHM
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>
#endif


/**
 * @class Slope a C. Fierz class ;-)
 */
class Slope {

 	public:
		Slope(const mio::Config& cfg);

		double prevailing_wind_dir;
		unsigned int nSlopes;
		unsigned int mainStation;  ///< main station, flat field or slope
		unsigned int sector;  ///< main station (0) or current slope sector (1:nSlopes)
		unsigned int first;   ///< first virtual slope station in computing sequence
		unsigned int luv;
		unsigned int lee;
		bool north, south;
		bool snow_erosion, mainStationDriftIndex;
		bool snow_redistribution, luvDriftIndex;

		unsigned int getSectorDir(const double& dir_or_expo) const;
		void setSlope(const unsigned int slope_sequence,
		              vector<SnowStation>& vecXdata, double& wind_dir);

 	private:
		double sector_width;  ///< width of slope sector: 360./std::max((unsigned)1, nSlopes-1) deg
};



/************************************************************
 * static section                                           *
 ************************************************************/

//Global variables in this file:
static string cfgfile = "io.ini";
static string mode = "RESEARCH";
static mio::Date dateBegin, dateEnd;
static vector<string> vecStationIDs;

/// @brief Main control parameters
struct MainControl {
		size_t nStep;        ///< Time step number
		size_t nAvg;    ///< Number of calculation time steps to average fluxes etc.
		size_t HzStep;  ///< Hazard step number (should be half of nStep in operational mode)
		bool TsDump;       ///< Flag for time series dump
		bool HzDump;       ///< Calculation of hazard information will be performed
		bool PrDump;       ///< Flag for profile dump
		bool XdataDump;    ///< Backup of Xdata will be performed
		bool sdbDump;      ///< Dump to data base if required in operational mode
		bool resFirstDump;  ///< Flag to dump initial state of snowpack
};

/************************************************************
 * non-static section                                       *
 ************************************************************/

Slope::Slope(const mio::Config& cfg)
		: prevailing_wind_dir(0.),
		  nSlopes(0),
		  mainStation(0),
		  sector(0),
		  first(1),
		  luv(0),
		  lee(0),
		  north(false),
		  south(false),
		  snow_erosion(false),
		  mainStationDriftIndex(false),
		  snow_redistribution(false),
		  luvDriftIndex(false),
		  sector_width(0) {
	cfg.getValue("NUMBER_SLOPES", "SnowpackAdvanced", nSlopes);
	cfg.getValue("SNOW_EROSION", "SnowpackAdvanced", snow_erosion);
	stringstream ss;
	ss << nSlopes;
	cfg.getValue("SNOW_REDISTRIBUTION", "SnowpackAdvanced", snow_redistribution);
	if (snow_redistribution && !(nSlopes > 1 && nSlopes % 2 == 1))
		throw mio::IOException(
		    "Please set NUMBER_SLOPES to 3, 5, 7, or 9 with SNOW_REDISTRIBUTION set! (nSlopes="
		        + ss.str() + ")",
		    AT);
	cfg.getValue("PREVAILING_WIND_DIR", "SnowpackAdvanced", prevailing_wind_dir,
	             mio::IOUtils::nothrow);
	sector_width = 360.
	    / static_cast<double>(std::max((unsigned) 1, nSlopes - 1));
}

/**
 * @brief Determine either direction of blowing wind or slope exposition.
 * NOTE that station slope.first always corresponds to the prevailing wind direction
 * @param dir_or_expo direction of wind or exposition
 **/
unsigned int Slope::getSectorDir(const double& dir_or_expo) const {
	double dir = dir_or_expo;
	if (dir > 360.)
		dir -= 360.;
	else if (dir < 0.)
		dir += 360.;
	const unsigned int sectorDir = (unsigned int) ((floor(
	    (dir + 0.5 * sector_width) / sector_width)) + 1);
	if (sectorDir >= nSlopes)
		return 1;
	else
		return sectorDir;
}

/**
 * @brief Set slope variables
 * @param slope_sequence computation sequence for slopes
 * @param vecXdata
 * @param wind_dir direction of wind
 **/
void Slope::setSlope(const unsigned int slope_sequence,
                     vector<SnowStation>& vecXdata, double& wind_dir) {
	mainStationDriftIndex = false;
	luvDriftIndex = false;
	switch (slope_sequence) {
		case 0:
			for (size_t kk = 0; kk < nSlopes; kk++) {
				vecXdata[kk].windward = false;
				vecXdata[kk].rho_hn = 0.;
				vecXdata[kk].hn = 0.;
			}
			if (nSlopes > 1) {
				luv = getSectorDir(wind_dir - prevailing_wind_dir);
				vecXdata[luv].windward = true;
				lee = (luv + nSlopes / 2) % (nSlopes - 1);
				if (lee == 0)
					lee = nSlopes - 1;
			} else {
				//requesting slope 0 of 0 expositions
				luv = lee = 0;
			}
			sector = mainStation;
			mainStationDriftIndex = ((nSlopes == 1) && snow_erosion);
			break;
		case 1:
			sector = luv;
			luvDriftIndex = snow_redistribution;
			break;
		default:
			sector++;
			if (sector == nSlopes)
				sector = 1;
	}
	north = (vecXdata[sector].meta.getSlopeAngle() > 0.
	    && vecXdata[sector].meta.getAzimuth() == 0.);
	south = (vecXdata[sector].meta.getSlopeAngle() > 0.
	    && vecXdata[sector].meta.getAzimuth() == 180.);
}


inline void Version() {
#ifdef _MSC_VER
	cout << "This version of Snowpack uses a BSD-licensed port of getopt for Visual C++. \n"
	<< "It therefore includes software developed by the University of "
	<< "California, Berkeley and its contributors." << endl;
#endif
	cout << "Snowpack version " << SN_VERSION << " compiled on " << __DATE__
	     << " " << __TIME__ << "\n" << "\tLibsnowpack "
	     << snowpack::getLibVersion() << "\n" << "\tMeteoIO "
	     << mio::getLibVersion() << endl;
}

inline void Usage(const string& programname) {
	Version();

	cout << "Usage: " << programname << endl
	     << "\t[-b, --begindate=YYYY-MM-DDTHH:MM] (e.g.:2007-08-11T09:00)\n"
	     << "\t[-e, --enddate=YYYY-MM-DDTHH:MM] (e.g.:2008-08-11T09:00 or NOW)\n"
	     << "\t[-c, --config=<ini file>] (e.g. io.ini)\n"
	     << "\t[-m, --mode=<operational or research>] (default: research)\n"
	     << "\t[-s, --stations=<comma delimited stationnames>] (e.g. DAV2,WFJ2)\n"
	     << "\t[-v, --version] Print the version number\n"
	     << "\t[-h, --help] Print help message and version information\n\n";
	cout
	    << "\tPlease not that the operational mode should only be used within SLF\n";

	cout << "Example: " << programname << " -c io.ini -e 1996-06-17T00:00\n\n";
}

inline void parseCmdLine(int argc, char **argv, string& begin_date_str,
                         string& end_date_str) {
	int longindex = 0, opt = -1;
	bool setEnd = false;

	struct option long_options[] = { { "begindate", required_argument, 0, 'b' }, {
	    "enddate", required_argument, 0, 'e' }, { "mode", required_argument, 0,
	    'm' }, { "config", required_argument, 0, 'c' }, { "stations",
	    required_argument, 0, 's' }, { "version", no_argument, 0, 'v' }, { "help",
	    no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

	if (argc == 1) {  //no arguments provided
		Usage(string(argv[0]));
		exit(1);
	}

	while ((opt = getopt_long(argc, argv, ":b:e:m:c:s:v:h", long_options,
	                          &longindex)) != -1) {
		switch (opt) {
			case 0:
				break;
			case 'b': {
				begin_date_str = string(optarg);  //we don't know yet the time zone, conversion will be done later
				break;
			}
			case 'e': {
				end_date_str = string(optarg);  //we don't know yet the time zone, conversion will be done later
				setEnd = true;
				break;
			}
			case 'm':
				mode = string(optarg);
				mio::IOUtils::toUpper(mode);
				if (!(mode == "RESEARCH" || mode == "OPERATIONAL")) {
					cerr << endl << "[E] Command line option '-" << char(opt)
					     << "' requires 'research' or 'operational' as operand\n";
					Usage(string(argv[0]));
					exit(1);
				}
				break;
			case 'c':
				cfgfile = string(optarg);
				break;
			case 's':
				mio::IOUtils::readLineToVec(string(optarg), vecStationIDs, ',');
				break;
			case ':':  //operand missing
				cerr << endl << "[E] Command line option '-" << char(opt)
				     << "' requires an operand\n";
				Usage(string(argv[0]));
				exit(1);
			case 'v':
				Version();
				exit(0);
			case 'h':
				Usage(string(argv[0]));
				exit(0);
			case '?':
				cerr << endl << "[E] Unknown argument detected\n";
				Usage(string(argv[0]));
				exit(1);
			default:
				cerr << endl << "[E] getopt returned character code " << opt << "\n";
				Usage(string(argv[0]));
				exit(1);
		}
	}

	if (!setEnd) {
		cerr << endl << "[E] You must specify an enddate for the simulation!\n";
		Usage(string(argv[0]));
		exit(1);
	}
}

inline void editMeteoData(mio::MeteoData& md, const string& variant,
                          const double& thresh_rain) {  //HACK: these should be handled by DataGenerators
	if (md(MeteoData::PSUM_PH) == IOUtils::nodata) {
		const double ta = md(MeteoData::TA);
		if (ta != IOUtils::nodata)
			md(MeteoData::PSUM_PH) = (ta >= IOUtils::C_TO_K(thresh_rain)) ? 1. : 0.;  //fallback: simple temp threshold
	}

	//Add the atmospheric emissivity as a parameter
	if (!md.param_exists("EA")) {
		md.addParameter("EA");
		md("EA") = SnLaws::AirEmissivity(md, variant);
	}
	// Snow stations without separate wind station use their own wind for local drifting and blowing snow
	if (!md.param_exists("VW_DRIFT"))
		md.addParameter("VW_DRIFT");
	if (md("VW_DRIFT") == mio::IOUtils::nodata)
		md("VW_DRIFT") = md(MeteoData::VW);
	if (!md.param_exists("DW_DRIFT"))
		md.addParameter("DW_DRIFT");
	if (md("DW_DRIFT") == mio::IOUtils::nodata)
		md("DW_DRIFT") = md(MeteoData::DW);
}

// Return true if snowpack can compute the next timestep, else false
inline bool validMeteoData(const mio::MeteoData& md, const string& StationName,
                           const string& variant,
                           const bool& enforce_snow_height,
                           const bool& advective_heat, const bool& soil_flux,
                           const unsigned int& nslopes) {
	bool miss_ta = false, miss_tsg = false, miss_rh = false, miss_precip = false,
	    miss_splitting = false, miss_hs = false;
	bool miss_rad = false, miss_ea = false, miss_wind = false, miss_drift = false,
	    miss_adv = false;

	if (md(MeteoData::TA) == mio::IOUtils::nodata)
		miss_ta = true;
	if (soil_flux == false && md(MeteoData::TSG) == mio::IOUtils::nodata)
		miss_tsg = true;
	if (md(MeteoData::RH) == mio::IOUtils::nodata)
		miss_rh = true;
	if ((variant != "ANTARCTICA")
	    && ((md(MeteoData::ISWR) == mio::IOUtils::nodata)
	        && (md(MeteoData::RSWR) == mio::IOUtils::nodata)))
		miss_rad = true;
	if (enforce_snow_height && (md(MeteoData::HS) == mio::IOUtils::nodata))
		miss_hs = true;
	if (!enforce_snow_height && (md(MeteoData::PSUM) == mio::IOUtils::nodata))
		miss_precip = true;
	if (!enforce_snow_height && (md(MeteoData::PSUM_PH) == mio::IOUtils::nodata))
		miss_splitting = true;
	if (md("EA") == mio::IOUtils::nodata)
		miss_ea = true;
	if (md(MeteoData::VW) == mio::IOUtils::nodata)
		miss_wind = true;
	if (nslopes > 1
	    && (md("DW_DRIFT") == mio::IOUtils::nodata
	        || md("VW_DRIFT") == mio::IOUtils::nodata))
		miss_drift = true;
	if (advective_heat && md("ADV_HEAT") == mio::IOUtils::nodata)
		miss_adv = true;

	if (miss_ta || miss_tsg || miss_rh || miss_rad || miss_precip
	    || miss_splitting || miss_hs || miss_ea || miss_wind || miss_drift
	    || miss_adv) {
		mio::Date now;
		now.setFromSys();
		cerr << "[E] [" << now.toString(mio::Date::ISO) << "] ";
		cerr << StationName << " missing { ";
		if (miss_ta)
			cerr << "TA ";
		if (miss_tsg)
			cerr << "TSG ";
		if (miss_rh)
			cerr << "RH ";
		if (miss_rad)
			cerr << "sw_radiation ";
		if (miss_hs)
			cerr << "HS ";
		if (miss_precip)
			cerr << "precipitation ";
		if (miss_splitting)
			cerr << "precip_splitting ";
		if (miss_ea)
			cerr << "lw_radiation ";
		if (miss_wind)
			cerr << "VW ";
		if (miss_drift)
			cerr << "drift ";
		if (miss_adv)
			cerr << "adv_heat ";
		cerr << "} on " << md.date.toString(mio::Date::ISO) << "\n";
		return false;
	}
	return true;
}

inline void copyMeteoData(const mio::MeteoData& md, CurrentMeteo& Mdata,
                          const double prevailing_wind_dir,
                          const double wind_scaling_factor) {
	Mdata.date = Date::rnd(md.date, 1);
	Mdata.ta = md(MeteoData::TA);
	Mdata.rh = md(MeteoData::RH);
	if (md.param_exists("RH_AVG"))
		Mdata.rh_avg = md("RH_AVG");
	Mdata.vw = md(MeteoData::VW);
	Mdata.dw = md(MeteoData::DW);
	Mdata.vw_max = md(MeteoData::VW_MAX);
	if (md.param_exists("VW_AVG"))
		Mdata.vw_avg = md("VW_AVG");

	Mdata.vw_drift = md("VW_DRIFT");
	if (Mdata.vw_drift != mio::IOUtils::nodata)
		Mdata.vw_drift *= wind_scaling_factor;
	Mdata.dw_drift = md("DW_DRIFT");
	if (Mdata.dw_drift == mio::IOUtils::nodata)
		Mdata.dw_drift = prevailing_wind_dir;

	Mdata.iswr = md(MeteoData::ISWR);
	Mdata.rswr = md(MeteoData::RSWR);

	Mdata.ea = md("EA");
	Mdata.tss = md(MeteoData::TSS);
	if (md.param_exists("TSS_A12H") && (md("TSS_A12H") != mio::IOUtils::nodata))
		Mdata.tss_a12h = md("TSS_A12H");
	else
		Mdata.tss_a12h = Constants::undefined;
	if (md.param_exists("TSS_A24H") && (md("TSS_A24H") != mio::IOUtils::nodata))
		Mdata.tss_a24h = md("TSS_A24H");
	else
		Mdata.tss_a24h = Constants::undefined;
	Mdata.ts0 = md(MeteoData::TSG);

	Mdata.psum_ph = md(MeteoData::PSUM_PH);
	Mdata.psum = md(MeteoData::PSUM);

	Mdata.hs = md(MeteoData::HS);
	if (md.param_exists("HS_A3H") && (md("HS_A3H") != mio::IOUtils::nodata))
		Mdata.hs_a3h = md("HS_A3H");
	else
		Mdata.hs_a3h = Constants::undefined;

	// Add measured new snow density if available
	if (md.param_exists("RHO_HN"))
		Mdata.rho_hn = md("RHO_HN");

	// Add geo_heat if available
	if (md.param_exists("GEO_HEAT"))
		Mdata.geo_heat = md("GEO_HEAT");
	else
		Mdata.geo_heat = mio::IOUtils::nodata;

	// Add advective heat (for permafrost) if available
	if (md.param_exists("ADV_HEAT"))
		Mdata.adv_heat = md("ADV_HEAT");
}

inline double getHS_last3hours(mio::IOManager &io,
                               const mio::Date& current_date) {
	std::vector<mio::MeteoData> MyMeteol3h;

	try {
		io.getMeteoData(current_date - 3.0 / 24.0, MyMeteol3h);  // meteo data with 3 h (left) lag
	} catch (...) {
		cerr << "[E] failed to read meteo data with 3 hours (left) lag\n";
		throw;
	}

	if (MyMeteol3h[0].param_exists("HS_A3H")
	    && (MyMeteol3h[0]("HS_A3H") != mio::IOUtils::nodata))
		return MyMeteol3h[0]("HS_A3H");
	else
		return Constants::undefined;
}

/**
 * @brief Make sure that both short wave fluxes get at least a "realistic" value but measured albedo only if both fluxes are measured
 * @note To be done only for flat field or single slope station
 * @param Mdata
 * @param Xdata
 * @param slope
 */
inline void setShortWave(CurrentMeteo& Mdata, const SnowStation& Xdata,
                         const bool& iswr_is_net) {
	if ((Mdata.iswr > 5.) && (Mdata.rswr > 3.) && !iswr_is_net)
		Mdata.mAlbedo = Mdata.rswr / Mdata.iswr;
	else
		Mdata.mAlbedo = Constants::undefined;

	const double cAlbedo = Xdata.Albedo;

	if (iswr_is_net) {
		const double netSW = Mdata.iswr;
		if (netSW == 0.) {  //this should only happen at night
			Mdata.iswr = 0.;
			Mdata.rswr = 0.;
			return;
		}
		Mdata.iswr = netSW / (1. - cAlbedo);
		Mdata.rswr = netSW / (1. / cAlbedo - 1.);
		return;
	}

	if (Mdata.iswr == mio::IOUtils::nodata)
		Mdata.iswr = Mdata.rswr / Xdata.Albedo;
	if (Mdata.rswr == mio::IOUtils::nodata)
		Mdata.rswr = Mdata.iswr * Xdata.Albedo;
}

/**
 * @brief determine which outputs need to be done for the current time step
 * @param mn_ctrl timestep control structure
 * @param step current time integration step
 * @param sno_step current step in the sno files (current sno profile)
 */

inline void getOutputControl(MainControl& mn_ctrl, const mio::Date& step,
                             const mio::Date& sno_step,
                             const double& calculation_step_length,
                             const double& tsstart, const double& tsdaysbetween,
                             const double& profstart,
                             const double& profdaysbetween,
                             const double& first_backup,
                             const double& backup_days_between) {
//HACK: put all tsstart, tsdaysbetween, etc in MainControl as well as current timestep
	const double Dstep = step.getJulian();
	const double Dsno_step = sno_step.getJulian();
	if (mn_ctrl.resFirstDump) {
		mn_ctrl.HzDump = false;
		mn_ctrl.TsDump = true;
		mn_ctrl.PrDump = true;
		mn_ctrl.resFirstDump = false;
	} else {
		// Hazard data, every half-hour
		mn_ctrl.HzDump = booleanTime(Dstep, 0.5 / 24., 0.0,
		                             calculation_step_length);
		// Time series (*.met)
		double bool_start = H_TO_D(tsstart);
		if (bool_start > 0.)
			bool_start += Dsno_step;
		mn_ctrl.TsDump = booleanTime(Dstep, tsdaysbetween, bool_start,
		                             calculation_step_length);
		// Profile (*.pro)
		bool_start = H_TO_D(profstart);
		if (bool_start > 0.)
			bool_start += Dsno_step;
		mn_ctrl.PrDump = booleanTime(Dstep, profdaysbetween, bool_start,
		                             calculation_step_length);

	}

	// Additional Xdata backup (*.<JulianDate>sno)
	const double bool_start = Dsno_step + first_backup;
	mn_ctrl.XdataDump = booleanTime(Dstep, backup_days_between, bool_start,
	                                calculation_step_length);
}

inline bool readSlopeMeta(mio::IOManager& io, SnowpackIO& snowpackio,
                          SnowpackConfig& cfg, const size_t& i_stn,
                          Slope& slope, mio::Date &current_date,
                          vector<SN_SNOWSOIL_DATA> &vecSSdata,
                          vector<SnowStation> &vecXdata, ZwischenData &sn_Zdata,
                          CurrentMeteo& Mdata, double &wind_scaling_factor,
                          double &time_count_deltaHS) {
	string snowfile;
	stringstream ss;
	ss << "SNOWFILE" << i_stn + 1;
	cfg.getValue(ss.str(), "Input", snowfile, mio::IOUtils::nothrow);

	//Read SSdata for every "slope" referred to as sector where sector 0 corresponds to the main station
	for (size_t sector = slope.mainStation; sector < slope.nSlopes; sector++) {
		try {
			if (sector == slope.mainStation) {
				if (snowfile.empty()) {
					snowfile = vecStationIDs[i_stn];
				} else {
					const size_t pos_dot = snowfile.rfind(".");
					const size_t pos_slash = snowfile.rfind("/");
					if (((pos_dot != string::npos) && (pos_dot > pos_slash))
					    || ((pos_dot != string::npos) && (pos_slash == string::npos)))  //so that the dot is not in a directory name
						snowfile.erase(pos_dot, snowfile.size() - pos_dot);
				}
				snowpackio.readSnowCover(snowfile, vecStationIDs[i_stn],
				                         vecSSdata[slope.mainStation], sn_Zdata, false);
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
				        "Reading snow cover data for station %s",
				        vecStationIDs[i_stn].c_str());
				// NOTE (Is it a HACK?) Reading station meta data provided in meteo data and prebuffering those data
				vector<mio::MeteoData> vectmpmd;
				if (current_date.isUndef())  //either force the start date or take it from the sno file
					current_date = Date::rnd(vecSSdata[slope.mainStation].profileDate, 1);
				else
					vecSSdata[sector].profileDate = current_date;
				io.getMeteoData(current_date, vectmpmd);
				if (vectmpmd.empty())
					throw mio::IOException(
					    "No data found for station " + vecStationIDs[i_stn] + " on "
					        + current_date.toString(mio::Date::ISO),
					    AT);
				Mdata.setMeasTempParameters(vectmpmd[i_stn]);
				vecSSdata[slope.mainStation].meta = mio::StationData::merge(
				    vectmpmd[i_stn].meta, vecSSdata[slope.mainStation].meta);
			} else {
				stringstream sec_snowfile;
				sec_snowfile << snowfile << sector;
				ss.str("");
				ss << vecSSdata[slope.mainStation].meta.getStationID() << sector;
				snowpackio.readSnowCover(sec_snowfile.str(), ss.str(),
				                         vecSSdata[sector], sn_Zdata, false);
				vecSSdata[sector].meta.position = vecSSdata[slope.mainStation].meta
				    .getPosition();
				vecSSdata[sector].meta.stationName = vecSSdata[slope.mainStation].meta
				    .getStationName();
			}
			vecXdata[sector].initialize(vecSSdata[sector], sector);  // Generate the corresponding Xdata
		} catch (const exception& e) {
			if (sector == slope.first) {
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
				        "No virtual slopes! Computation for main station %s only!",
				        vecStationIDs[i_stn].c_str());
				slope.nSlopes = 1;
				if ((mode == "OPERATIONAL")
				    && (vecSSdata[slope.mainStation].meta.getSlopeAngle()
				        > Constants::min_slope_angle)) {
					cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "true");
				}
				break;
			} else {
				cout << e.what();
				throw;
			}
		}

		// Operational mode ONLY: Pass ... wind_factor, snow depth discrepancy time counter
		if ((sector == slope.mainStation) && (mode == "OPERATIONAL")) {
			wind_scaling_factor = vecSSdata[slope.mainStation].WindScalingFactor;
			time_count_deltaHS = vecSSdata[slope.mainStation].TimeCountDeltaHS;
		}
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
	        "Finished initializing station %s", vecStationIDs[i_stn].c_str());

	//CHECK date inconsistencies between sno files
	bool dates_consistent(true);
	for (size_t sector = slope.first; sector < slope.nSlopes; sector++) {
		if (vecSSdata[sector].profileDate
		    != vecSSdata[slope.mainStation].profileDate) {
			prn_msg(
			    __FILE__,
			    __LINE__,
			    "err",
			    mio::Date(),
			    "%s : Date of profile on virtual slope %d inconsistent with flat field",
			    vecStationIDs[i_stn].c_str(), sector);
			dates_consistent = false;

		}
	}
	if (!dates_consistent)
		return false;  //go to next station

	// Do not go ahead if starting time is larger than maxtime!
	if (vecSSdata[slope.mainStation].profileDate > dateEnd) {
		prn_msg(__FILE__, __LINE__, "err", mio::Date(),
		        "%s : Starting time (%.5lf) larger than end time(%.5lf)",
		        vecStationIDs[i_stn].c_str(),
		        vecSSdata[slope.mainStation].profileDate.getJulian(),
		        dateEnd.getJulian());
		return false;  //goto next station
	}

	return true;
}

inline void addSpecialKeys(SnowpackConfig &cfg) {
	const string variant = cfg.get("VARIANT", "SnowpackAdvanced");

	// Add keys to perform running mean in Antarctic variant
	if (variant == "ANTARCTICA") {
		cfg.addKey("VW_AVG::COPY", "Input", "VW");
		cfg.addKey("RH_AVG::COPY", "Input", "RH");

		cfg.addKey("VW_AVG::filter1", "Filters", "AGGREGATE");
		cfg.addKey("VW_AVG::arg1::type", "Filters", "MEAN");
		cfg.addKey("VW_AVG::arg1::soft", "Filters", "true");
		cfg.addKey("VW_AVG::arg1::min_pts", "Filters", "101");
		cfg.addKey("VW_AVG::arg1::min_span", "Filters", "360000");
		cfg.addKey("RH_AVG::filter1", "Filters", "AGGREGATE");
		cfg.addKey("RH_AVG::arg1::type", "Filters", "MEAN");
		cfg.addKey("RH_AVG::arg1::soft", "Filters", "true");
		cfg.addKey("RH_AVG::arg1::min_pts", "Filters", "101");
		cfg.addKey("RH_AVG::arg1::min_span", "Filters", "360000");
	}

	const std::string tst_sw_mode = cfg.get("SW_MODE", "Snowpack");  // Test settings for SW_MODE
	if (tst_sw_mode == "BOTH") {  //HACK: this is only for INP!
		// Make sure there is not only one of ISWR and RSWR available
		bool iswr_inp = true, rswr_inp = true;
		cfg.getValue("ISWR_INP", "Input", iswr_inp, IOUtils::nothrow);
		cfg.getValue("RSWR_INP", "Input", rswr_inp, IOUtils::nothrow);
		if (!(iswr_inp && rswr_inp)) {
			cerr
			    << "[E] SW_MODE = "
			    << tst_sw_mode
			    << ": Please set both ISWR_INP and RSWR_INP to true in [Input]-section of io.ini!\n";
			exit(1);
		}
	}

	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	bool detect_grass = cfg.get("DETECT_GRASS", "SnowpackAdvanced");
	if (mode == "OPERATIONAL") {
		cfg.addKey("RESEARCH", "SnowpackAdvanced", "false");
		cfg.addKey("AVGSUM_TIME_SERIES", "Output", "false");
		if (useCanopyModel) {
			throw mio::IOException("Please don't set CANOPY to 1 in OPERATIONAL mode",
			                       AT);
		}
		if (!detect_grass) {
			cfg.addKey("DETECT_GRASS", "SnowpackAdvanced", "true");
			detect_grass = true;
		}
	}

	if (detect_grass) {
		// we need various average values of tss and hs, all for "past" windows (left)
		// Require at least one value per 3 hours
		cfg.addKey("TSS_A24H::COPY", "Input", "TSS");
		cfg.addKey("TSS_A24H::filter1", "Filters", "AGGREGATE");
		cfg.addKey("TSS_A24H::arg1::type", "Filters", "MEAN");
		cfg.addKey("TSS_A24H::arg1::centering", "Filters", "left");
		cfg.addKey("TSS_A24H::arg1::min_pts", "Filters", "48");  //TODO change # data required to 4
		cfg.addKey("TSS_A24H::arg1::min_span", "Filters", "86340");

		cfg.addKey("TSS_A12H::COPY", "Input", "TSS");
		cfg.addKey("TSS_A12H::filter1", "Filters", "AGGREGATE");
		cfg.addKey("TSS_A12H::arg1::type", "Filters", "MEAN");
		cfg.addKey("TSS_A12H::arg1::centering", "Filters", "left");
		cfg.addKey("TSS_A12H::arg1::min_pts", "Filters", "24");  //TODO change # data required to 2
		cfg.addKey("TSS_A12H::arg1::min_span", "Filters", "43140");

		cfg.addKey("HS_A3H::COPY", "Input", "HS");
		cfg.addKey("HS_A3H::filter1", "Filters", "AGGREGATE");
		cfg.addKey("HS_A3H::arg1::type", "Filters", "MEAN");
		cfg.addKey("HS_A3H::arg1::centering", "Filters", "left");
		cfg.addKey("HS_A3H::arg1::min_pts", "Filters", "6");  //TODO change # data required to 1
		cfg.addKey("HS_A3H::arg1::min_span", "Filters", "10740");
	}

	//warn the user if the precipitation miss proper re-accumulation
	const bool HS_driven = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack");
	if (mode != "OPERATIONAL" && !HS_driven) {
		const bool psum_key_exists = cfg.keyExists("PSUM::resample",
		                                           "Interpolations1D");
		const std::string psum_resampling =
		    (psum_key_exists) ?
		        IOUtils::strToUpper(cfg.get("PSUM::resample", "Interpolations1D")) :
		        "LINEAR";
		if (psum_resampling != "ACCUMULATE") {
			std::cerr
			    << "[W] The precipitation should be re-accumulated over CALCULATION_STEP_LENGTH, not doing it is most probably an error!\n";
		} else {
			const double psum_accumulate = cfg.get("PSUM::accumulate::period",
			                                       "Interpolations1D");
			const double sn_step_length = cfg.get("CALCULATION_STEP_LENGTH",
			                                      "Snowpack");
			if (sn_step_length * 60. != psum_accumulate)
				std::cerr
				    << "[W] The precipitation should be re-accumulated over CALCULATION_STEP_LENGTH (currently, over "
				    << psum_accumulate << "s)\n";
		}
	}
}

inline void writeForcing(Date d1, const Date& d2, const double& Tstep,
                         IOManager &io) {
	std::vector<std::vector<MeteoData> > vecMeteo;
	prn_msg(__FILE__, __LINE__, "msg", mio::Date(),
	        "Reading and writing out forcing data...");

	const std::string experiment = io.getConfig().get("EXPERIMENT", "Output");
	std::map<std::string, size_t> mapIDs;  //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo;  //we need some intermediate storage, for storing data sets for 1 timestep

	for (; d1 <= d2; d1 += Tstep) {  //time loop
		io.getMeteoData(d1, Meteo);  //read 1 timestep at once, forcing resampling to the timestep
		for (size_t ii = 0; ii < Meteo.size(); ii++) {
			const std::string stationID(Meteo[ii].meta.stationID);
			if (mapIDs.count(stationID) == 0) {  //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[stationID] = ii;
				if (ii >= vecMeteo.size())
					vecMeteo.push_back(std::vector<MeteoData>());  //allocate a new station
			}
			Meteo[ii].meta.stationID = Meteo[ii].meta.stationID + "_" + experiment
			    + "_forcing";
			vecMeteo[mapIDs[stationID]].push_back(Meteo[ii]);  //fill the data manually into the vector of vectors
		}
	}
	io.writeMeteoData(vecMeteo);

	prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Forcing data written out");
}

inline void printStartInfo(const SnowpackConfig& cfg, const std::string& name) {
	const bool useSoilLayers = cfg.get("SNP_SOIL", "Snowpack");
	if (useSoilLayers) {
		bool soil_flux = false;
		cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux);
		prn_msg(__FILE__, __LINE__, "msg", mio::Date(),
		        "Start SNOWPACK w/ soil layers in %s mode", mode.c_str());
	} else {
		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Start SNOWPACK in %s mode",
		        mode.c_str());
	}

	const string variant = cfg.get("VARIANT", "SnowpackAdvanced");
	if (variant != "DEFAULT") {
		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Variant is '%s'",
		        variant.c_str());
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "%s compiled on %s at %s",
	        name.c_str(), __DATE__, __TIME__);

	if (mode != "OPERATIONAL") {
		const string experiment = cfg.get("EXPERIMENT", "Output");
		const string outpath = cfg.get("METEOPATH", "Output");
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Experiment : %s",
		        experiment.c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Output dir : %s",
		        outpath.c_str());
	}
}

// SNOWPACK MAIN **************************************************************
inline void real_main(int argc, char *argv[]) {
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );  //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	//parse the command line arguments
	std::string begin_date_str, end_date_str;
	parseCmdLine(argc, argv, begin_date_str, end_date_str);

	mio::Timer meteoRead_timer;
	mio::Timer run_timer;
	run_timer.start();

	MainControl mn_ctrl;  //Time step control parameters

	SnowpackConfig cfg(cfgfile);
	addSpecialKeys(cfg);

	const double i_time_zone = cfg.get("TIME_ZONE", "Input");  //get user provided input time_zone
	if (!begin_date_str.empty()) {
		mio::IOUtils::convertString(dateBegin, begin_date_str, i_time_zone);
	}
	if (end_date_str == "NOW") {  //interpret user provided end date
		dateEnd.setFromSys();
		dateEnd.setTimeZone(i_time_zone);
		dateEnd.rnd(1800, mio::Date::DOWN);
	} else {
		mio::IOUtils::convertString(dateEnd, end_date_str, i_time_zone);
	}

	const string variant = cfg.get("VARIANT", "SnowpackAdvanced");
	const string experiment = cfg.get("EXPERIMENT", "Output");
	const string outpath = cfg.get("METEOPATH", "Output");
	const bool useSoilLayers = cfg.get("SNP_SOIL", "Snowpack");
	const bool useCanopyModel = cfg.get("CANOPY", "Snowpack");
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH",
	                                               "Snowpack");

	int nSolutes = Constants::iundefined;
	cfg.getValue("NUMBER_OF_SOLUTES", "Input", nSolutes, mio::IOUtils::nothrow);
	if (nSolutes > 0)
		SnowStation::number_of_solutes = static_cast<short unsigned int>(nSolutes);

	//Interval between profile backups (*.sno\<JulianDate\>) (d)
	double backup_days_between = 400.;
	cfg.getValue("BACKUP_DAYS_BETWEEN", "Output", backup_days_between,
	             mio::IOUtils::nothrow);

	//First additional profile backup (*.sno\<JulianDate\>) since start of simulation (d)
	double first_backup = 0.;
	cfg.getValue("FIRST_BACKUP", "Output", first_backup, mio::IOUtils::nothrow);

	const double profstart = cfg.get("PROF_START", "Output");
	const double profdaysbetween = cfg.get("PROF_DAYS_BETWEEN", "Output");
	const double tsstart = cfg.get("TS_START", "Output");
	const double tsdaysbetween = cfg.get("TS_DAYS_BETWEEN", "Output");

	const double thresh_rain = cfg.get("THRESH_RAIN", "SnowpackAdvanced");  //Rain only for air temperatures warmer than threshold (degC)
	const bool advective_heat = cfg.get("ADVECTIVE_HEAT", "SnowpackAdvanced");
	const bool soil_flux = (useSoilLayers) ? cfg.get("SOIL_FLUX", "Snowpack") : false;

	//If the user provides the stationIDs - operational use case
	if (!vecStationIDs.empty()) {  //operational use case: stationIDs provided on the command line
		for (size_t i_stn = 0; i_stn < vecStationIDs.size(); i_stn++) {
			stringstream ss;
			ss << "STATION" << i_stn + 1;
			cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]);
		}
	}

	SnowpackIO snowpackio(cfg);
	mio::IOManager io(cfg);
	io.setMinBufferRequirements(IOUtils::nodata, 1.1);  //we require the buffer to contain at least 1.1 day before the current point

	if (vecStationIDs.empty()) {  //research use case: stationIDs provided by the available input files
		vector<StationData> accessible_stations;
		io.getStationData(dateEnd, accessible_stations);  //we are retrieving meta information from MeteoIO
		for (size_t ii = 0; ii < accessible_stations.size(); ii++) {
			vecStationIDs.push_back(accessible_stations[ii].getStationID());  //HACK: accessible_stations should be directly used
		}
	}

	//now, let's start!
	printStartInfo(cfg, string(argv[0]));

	// START LOOP OVER ALL STATIONS
	bool write_forcing = cfg.get("WRITE_PROCESSED_METEO", "Output");  //it will be set to false once it has been done
	for (size_t i_stn = 0; i_stn < vecStationIDs.size(); i_stn++) {
		cout << endl;
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Run on meteo station %s",
		        vecStationIDs[i_stn].c_str());
		run_timer.reset();
		meteoRead_timer.reset();

		Slope slope(cfg);

		// Used to scale wind for blowing and drifting snowpack (from statistical analysis)
		double wind_scaling_factor = cfg.get("WIND_SCALING_FACTOR",
		                                     "SnowpackAdvanced");

		// Control of time window: used for adapting diverging snow depth in operational mode
		double time_count_deltaHS = 0.;

		// Snowpack data (input/output)
		ZwischenData sn_Zdata;  // "Memory"-data, required for every operational station
		vector<SN_SNOWSOIL_DATA> vecSSdata(slope.nSlopes,
		                                   SN_SNOWSOIL_DATA(/*number_of_solutes*/));
		vector<SnowStation> vecXdata;
		for (size_t ii = 0; ii < slope.nSlopes; ii++)  //fill vecXdata with *different* SnowStation objects
			vecXdata.push_back(
			    SnowStation(useCanopyModel, useSoilLayers,
			                (variant == "SEAICE")/*, number_of_solutes*/));

		// Create meteo data object to hold interpolated current time steps
		CurrentMeteo Mdata(cfg);

		mio::Date current_date(dateBegin);
		meteoRead_timer.start();
		if (mode == "OPERATIONAL")
			cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "false");
		const bool read_slope_status = readSlopeMeta(io, snowpackio, cfg, i_stn,
		                                             slope, current_date, vecSSdata,
		                                             vecXdata, sn_Zdata, Mdata,
		                                             wind_scaling_factor,
		                                             time_count_deltaHS);
		meteoRead_timer.stop();
		if (!read_slope_status)
			continue;  //something went wrong, move to the next station

		memset(&mn_ctrl, 0, sizeof(MainControl));
		if (mode == "RESEARCH") {
			mn_ctrl.resFirstDump = true;  //HACK to dump the initial state in research mode
			deleteOldOutputFiles(outpath, experiment, vecStationIDs[i_stn],
			                     slope.nSlopes, snowpackio.getExtensions());
			cfg.write(
			    outpath + "/" + vecStationIDs[i_stn] + "_" + experiment + ".ini");  //output config
			current_date -= calculation_step_length / (24. * 60.);
		} else {
			const string db_name = cfg.get("DBNAME", "Output", "");
			if (db_name == "sdbo" || db_name == "sdbt")
				mn_ctrl.sdbDump = true;
		}

		SunObject sun(vecSSdata[slope.mainStation].meta.position.getLat(),
		              vecSSdata[slope.mainStation].meta.position.getLon(),
		              vecSSdata[slope.mainStation].meta.position.getAltitude());
		sun.setElevationThresh(0.6);
		vector<ProcessDat> qr_Hdata;     //Hazard data for t=0...tn
		vector<ProcessInd> qr_Hdata_ind;  //Hazard data Index for t=0...tn
		const double duration = (dateEnd.getJulian() - current_date.getJulian()
		    + 0.5 / 24) * 24 * 3600;  //HACK: why is it computed this way?
		Hazard hazard(cfg, duration);
		hazard.initializeHazard(sn_Zdata.drift24,
		                        vecXdata.at(0).meta.getSlopeAngle(), qr_Hdata,
		                        qr_Hdata_ind);

		double meteo_step_length = -1.;
		const bool enforce_snow_height = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS",
		                                         "Snowpack");

		//from current_date to dateEnd, if necessary write out meteo forcing
		if (write_forcing == true) {
			writeForcing(current_date, dateEnd, calculation_step_length / 1440, io);
			write_forcing = false;  //no need to call it again for the other stations
		}

		// START TIME INTEGRATION LOOP
		do {
			current_date += calculation_step_length / 1440;
			mn_ctrl.nStep++;
			mn_ctrl.nAvg++;

			// Get meteo data
			vector<mio::MeteoData> vecMyMeteo;
			meteoRead_timer.start();
			io.getMeteoData(current_date, vecMyMeteo);
			if (meteo_step_length < 0.) {
				std::stringstream ss2;
				meteo_step_length = io.getAvgSamplingRate();
				ss2 << meteo_step_length;
				cfg.addKey("METEO_STEP_LENGTH", "Snowpack", ss2.str());
			}
			meteoRead_timer.stop();
			editMeteoData(vecMyMeteo[i_stn], variant, thresh_rain);
			if (!validMeteoData(vecMyMeteo[i_stn], vecStationIDs[i_stn], variant,
			                    enforce_snow_height, advective_heat, soil_flux,
			                    slope.nSlopes)) {
				prn_msg(__FILE__, __LINE__, "msg-", current_date,
				        "No valid data for station %s on [%s]",
				        vecStationIDs[i_stn].c_str(),
				        current_date.toString(mio::Date::ISO).c_str());
				current_date -= calculation_step_length / 1440;
				break;
			}

			//determine which outputs will have to be done
			getOutputControl(mn_ctrl, current_date,
			                 vecSSdata[slope.mainStation].profileDate,
			                 calculation_step_length, tsstart, tsdaysbetween,
			                 profstart, profdaysbetween, first_backup,
			                 backup_days_between);

			//Radiation data
			sun.setDate(current_date.getJulian(), current_date.getTimeZone());


			// START LOOP OVER ASPECTS
			for (unsigned int slope_sequence = 0; slope_sequence < slope.nSlopes;
			    slope_sequence++) {

				SnowpackConfig tmpcfg(cfg);

				//fill Snowpack internal structure with forcing data
				copyMeteoData(vecMyMeteo[i_stn], Mdata, slope.prevailing_wind_dir,
				              wind_scaling_factor);
				Mdata.copySnowTemperatures(vecMyMeteo[i_stn], slope_sequence);
				Mdata.copySolutes(vecMyMeteo[i_stn], SnowStation::number_of_solutes);
				slope.setSlope(slope_sequence, vecXdata, Mdata.dw_drift);

				// SNOWPACK model (Temperature and Settlement computations)
				Snowpack snowpack(tmpcfg);  //the snowpack model to use

				// TODO
				// Read-in top snowpack layers and ensure they are in sync with the meteo data

				// Snow albedo comparison
				cout
				    << snowpack.getParameterizedAlbedo(vecXdata[slope.sector], Mdata)
				        - snowpack.getModelAlbedo(vecXdata[slope.sector], Mdata);  //either parametrized or measured

			}  //end loop on slopes

		} while ((dateEnd.getJulian() - current_date.getJulian())
		    > calculation_step_length / (2. * 1440));

	}


}

int main(int argc, char *argv[]) {
	//try {
	real_main(argc, argv);
	/*} catch (const std::exception &e) {
	 std::cerr << e.what() << endl;
	 throw;
	 }*/

	return EXIT_SUCCESS;
}
