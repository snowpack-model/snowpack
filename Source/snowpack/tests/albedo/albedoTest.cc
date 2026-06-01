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

/************************************************************
 * static section                                           *
 ************************************************************/

//Global variables in this file:
static string cfgfile = "io.ini";
static string mode = "RESEARCH";
static bool restart = false;
static mio::Date dateBegin, dateEnd;
static vector<string> vecStationIDs;

/************************************************************
 * non-static section                                       *
 ************************************************************/

inline void Version()
{
#ifdef _MSC_VER
	cout << "This version of Snowpack uses a BSD-licensed port of getopt for Visual C++. \n"
		<< "It therefore includes software developed by the University of "
		<< "California, Berkeley and its contributors." << endl;
#endif
	cout << "Snowpack version " << SN_VERSION << " compiled on " << __DATE__ << " " << __TIME__ << "\n"
		<< "\tLibsnowpack " << snowpack::getLibVersion() << "\n"
		<< "\tMeteoIO " << mio::getLibVersion() << endl;
}

inline void Usage(const string& programname)
{
	Version();

	cout << "Usage: " << programname << endl
		<< "\t[-b, --begindate=YYYY-MM-DDTHH:MM] (e.g.:2007-08-11T09:00)\n"
		<< "\t[-e, --enddate=YYYY-MM-DDTHH:MM] (e.g.:2008-08-11T09:00 or NOW)\n"
		<< "\t[-c, --config=<ini file>] (e.g. io.ini)\n"
		<< "\t[-m, --mode=<operational or research>] (default: research)\n"
		<< "\t[-r, --restart (skip first time step, only in research mode)\n"
		<< "\t[-s, --stations=<comma delimited stationnames>] (e.g. DAV2,WFJ2)\n"
		<< "\t[-v, --version] Print the version number\n"
		<< "\t[-h, --help] Print help message and version information\n\n";
	cout << "\tPlease note that the operational mode should only be used within SLF\n";

	cout << "Example: " << programname << " -c io.ini -e 1996-06-17T00:00\n\n";
}

inline void printStartInfo(const SnowpackConfig& cfg, const std::string& name)
{
	const bool useSoilLayers = cfg.get("SNP_SOIL", "Snowpack");
	if (useSoilLayers) {
		bool soil_flux = false;
		cfg.getValue("SOIL_FLUX", "Snowpack", soil_flux);
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK w/ soil layers in %s mode", mode.c_str());
	} else {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Start SNOWPACK in %s mode", mode.c_str());
	}

	const std::string variant = cfg.get("VARIANT", "SnowpackAdvanced");
	if (variant != "DEFAULT") {
		prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Variant is '%s'", variant.c_str());
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
	        "%s compiled on %s at %s", name.c_str(), __DATE__, __TIME__);

	if (mode != "OPERATIONAL") {
		const std::string experiment = cfg.get("EXPERIMENT", "Output");
		const std::string outpath = cfg.get("METEOPATH", "Output");
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Experiment : %s", experiment.c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Output dir : %s", outpath.c_str());
	}
}

inline void parseCmdLine(int argc, char **argv, string& begin_date_str, string& end_date_str)
{
	int longindex=0, opt=-1;
	bool setEnd = false;

	struct option long_options[] =
	{
		{"begindate", required_argument, nullptr, 'b'},
		{"enddate", required_argument, nullptr, 'e'},
		{"mode", required_argument, nullptr, 'm'},
		{"restart", no_argument, nullptr, 'r'},
		{"config", required_argument, nullptr, 'c'},
		{"stations", required_argument, nullptr, 's'},
		{"version", no_argument, nullptr, 'v'},
		{"help", no_argument, nullptr, 'h'},
		{nullptr, 0, nullptr, 0}
	};

	if (argc==1) { //no arguments provided
		Usage(string(argv[0]));
		exit(1);
	}

	while ((opt=getopt_long( argc, argv, ":b:e:m:rc:s:v:h", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0:
			break;
		case 'b': {
			begin_date_str=string(optarg); //we don't know yet the time zone, conversion will be done later
			break;
		}
		case 'e': {
			end_date_str=string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		}
		case 'm':
			mode = string(optarg);
			mio::IOUtils::toUpper(mode);
			if (!(mode == "RESEARCH" || mode == "OPERATIONAL")) {
				cerr << endl << "[E] Command line option '-" << char(opt) << "' requires 'research' or 'operational' as operand\n";
				Usage(string(argv[0]));
				exit(1);
			}
			break;
		case 'r':
			restart = true;
			break;
		case 'c':
			cfgfile = string(optarg);
			break;
		case 's':
			mio::IOUtils::readLineToVec(string(optarg), vecStationIDs, ',');
			break;
		case ':': //operand missing
			cerr << endl << "[E] Command line option '-" << char(opt) << "' requires an operand\n";
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
			cerr << endl << "[E] getopt returned character code " <<  opt << "\n";
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

inline bool readSlopeMeta(mio::IOManager& io, SnowpackIO& snowpackio, SnowpackConfig& cfg, const size_t& i_stn,
                   Slope& slope, mio::Date &current_date, vector<SN_SNOWSOIL_DATA> &vecSSdata,
                   vector<SnowStation> &vecXdata, ZwischenData &sn_Zdata, CurrentMeteo& Mdata)
{
	std::string snowfile;
	stringstream ss;
	ss << "SNOWFILE" << i_stn+1;
	cfg.getValue(ss.str(), "Input", snowfile, mio::IOUtils::nothrow);
	const bool slope_from_sno = cfg.get("SLOPE_FROM_SNO", "Input", true);

	//Read SSdata for every "slope" referred to as sector where sector 0 corresponds to the main station
	for (size_t sector=slope.mainStation; sector<slope.nSlopes; sector++) {
		try {
			if (sector == slope.mainStation) {
				if (snowfile.empty()) {
					snowfile = vecStationIDs[i_stn];
				} else {
					const size_t pos_dot = snowfile.rfind(".");
					const size_t pos_slash = snowfile.rfind("/");
					if (((pos_dot != string::npos) && (pos_dot > pos_slash)) ||
						((pos_dot != string::npos) && (pos_slash == string::npos))) //so that the dot is not in a directory name
						snowfile.erase(pos_dot, snowfile.size()-pos_dot);
				}
				snowpackio.readSnowCover(snowfile, vecStationIDs[i_stn], vecSSdata[slope.mainStation], sn_Zdata, (vecXdata[sector].Seaice!=NULL));
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Reading snow cover data for station %s",
				        vecStationIDs[i_stn].c_str());
				// Reading station meta data provided in meteo data and prebuffering those data
				std::vector<mio::MeteoData> vectmpmd;
				if (current_date.isUndef()) //either force the start date or take it from the sno file
					current_date = Date::rnd(vecSSdata[slope.mainStation].profileDate, 1);
				else
					vecSSdata[sector].profileDate = current_date;
				io.getMeteoData(current_date, vectmpmd);
				if (vectmpmd.empty())
					throw mio::IOException("No data found for station " + vecStationIDs[i_stn] + " on "
					                       + current_date.toString(mio::Date::ISO), AT);
				Mdata.setMeasTempParameters(vectmpmd[i_stn]);

				//either get the slope metadata from the sno file or from the meteo data
				if (slope_from_sno) { //position from the meteo forcings, slope and name from the sno file
					vecSSdata[slope.mainStation].meta.position = vectmpmd[i_stn].meta.position;
				} else { //all metadata from the meteo forcings
					vecSSdata[slope.mainStation].meta = vectmpmd[i_stn].meta;
					if (vecSSdata[slope.mainStation].meta.getSlopeAngle()==mio::IOUtils::nodata || vecSSdata[slope.mainStation].meta.getAzimuth()==mio::IOUtils::nodata)
						throw mio::NoDataException("SLOPE_FROM_SNO has been set to false, but slope information is incomplete in the meteorological forcings for station " + vecStationIDs[i_stn], AT);
				}
			} else {
				std::stringstream sec_snowfile;
				sec_snowfile << snowfile << sector;
				ss.str("");
				ss << vecSSdata[slope.mainStation].meta.getStationID() << sector;
				snowpackio.readSnowCover(sec_snowfile.str(), ss.str(), vecSSdata[sector], sn_Zdata, (vecXdata[sector].Seaice!=NULL));
				vecSSdata[sector].meta.position = vecSSdata[slope.mainStation].meta.getPosition();
				vecSSdata[sector].meta.stationName = vecSSdata[slope.mainStation].meta.getStationName();
				if (!current_date.isUndef()) vecSSdata[sector].profileDate = current_date; //this should have been set when processing the main station
			}
			vecXdata[sector].initialize(vecSSdata[sector], sector); // Generate the corresponding Xdata
		} catch (const exception& e) {
			if (sector == slope.first) {
				prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
				        "No virtual slopes! Computation for main station %s only!", vecStationIDs[i_stn].c_str());
				slope.nSlopes = 1;
				if ((mode == "OPERATIONAL")
					&& (vecSSdata[slope.mainStation].meta.getSlopeAngle() > Constants::min_slope_angle)) {
					cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "true");
				}
				break;
			} else {
				cout << e.what();
				throw;
			}
		}
	}
	prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Finished initializing station %s", vecStationIDs[i_stn].c_str());

	//CHECK date inconsistencies between sno files
	bool dates_consistent(true);
	for (size_t sector=slope.first; sector<slope.nSlopes; sector++) {
		if (vecSSdata[sector].profileDate != vecSSdata[slope.mainStation].profileDate) {
			prn_msg(__FILE__, __LINE__, "err", mio::Date(),
				"%s : Date of profile on virtual slope %d inconsistent with flat field", vecStationIDs[i_stn].c_str(), sector);
			dates_consistent = false;

		}
	}
	if (!dates_consistent) return false; //go to next station

	// Do not go ahead if starting time is larger than maxtime!
	if (vecSSdata[slope.mainStation].profileDate > dateEnd) {
		prn_msg(__FILE__, __LINE__, "err", mio::Date(),
			"%s : Starting time (%.5lf) larger than end time(%.5lf)",
			vecStationIDs[i_stn].c_str(), vecSSdata[slope.mainStation].profileDate.getJulian(), dateEnd.getJulian());
		return false; //goto next station
	}

	return true;
}

inline void writeForcing(Date d1, const Date& d2, const double& Tstep, IOManager &io)
{
	std::vector< std::vector<MeteoData> > vecMeteo;
	prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Reading and writing out forcing data...");

	const std::string experiment = io.getConfig().get("EXPERIMENT", "Output");
	std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep

	for(; d1<=d2; d1+=Tstep) { //time loop
		io.getMeteoData(d1, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		for(size_t ii=0; ii<Meteo.size(); ii++) {
			const std::string stationID( Meteo[ii].meta.stationID );
			if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[ stationID ] = ii;
				if (ii>=vecMeteo.size()) vecMeteo.push_back( std::vector<MeteoData>() ); //allocate a new station
			}
			Meteo[ii].meta.stationID = Meteo[ii].meta.stationID + "_" + experiment + "_forcing";
			vecMeteo[ mapIDs[stationID] ].push_back(Meteo[ii]); //fill the data manually into the vector of vectors
		}
	}
	io.writeMeteoData(vecMeteo);

	prn_msg(__FILE__, __LINE__, "msg",  mio::Date(), "Forcing data written out");
}

// SNOWPACK MAIN **************************************************************
inline void real_main (int argc, char *argv[])
{
	setbuf(stdout, NULL); //always flush stdout
	setbuf(stderr, NULL); //always flush stderr
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	//parse the command line arguments
	std::string begin_date_str, end_date_str;
	parseCmdLine(argc, argv, begin_date_str, end_date_str);
	SnowpackConfig cfg(cfgfile);
	cfg.checkStandaloneKeys(mode); //perform some consistency checks on the configuration
	cfg.addStandaloneKeys(mode); //add special keys to run Snowpack standalone

	mio::Timer meteoRead_timer;
	mio::Timer run_timer;
	run_timer.start();
	time_t nowSRT = time(NULL);
	MainControl mn_ctrl(cfg); //Time step control parameters

	const double i_time_zone = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	if (!begin_date_str.empty()) {
		mio::IOUtils::convertString(dateBegin, begin_date_str, i_time_zone);
	}
	if (end_date_str == "NOW") { //interpret user provided end date
		dateEnd.setFromSys();
		dateEnd.setTimeZone(i_time_zone);
		dateEnd.rnd(1800, mio::Date::DOWN);
	} else {
		mio::IOUtils::convertString(dateEnd, end_date_str, i_time_zone);
	}

	const std::string variant = cfg.get("VARIANT", "SnowpackAdvanced");
	const std::string experiment = cfg.get("EXPERIMENT", "Output");
	const std::string outpath = cfg.get("METEOPATH", "Output");

	int nSolutes = Constants::iundefined;
	cfg.getValue("NUMBER_OF_SOLUTES", "Input", nSolutes, mio::IOUtils::nothrow);
	if (nSolutes > 0) SnowStation::number_of_solutes = static_cast<short unsigned int>(nSolutes);

	const bool prn_check = cfg.get("INFLATE_INFO", "Snowpack", false);
	const double inflate_tracking_period = cfg.get("INFLATE_TRACKING_PERIOD", "Snowpack", 1.); // in days
	const bool grooming = cfg.get("SNOW_GROOMING", "TechSnow");
	const bool classify_profile = cfg.get("CLASSIFY_PROFILE", "Output");
	const bool profwrite = cfg.get("PROF_WRITE", "Output");
	const bool tswrite = cfg.get("TS_WRITE", "Output");
	const bool snow_write = cfg.get("SNOW_WRITE", "Output");

	const bool precip_rates = cfg.get("PRECIP_RATES", "Output");
	const bool avgsum_time_series = cfg.get("AVGSUM_TIME_SERIES", "Output");
	const bool cumsum_mass = cfg.get("CUMSUM_MASS", "Output");

	//If the user provides the stationIDs - operational use case
	if (!vecStationIDs.empty()) { //operational use case: stationIDs provided on the command line
		for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
			stringstream ss;
			ss << "METEOFILE" << i_stn+1;
			cfg.addKey(ss.str(), "Input", vecStationIDs[i_stn]);
		}
	}

	SnowpackIO snowpackio(cfg);
	mio::IOManager io(cfg);
	io.setMinBufferRequirements(IOUtils::nodata, 1.1); //we require the buffer to contain at least 1.1 day before the current point

	if (vecStationIDs.empty()) { //research use case: stationIDs provided by the available input files
		vector<StationData> accessible_stations;
		io.getStationData(dateEnd, accessible_stations); //we are retrieving meta information from MeteoIO
		for (size_t ii=0; ii<accessible_stations.size(); ii++) {
			vecStationIDs.push_back( accessible_stations[ii].getStationID() ); //HACK: accessible_stations should be directly used
		}
	}

	//now, let's start!
	printStartInfo(cfg, string(argv[0]));

	// START LOOP OVER ALL STATIONS
	bool write_forcing = cfg.get("WRITE_PROCESSED_METEO", "Output"); //it will be set to false once it has been done
	for (size_t i_stn=0; i_stn<vecStationIDs.size(); i_stn++) {
		Meteo meteo(cfg);
		cout << endl;
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Run on meteo station %s", vecStationIDs[i_stn].c_str());
		run_timer.reset();
		meteoRead_timer.reset();

		Slope slope(cfg);
		Cumsum cumsum(slope.nSlopes);

		double lw_in = Constants::undefined;    // Storage for LWin from flat field energy balance

		// Snowpack data (input/output)
		ZwischenData sn_Zdata;   // "Memory"-data, required for every operational station
		vector<SN_SNOWSOIL_DATA> vecSSdata(slope.nSlopes, SN_SNOWSOIL_DATA(/*number_of_solutes*/));
		vector<SnowStation> vecXdata;
		for (size_t ii=0; ii<slope.nSlopes; ii++) { //fill vecXdata with *different* SnowStation objects
			vecXdata.push_back( SnowStation(cfg, false /*Is A3d?*/, (mode == "OPERATIONAL")) );
			if (vecXdata.back().Seaice != NULL) vecXdata[ii].Seaice->ConfigSeaIce(cfg);
		}

		// Create meteo data object to hold interpolated current time steps
		CurrentMeteo Mdata(cfg);
		// To collect surface exchange data for output
		SurfaceFluxes surfFluxes/*(number_of_solutes)*/;
		// Boundary condition (fluxes)
		BoundCond sn_Bdata;

		mio::Date current_date( dateBegin );
		meteoRead_timer.start();
		if (mode == "OPERATIONAL")
			cfg.addKey("PERP_TO_SLOPE", "SnowpackAdvanced", "false");
		const bool read_slope_status = readSlopeMeta(io, snowpackio, cfg, i_stn, slope, current_date, vecSSdata, vecXdata, sn_Zdata, Mdata);
		meteoRead_timer.stop();
		if (!read_slope_status) continue; //something went wrong, move to the next station

		mn_ctrl.reset(); //values read from cfg are kept, others set back to 0 / false
		if (mode == "RESEARCH") {
			if (!restart) {
				mn_ctrl.resFirstDump = true;  //HACK to dump the initial state in research mode
				current_date -= mn_ctrl.calculation_step_length/(24.*60.); //Do a first time step to fill all output fields
			} else {
				mn_ctrl.resFirstDump = false; //No initial state dump when doing a restart
			}
			deleteOldOutputFiles(outpath, experiment, vecStationIDs[i_stn], slope.nSlopes, snowpackio.getExtensions());
			cfg.write(outpath + "/" + vecStationIDs[i_stn] + "_" + experiment + ".ini"); //output config
		} else {
			const std::string db_name = cfg.get("DBNAME", "Output", "");
			if (db_name == "sdbo" || db_name == "sdbt")
				mn_ctrl.sdbDump = true;
		}

		SunObject sun(vecSSdata[slope.mainStation].meta.position.getLat(), vecSSdata[slope.mainStation].meta.position.getLon(), vecSSdata[slope.mainStation].meta.position.getAltitude());
		sun.setElevationThresh(0.6);
		vector<ProcessDat> qr_Hdata;     //Hazard data for t=0...tn
		vector<ProcessInd> qr_Hdata_ind; //Hazard data Index for t=0...tn
		const double duration = (dateEnd.getJulian() - current_date.getJulian() + 0.5/24)*24*3600; //HACK: why is it computed this way?
		Hazard hazard(cfg, duration);
		hazard.initializeHazard(sn_Zdata.drift24, vecXdata.at(0).meta.getSlopeAngle(), qr_Hdata, qr_Hdata_ind);

		prn_msg(__FILE__, __LINE__, "msg", mio::Date(), "Start simulation for %s on %s",
			vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO_TZ).c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "End date specified by user: %s",
		        dateEnd.toString(mio::Date::ISO_TZ).c_str());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Integration step length: %f min",
		        mn_ctrl.calculation_step_length);

		bool computed_one_timestep = false;
		double meteo_step_length = -1.;
		const bool enforce_snow_height = cfg.get("ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack");

		//from current_date to dateEnd, if necessary write out meteo forcing
		if (write_forcing==true) {
			writeForcing(current_date, dateEnd, mn_ctrl.calculation_step_length/1440, io);
			write_forcing = false; //no need to call it again for the other stations
		}

		// START TIME INTEGRATION LOOP
		do {
			current_date += mn_ctrl.calculation_step_length/1440; //converting calculation_step_length in days
			mn_ctrl.nStep++;
			mn_ctrl.nAvg++;

			// Get meteo data
			vector<mio::MeteoData> vecMyMeteo;
			meteoRead_timer.start();
			io.getMeteoData(current_date, vecMyMeteo);
			if (vecMyMeteo.empty()) {
				prn_msg(__FILE__, __LINE__, "msg-", current_date, "No forcing data provided for [%s]",
				        current_date.toString(mio::Date::ISO).c_str());
				current_date -= mn_ctrl.calculation_step_length/1440;
				break;
			}
			if(meteo_step_length<0.) {
				std::stringstream ss2;
				meteo_step_length = io.getAvgSamplingRate();
				ss2 << meteo_step_length;
				cfg.addKey("METEO_STEP_LENGTH", "Snowpack", ss2.str());
			}
			meteoRead_timer.stop();
			Mdata.editMeteoData(vecMyMeteo[i_stn]);
			if (!Mdata.validMeteoData(vecMyMeteo[i_stn], vecStationIDs[i_stn], slope.nSlopes)) {
				prn_msg(__FILE__, __LINE__, "msg-", current_date, "No valid data for station %s on [%s]",
				        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str());
				current_date -= mn_ctrl.calculation_step_length/1440;
				break;
			}

			//determine which outputs will have to be done
			mn_ctrl.getOutputControl(current_date, vecSSdata[slope.mainStation].profileDate);

			//Radiation data
			sun.setDate(current_date.getJulian(), current_date.getTimeZone());
			const double hs_a3hl6 = Mdata.getHS_last3hours(io, current_date);

			// START LOOP OVER ASPECTS
			for (unsigned int slope_sequence=0; slope_sequence<slope.nSlopes; slope_sequence++) {
				double tot_mass_in = 0.; // To check mass balance over one CALCULATION_STEP_LENGTH if MASS_BALANCE is set
				SnowpackConfig tmpcfg(cfg);

				//fill Snowpack internal structure with forcing data
				Mdata.setMeteoData(vecMyMeteo[i_stn], slope.prevailing_wind_dir);
				Mdata.copySnowTemperatures(vecMyMeteo[i_stn], slope_sequence);
				Mdata.copySolutes(vecMyMeteo[i_stn], SnowStation::number_of_solutes);
				slope.setSlope(slope_sequence, vecXdata, Mdata.dw_drift);
				Mdata.dataForCurrentTimeStep(surfFluxes, vecXdata, slope, tmpcfg,
                                       sun, cumsum.precip, lw_in, hs_a3hl6,
                                       tot_mass_in);

				// Find the Wind Profile Parameters, w/ or w/o canopy; take care of canopy
				meteo.compMeteo(Mdata, vecXdata[slope.sector], true);

				// Store ea from main slope to be used on the virtual slopes (in case it was modified in copyMeteoData or dataForCurrentTimeStep)
				if (slope.sector == slope.mainStation) vecMyMeteo[i_stn]("EA") = Mdata.ea;

				// Notify user every fifteen days of date being processed
				const double notify_start = floor(vecSSdata[slope.mainStation].profileDate.getJulian()) + 15.5;
				if ((mode == "RESEARCH") && (slope.sector == slope.mainStation)
				        && booleanTime(current_date.getJulian(), 15., notify_start, mn_ctrl.calculation_step_length)) {
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					            "Station %s (%d slope(s)): advanced to %s station time",
					                vecSSdata[slope.mainStation].meta.stationID.c_str(), slope.nSlopes,
					                    current_date.toString(mio::Date::DIN).c_str());
				}

				// SNOWPACK model (Temperature and Settlement computations)
				Snowpack snowpack(tmpcfg); //the snowpack model to use
				Stability stability(tmpcfg, classify_profile);
				snowpack.runSnowpackModel(Mdata, vecXdata[slope.sector], cumsum.precip, sn_Bdata, surfFluxes);

				if (grooming)
					snowpack.snowPreparation(current_date, vecXdata[slope.sector] );

				stability.checkStability(Mdata, vecXdata[slope.sector]);

				/***** OUTPUT SECTION *****/
				surfFluxes.collectSurfaceFluxes(sn_Bdata, vecXdata[slope.sector], Mdata);
				if (slope.sector == slope.mainStation) { // main station only (usually flat field)
					// Calculate consistent lw_in for virtual slopes
					if ( vecXdata[slope.mainStation].getNumberOfElements() > 0 ) {
						double k_eff, gradT;
						k_eff =
						    vecXdata[slope.mainStation].Edata[vecXdata[slope.mainStation].getNumberOfElements()-1].k[TEMPERATURE];
						gradT =
						    vecXdata[slope.mainStation].Edata[vecXdata[slope.mainStation].getNumberOfElements()-1].gradT;
						lw_in = k_eff*gradT + sn_Bdata.lw_out - sn_Bdata.qs - sn_Bdata.ql - sn_Bdata.qr;
					} else {
						lw_in = Constants::undefined;
					}
					// Deal with new snow densities
					if (vecXdata[slope.mainStation].hn > 0.) {
						surfFluxes.cRho_hn = vecXdata[slope.mainStation].rho_hn;
						surfFluxes.mRho_hn = Mdata.rho_hn;
					}
					if (slope.snow_erosion) {
						// Update drifting snow index (VI24),
						//   from erosion at the main station only if no virtual slopes are available
						if (slope.mainStationDriftIndex)
							cumulate(cumsum.drift, surfFluxes.drift);
						// Update erosion mass from main station
						// NOTE cumsum.erosion[] will be positive in case of real erosion at any time during the output time step
						if (vecXdata[slope.mainStation].ErosionMass > Constants::eps) {
							// Real erosion
							if (cumsum.erosion[slope.mainStation] > Constants::eps)
								cumsum.erosion[slope.mainStation] += vecXdata[slope.mainStation].ErosionMass;
							else
								cumsum.erosion[slope.mainStation] = vecXdata[slope.mainStation].ErosionMass;
						} else {
							// Potential erosion at main station only
							if (cumsum.erosion[slope.mainStation] < -Constants::eps)
								cumsum.erosion[slope.mainStation] -= surfFluxes.mass[SurfaceFluxes::MS_WIND];
							else if (!(cumsum.erosion[slope.mainStation] > Constants::eps))
								cumsum.erosion[slope.mainStation] = -surfFluxes.mass[SurfaceFluxes::MS_WIND];
						}
					}

					const size_t i_hz = mn_ctrl.HzStep;
					if (mode == "OPERATIONAL" && !cumsum_mass) { // Cumulate flat field runoff in operational mode
						qr_Hdata.at(i_hz).runoff += surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];
						cumsum.runoff += surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF];
					}

					//check if inflate-deflate is required and perform it if necessary
					if (enforce_snow_height && vecXdata[slope.mainStation].allow_inflate) {
						if (vecXdata[slope.mainStation].needDeflateInflate(mn_ctrl.sn_dt, inflate_tracking_period)) {
							vecXdata[slope.mainStation].deflateInflate(Mdata, qr_Hdata.at(i_hz).dhs_corr, qr_Hdata.at(i_hz).mass_corr, prn_check);
							vecXdata[slope.mainStation].TimeCountDeltaHS = 0.; //Time counter tracking erroneous settlement in operational mode
						}
					}

					if (mn_ctrl.HzDump) { // Save hazard data ...
						qr_Hdata.at(i_hz).stat_abbrev = vecStationIDs[i_stn];
						if (mode == "OPERATIONAL") {
							qr_Hdata.at(i_hz).loc_for_snow = (unsigned char)vecStationIDs[i_stn][vecStationIDs[i_stn].length()-1];
							//TODO: WHAT SHOULD WE SET HERE? wstat_abk (not existing yet in DB) and wstao_nr, of course;-)
							qr_Hdata_ind.at(i_hz).loc_for_wind = -1;
						} else {
							qr_Hdata.at(i_hz).loc_for_snow = 2;
							qr_Hdata.at(i_hz).loc_for_wind = 1;
						}
						hazard.getHazardDataMainStation(qr_Hdata.at(i_hz), qr_Hdata_ind.at(i_hz),
						                                sn_Zdata, cumsum.drift, slope.mainStationDriftIndex,
						                                vecXdata[slope.mainStation], Mdata, surfFluxes);
						if (slope.nSlopes==1) { //only one slope, so set lwi_N and lwi_S to the same value
							const double lwi = vecXdata[slope.mainStation].getLiquidWaterIndex();
							if ((lwi < -Constants::eps) || (lwi >= 10.))
								qr_Hdata_ind.at(i_hz).lwi_N = qr_Hdata_ind.at(i_hz).lwi_S = false;
							qr_Hdata.at(i_hz).lwi_N = lwi;
							qr_Hdata.at(i_hz).lwi_S = lwi;
						}
						mn_ctrl.HzStep++;
						if (slope.mainStationDriftIndex)
							cumsum.drift = 0.;
						surfFluxes.hoar = 0.;
						// Inflate/deflate sums
						cumsum.dhs_corr += qr_Hdata.at(i_hz).dhs_corr;
						cumsum.mass_corr += qr_Hdata.at(i_hz).mass_corr;
					}
					// New snow water equivalent (kg m-2), rain was dealt with in Watertransport.cc
					surfFluxes.mass[SurfaceFluxes::MS_HNW] += vecXdata[slope.mainStation].hn
					                                              * vecXdata[slope.mainStation].rho_hn;
					if (!avgsum_time_series) { // Sum up precipitations
						cumsum.rain += surfFluxes.mass[SurfaceFluxes::MS_RAIN];
						cumsum.snow += surfFluxes.mass[SurfaceFluxes::MS_HNW];
					}
				} else {
					const size_t i_hz = (mn_ctrl.HzStep > 0) ? mn_ctrl.HzStep-1 : 0;
					if (slope.luvDriftIndex) {
						// Update drifting snow index (VI24),
						// considering only snow eroded from the windward slope
						cumulate(cumsum.drift, surfFluxes.drift);
					}
					if (mn_ctrl.HzDump) {
						// NOTE qr_Hdata was first saved at the end of the mainStation simulation, at which time the drift index could not be dumped!
						hazard.getHazardDataSlope(qr_Hdata.at(i_hz), qr_Hdata_ind.at(i_hz),
						                          sn_Zdata.drift24, cumsum.drift, vecXdata[slope.sector],
						                          slope.luvDriftIndex, slope.north, slope.south);
						if(slope.luvDriftIndex) cumsum.drift = 0.;
					}

					// Update erosion mass from windward virtual slope
					cumsum.erosion[slope.sector] += vecXdata[slope.sector].ErosionMass;
				}

				// TIME SERIES (*.met)
				if (tswrite && mn_ctrl.TsDump) {
					// Average fluxes
					if (avgsum_time_series) {
						averageFluxTimeSeries(mn_ctrl.nAvg, vecXdata[slope.sector].useCanopyModel, surfFluxes, vecXdata[slope.sector]);
					} else {
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] = cumsum.rain;
						surfFluxes.mass[SurfaceFluxes::MS_HNW] = cumsum.snow;
						// Add eroded snow from luv to precipitations on lee slope
						if (slope.sector == slope.lee && cumsum.erosion[slope.luv] > Constants::eps)
							surfFluxes.mass[SurfaceFluxes::MS_HNW] += cumsum.erosion[slope.luv] / vecXdata[slope.luv].cos_sl;
					}

					if (precip_rates) { // Precip rates in kg m-2 h-1
						surfFluxes.mass[SurfaceFluxes::MS_RAIN] /= static_cast<double>(mn_ctrl.nAvg)*M_TO_H(mn_ctrl.calculation_step_length);
						surfFluxes.mass[SurfaceFluxes::MS_HNW] /= static_cast<double>(mn_ctrl.nAvg)*M_TO_H(mn_ctrl.calculation_step_length);
						if ((mode == "OPERATIONAL") && (!cumsum_mass)) {
							surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] = cumsum.runoff;
							surfFluxes.mass[SurfaceFluxes::MS_SNOWPACK_RUNOFF] /= static_cast<double>(mn_ctrl.nAvg)*M_TO_H(mn_ctrl.calculation_step_length);
							cumsum.runoff = 0.;
						}
					}

					// Erosion mass rate in kg m-2 h-1
					surfFluxes.mass[SurfaceFluxes::MS_WIND] = cumsum.erosion[slope.sector];
					surfFluxes.mass[SurfaceFluxes::MS_WIND] /= static_cast<double>(mn_ctrl.nAvg)*M_TO_H(mn_ctrl.calculation_step_length);

					// Dump
					const size_t i_hz = (mn_ctrl.HzStep > 0) ? mn_ctrl.HzStep - 1 : 0;
					size_t i_hz0 = (mn_ctrl.HzStep > 1) ? mn_ctrl.HzStep - 2 : 0;
					if (slope.mainStationDriftIndex)
						i_hz0 = i_hz;
					const double wind_trans24 = (slope.sector == slope.mainStation) ? qr_Hdata.at(i_hz0).wind_trans24 : qr_Hdata.at(i_hz).wind_trans24;
					ProcessDat tmpHdata = qr_Hdata.at(i_hz);			// Temporary ProcessDat object for output
					tmpHdata.dhs_corr = cumsum.dhs_corr; cumsum.dhs_corr = 0.;	// overwrite the inflate/deflate variables
					tmpHdata.mass_corr = cumsum.mass_corr; cumsum.mass_corr = 0.;
					snowpackio.writeTimeSeries(vecXdata[slope.sector], surfFluxes, Mdata,
					                           tmpHdata, wind_trans24);

					if (avgsum_time_series) {
						surfFluxes.reset(cumsum_mass);
						if (vecXdata[slope.sector].useCanopyModel) vecXdata[slope.sector].Cdata.reset(cumsum_mass);
					}
					surfFluxes.cRho_hn = Constants::undefined;
					surfFluxes.mRho_hn = Constants::undefined;
					// reset cumulative variables
					if (slope_sequence == slope.nSlopes-1) {
						cumsum.erosion.assign(cumsum.erosion.size(), 0.);
						cumsum.rain = cumsum.snow = 0.;
						mn_ctrl.nAvg = 0;
					}
				}

				// SNOW PROFILES ...
				// ... for visualization (*.pro), etc. (*.prf)
				if (profwrite && mn_ctrl.PrDump)
					snowpackio.writeProfile(current_date, vecXdata[slope.sector]);

				// ... backup Xdata (*.sno<JulianDate>)
				if (mn_ctrl.XdataDump) {
					std::stringstream ss;
					ss << vecStationIDs[i_stn];
					if (slope.sector != slope.mainStation) ss << slope.sector;
					snowpackio.writeSnowCover(current_date, vecXdata[slope.sector], sn_Zdata, true);
					prn_msg(__FILE__, __LINE__, "msg", current_date,
					        "Backup Xdata dumped for station %s [%.2f days, step %d]", ss.str().c_str(),
					        (current_date.getJulian()
					            - (vecSSdata[slope.mainStation].profileDate.getJulian() + 0.5/24)),
					        mn_ctrl.nStep);
				}

				// check mass balance if AVGSUM_TIME_SERIES is not set (screen output only)
				if (!avgsum_time_series) {
					const bool mass_balance = cfg.get("MASS_BALANCE", "SnowpackAdvanced");
					if (mass_balance) {
						if (massBalanceCheck(vecXdata[slope.sector], surfFluxes, tot_mass_in) == false)
							prn_msg(__FILE__, __LINE__, "msg+", current_date, "Mass error at end of time step!");
					}
				}

				if (tswrite && mn_ctrl.TsDump) vecXdata[slope.sector].resetSlopeParFlux();

				// Snow albedo comparison
				cout
				    << snowpack.getParameterizedAlbedo(vecXdata[slope.sector], Mdata)
				        - snowpack.getModelAlbedo(vecXdata[slope.sector], Mdata);  //either parametrized or measured

			} //end loop on slopes
			computed_one_timestep = true;
		} while ((dateEnd.getJulian() - current_date.getJulian()) > mn_ctrl.calculation_step_length/(2.*1440));
		//end loop on timesteps

		// If the simulation run for at least one time step,
		//   dump the PROFILEs (Xdata) for every station referred to as sector where sector 0 corresponds to the main station
		if (computed_one_timestep && snow_write) {
			for (size_t sector=slope.mainStation; sector<slope.nSlopes; sector++) {
				snowpackio.writeSnowCover(current_date, vecXdata[sector], sn_Zdata);
				if (sector == slope.mainStation) {
					prn_msg(__FILE__, __LINE__, "msg", mio::Date(),
					        "Writing data to sno file(s) for %s (station %s) on %s",
					        vecSSdata[slope.mainStation].meta.getStationName().c_str(),
					        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str());
				}
			}
			// Dump time series to snowpack.ams_pmod@SDBx (hazard data)
			if (mn_ctrl.sdbDump) {
				mio::Timer sdbDump_timer;
				sdbDump_timer.reset();
				sdbDump_timer.start();
				if (snowpackio.writeHazardData(vecStationIDs[i_stn], qr_Hdata, qr_Hdata_ind, mn_ctrl.HzStep)) {
					sdbDump_timer.stop();
					prn_msg(__FILE__, __LINE__, "msg-", mio::Date(),
					        "Finished writing Hdata to SDB for station %s on %s (%lf s)",
					        vecStationIDs[i_stn].c_str(), current_date.toString(mio::Date::ISO).c_str(), sdbDump_timer.getElapsed());
				}
			}
		}
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Total time to read meteo data : %lf s",
		        meteoRead_timer.getElapsed());
		prn_msg(__FILE__, __LINE__, "msg-", mio::Date(), "Runtime for station %s: %lf s",
		        vecStationIDs[i_stn].c_str(), run_timer.getElapsed());
	}

	time_t nowEND=time(NULL);
	cout << endl;
	cout << "[i] []                 STARTED  running SLF " << mode << " Snowpack Model on " << ctime(&nowSRT);
	if (mode == "OPERATIONAL"){
		cout << "                       ===========================================================================" << endl;
	} else {
		cout << "                       ========================================================================" << endl;
	}
	cout << "                       FINISHED running SLF " << mode << " Snowpack Model on " << ctime(&nowEND) << endl;
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
