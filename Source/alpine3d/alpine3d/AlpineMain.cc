/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef DEBUG_ARITHM
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif
	#include <fenv.h>
#endif

#include <exception>
#include <iostream>
#include <csignal>
#include <unistd.h>
#include <getopt.h>
#include <ctime>

#include <meteoio/MeteoIO.h>
#include <snowpack/libsnowpack.h>

#include <alpine3d/AlpineMain.h>
#include <alpine3d/AlpineControl.h>
#include <alpine3d/MPIControl.h>

#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/DataAssimilation.h>

using namespace std;
using namespace mio;

//Global variables for AlpineMain
static int steps = 0; //how many hours to simulate
static Date startdate;

static bool enable_eb=false, enable_drift=false, enable_runoff=false, enable_da=false, nocompute=false, restart=false;
static int npsnowpack=1, npebalance=1;

inline void Usage(char *prog)
{
	cout << "Usage: " << prog << "\n";
	cout << "\t-a, --startdate=YYYY-MM-DDTHH:MM (e.g.:2008-08-11T09:00)\n";
	cout << "\t-z, --enddate=YYYY-MM-DDTHH:MM (e.g.:2008-08-11T09:00) OR ";
	cout << "-n, --steps=<number of timesteps>\n";
	cout << "\t[--enable-eb] (default:off)\n";
	cout << "\t[--enable-runoff] (default:off)\n";
	cout << "\t[--enable-drift] (default:off)\n";
	cout << "\t[--enable-da] (default:off)\n";
	cout << "\t[--no-compute] Don't compute the simulation, only parse and validate the inputs (default:off)\n";
	cout << "\t[--restart] Restart simulation from previous .sno files (default:off)\n";
	cout << "\t[-p, --np-snowpack=<number of processors for SnowPack>]\n";
	cout << "\t[-b, --np-ebalance=<number of processors for Ebalance>]\n";
	cout << "\t[-i, --iofile=<MeteoIO ini file> (default:./io.ini)]\n";
	cout << "\t[-h, --help] Print help message and version information\n";

	exit(1);
}

inline void parseCmdLine(int argc, char **argv, Config &cfg)
{
	std::string iofile( "io.ini" );
	int eeb=0, edr=0, ero=0, eda=0, nco=0, rst=0;
	int longindex=0, opt = -1;
	bool setStart=false, setEnd=false, setSteps=false;
	std::string start_date_input, end_date_input;

	const struct option long_options[] =
	{
		/* These options set a flag. */
		{"enable-eb",		no_argument,&eeb, 1},
		{"enable-drift",	no_argument,&edr, 1},
		{"enable-runoff",	no_argument,&ero, 1},
		{"enable-da",		no_argument,&eda, 1},
		{"no-compute",		no_argument,&nco, 1},
		{"restart",		no_argument,&rst, 1},

		/* These options don't set a flag.
			We distinguish them by their indices. */
		{"np-snowpack",		required_argument, 0, 'p'},
		{"np-ebalance",		required_argument, 0, 'b'},
		{"iofile",		required_argument, 0, 'i'},
		{"startdate",		required_argument, 0, 'a'},
		{"enddate",		required_argument, 0, 'z'},
		{"steps",		required_argument, 0, 'n'},
		{"help",		no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	while ((opt=getopt_long( argc, argv, ":d:p:r:s:e:i:f:a:z:n:h", long_options, &longindex)) != -1) {
		switch (opt){
		case 0:
			break;
		case 'b':
			IOUtils::convertString(npebalance, string(optarg));
			if (npebalance < 1)
				throw IOException("Number of processors for Ebalance computation is not valid!", AT);
			break;
		case 'p':
			IOUtils::convertString(npsnowpack, string(optarg));
			if (npsnowpack < 1)
				throw IOException("Number of processors for SnowPack computation is not valid!", AT);
			break;
		case 'i':
			iofile = string(optarg);
			break;
		case 'a':
			start_date_input = string(optarg);
			setStart = true;
			break;
		case 'z':
			end_date_input = string(optarg);
			setEnd = true;
			break;
		case 'n':
			IOUtils::convertString(steps, string(optarg));
			setSteps = true;
			break;

		case ':': //operand missing
			cout << "Command line parameter '-" << char(optopt) << "' requires an operand";
			break;
		case 'h':
			Usage(argv[0]);
			break;
		case '?':
			cerr << "\nERROR: Unknown argument detected\n";
			Usage(argv[0]);
			break;
		default:
			cerr << "\nERROR: getopt returned character code " <<  opt << "\n";
			Usage(argv[0]);
		}
	}

	//check that the minimum flags have been provided by the user
	if (!((setStart && setSteps) || (setStart && setEnd)) && MPIControl::instance().master()) {
		cout << "\nERROR: You must at least specify the 'startdate' and the 'steps' parameters"
		     << " or the 'startdate' and the 'enddate' parameters\n\n";
		Usage(argv[0]);
	}

	enable_eb	= (eeb!=0);
	enable_drift	= (edr!=0);
	enable_runoff	= (ero!=0);
	enable_da	= (eda!=0);
	nocompute	= (nco!=0);
	restart		= (rst!=0);

	//now, we read the config file and set the start and end dates
	//only the master requires the io.ini, the rest receives a the cfg object:
	const bool isMaster = MPIControl::instance().master(); // Check if this process is the master (always true for non-parallel mode)
	if (isMaster) {
		cfg.addFile(iofile);
	}
	MPIControl::instance().broadcast(cfg);
	
	const double TZ = cfg.get("TIME_ZONE", "INPUT");
	startdate.setTimeZone(TZ);
	Date enddate;
	enddate.setTimeZone(TZ);
	IOUtils::convertString(startdate, start_date_input, TZ);
	if (end_date_input == "NOW") { //interpret user provided end date
		enddate.setFromSys();
		enddate.setTimeZone(TZ);
		enddate.rnd(1800, Date::DOWN);
	} else {
		IOUtils::convertString(enddate, end_date_input, TZ);
	}

	//check that the start and end dates make sense
	if (setStart && setEnd) {
		steps = (int)floor((enddate.getJulian(true) - startdate.getJulian(true))*24.+0.01) + 1; //to be on the safe side when rounding
		if (steps <= 0) {
		  throw InvalidArgumentException("startdate="+startdate.toString(Date::ISO)+" > enddate="+enddate.toString(Date::ISO),AT);
		}
	}
	
	//SNOWPARAM: Check for time step consistency
	const double sp_dt = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	if ( dt_main != M_TO_S(sp_dt) ) {
		if (isMaster) cout << "[I] Alpine3D's dt_main=" << S_TO_M(dt_main) << " [min], Snowpack's CALCULATION_STEP_LENGTH=" << sp_dt << " [min] (in io.ini)\t";

		if ( (int)(dt_main)%(int)(M_TO_S(sp_dt)) == 0) {
			if (isMaster) cout << "--> " << (int)(dt_main/M_TO_S(sp_dt)) << " Snowpack simulations per dt_main\n";
		} else {
			if (isMaster) cout << "Time step inconsistency between Alpine3D's dt_main=" << S_TO_M(dt_main) << " [min], Snowpack's CALCULATION_STEP_LENGTH=" << sp_dt << " [min] (in io.ini). Please set Snowpack's time step to an integral fraction of dt_main!\n";
			exit(1);
		}
	}
}

inline void setStaticData(const Config &cfg, IOManager& io, DEMObject &dem, Grid2DObject &landuse, std::vector<Coords> &vec_pts)
{
	const bool isMaster = MPIControl::instance().master(); // Check if this process is the master (always true for non-parallel mode)

	if (isMaster) { //Reading DEM, LANDUSE, special PTS on master process only
		const std::string slope_algorithm = cfg.get("SLOPE_ALGORITHM", "Input", "CORRIPIO");
		dem.setDefaultAlgorithm(slope_algorithm);
		dem.setUpdatePpt((DEMObject::update_type)(DEMObject::SLOPE | DEMObject::NORMAL | DEMObject::CURVATURE));

		io.readDEM(dem);

		//cleanup points that might have elevation but no slope: this makes for much easier checks!
		dem.sanitize();
		io.readLanduse(landuse);
		//TODO: generate a map of bool saying for each pixel if it is a valid pixel or not
		//-> if (isValid(ix,iy)==true) -> calculate (eb, sn, etc)
		cerr << "[i] Read 2D Grid DEM: " << dem.getNx() << "x" << dem.getNy() << "\n";
		cerr << "[i] Read 2D Grid Landuse: " << landuse.getNx() << "x" << landuse.getNy() << "\n";
		if (!landuse.isSameGeolocalization(dem))
			throw IOException("The landuse and the provided DEM don't have the same geolocalization!", AT);

		if (cfg.keyExists("POI", "Input")) {
			io.readPOI(vec_pts);
			if (!dem.gridify(vec_pts, true)) { //keep invalid points
				if (isMaster) cerr << "[E] Some POI are invalid or outside the DEM:\n";
				for (size_t ii=0; ii<vec_pts.size(); ii++) 
					if (!vec_pts[ii].indexIsValid() && isMaster) 
						std::cout  << "[E] Point " << ii << "\t" << vec_pts[ii].toString(Coords::CARTESIAN) << "\n";
				throw InvalidArgumentException("Invalid POI, please check in the logs", AT);
			} else if (isMaster)
				std::cout << "[i] Using " << vec_pts.size() << " POI\n";
		}

		bool local_coords = false;
		cfg.getValue("COMPUTE_IN_LOCAL_COORDS", "Input", local_coords, IOUtils::nothrow);
		if (local_coords) {
			dem.llcorner.setProj("LOCAL","");
			dem.llcorner.setLocalRef(dem.llcorner.getLat(), dem.llcorner.getLon());
		}
		landuse.llcorner.copyProj(dem.llcorner);
	}

	MPIControl::instance().broadcast(dem);
	MPIControl::instance().broadcast(landuse);
	MPIControl::instance().broadcast(vec_pts);
}

inline void setModules(const Config &cfg, IOManager& io, const DEMObject &dem, const Grid2DObject &landuse, const std::vector<Coords> &vec_pts, SnowDriftA3D*& drift, EnergyBalance*& eb, SnowpackInterface*& snowpack, DataAssimilation*& da, Runoff*& runoff)
{
	const bool isMaster = MPIControl::instance().master(); // Check if this process is the master (always true for non-parallel mode)
	
	//EBALANCE
	if (enable_eb && !nocompute) {
		try {
			eb = new EnergyBalance(npebalance, cfg, dem);
		} catch(std::exception& e) {
			std::cout << "[E] Exception in EnergyBalance constructor\n";
			cout << e.what() << endl;
			throw;
		}
	}
	
	//SNOWDRIFT
	if (enable_drift && !nocompute && isMaster) {
		try {
			drift = new SnowDriftA3D(dem, cfg);
		} catch(std::exception& e) {
			std::cout << "[E] Exception in SnowDriftA3D constructor\n";
			cout << e.what() << endl;
			throw;
		}
	}

	//DATA ASSIMILATION
	if (enable_da && !nocompute && isMaster) {
		da = new DataAssimilation(io);
	}
	
	//RUNOFF
	if (enable_runoff) {
		SnowpackConfig sn_cfg( cfg ); //so we also get the default value of "THRESH_RAIN"
		const double thresh_rain = sn_cfg.get("THRESH_RAIN", "SnowpackAdvanced");
		runoff = new Runoff (cfg, dem, IOUtils::C_TO_K(thresh_rain));
	}

	//SNOWPACK
	const bool glacier_katabatic_flow = cfg.get("GLACIER_KATABATIC_FLOW", "Snowpack", false);
	if (!nocompute || glacier_katabatic_flow){
		try {
			std::string grids_requirements;
			if (enable_eb) grids_requirements = grids_requirements+" "+eb->getGridsRequirements();
			if (enable_drift) grids_requirements = grids_requirements+" "+drift->getGridsRequirements();
			if (enable_da) grids_requirements = grids_requirements+" "+da->getGridsRequirements();
			if (enable_runoff) grids_requirements = grids_requirements+" "+runoff->getGridsRequirements();
			snowpack = new SnowpackInterface(cfg, npsnowpack, dem, landuse, vec_pts, startdate, grids_requirements, restart);
		} catch(std::exception& e) {
			std::cout << "[E] Exception in SnowpackInterface constructor\n";
			cout << e.what() << endl;
			throw;
		}
	}

	//set callback pointers
	if (!nocompute) {
		if (drift) {
			snowpack->setSnowDrift(*drift);
			drift->setSnowPack(*snowpack);
			if (eb) drift->setEnergyBalance(*eb);
			if (isMaster) cout << "[i] Snowpack and Snowdrift interfaces exchanged\n";
		}

		if (eb) {
			snowpack->setEnergyBalance(*eb);
			eb->setSnowPack(*snowpack);
			if (isMaster) cout << "[i] Snowpack and Ebalance interfaces exchanged\n";
		}

		if (da) {
			snowpack->setDataAssimilation(*da);
			da->setSnowPack(*snowpack);
			if (isMaster) cout << "[i] Snowpack and Ebalance interfaces exchanged\n";
		}
		
		if (runoff) {
			snowpack->setRunoff(*runoff);
			runoff->setSnowPack(*snowpack);
			if (isMaster) cout << "[i] Snowpack and Ebalance interfaces exchanged\n";
		}
	}
}

inline void cleanDestroyAll(SnowDriftA3D*& drift, EnergyBalance*& eb, SnowpackInterface*& snowpack, DataAssimilation*& da, Runoff*& runoff)
{
	if (da) da->Destroy();
	if (eb) eb->Destroy();
	if (drift) drift->Destroy();
	
	if (da) delete da;
	if (eb) delete eb;
	if (drift) delete drift;
	if (runoff) delete runoff;
	if (snowpack) delete snowpack;
}

inline void start_message(int argc, char **argv) 
{
	MPIControl& mpicontrol = MPIControl::instance();
	
	cout << argv[0] << " " <<  A3D_VERSION << " compiled on " << __DATE__ << " " << __TIME__ << "\n";
	cout << "\tLibsnowpack " << snowpack::getLibVersion() << "\n";
	cout << "\tMeteoIO " << mio::getLibVersion() << "\n";
	if (argc==1) Usage(argv[0]);
	
	const size_t n_process = mpicontrol.size();
	const size_t n_threads = mpicontrol.max_threads();
	cout << "\nRunning as " << n_process << " process";
	if (n_process>1)  cout << "es";
	cout << " and " << n_threads << " thread";
	if (n_threads>1)  cout << "s";
	cout << " with command line:";
	for (int args=0; args<argc; args++) 
		cout << " " << argv[args];
	cout << endl;

}

inline void real_main(int argc, char **argv)
{
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
	MPIControl& mpicontrol = MPIControl::instance();
	mpicontrol.barrier(); //make sure all nodes are initialized
	if (mpicontrol.master()) {
		start_message(argc, argv);
	}

	Config cfg;
	parseCmdLine(argc, argv, cfg); //pass command line arguments, set config file for cfg
	IOManager io(cfg);
	const time_t start = time(NULL);

	//Try doing a struct "modules" with pointers on modules -> easier to pass?
	SnowDriftA3D *drift(NULL);
	EnergyBalance *eb(NULL);
	SnowpackInterface *snowpack(NULL);
	DataAssimilation* da(NULL);
	Runoff* runoff(NULL);

	DEMObject dem;
	Grid2DObject landuse;
	std::vector<Coords> vec_pts;
	
	//make sure that if the sno files are not written out, they will still be available for the POI
	const std::string meteo_outpath = cfg.get("METEOPATH", "Output");
	const bool snow_write = cfg.get("SNOW_WRITE", "Output");
	if (!snow_write) {
		cfg.addKey("SNOWPATH", "Output", meteo_outpath);
	}
	cfg.write(meteo_outpath + "/io.ini"); //backup the ini file

	try { //main integration loop
		setStaticData(cfg, io, dem, landuse, vec_pts);
		setModules(cfg, io, dem, landuse, vec_pts, drift, eb, snowpack, da, runoff);
		AlpineControl control(snowpack, drift, eb, da, runoff, cfg, dem);
		control.setNoCompute(nocompute);
		control.Run(startdate, steps);
	} catch (std::exception& e) {
		cerr << "[E] exception thrown: " << e.what() << endl;
		cleanDestroyAll(drift, eb, snowpack, da, runoff);
		exit(1);
	} catch (...) {
		cout << "[E] exception thrown!" << endl;
		cleanDestroyAll(drift, eb, snowpack, da, runoff);
		exit(1);
	}

	cleanDestroyAll(drift, eb, snowpack, da, runoff);
	const time_t  end = time(NULL);
	if (mpicontrol.master()) {
		printf("\n");
		printf(" STARTED  Alpine3D Model on %s", ctime(&start));
		printf(" ===========================================================================\n");
		printf(" FINISHED Alpine3D Model on %s", ctime(&end) );
		printf(" ===========================================================================\n");
	}
}

int main(int argc, char *argv[]) {
	try {
		real_main(argc, argv);
	} catch (const std::exception &e) {
		std::cerr << e.what() << endl;
		exit(1);
	}

	return 0;
}
