/*
 *  meteoio_timeseries
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <cstdio>
#include <string.h>
#include <map>
#include <vector>
#include <meteoio/MeteoIO.h>

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

#ifdef DEBUG_ARITHM
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif
	#ifndef __USE_GNU
		#define __USE_GNU
	#endif
	#include <fenv.h>
#endif

using namespace mio; //The MeteoIO namespace is called mio

//Global variables in this file:
static std::string cfgfile( "io.ini" );
static mio::Date dateBegin, dateEnd;
static double samplingRate = 60.;

inline void Version()
{
#ifdef _MSC_VER
	std::cout << "This version of meteoio_timeseries uses a BSD-licensed port of getopt for Visual C++. \n"
		<< "It therefore includes software developed by the University of "
		<< "California, Berkeley and its contributors." << std::endl;
#endif
	std::cout << "MeteoIO version " << mio::getLibVersion() << std::endl;
}

inline void Usage(const std::string& programname)
{
	Version();

	std::cout << "Usage: " << programname << std::endl
		<< "\t[-b, --begindate=YYYY-MM-DDTHH:MM] (e.g.:2007-08-11T09:00)\n"
		<< "\t[-e, --enddate=YYYY-MM-DDTHH:MM] (e.g.:2008-08-11T09:00 or NOW)\n"
		<< "\t[-c, --config=<ini file>] (e.g. io.ini)\n"
		<< "\t[-s, --sampling-rate=<sampling rate in minutes>] (e.g. 60)\n"
		<< "\t[-p, --progress] Show progress\n"
		<< "\t[-v, --version] Print the version number\n"
		<< "\t[-h, --help] Print help message and version information\n\n";

	std::cout << "Example: " << programname << " -c io.ini -b 1996-06-17T00:00 -e NOW\n\n";
}

inline void parseCmdLine(int argc, char **argv, std::string& begin_date_str, std::string& end_date_str, bool& showProgress)
{
	int longindex=0, opt=-1;
	bool setEnd = false;

	struct option long_options[] =
	{
		{"begindate", required_argument, 0, 'b'},
		{"enddate", required_argument, 0, 'e'},
		{"config", required_argument, 0, 'c'},
		{"sampling-rate", required_argument, 0, 's'},
		{"progress", no_argument, 0, 'p'},
		{"version", no_argument, 0, 'v'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	if (argc==1) { //no arguments provided
		Usage(std::string(argv[0]));
		exit(1);
	}

	while ((opt=getopt_long( argc, argv, ":b:e:c:s:pvh", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0:
			break;
		case 'b': {
			begin_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			break;
		}
		case 'e': {
			end_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		}
		case 'c':
			cfgfile = std::string(optarg);
			break;
		case 's':
			mio::IOUtils::convertString(samplingRate, std::string(optarg));
			break;
		case ':': //operand missing
			std::cerr << std::endl << "[E] Command line option '-" << char(opt) << "' requires an operand\n";
			Usage(std::string(argv[0]));
			exit(1);
		case 'p':
			showProgress = true;
			break;
		case 'v':
			Version();
			exit(0);
		case 'h':
			Usage(std::string(argv[0]));
			exit(0);
		case '?':
			std::cerr << std::endl << "[E] Unknown argument detected\n";
			Usage(std::string(argv[0]));
			exit(1);
		default:
			std::cerr << std::endl << "[E] getopt returned character code " <<  opt << "\n";
			Usage(std::string(argv[0]));
			exit(1);
		}
	}

	if (!setEnd) {
		std::cerr << std::endl << "[E] You must specify an enddate!\n";
		Usage(std::string(argv[0]));
		exit(1);
	}
}

static void real_main(int argc, char* argv[])
{
	bool showProgress = false;
	std::string begin_date_str, end_date_str;
	parseCmdLine(argc, argv, begin_date_str, end_date_str, showProgress);
	
	Config cfg(cfgfile);
	const double TZ = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	if (!begin_date_str.empty()) {
		mio::IOUtils::convertString(dateBegin, begin_date_str, TZ);
	}
	if (end_date_str == "NOW") { //interpret user provided end date
		dateEnd.setFromSys();
		dateEnd.setTimeZone(TZ);
		dateEnd.rnd(10, mio::Date::DOWN); //rounding 10' down
	} else {
		mio::IOUtils::convertString(dateEnd, end_date_str, TZ);
	}

	samplingRate /= 24.*60; //convert to sampling rate in days
	
	IOManager io(cfg);
	std::cout << "Powered by MeteoIO " << getLibVersion() << "\n";
	std::cout << "Reading data from " << dateBegin.toString(Date::ISO) << " to " << dateEnd.toString(Date::ISO) << "\n";

	Timer timer;
	timer.start();

	std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	std::vector< std::vector<MeteoData> > vecMeteo; //so we can keep and output the data that has been read

	size_t insert_position = 0;
	for (Date d=dateBegin; d<=dateEnd; d+=samplingRate) { //time loop
		if(showProgress) std::cout << d.toString(Date::ISO) << "\n";
		io.getMeteoData(d, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		for(size_t ii=0; ii<Meteo.size(); ii++) { //loop over all stations
			if (Meteo[ii].isNodata()) continue;
			const std::string stationID( Meteo[ii].meta.stationID );
			if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[ stationID ] = insert_position++;
				vecMeteo.push_back( std::vector<MeteoData>() ); //allocating the new station
				const size_t nr_samples = static_cast<size_t>(Optim::ceil( (dateEnd.getJulian() - d.getJulian()) / samplingRate ) + 1);
				vecMeteo[ mapIDs[stationID] ].reserve( nr_samples ); //to avoid memory re-allocations with push_back()
			}
			vecMeteo[ mapIDs[stationID] ].push_back(Meteo[ii]); //fill the data manually into the vector of vectors
		}
	}

	//In both case, we write the data out
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecMeteo);

	timer.stop();
	std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
}

int main(int argc, char** argv)
{
	setbuf(stdout, NULL); //always flush stdout
	setbuf(stderr, NULL); //always flush stderr
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	
	try {
		real_main(argc, argv);
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	return 0;
}
