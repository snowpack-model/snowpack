// SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <csignal>
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

/**
 * @page mio_standalone Standalone usage
 * @section mio_pc Data processing on a personal computer
 * It is possible to use MeteoIO as a standalone tool instead of as a library embedded within
 * a numerical model. For this purpose, a program is provided (compiled by default) that
 * processes data and generates timeseries, either on the command line or graphically
 * when used through the Inishell GUI.
 * 
 * @subsection meteoio_timeseries The meteoio_timeseries program
 * In the "bin" sub-directory of MeteoIO, you can find the meteoio_timeseries program. It reads
 * what to do on the data from a provided ini file and produces timeseries as output. It can take
 * several arguments to control the data processing:
 *    * the begin date, ISO formatted with the "-b" option;
 *    * the end date, ISO formatted with the "-e" option;
 *    * alternatively, a duration in days with the "-d" option;
 *    * the configuration file to use, with the "-c" option;
 *    * the output sampling rate in minutes with the "-s" option;
 *    * a progressive writing out of the data by specifying how many timesteps to buffer before writing
 * with the "-o" option;
 *    * some progress indicator with the "-p" option;
 *    * the possibility to generate spatially distributed data (gridded data) with the "--2d" option.
 * 
 * For example, in order to output data from the 1st of September 2004 at noon until the 3rd of April 2008 with 
 * half-hourly values, using the ini file "io_myStation.ini" in the "cfgfiles" sub-directory:
 * @code
 * meteoio_timeseries -c cfgfiles/io_myStation.ini -b 2004-09-01T12:00 -e 2008-04-03 -s 30
 * @endcode
 * Or to generate hourly air temperature grids (in this case, you need to define the spatial interpolations 
 * in the [Interpolations2D] section as well as list which meteorological parameters to interpolate in the [Output] 
 * section with the key "GRID2D_OUTPUT_PARAMS". Please also note that computing gridded data is much more 
 * time consuming than point-wise timeseries!):
 * @code
 * meteoio_timeseries -c cfgfiles/io_myStation.ini -b 2004-09-01T12:00 -e 2008-04-03 -s 60 --2d
 * @endcode
 * 
 * @subsection Inishell The Inishell Graphical User Interface
 * The <a href="https://code.wsl.ch/snow-models/inishell">Inishell</a> software is a tool that automatically and semantically 
 * generates GUIs for scientific models, including for MeteoIO. It makes writing and editing ini files much easier, 
 * helps reduce syntax errors (through input validation) and gives a good overview of the available options. It also 
 * often directly links a given configuration key to its online documentation. When reading an existing ini file,
 * all keys that are not recognized (as well as comments) are preserved. And it can directly run meteoio_timeseries while offering
 * a GUI to configure its options. It is therefore the recommended way of configuring MeteoIO.
 * 
 * Finally, as Inishell generates ini files, it is possible to use it to configure MeteoIO and then run the processing from the
 * command line on another computer (such as a server).
 * 
 * \image html Inishell.png "The Inishell GUI"
 * \image latex Inishell.eps "The Inishell GUI" width=0.9\textwidth
 * 
 * For more information, see (Bavay, M., Reisecker, M., Egger, T., and Korhammer, D.: 
 * <i>Inishell 2.0: semantically driven automatic GUI generation for scientific models</i>, Geosci. Model Dev., 15, 
 * 365–378, doi: <a href="https://doi.org/10.5194/gmd-15-365-2022">10.5194/gmd-15-365-2022</a>, 2022.)
 *
 * @section mio_server Data processing services on a server
 * Some features have been added to the meteoio_timeseries program in order to make it easier to manage when providing
 * data processing services on a server. So it is generally a good idea to first get familiar with meteoio_timeseries
 * on a personal computer, get familiar with Inishell to configure MeteoIO and then add the few tricks that help when running
 * MeteoIO on a server.
 * 
 * @subsection server_workflow Recommended workflow
 * It is highly recommended to create and fine tune the ini files for MeteoIO on a workstation first, relying on Inishell
 * and manually running meteoio_timeseries through Inishell. This should help correct errors and optimize the processing
 * with little efforts. It is also often a good idea to output the results in the smet format so the outputs can easily 
 * be visualized (for example with <a href="https://run.niviz.org">niViz</a>) and examined with a text editor if necessary
 * before moving to binary formats (for example, smet and netcdf both use the same ACDD metadata so even the metadata can be
 * fine tuned this way).
 * 
 * Once the processing works satisfactorily, move the ini file(s) in place to the server.
 * 
 * @subsection keeping_things_simple Keeping things simple
 * In order to leave less room for failures, it is better to keep the general workflow around MeteoIO quite simple. Several constructs
 * in MeteoIO support this effort by reducing the complexity of scripts:
 *    * any key can be \ref config_special_syntax "dynamically generated from an environment variable", so there less need to generate ini files from scripts;
 *    * keys can be defined by \ref config_special_syntax "arithmetic expressions or as references to another key";
 *    * if some ini files must still be programatically generated, please consider using the \ref config_import "IMPORT_BEFORE / IMPORT_AFTER features": 
 * this allows to keep a standard ini file and dynamically overwrite or extend it with programatically generated small ini files fragments. You can either
 * always generate an ini file that first import the base file and then adds a few keys or do it the other way around. Please note that you
 * have an unlimited number and levels of imports at your disposal!
 *    * Using the \ref config_import "IMPORT_BEFORE / IMPORT_AFTER features" can be used to define different data processing levels: for example a base
 * ini file contains everything to read the raw data and standardize it (including renaming meteorological parameters when necessary). This base
 * ini file is imported by another ini file that excludes some time periods when it is known that some sensors where malfunctioning and it is
 * itself imported by another ini file that applies some filters to the data. A last ini file imports it and defines resampling, data generators, etc
 * in order to deliver an output suitable for running numerical models. Thus generating a dataset of any given processing level is achieved by running
 * meteoio_timeseries on one of these ini files and no configuration work is duplicated.
 *    * several keys (such as the line number exclusions in the CSV plugin) are redundant in the sense that they convey the same information but
 * let you provide the said information in the easiest and most logical way for your case. For example, you can provide the lines to keep or
 * on the opposite the lines to reject.
 * 
 * @subsection dealing_with_the_unexpected Dealing with the unexpected
 * Unfortunately, bugs can always happen. Some might be MeteoIO's fault, some might be another component (such as the file system, 
 * a runaway process, etc). In order to get something as robust as possible, a few tricks are at your disposal:
 *    * the meteoio_timeseries program has a "timeout" option that kills the program after the expiration of the timeout if
 * it is still running. This could prevent processes from getting stuck waiting for I/O (for example from a network
 * drive that is not responding anymore) or to prevent a client from doing a Denial Of Service by requesting a very
 * long processing (although MeteoIO is usually very fast, the right combination of data time coverage and sampling rate
 * could achieve this).
 *    * the meteoio_timeseries program has another option to control the duration of its execution: it catches the SIGTERM 
 * signal in order to print a stack trace before exiting (on supported platforms). The goal is that before starting a 
 * new meteoio_timeseries run, you can search for meteoio_timeseries processes and kill them with SIGTERM while collecting 
 * stack traces that would help explain why they were still running (it could be that you've sent a new job too early, 
 * it could be that a running process was waiting for some I/O for example).
 *    * it is possible to compile MeteoIO with the "DEBUG_ARITHM" option (on supported platforms). In this mode, 
 * the processor is configured to throw an exception in case of arithmetic exception. If you configure the execution
 * environment to allow core dumps, this would let you open the core dump in a debugger and come directly where some
 * very wrong arithmetic was attempted so you can fix it. This is very useful if you develop your own filters or 
 * data generators as it is better to stop the processing and fix the problem than compute junk.
 *    * it is also possible to compile MeteoIO with the "LEAKS_CHECK" option. This runs slower than a normal binary must
 * still much faster than in an environment such as <a href="https://valgrind.org/">valgrind</a>, so you can let it run
 * through real, daily operational service and be informed of any memory leaks. This is very useful if you develop your own
 * enhancements to MeteoIO (or to check what the MeteoIO developers did!).
 */

//Global variables in this file:
static std::string cfgfile( "io.ini" );
static double samplingRate = IOUtils::nodata;
static size_t outputBufferSize = 0;
static unsigned int timeout_secs = 0;

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
		<< "\t[-b YYYY-MM-DDTHH:MM] Begin date, ISO formatted (e.g.:2007-08-11T09:00)\n"
		<< "\t[-e YYYY-MM-DDTHH:MM] End date, ISO formatted (e.g.:2008-08-11T09:00 or NOW)\n"
		<< "\t[-d <days>] Duration from the begin date, in days (e.g.: 30)\n"
		<< "\t[-c <ini file>] Configuration file to use (e.g. io.ini)\n"
		<< "\t[-s <minutes>] Output sampling rate, in minutes (default: 60)\n"
		<< "\t[-o <number of timesteps>] Output buffer size, in number of timesteps (e.g. 24, requires APPEND mode enabled in output plugin)\n"
		<< "\t[-p] Show progress\n"
		<< "\t[-t <seconds>] Kill the process after that many seconds if still running\n"
		<< "\t[--2d] Generate spatially interpolated data (2D grids) instead of point-wise timeseries (default: false)\n"
		<< "\t[-v] Print the version number\n"
		<< "\t[-h] Print help message and version information\n\n";

	std::cout << "If the key DATA_QA_LOGS in the [General] section is set to true, you can also provide a list of\n";
	std::cout << "meteo parameters to check for valid values in the QA_CHECK_MISSING key of the [General] section.\n";
	std::cout << "You can define the output sampling rate in the SAMPLING_RATE_MIN key of the [Output] section\n";
	std::cout << "but the command line parameter has priority.\n";
	std::cout << "Example use:\n\t" << programname << " -c io.ini -b 1996-06-17T00:00 -e NOW\n\n";
}

inline Date getDate(const std::string& date_str, const double& TZ)
{
	Date parsedDate;
	if (date_str == "NOW") { //interpret user provided start date
		parsedDate.setFromSys();
		parsedDate.setTimeZone(TZ);
		parsedDate.rnd(10, mio::Date::DOWN); //rounding 10' down
	} else {
		mio::IOUtils::convertString(parsedDate, date_str, TZ);
	}
	
	return parsedDate;
}

inline void parseCmdLine(int argc, char **argv, Config &cfg, Date& begin_date, Date& end_date, bool& showProgress, bool& generate2D)
{
	std::string begin_date_str, end_date_str;
	double duration = IOUtils::nodata;
	int longindex=0, opt=-1;
	int output2D = (generate2D)? 0 : 1;
	bool setStart = false, setEnd = false, setDuration = false;
	
	if (argc==1) { //no arguments provided
		Usage(std::string(argv[0]));
		exit(1);
	}

	const struct option long_options[] =
	{
		{"2d", no_argument, &output2D, 0},
		{nullptr, 0, nullptr, 0}
	};

	while ((opt=getopt_long( argc, argv, "b:e:d:c:s:o:t:pvh", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0: // getopt_long() set a variable
			break;
		case 'b':
			begin_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			setStart = true;
			break;
		case 'e':
			end_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		case 'd':
			if (!mio::IOUtils::convertString(duration, std::string(optarg)))
				throw ConversionFailedException("Could not parse the duration argument '"+std::string(optarg)+"'", AT);
			setDuration = true;
			break;
		case 'c':
			cfgfile = std::string(optarg);
			break;
		case 's':
			if (!mio::IOUtils::convertString(samplingRate, std::string(optarg)))
				throw ConversionFailedException("Could not parse the sampling rate argument '"+std::string(optarg)+"'", AT);
			break;
		case 'o':
			if (!mio::IOUtils::convertString(outputBufferSize, std::string(optarg)))
				throw ConversionFailedException("Could not parse the output-buffer argument '"+std::string(optarg)+"'", AT);
			break;
		case ':': //operand missing, but it seems that this check does not work properly
			std::cerr << std::endl << "[E] Command line option '-" << char(opt) << "' requires an operand\n";
			Usage(std::string(argv[0]));
			exit(1);
		case 'p': 
			showProgress = true;
			break;
		case 't': 
			if (!mio::IOUtils::convertString(timeout_secs, std::string(optarg)))
				throw ConversionFailedException("Could not parse the timeout argument '"+std::string(optarg)+"'", AT);
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
	
	//synchronize with our boolean
	generate2D = (output2D==0);

	const bool validDateRange = (setStart && setEnd && !setDuration) || (setStart && !setEnd && setDuration) || (!setStart && setEnd && setDuration);
	if (!validDateRange) {
		std::cerr << std::endl << "[E] You must specify either {startdate and enddate}, or {startdate and duration} or {enddate and duration}!\n\n";
		Usage(std::string(argv[0]));
		exit(1);
	}
	
	cfg.addFile(cfgfile);
	const double TZ = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	
	//the date range specification has been validated above
	if (!begin_date_str.empty()) begin_date = getDate( begin_date_str, TZ );
	if (!end_date_str.empty()) 
		end_date = getDate( end_date_str, TZ );
	else
		end_date = begin_date + duration;
	if (begin_date.isUndef())
		begin_date = end_date - duration;
	
	//we don't overwrite command line options if set
	if (samplingRate==IOUtils::nodata)
		samplingRate = cfg.get("SAMPLING_RATE_MIN", "Output", 60.);
	samplingRate /= 24.*60; //convert to sampling rate in days
	if (samplingRate<=0)
		throw InvalidArgumentException("The sampling rate argument must be > 0! (check both on the command line and as configuration key)", AT);
}

static void signal_handler( int signal_num ) 
{
	throw IOException("Aborting after receiving signal "+IOUtils::toString(signal_num), AT);
}

static void signals_catching(void) 
{
#ifdef _WIN32
	typedef void(*SignalHandlerPointer)(int);
	SignalHandlerPointer previousHandler = signal(SIGTERM, signal_handler);
#else
	struct sigaction catch_signal;
	catch_signal.sa_handler = signal_handler;
	sigemptyset(&catch_signal.sa_mask); // We don't want to block any other signals
	catch_signal.sa_flags = 0;
	
	sigaction(SIGTERM, &catch_signal, nullptr);
	//sigaction(SIGHUP, &catch_signal, nullptr);
#endif
}

static void validMeteoData(const std::vector<std::string>& enforce_variables, const mio::MeteoData& md)
{
	const std::string msg_head( "[DATA_QA] Missing "+md.meta.getStationID()+"::" );

	for (const std::string& var : enforce_variables) {
		if (md(var) == mio::IOUtils::nodata)
			std::cout << msg_head << var << " " << md.date.toString(mio::Date::ISO) << " [" << md.date.toString(mio::Date::ISO_WEEK) << "]\n";
	}
}

static void generatePointTimeseries(const Config& cfg, const Date& dateBegin, const Date& dateEnd, const bool& showProgress)
{
	std::cout << "Reading data from " << dateBegin.toString(Date::ISO) << " to " << dateEnd.toString(Date::ISO) << "\n";
	
	IOManager io(cfg);
	const bool data_qa = cfg.get("DATA_QA_LOGS", "General", false);
	std::vector<std::string> enforce_variables;
	if (data_qa) cfg.getValue("QA_CHECK_MISSING", "General", enforce_variables, IOUtils::nothrow);

	std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	std::vector< std::vector<MeteoData> > vecMeteo; //so we can keep and output the data that has been read

	size_t insert_position = 0;
	size_t count = 0;
	for (Date d=dateBegin; d<=dateEnd; d+=samplingRate) { //time loop
		if (showProgress) std::cout << d.toString(Date::ISO) << "\n";
		count++;
		io.getMeteoData(d, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		
		for (const MeteoData& md : Meteo) {
			if (data_qa) validMeteoData( enforce_variables, md ); //check that we have everything we need
			if (md.isNodata()) continue;
			
			const std::string stationID( md.meta.stationID );
			if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[ stationID ] = insert_position++;
				vecMeteo.push_back( std::vector<MeteoData>() ); //allocating the new station
				const size_t nr_samples = static_cast<size_t>(Optim::ceil( (dateEnd.getJulian() - d.getJulian()) / samplingRate ) + 1);
				const size_t nr_samples_buffered = (outputBufferSize > 0) ? (outputBufferSize) : (nr_samples);
				vecMeteo[ mapIDs[stationID] ].reserve( std::min(nr_samples, nr_samples_buffered) ); //to avoid memory re-allocations with push_back()
			}
			vecMeteo[ mapIDs[stationID] ].push_back(md); //fill the data manually into the vector of vectors
		}
		
		if (outputBufferSize > 0 && count%outputBufferSize == 0) {	// Check for buffered output
			std::cout << "Writing output data and clearing buffer" << std::endl;
			//handle extra parameters changing over time
			for (auto& station_data : vecMeteo) MeteoData::unifyMeteoData( station_data );
			io.writeMeteoData(vecMeteo);
			for (const MeteoData& md : Meteo) {
				const std::string stationID( md.meta.stationID );
				auto it = mapIDs.find(stationID);
				if (it != mapIDs.end()) {
					vecMeteo[it->second].clear();
					vecMeteo[it->second].shrink_to_fit();
				}
			}
		}
	}

	if (data_qa) {
		std::map<std::string, size_t>::const_iterator it_stats;
		for (it_stats=mapIDs.begin(); it_stats != mapIDs.end(); ++it_stats) std::cout << "[DATA_QA] Processing " << it_stats->first << "\n";
	}
	
	//In any case, we write the data out
	std::cout << "Writing output data" << std::endl;
	//handle extra parameters changing over time
	for(size_t ii=0; ii<vecMeteo.size(); ii++) MeteoData::unifyMeteoData( vecMeteo[ii] );
	
	io.writeMeteoData(vecMeteo);
	std::cout << "Number of timesteps: " << count << "\n";
}

static void generate2DTimeseries(const Config& cfg, const Date& dateBegin, const Date& dateEnd, const bool& showProgress)
{
	std::cout << "Reading data from " << dateBegin.toString(Date::ISO) << " to " << dateEnd.toString(Date::ISO) << "\n";
	
	//reading which parameters the user wants
	std::vector< std::string > output_variables;
	cfg.getValue("GRID2D_OUTPUT_PARAMS", "output", output_variables, IOUtils::nothrow);
	std::map< MeteoData::Parameters, MeteoGrids::Parameters > std_variables;
	std::set< std::string > extra_variables;
	for (auto& var : output_variables) {
		const size_t var_idx = MeteoData::getStaticParameterIndex( var );
		if ( var_idx != IOUtils::npos ) {
			std_variables[ static_cast<MeteoData::Parameters>(var_idx) ] = MeteoData::findGridParam( static_cast<MeteoData::Parameters>(var_idx) );
		} else {
			extra_variables.insert( var );
		}
	}
	
	IOManager io(cfg);
	
	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	dem.setUpdatePpt((DEMObject::update_type)(DEMObject::SLOPE | DEMObject::NORMAL | DEMObject::CURVATURE));
	io.readDEM(dem);
	Grid2DObject grid;
	
	//running the time loop
	size_t count = 0;
	for (Date d=dateBegin; d<=dateEnd; d+=samplingRate) { //time loop
		if (showProgress) std::cout << d.toString(Date::ISO) << "\n";
		count++;
		
		//for parameters listed in MeteoData::Parameters
		for (auto& var : std_variables) {
			io.getMeteoData(d, dem, var.first, grid);
			io.write2DGrid(grid, var.second, d);
		}
		
		for (auto& var : extra_variables) {
			io.getMeteoData(d, dem, var, grid);
			const std::string argument( var+"@"+d.toString(Date::ISO) );
			io.write2DGrid(grid, argument);
		}
	}
	
	std::cout << "Number of timesteps: " << count << "\n";
}

int main(int argc, char** argv)
{
	setbuf(stdout, nullptr); //always flush stdout
	setbuf(stderr, nullptr); //always flush stderr
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	signals_catching(); //trigger call stack print in case of SIGTERM
	
	bool generate2D = false;
	bool showProgress = false;
	Config cfg;
	Date dateBegin, dateEnd;
	parseCmdLine(argc, argv, cfg, dateBegin, dateEnd, showProgress, generate2D);
	if (timeout_secs!=0) WatchDog watchdog(timeout_secs); //set to kill itself after that many seconds
	std::cout << "Powered by MeteoIO " << getLibVersion() << "\n";
	
	Timer timer;
	timer.start();
	
	try {
		if (!generate2D) {
			generatePointTimeseries( cfg, dateBegin, dateEnd, showProgress );
		} else {
			generate2DTimeseries( cfg, dateBegin, dateEnd, showProgress );
		}
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	
	timer.stop();
	std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
	return 0;
}
