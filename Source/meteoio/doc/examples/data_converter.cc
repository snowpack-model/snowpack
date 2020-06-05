/*This example takes two ISO-formatted dates and a sampling rate (in minutes) on the command line
* for example ./data_converter 2008-12-01T00:00:00 2008-12-31T23:00 60
* It will retrieve the data for this time interval and write it out once per 60 minutes as specified
* in the io.ini configuration
*/
#include <iostream>
#include <cstdio>
#include <string.h>
#include <map>
#include <vector>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

static void real_main(int argc, char** argv)
{
	setbuf(stdout, NULL); //always flush stdout
	const std::string config_filename = (argc==5)? argv[4] : "io.ini";
	Config cfg(config_filename);
	const double TZ = cfg.get("TIME_ZONE", "Input");

	double Tstep;
	IOUtils::convertString(Tstep, argv[3]);
	Tstep /= 24.*60; //convert to sampling rate in days

	Date d_start, d_end;
	IOUtils::convertString(d_start, argv[1], TZ);
	IOUtils::convertString(d_end, argv[2], TZ);

	IOManager io(cfg);
	std::cout << "Powered by MeteoIO " << getLibVersion() << "\n";
	std::cout << "Reading data from " << d_start.toString(Date::ISO) << " to " << d_end.toString(Date::ISO) << "\n";

	Timer timer;
	timer.start();

	std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	std::vector< std::vector<MeteoData> > vecMeteo; //so we can keep and output the data that has been read

	size_t insert_position = 0;
	for (Date d=d_start; d<=d_end; d+=Tstep) { //time loop
		io.getMeteoData(d, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		for(size_t ii=0; ii<Meteo.size(); ii++) { //loop over all stations
			if (Meteo[ii].isNodata()) continue;
			const std::string stationID( Meteo[ii].meta.stationID );
			if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[ stationID ] = insert_position++;
				vecMeteo.push_back( std::vector<MeteoData>() ); //allocating the new station
				const size_t nr_samples = static_cast<size_t>(Optim::ceil( (d_end.getJulian() - d.getJulian()) / Tstep ) + 1);
				vecMeteo[ mapIDs[stationID] ].reserve( nr_samples ); //to avoid memory re-allocations with push_back()
			}
			vecMeteo[ mapIDs[stationID] ].push_back(Meteo[ii]); //fill the data manually into the vector of vectors
		}
	}

	//io.getMeteoData(d_start, d_end, vecMeteo); //This would be the call that does NOT resample the data, instead of the above "for" loop

	//In both case, we write the data out
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecMeteo);

	timer.stop();
	std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
}

int main(int argc, char** argv)
{
	if(argc<4 || argc>5) {
		std::cout << "Invalid number of arguments! Please provide a date range and a sampling rate (in minutes)\n";
		exit(0);
	}

	try {
		real_main(argc, argv);
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	return 0;
}
