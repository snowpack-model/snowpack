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
#include <string.h>
#include <map>
#include <vector>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

static std::vector< std::vector<MeteoData> > read_data(const std::string& io_file, const bool& write_data)
{
	Config cfg(io_file);
	const double TZ = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone

	const Date dateBegin(2008, 12, 1, 0, 0, TZ);
	const Date dateEnd(2009, 1, 31, 23, 0, TZ);
	const double samplingRate = 1./24.;
	
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

	//write the data out
	if (write_data) {
		std::cout << "Writing output data" << std::endl;
		io.writeMeteoData(vecMeteo);
	}

	timer.stop();
	std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
	
	return vecMeteo;
}

static void compare_data(const std::vector<MeteoData> &ref, const std::vector<MeteoData> &test)
{
	if (test == ref) return;
	
	if (test.size() != ref.size()) {
		const size_t test_size = test.size();
		const size_t ref_size = ref.size();
		std::cout << "Station '" << ref.front().meta.getStationID() << "' does not has the same length between REF (" << ref_size << ") and TEST (" << test_size << ")\n";
		if (test_size!=0 && ref_size!=0)
			std::cout << "Station '" << ref.front().meta.getStationID() << "' covers " << ref.front().date.toString(Date::ISO) << " - " << ref.back().date.toString(Date::ISO) << " versus " << test.front().date.toString(Date::ISO) << " - " << test.back().date.toString(Date::ISO) << "\n";
		return;
	}
	
	std::set< std::string > diff_params;
	size_t diff_nr_params = 0;
	for (size_t ii=0; ii<ref.size(); ii++) {
		if (test[ii]==ref[ii]) continue;
		
		if (test[ii].getNrOfParameters() != ref[ii].getNrOfParameters()) {
			diff_nr_params++;
			continue;
		}
		
		for (size_t paridx=0; paridx<ref[ii].getNrOfParameters(); paridx++) {
			if (test[ii].date != ref[ii].date)
				std::cout << "Mis-matched dates at REF=" << ref[ii].date.toString(Date::ISO) << "\n";
			
			if (test[ii].meta != ref[ii].meta)
				std::cout << "Mis-matched metadata at " << ref[ii].date.toString(Date::ISO) << "\n";
			
			if (test[ii](paridx) != ref[ii](paridx)) {
				std::cout << "At " << ref[ii].date.toString(Date::ISO) << " " << ref[ii](paridx) << " != " << test[ii](paridx) << "\n";
				diff_params.insert( ref[ii].getNameForParameter(paridx) );
			}
		}
	}
	
	std::cout << "Station '" << ref.front().meta.getStationID() << "' has differences between REF and TEST on the parameters: ";
	for (auto paramName : diff_params) std::cout << " " << paramName;
	std::cout << "\n";
}

static bool compare_data(const std::vector< std::vector<MeteoData> > &ref, const std::vector< std::vector<MeteoData> > &test)
{
	if (test==ref) return true;
	if (test.size()==0) return true; //we enforce ref.size()==test.size() later
	
	//check number of stations and order
	std::vector<std::string> vecRefIds;
	for (size_t ii=0; ii<ref.size(); ii++) vecRefIds.push_back( ref[ii].front().meta.getStationID() );
	std::vector<std::string> vecTestIds;
	for (size_t ii=0; ii<test.size(); ii++) vecTestIds.push_back( test[ii].front().meta.getStationID() );
	if (vecTestIds != vecRefIds) {
		if (vecTestIds.size() != vecRefIds.size())
			std::cout << "Not the same number of stations in REF and TEST\n";
		else
			std::cout << "Not the same content / order for the stations in REF and TEST\n";
		
		std::cout << "Ref stations:";
		for (size_t ii=0; ii<ref.size(); ii++) std::cout << " " << ref[ii].front().meta.getStationID();
		std::cout << "\n";
		std::cout << "Test stations:";
		for (size_t ii=0; ii<test.size(); ii++) std::cout << " " << test[ii].front().meta.getStationID();
		std::cout << "\n";
		
		return false;
	}
	
	//find which stations have differences
	for (size_t ii=0; ii<ref.size(); ii++) { //test and ref have the same size and stations' ordering
		compare_data(ref[ii], test[ii]);
	}
	
	return false;
}

int main()
{
	setbuf(stdout, nullptr); //always flush stdout
	setbuf(stderr, nullptr); //always flush stderr

	try {
		//the boolean flag can be set to TRUE in order to write the data out
		std::vector< std::vector<MeteoData> > vecMeteoRef( read_data("io_ref_data.ini", false) );
		std::vector< std::vector<MeteoData> > vecMeteoTest( read_data("io.ini", false) );
		
		const bool status = compare_data(vecMeteoRef, vecMeteoTest);
		if (status) {
			std::cout << "All OK\n";
			exit(0);
		} else {
			std::cout << "REF and TEST datasets differ!!\n";
			exit(1);
		}
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	return 0;
}
