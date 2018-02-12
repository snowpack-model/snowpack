/* This is a test program to run single time steps simulations of the katabatic flow temperature correction
 * compile with something like:
 * gcc glaciers.cc -g -I ~/usr/include/ -o glaciers -lmeteoio -lalpine3d -L ~/usr/lib/ -lstdc++
 */

#include <meteoio/MeteoIO.h>
#include <alpine3d/Alpine3D.h>

using namespace mio;

const Grid2DObject computeGlacierMap(const Grid2DObject& lus)
{
	Grid2DObject map(lus); //only glaciated pixels should be IOUtils::nodata
	for (size_t ii=0; ii<lus.size(); ii++) {
		if (map(ii)==IOUtils::nodata) {
			map(ii) = 0;
			continue;
		}
		const int lus_code = ((int)(map(ii) + 0.0001) - 10000) / 100;
		if (lus_code==14)
			map(ii) = IOUtils::nodata;
	}
	
	return map;
}

void setPOI(const DEMObject& dem, std::vector<Coords>& pts, std::vector<StationData> &vecStation)
{
	vecStation.clear();
	
	for(unsigned int i=0; i<pts.size(); i++) {
		dem.gridify(pts[i]);
		StationData station;
		std::stringstream sname, sid;
		sname << "Point_" << pts[i].getGridI() << "_" << pts[i].getGridJ();
		station.stationName = sname.str();
		sid << "pt_" << i+1;
		station.stationID = sid.str();
		station.position = pts[i];
		const double altitude = dem.grid2D( pts[i].getGridI(), pts[i].getGridJ() );
		station.position.setAltitude(altitude, false);
		vecStation.push_back(station);
	}
}

int main(int /*argc*/, char** argv) {
	Config cfg( "io.ini" );
	cfg.addKey("GRID2DPATH", "Output", "./");
	IOManager io( cfg );
	
	DEMObject dem;
	dem.setUpdatePpt(DEMObject::SLOPE);
	io.readDEM(dem);
	
	//read and initialize the POIs
	std::vector<StationData> vecStation;
	std::vector<Coords> pts;
	io.readPOI(pts);
	setPOI(dem, pts, vecStation);
	
	Grid2DObject lus;
	io.read2DGrid(lus, "arolla_rad.lus");
	
	Glaciers glacier(cfg, dem);
	glacier.setGlacierMap( computeGlacierMap(lus) );
	/*Grid2DObject alt, dist;
	glacier.getGrids(alt, dist);
	io.write2DGrid(computeGlacierMap(lus), "lus.asc");
	io.write2DGrid(alt, "alt.asc");
	io.write2DGrid(dist, "dist.asc");
	exit(0);*/
	
	const double TZ = cfg.get("TIME_ZONE", "Input");
	Date d1, d2;
	IOUtils::convertString(d1, argv[1], TZ);
	IOUtils::convertString(d2, argv[2], TZ);
	const double Tstep = 1./24.;
	
	//loop over the timesteps and fill the POIs
	std::vector< std::vector<MeteoData> > vecvecMeteo;
	vecvecMeteo.insert(vecvecMeteo.begin(), vecStation.size(), std::vector<MeteoData>()); //allocation for the vectors
	
	for(; d1<=d2; d1+=Tstep) { //time loop
		Grid2DObject ta;
		io.getMeteoData(d1, dem, MeteoData::TA, ta);
		Grid2DObject ta_corr(ta);
		glacier.correctTemperatures(ta_corr);
		
		for(unsigned int ii=0; ii<pts.size(); ii++) { //data extraction for 1 station at a time
			MeteoData meteo(d1, vecStation[ii]);
			meteo(MeteoData::TA) = ta.grid2D( pts[ii].getGridI(), pts[ii].getGridJ());
			meteo.addParameter("TA_corr");
			meteo("TA_corr") = ta_corr.grid2D( pts[ii].getGridI(), pts[ii].getGridJ());
			vecvecMeteo[ii].push_back( meteo );
		}
	}

	//write out the POIs' timeseries
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecvecMeteo);
}







