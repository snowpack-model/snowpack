#include <meteoio/MeteoIO.h>
// compile with something like:
// gcc pt_extract.cc -I ~/usr/include/ -o pt_extract -lmeteoio -L ~/usr/lib/ -lstdc++
//Then run it without any arguments from the path where the hard-coded GRID2DPATH, etc work
//Since it creates internally its own ini file, there is no need to have one on the disk.

using namespace std;
using namespace mio; //The MeteoIO namespace is called mio

//edit this function to configure pt_extract
void setINI(Config& cfg) {
	const string pts_file = "./input/surface-grids/arolla.poi";
	const string dem_file = "./input/surface-grids/arolla.dem";

	cfg.addKey("GRID2D", "Input", "ARC");
	cfg.addKey("GRID2DPATH", "Input", "./output/grids");
	cfg.addKey("DEM", "Input", "ARC");
	cfg.addKey("DEMFILE", "Input", dem_file);
	cfg.addKey("COORDSYS", "Input", "CH1903");
	cfg.addKey("TZ", "Input", "1");
	cfg.addKey("POI", "Input", "SMET");
	cfg.addKey("POIFILE", "Input", pts_file);
	cfg.addKey("METEO", "Output", "SMET");
	cfg.addKey("METEOPATH", "Output", "./output");
	cfg.addKey("TZ", "Output", "1");
	cfg.addKey("COORDSYS", "Output", "CH1903");
}

///////////////////////////// no need to edit anything here below ////////////////////////////

Date getDate(const string& name, Config& cfg) {
	double TZ;
	cfg.getValue("TZ", "Input", TZ);
	Date date;
	if( IOUtils::convertString(date, name, TZ)==false ) {
		stringstream ss;
		ss << "Can not parse date in \"" << name << "\"";
		throw ConversionFailedException(ss.str(), AT);
	}
	return date;
}

void fillMeteo(MeteoData& meteo, const Grid2DObject& grid, const string& fname, const Coords& pt) {
	size_t beg = fname.find_last_of(".");
	const string ext = fname.substr(beg+1, fname.length()-beg);

	if(ext.compare("sdp")==0) meteo(MeteoData::HS) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("tss")==0) meteo(MeteoData::TSS) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("ta")==0) meteo(MeteoData::TA) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("lwr")==0) meteo(MeteoData::ILWR) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("swr")==0) meteo(MeteoData::ISWR) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("soil")==0) meteo(MeteoData::TSG) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("rh")==0) meteo(MeteoData::RH) = grid.grid2D( pt.getGridI(), pt.getGridJ());
	if(ext.compare("psum")==0) meteo(MeteoData::PSUM) = grid.grid2D( pt.getGridI(), pt.getGridJ());
}

int main(int /*argc*/, char** /*argv*/) {
	Config cfg;
	setINI(cfg);

	//create iomanager
	IOManager io(cfg);
	Grid2DObject grid;

	//get all files from directory
	const string path = cfg.get("GRID2DPATH", "Input");
	const string pattern = "20";
	list<string> dirlist;
	FileUtils::readDirectory(path, dirlist, pattern);
	dirlist.sort();

	//get special points and create meteo vectors
	vector<Coords> pts;
	io.readPOI(pts);
	DEMObject dem;
	io.readDEM(dem);
	io.read2DGrid(grid, *dirlist.begin());
	vector< vector<MeteoData> > vecvecMeteo;
	vector<StationData> vecStation;
	for(unsigned int i=0; i<pts.size(); i++) {
		grid.gridify(pts[i]);
		StationData station;
		stringstream sname, sid;
		sname << "Point_" << pts[i].getGridI() << "_" << pts[i].getGridJ();
		station.stationName = sname.str();
		sid << "pt_" << i+1;
		station.stationID = sid.str();
		station.position = pts[i];
		const double altitude = dem.grid2D( pts[i].getGridI(), pts[i].getGridJ() );
		station.position.setAltitude(altitude, false);
		vecStation.push_back(station);
	}
	vecvecMeteo.insert(vecvecMeteo.begin(), vecStation.size(), vector<MeteoData>()); //allocation for the vectors

	vector<Date> vecDate;
	list<string>::iterator itr;
	for( itr = dirlist.begin(); itr != dirlist.end(); itr++ ) {
		io.read2DGrid(grid, *itr); //itr does not contain the path
		const Date date = getDate(*itr, cfg);

		for(unsigned int i=0; i<pts.size(); i++) { //data extraction for 1 station at a time
			vector<MeteoData>& vecMeteo = vecvecMeteo[i];
			bool found=false;
			unsigned int j=0;
			for(; j<vecMeteo.size(); j++) {
				if(vecMeteo[j].date==date) {
					found=true;
					break;
				}
			}

			if(found==false) {
				MeteoData meteo;
				meteo.meta = vecStation[i];
				meteo.date = date;
				vecMeteo.push_back(meteo);
			}

			fillMeteo(vecMeteo[j], grid, *itr, pts[i]);
		}
	}

	const string outpath = cfg.get("METEOPATH", "Output");
	io.writeMeteoData(vecvecMeteo, outpath);

	return 0;
}
