#include <meteoio/MeteoIO.h>

using namespace std;
using namespace mio;

int main() {
	Coords point1("CH1903","");
	point1.setXY(785425. , 191124., 1400.);

	//First, we check the matching lat/long
	cout << "CH1903 to Lat/long\n";
	if(!IOUtils::checkEpsilonEquality(point1.getLat(), 46.8453968543, IOUtils::lat_epsilon) || !IOUtils::checkEpsilonEquality(point1.getLon(), 9.87000605412, IOUtils::lon_epsilon)) {
		cerr << setprecision(12) << "calculated lat=" << point1.getLat() << " calculated lon=" << point1.getLon() << "\n";
		exit(1);
	}

	//Now we set up a local coordinate system, using as (0,0) the old observatory of Berne
	//We use the Vincenty algorithm for the distance calculations (more precise but more intensive than the default)
	point1.setProj("LOCAL","(46°57'08.66\",7°26'22.50\")");
	point1.setDistances(Coords::GEO_VINCENTY);

	//we print the new (x,y) coordinates
	//we add (600000,200000) because in CH1903, the old observatory of Berne is (600000,200000)
	cout << "LOCAL to Easting/Northing\n";
	const double local_x=point1.getEasting()+600000., local_y=point1.getNorthing()+200000.;
	if(!IOUtils::checkEpsilonEquality(local_x, 785351.221299, 1e-4) || !IOUtils::checkEpsilonEquality(local_y, 190975.590555, 1e-4)) {
		cerr << setprecision(12) << "calculated Easting=" << local_x << " calculated northing=" << local_y << "\n";
		exit(1);
	}

	//we want again CH1903 coodinates -> we set up the projection to CH1903
	cout << "Lat/long to CH1903\n";
	point1.setProj("CH1903","");
	const double ch_x=point1.getEasting(), ch_y=point1.getNorthing();
	if(!IOUtils::checkEpsilonEquality(ch_x, 785425.055232, 1e-4) || !IOUtils::checkEpsilonEquality(ch_y, 191123.808706, 1e-4)) {
		cerr << setprecision(12) << "calculated Easting=" << ch_x << " calculated northing=" << ch_y << "\n";
		exit(1);
	}

	//A nice finishing touch: we print nicely formatted lat/lon
	cout << "Pretty printing\n";
	stringstream ss;
	ss << CoordsAlgorithms::printLatLon(point1.getLat(), point1.getLon());
	if(ss.str() != string("(46°50'43.428675\" , 9°52'12.021795\")")) {
		cerr << "Pretty printing failed: " << ss.str() << "\n";
		exit(1);
	}

	//Now to UTM
	cout << "Lat/long to UTM\n";
	point1.setProj("UTM","32T");
	const double utm_x=point1.getEasting(), utm_y=point1.getNorthing();
	if(!IOUtils::checkEpsilonEquality(utm_x, 566333.101179, 1e-4) || !IOUtils::checkEpsilonEquality(utm_y, 5188351.26369, 1e-4)) {
		cerr << setprecision(12) << "calculated Easting=" << utm_x << " calculated northing=" << utm_y << "\n";
		exit(1);
	}

	//Now to UTM
	cout << "UTM to Lat/long\n";
	point1.setProj("UTM","32T");
	point1.setXY(566333.101179 , 5188351.26369, 1400.);
	const double utm_lat=point1.getLat(), utm_lon=point1.getLon();
	if(!IOUtils::checkEpsilonEquality(utm_lat, 46.8453968533, 1e-4) || !IOUtils::checkEpsilonEquality(utm_lon, 9.87000605532, 1e-4)) {
		cerr << setprecision(12) << "calculated lat=" << utm_lat << " calculated lon=" << utm_lon << "\n";
		exit(1);
	}

	//Now to UPS
	cout << "Lat/long to UPS\n";
	point1.setProj("UPS","N");
	point1.setLatLon(84.2872338889, -132.247989167, 0.);
	const double ups_x=point1.getEasting(), ups_y=point1.getNorthing();
	if(!IOUtils::checkEpsilonEquality(ups_x, 1530125.78038, 1e-4) || !IOUtils::checkEpsilonEquality(ups_y, 2426773.59547, 1e-4)) {
		cerr << setprecision(12) << "calculated Easting=" << ups_x << " calculated northing=" << ups_y << "\n";
		exit(1);
	}

	//Now to UPS
	cout << "UPS to Lat/long\n";
	point1.setProj("UPS","S");
	point1.setXY(2500000., 1500000., 0.);
	const double ups_lat=point1.getLat(), ups_lon=point1.getLon();
	if(!IOUtils::checkEpsilonEquality(ups_lat, -83.6373175, 1e-4) || !IOUtils::checkEpsilonEquality(ups_lon, 135., 1e-4)) {
		cerr << setprecision(12) << "calculated lat=" << ups_lat << " calculated lon=" << ups_lon << "\n";
		exit(1);
	}

	cout << "Pretty printing of negative latitudes\n";
	stringstream ss2;
	ss2 << point1.toString(Coords::LATLON);
	if(ss2.str() != string("(-83.637318 , 135.000000)")) {
		cerr << "Pretty printing failed: " << ss2.str() << "\n";
		exit(1);
	}

	return 0;
}
