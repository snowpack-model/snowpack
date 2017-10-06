#include <meteoio/MeteoIO.h>
#include <cstdio>

using namespace std;
using namespace mio; //The MeteoIO namespace is called mio

int main() {
	//we create one point in the CH1903 coordinate system
	//and we attribute it a set of coordinates (this is the Laret lake in Davos)
	Coords point1("CH1903","");
	point1.setXY(785425. , 191124., 1400.);

	//First, we check the matching lat/long
	printf("CH1903 to Lat/long\n");
	printf("\t(%g , %g) -> (%g , %g)\n", point1.getEasting(), point1.getNorthing(), point1.getLat(), point1.getLon());

	//Now we set up a local coordinate system, using as (0,0) the old observatory of Berne
	//We use the Vincenty algorithm for the distance calculations (more precise but more intensive than the default)
	point1.setProj("LOCAL","(46°57'08.66\",7°26'22.50\")");
	point1.setDistances(Coords::GEO_VINCENTY);

	//we print the new (x,y) coordinates
	//we add (600000,200000) because in CH1903, the old observatory of Berne is (600000,200000)
	printf("Lat/long to Local\n\t(%g , %g) -> (%g , %g)\n", point1.getLat(), point1.getLon(), point1.getEasting()+600000., point1.getNorthing()+200000.);

	//we want again CH1903 coodinates -> we set up the projection to CH1903
	point1.setProj("CH1903","");
	//we print the (lat,long) and (x,y) coordinates once more to check that they still match what we had a the begining
	printf("Lat/long to CH1903\n\t(%g , %g) -> (%g , %g)\n", point1.getLat(), point1.getLon(), point1.getEasting(), point1.getNorthing());

	//A nice finishing touch: we print nicely formatted lat/lon
	cout << "Pretty printing\n";
	cout << "\t(" << point1.getLat() << " , " << point1.getLon() << ") = " << point1.toString(Coords::LATLON) << endl;

	return 0;
}
