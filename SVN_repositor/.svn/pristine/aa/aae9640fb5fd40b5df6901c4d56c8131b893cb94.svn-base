#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio;

//This is the a basic example of time manipulation
//provide date as ISO formatted, for example 2008-12-01T15:35:00
int main(int argc, char** argv) {
	(void)argc;
	const double TZ=+1.;

	Date now;
	now.setFromSys();
	//now.setTimeZone(TZ);
	std::cout << "now=" << now.toString() << "\n";
	now.rnd(1800, Date::DOWN);
	std::cout << "Rounded now=" << now.toString() << "\n";

	Date d1;
	IOUtils::convertString(d1,argv[1], 0);
	std::cout << "In timezone GMT+0:\n";
	std::cout << d1.toString() << "\n";

	std::cout << "In timezone GMT" << std::showpos << TZ << std::noshowpos << ":\n";
	d1.setTimeZone(TZ,false);
	std::cout << d1.toString() << "\n";

	std::cout << "Same, directly read in timezone GMT" << std::showpos << TZ << std::noshowpos << ":\n";
	d1.setTimeZone(TZ,false);
	IOUtils::convertString(d1,argv[1], TZ);
	std::cout << d1.toString() << "\n";

	std::cout << "And swapped back to timezone GMT+0:\n";
	d1.setTimeZone(0.,false);
	std::cout << d1.toString() << "\n";

	return 0;
}
