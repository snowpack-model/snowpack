#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

// ----- Constants ----
static const int minutes_per_day = 60 * 24;

mio::SunObject Sun(46.77181, 9.86820, 2192.); //Stillberg station
static const double slope_azi=38., slope_elev=35.;   //Stillberg station
static const double TZ=1.;
static const double TA = 273.15+11., RH = 0.5, mean_albedo = 0.5;
static const std::string f_output("curr_output.txt");

// ---- INPUT Data to controll----
// if you add an date here, don't forget to add it to list with push_pack in main !!!!
static const mio::Date d1(2000, 12, 16, 01, 0, 0, TZ); // Winter, good Weather
static const mio::Date d2(2000, 12, 26, 01, 0, 0, TZ); // Winter, "bad" Weather
static const mio::Date d3(2001, 05, 29, 01, 0, 0, TZ); // Summer, good Weather
static const mio::Date d4(2001, 05, 26, 01, 0, 0, TZ); // Summer, "bad" Weather
static const mio::Date d5(2001, 04, 14, 01, 0, 0, TZ); // Higher income then possible
static const mio::Date d6(2001, 05, 22, 01, 0, 0, TZ); // double peaks

static const double iswr_ref []= {-1000., -100., -10., -1., -0.1, -0.01, 0, 0.01, 0.1, 1., 10., 100., 1000};

// write out ref file
bool writeSun24h(ofstream& os, const mio::Date start_date, const double i_iswr_ref) 
{
	for (mio::Date date(start_date); date <= (start_date+1.); date+=(10./minutes_per_day)) { //every 10 minutes
		os << date.toString(Date::ISO) << "\t";
		os << std::setprecision(10);

		Sun.setDate(date.getJulian(), date.getTimeZone()); //local julian date and timezone
		Sun.calculateRadiation(TA, RH, mean_albedo);

		//radiation in the beam'
		double b_toa, b_direct, b_diffuse, md_beam;
		Sun.getBeamRadiation(b_toa, b_direct, b_diffuse);
		md_beam=Sun.getSplitting(b_toa,i_iswr_ref);

		os << b_toa << "\t" << b_direct << "\t" << b_diffuse << "\t";
		os << md_beam << "\t";

		//radiation on the horizontal
		double h_toa, h_direct, h_diffuse, md_horizontal;
		Sun.getHorizontalRadiation(h_toa, h_direct, h_diffuse);
		md_horizontal=Sun.getSplitting(h_toa,i_iswr_ref);

		os << h_toa << "\t" << h_direct << "\t" << h_diffuse << "\t";
		os << md_horizontal << "\t";

		//radiation on the slope
		double s_toa, s_direct, s_diffuse, md_slope;
		Sun.getSlopeRadiation(slope_azi, slope_elev, s_toa, s_direct, s_diffuse);
		md_slope=Sun.getSplitting(s_toa,i_iswr_ref);

		os << s_toa << "\t" << s_direct << "\t" << s_diffuse << "\t";
		os << md_slope << "\t";

		//other sun stuff
		double solar_azi, solar_elev, eccentricity;
		double sunrise, sunset, tmp_daylight;
		Sun.position.getHorizontalCoordinates(solar_azi, solar_elev, eccentricity);
		Sun.position.getDaylight(sunrise, sunset, tmp_daylight);

		os << Sun.getElevationThresh() << "\t";
		os << Sun.position.getAngleOfIncidence(slope_azi,slope_elev) << "\t";
		os << eccentricity << "\t";
		//os << Hour Angle << "\t";
		os << solar_azi << "\t" << solar_elev << "\t";
		os << Date::printFractionalDay(sunrise)<< "\t";
		os << Date::printFractionalDay(sunset)<< "\t";
		os << Date::printFractionalDay(tmp_daylight/minutes_per_day)<< "\t";

		// end line
		os << endl;
	}
	return true;
}


// print out header to know which line is which
void printHeader(ostream& os){
	os << "# date" << "\t" ;
	os << "r. beam total" << "\t"<< "r. beam direct" << "\t"<< "r. beam diffues" << "\t"<< "r. beam splitting" << "\t";
	os << "r. horizontal total" << "\t"<< "r. horizontal direct" << "\t"<< "r. horizontal diffuse" << "\t"<< "r. horizontal splitting" << "\t";
	os << "r. slope total" << "\t"<< "r. slope direct" << "\t"<< "r. slope diffuse" << "\t"<< "r. slope splitting" << "\t";
	os << "Elevation Thresh" << "\t"<< "Eccentricity " << "\t"<< "solar Azimut" << "\t"<< "solar elevation" << "\t"<< "Sunrise Time" << "\t"<< "Sunset Time" << "\t"<< "Sunglight" << "\t";
	os << endl;
}

//Test if sun simulation at Stilberg station at different dates and different iswr_ref
int main() 
{
	// ----- Cenerate list for loops --------
	cout << " --- Init Variables \n";

	std::list<mio::Date> date;
	date.push_back(d1);
	date.push_back(d2);
	date.push_back(d3);
	date.push_back(d4);
	date.push_back(d5);
	date.push_back(d6);
	
	std::list<double> iswr(iswr_ref, iswr_ref + sizeof(iswr_ref) / sizeof(double));

	// ----- Write reference file ------
	std::cout << " --- Start writing Output file \"" << f_output << "\"\n";
	std::ofstream ofs(f_output.c_str(), ofstream::out | ofstream::trunc);

	printHeader(ofs);
	for (list<mio::Date>::iterator it_date = date.begin(); it_date != date.end(); it_date++) {
		std::cout << " -- read information for date : " << (*it_date).toString(Date::ISO) << "\n -- Processing reference iswr [";
		for (list<double>::iterator it_iswr = iswr.begin(); it_iswr != iswr.end(); it_iswr++) {
			cout << " " << *it_iswr;
			if(!writeSun24h(ofs, *it_date, *it_iswr)){
				exit(1);
			}
		}
		cout << " ]" << endl;
	}
	ofs.close();

	return 0;
}
