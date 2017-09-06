#include <stdio.h>
#include <stdlib.h>
#include <meteoio/MeteoIO.h>

int main(int /*argc*/, char** argv) {
	mio::SunObject Sun(46.77181, 9.86820, 2192.); //Stillberg station
	const double slope_azi=38., slope_elev=35.;   //Stillberg station
	const double TZ=1.;
	const double TA = 273.15+11., RH = 0.5, mean_albedo = 0.5;

	mio::Date d1;
	d1.setTimeZone(TZ);
	mio::IOUtils::convertString(d1,argv[1], TZ); //get date and time
	Sun.setDate(d1.getJulian(), d1.getTimeZone()); //local julian date and timezone

	double iswr_ref;
	mio::IOUtils::convertString(iswr_ref,argv[2]); //get measured global incoming radiation

	Sun.calculateRadiation(TA, RH, mean_albedo);
	std::cout << Sun.toString();

	//radiation in the beam
	double B_toa, B_direct, B_diffuse;
	Sun.getBeamRadiation(B_toa, B_direct, B_diffuse);
	const double Md_beam=Sun.getSplitting(B_toa,iswr_ref);

	//radiation on the horizontal
	double H_toa, H_direct, H_diffuse;
	Sun.getHorizontalRadiation(H_toa, H_direct, H_diffuse);
	const double Md_horizontal=Sun.getSplitting(H_toa,iswr_ref);

	//radiation on the slope
	double S_toa, S_direct, S_diffuse;
	Sun.getSlopeRadiation(slope_azi, slope_elev, S_toa, S_direct, S_diffuse);
	const double Md_slope=Sun.getSplitting(S_toa,iswr_ref);

	//outputing the results for beam, horizontal and slope
	const unsigned int colw=10;
	std::cout << std::setw(colw) << "Slope" << std::fixed << std::setw(colw) << std::setprecision(1) << S_toa;
	std::cout << std::fixed << std::setw(colw) << std::setprecision(1) << S_direct;
	std::cout << std::fixed << std::setw(colw) << std::setprecision(1) << S_diffuse;
	std::cout << std::fixed << std::setw(colw) << std::setprecision(1) << S_direct+S_diffuse << "\n\n";
	printf("Slope: azi=%g angle=%g ",slope_azi,slope_elev);
	printf("Angle of incidence=%g\n",Sun.position.getAngleOfIncidence(slope_azi,slope_elev));
	printf("Beam Splitting coefficient=%g\n",Md_beam);
	printf("Horizontal Splitting coefficient=%g\n",Md_horizontal);
	printf("Slope Splitting coefficient=%g\n",Md_slope);

	return 0;
}
