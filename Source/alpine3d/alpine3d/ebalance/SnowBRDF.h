#ifndef SNOWBRDF_H
#define SNOWBRDF_H

#include <meteoio/MeteoIO.h>

class SnowBRDF{

	public:
		SnowBRDF(){};
		SnowBRDF(const mio::Config& cfg);
		~SnowBRDF();

		double get_RF(double cth_i, double cphi, double cth_v);
		bool albedo_inc=false;

	private:

		void initialize(const mio::Config& cfg);
		double BRDF_data[33][120][60];
};

#endif
