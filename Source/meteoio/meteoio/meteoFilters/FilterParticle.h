#ifndef FILTERPARTICLE_H
#define FILTERPARTICLE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <meteoio/Eigen/Dense>

#include <inttypes.h> //for interaction with RNG
#include <string>
#include <vector>

namespace mio {

class FilterParticle : public ProcessingBlock {
	public:
		FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void fill_state(Eigen::VectorXd& xx, const unsigned int& param, const std::vector<MeteoData>& ivec,
		        RandomNumberGenerator& RNGU);
		void readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out);
		void vecMeteoToEigen(const std::vector<MeteoData>& vec, Eigen::VectorXd& eig, const unsigned int& param);

		typedef enum PF_FILTER_ALGORITHM {
			SIR, //only SIR is implemented so far
			SIS
		} pf_filter_algorithm;
		typedef enum PF_RESAMPLE_ALGORITHM {
			SYSTEMATIC, //only SYSTEMATIC implemented so far
			EPANECHNIKOV
		} pf_resample_algorithm;

		pf_filter_algorithm filter_alg;
		pf_resample_algorithm resample_alg; //resampling of particle paths

		unsigned int NN; //number of particles
		bool path_resampling; //has nothing to do with temporal or spatial meteo resampling

		std::string model_appendix; //naming convention for modeled data in the meteo set
		std::string model_expression; //model formula (as opposed to file input)
		double model_x0; //initial state at T=0

		struct rng_settings {
			RandomNumberGenerator::RNG_TYPE algorithm;
			RandomNumberGenerator::RNG_DISTR distribution;
			std::vector<double> parameters;
			std::vector<uint64_t> seed;
			rng_settings() : algorithm(RandomNumberGenerator::RNG_XOR),
			                 distribution(RandomNumberGenerator::RNG_GAUSS),
			                 parameters(),
			                 seed() {}
		};
		struct rng_settings rng_model; //noise generator for model function
		struct rng_settings rng_prior; //prior pdf generator

		bool add_model_noise; //when the model function is calculated, add noise?
		bool be_verbose; //output warnings/info?

};

} //namespace

#endif //FILTERPARTICLE_H
