/* DO NOT USE this file until you see the usual header here */

#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

#include <ctime> //for time seed
#include <inttypes.h> //for uint64_t etc. (ULL guarantees at least 64 bits and should be ok)
#include <string> //for toString()
#include <vector> //for internal states

namespace mio {

class RngCore {
	public:
		bool hardware_seed_success; //store if hardware seed went as planned (Windows?) 
	
		RngCore();	
		virtual ~RngCore();
	
		virtual uint64_t int64() = 0;
		virtual uint32_t int32() = 0;
		virtual void getState(std::vector<uint64_t>& ovec_seed) const = 0;
		virtual void setState(const std::vector<uint64_t>& ivec_seed) = 0;
		//hardware or time seed; everyone may retrieve those from our RNG from outside:
		bool getUniqueSeed(uint64_t& store) const;

	protected: //some lower level functions
		uint64_t combine32to64(const uint32_t& low, const uint32_t& high) const;
		double doubFromInt(const uint64_t& rn) const;
		double trueDoub(); //[0, 1]

	private:
		bool getEntropy(uint64_t& store) const; //hardware seed
		uint64_t timeMixer(const time_t& tt, const clock_t& cc) const;
		uint32_t hash(const uint32_t& nn) const;
		size_t countLeadingZeros(const uint64_t& nn) const;
};

class RandomNumberGenerator : private RngCore {
	public:
		enum RNG_TYPE //computation method
		{
			RNG_XOR,
			RNG_PCG,
			RNG_MTW
		};
		
//CUSTOM_DIST step 1/6: Give your distribution a name in this enum
		enum RNG_DISTR //desired distribution, only used for doubles!
		{
			RNG_UNIFORM,
			RNG_GAUSS,
			RNG_NORMAL, //= RNG_GAUSS
			RNG_GAMMA
		};
		enum RNG_BOUND //return uniform double respecting these boundaries
		{
			RNG_AINCBINC, //[0, 1]
			RNG_AINCBEXC, //[0, 1)
			RNG_AEXCBINC, //(0, 1]
			RNG_AEXCBEXC  //(0, 1)
		};

		//distribution_params's default is an empty vector, which means choose default params
		RandomNumberGenerator(const RNG_TYPE& type = RNG_XOR, const RNG_DISTR& distribution = RNG_UNIFORM,
    		    const std::vector<double>& distribution_params = std::vector<double>());
		RandomNumberGenerator(const RandomNumberGenerator& rng);
		virtual ~RandomNumberGenerator();

		RandomNumberGenerator& operator=(const RandomNumberGenerator& rng);

		uint64_t int64();
		uint32_t int32();
		double doub(); //we keep this separate for speed
		double doub(const RNG_BOUND& bounds, const bool& true_double = false);
		double draw(); //alias for uniform double

		double pdf(const double& xx); //probability density function
		double cdf(const double& xx); //cumulative probability distribution function

		uint64_t range64(const uint64_t& aa, const uint64_t& bb); //[a, b]
		uint32_t range32(const uint32_t& aa, const uint32_t& bb); //[a, b]
		bool trueRange32(const uint32_t& aa, const uint32_t& bb, uint32_t& result,
		    const unsigned int& nmax = 1e6); //[a, b]
		
		void getState(std::vector<uint64_t>& ovec_seed) const;
		void setState(const std::vector<uint64_t>& ivec_seed);

		RNG_DISTR getDistribution(std::vector<double>& vec_params) const;
		void setDistribution(const RNG_DISTR& distribution, const std::vector<double>& vec_params =
		    std::vector<double>()); //construct empty vector as default
		double getDistributionParameter(const std::string& param_name) const;
		void setDistributionParameter(const std::string& param_name, const double& param_val);
		
		bool getHardwareSeedSuccess() const;
		bool getUniqueSeed(uint64_t& store) const; //allow for outside calls to the seeding function
		std::string toString();

	private:
		RngCore* rng_core; //generator algorithm
		RNG_TYPE rng_type; //for output only so far
		RNG_DISTR rng_distribution;
		std::vector<double> DistributionParameters; //anything needed by the distributions can be stored here
		
		bool rng_muller_generate; //bookkeeping Box-Muller transform
		double rng_muller_z1; //cache
		//(tradeoff between readability with the vector and speed with globals)
		
		double (RandomNumberGenerator::*doubFunc)(); //double random numbers algorithm for distribution
		double (RandomNumberGenerator::*pdfFunc)(const double& xx) const; //probability density function
		double (RandomNumberGenerator::*cdfFunc)(const double& xx) const; //cumulative distribution function
		
//CUSTOM_DIST step 2/6: Add your distribution function, its pdf and cdf here, matching exactly this type: 
		double doubUniform();
		double pdfUniform(const double& xx) const;
		double cdfUniform(const double& xx) const;
		double doubGauss(); //=normal
		double pdfGauss(const double& xx) const;
		double cdfGauss(const double& xx) const;
};

class RngXor : public RngCore {
	public: //new generators must provide these
		RngXor();
		uint64_t int64();
		uint32_t int32();
		void getState(std::vector<uint64_t>& ovec_seed) const;
		void setState(const std::vector<uint64_t>& ivec_seed);

	private:
		uint64_t state;
		uint64_t uu, vv, ww;

		bool initAllStates();
};

class RngPcg : public RngCore {
	public:
		RngPcg();
		uint64_t int64();
		uint32_t int32( );
		void getState(std::vector<uint64_t>& ovec_seed) const;
		void setState(const std::vector<uint64_t>& ivec_seed);

	private:
		uint64_t state;
		uint64_t inc;

		bool initAllStates();
};

class RngFactory { //factory for the generator algorithm
	public: //(so that memory dedicated to the states lives only as long as the RNG)
		static RngCore* getCore(const RandomNumberGenerator::RNG_TYPE& algorithm);
};

} //namespace

#endif

