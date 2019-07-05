#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/tinyexpr.h>

#include <sstream> //for readLineToVec

namespace mio {

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
          model_appendix("_MOD"), model_expression(""), model_x0(IOUtils::nodata),
          rng_model(), rng_prior(), add_model_noise(true),
          be_verbose(true)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void FilterParticle::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{
	const size_t TT = ivec.size(); //number of time steps

	//init random number generators:
	RandomNumberGenerator RNGU(rng_model.algorithm, rng_model.distribution, rng_model.parameters);
	if (!rng_model.seed.empty()) //seed integers given in ini file
		RNGU.setState(rng_model.seed);
	RandomNumberGenerator RNG0(rng_prior.algorithm, rng_prior.distribution, rng_prior.parameters);
	if (!rng_prior.seed.empty())
		RNGU.setState(rng_prior.seed);

	Eigen::VectorXd xx(TT); //state / model
	Eigen::VectorXd zz(TT); //observation
	vecMeteoToEigen(ivec, zz, param);
	fill_state(xx, param, ivec, RNGU); //get state/model values

	ovec = ivec;
}

void FilterParticle::fill_state(Eigen::VectorXd& xx, const unsigned int& param, const std::vector<MeteoData>& ivec,
        RandomNumberGenerator& RNGU)
{
	const size_t TT = ivec.size(); //number of time steps

	const bool has_model_expression = !model_expression.empty(); //user has entered a model equation, valid or not

	//check if model data exists in meteo data as per our naming convention:
	const std::string param_name = ivec.front().getNameForParameter(param);
	const size_t param_mod = ivec.front().getParameterIndex(param_name + model_appendix);
	const bool has_param_mod = (param_mod != IOUtils::npos); //user has supplied model data alongside the meteo data

	if (has_model_expression && has_param_mod)
		if (be_verbose) std::cerr << "[W] Both model data and model function are supplied for the particle filter. Ignoring model function.\n";

	if (has_param_mod) {

		for (size_t kk = 1; kk < TT; ++kk) {
			const double val_model = ivec[kk](param_mod);
			if (val_model != IOUtils::nodata && ivec[kk](param) != IOUtils::nodata) { //TODO: disallow nodata or skip in core
				xx[kk] = val_model; //TODO: disable range check after development
			} else
				throw NoDataException("No model data at " + ivec[kk].date.toString(Date::ISO, false) + ". Please refine your model output to match the timestamps or enable meteo resampling.");
		}

	} else if (has_model_expression) {

		double te_kk, te_x_km1; //2 substitutions available: index and state value at previous time step
		te_variable te_vars[] = {{"kk", &te_kk, 0, 0}, {"x_km1", &te_x_km1, 0, 0}}; //read: "x_(k-1)"

		int te_err;
		te_expr *expr = te_compile(model_expression.c_str(), te_vars, 2, &te_err); //ready the lazy equation with variables including syntax check

		if (expr) {
			if (model_x0 == IOUtils::nodata) {
				if (be_verbose) std::cerr << "[W] No initial state value x_0 provided for particle filter; using 1st measurement.\n";
				model_x0 = ivec.front()(param);
			}
			xx[0] = model_x0;

			for (size_t kk = 1; kk < TT; ++kk) { //fill rest of vector
				te_kk = (double)kk;
				te_x_km1 = xx[kk-1];
				const double res = te_eval(expr); //evaluate expression with current substitution values
				xx[kk] = res;
			}
			te_free(expr);

	    } else {
	    	throw InvalidFormatException("Arithmetic expression '" + model_expression +
	    	    "'could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err), AT);
	    } //endif expr

	} else {
		throw InvalidArgumentException("No modeled data found in meteo data (expected '"
		    + param_name + model_appendix + "' for the particle filter). No model function found either.");
	} //endif model_expression available

	if (add_model_noise) { //user asks to add noise to model function
		for (size_t kk = 0; kk < TT; ++kk)
			xx[kk] += RNGU.draw(); //add process noise uu
		std::cout << RNGU.toString();
	} //endif noise

}

void FilterParticle::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::" + block_name);

	std::string rng_model_algorithm, rng_model_distribution, rng_model_parameters;
	bool has_prior(false);

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {

		/*** FILTER settings ***/
		if (vecArgs[ii].first == "NO_OF_PARTICLES") {
			IOUtils::parseArg(vecArgs[ii], where, NN);
		} else if (vecArgs[ii].first == "PATH_RESAMPLING") {
			IOUtils::parseArg(vecArgs[ii], where, path_resampling);

		/*** MODEL FUNCTION settings ***/
		} else if (vecArgs[ii].first == "MODEL_APPENDIX") {
			model_appendix = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "MODEL_FUNCTION") {
			model_expression = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "MODEL_X0") {
			IOUtils::parseArg(vecArgs[ii], where, model_x0);
		} else if (vecArgs[ii].first == "MODEL_ADD_NOISE") {
			IOUtils::parseArg(vecArgs[ii], where, add_model_noise);

		/*** MODEL RNG settings ***/
		} else if (vecArgs[ii].first == "MODEL_RNG_ALGORITHM") { //everything RNG will be defaulted if not provided - cf. RNG doc
			rng_model.algorithm = RandomNumberGenerator::strToRngtype(vecArgs[ii].second);
		} else if (vecArgs[ii].first == "MODEL_RNG_DISTRIBUTION") {
			rng_model.distribution = RandomNumberGenerator::strToRngdistr(vecArgs[ii].second); //convert from int to enum RNG_DISTR
		} else if (vecArgs[ii].first == "MODEL_RNG_PARAMETERS") {
			IOUtils::readLineToVec(vecArgs[ii].second, rng_model.parameters);
		} else if (vecArgs[ii].first == "MODEL_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, rng_model.seed);

		/*** PRIOR PDF RNG settings ***/
		} else if (vecArgs[ii].first == "PRIOR_RNG_ALGORITHM") { //everything RNG will be defaulted if not provided - cf. RNG doc
			rng_prior.algorithm = RandomNumberGenerator::strToRngtype(vecArgs[ii].second);
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_DISTRIBUTION") {
			rng_prior.distribution = RandomNumberGenerator::strToRngdistr(vecArgs[ii].second); //convert from int to enum RNG_DISTR
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_PARAMETERS") {
			IOUtils::readLineToVec(vecArgs[ii].second, rng_prior.parameters);
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, rng_prior.seed);
			has_prior = true;

		/*** MISC settings ***/
		} else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		}

		if (!has_prior) //not one prior pdf RNG setting -> pick importance density
			rng_prior = rng_model; //if there's a least one, then r_p deviates from the default constructor, not from r_m

	} //endfor ii
}

void FilterParticle::vecMeteoToEigen(const std::vector<MeteoData>& vec, Eigen::VectorXd& eig, const unsigned int& param)
{ //naive copy of 1 meteo parameter in STL vector to Eigen vector
	eig.resize(vec.size());
	for (int i = 0; i < vec.size(); ++i)
		eig[i] = vec[i](param);
}

void FilterParticle::readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out)
{ //uint64 type version
	vec_out.clear();
	std::istringstream iss(line_in);
	uint64_t val;
	while (!iss.eof()) {
		iss >> std::skipws >> val;
		if (iss.fail())
			throw InvalidFormatException("Unable to parse process noise seed integers for particle filter.", AT);
		vec_out.push_back(val);
	}
}

} //namespace
