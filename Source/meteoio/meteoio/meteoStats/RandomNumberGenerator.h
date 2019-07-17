/***********************************************************************************/
/*  Copyright 2018 Michael Reisecker and work cited in documentation and source    */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

#include <ctime> //for time seed
#include <inttypes.h> //for uint64_t etc. (ULL guarantees at least 64 bits and should be ok)
#include <string> //for toString()
#include <vector> //for internal states

namespace mio {

/**
 * @class RandomNumberGenerator
 * @section rng_rng Random Number Generator
 * Random integer:
 * @code
 * mio::RandomNumberGenerator RNG;
 * int nn = RNG.int64();
 * @endcode
 * Random double with Gaussian distribution:
 * @code
 * RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAUSS)
 * double rr = RNG.doub();
 * @endcode
 *
 * @section rng_purpose Purpose of this class
 * We offer two inherently 32 bit generators, and an inherently 64 bit generator
 * (although all three need 64 bits space), aswell as some convenience methods.
 *
 * The goal is to have a generator suite that satisfies all needs for statistical filters / Monte Carlo methods (and not
 * more), especially when working within MeteoIO. In a way, statistical filters are what ultimately justify this class,
 * and therefore it is meant to be tailored to their needs (and be ANSI C).
 *
 * So, if you are currently using this (cf. \ref rng_appendix_A):
 * @code
 *     srand( time(NULL) );
 *     return rand() % range;
 *     return rand() / double(RAND_MAX + 1);
 * @endcode
 * then switch to MeteoIO's RNG. If however you rely heavily on the best quality random numbers, maybe even crypto-secure,
 * there are some links to dedicated libraries in the \ref rng_bibliography.
 * Apart from the generators and distributions, this class aims to take away all the small steps that are often quickly
 * deemed good enough, i.e. \ref rng_algorithms "generator choice", \ref rng_seeding "seeding", \ref rng_seeding "saving states",
 * \ref rng_range "range calculations", ...
 * @note There is an example program exercising most of the RNG's features in the `/doc/examples` folder.
 *
 * @section rng_overview Overview
 * What it can already do:
 *  - produce quality 64 bit, 32 bit and double random numbers with one simple call
 *  - doubles with different probability distributions
 *  - some probability density functions and cumulative distribution functions
 *  - make use of quality hardware and time seeds
 *  - fast downscaling of random numbers to a range
 *  - true floating point random numbers without rounding
 *  - can be resumed from a saved state
 *  - sidesteps some widespread misuse of quick & dirty solutions
 *  - sidesteps some issues with the insidious standard library
 *  - offers a ready-to-use interface for implementing new distributions (or even generators)
 *  - passes statistical tests
 *  - very good benchmarks for the generator cores
 *
 * What's left to do:
 *  - some distributions
 *  - Monte Carlo sampling template for arbitrary distribution functions
 *
 * @section rng_integers Random integers
 * To draw random integers, you can use either the int32() or int64() function call to receive
 * 32 or 64 bit pseudo-random values, respectively:
 * @code
 *  uint32_t rn = RNG.int32();
 *  uint64_t rm = RNG.int64();
 * @endcode
 * @note The generators guarantee exactly 32 or 64 random bits, so the appropriate
 * types defined by `inttypes.h` are suited best. However, if you can live with
 * `conversion` compiler warnings then any integer type will do. Of course, you
 * can also simply use the type uintNN_t is mapped to on your machine. In short,
 * `int rn = RNG.int32()` will work for a quick try.
 *
 * @section rng_doubles Double values
 *
 * Usually, you will simply get a double value like this:
 * @code
 * double rr = RNG.doub();
 * @endcode
 *
 * You can call the doub() function with an `RNG_BOUND` argument inkluding or excluding 0 and 1 (cf. the enum below).
 * This can only be done for the uniform distribution, where it's clear what the borders are.
 * @code
 *     double rr = RNG.doub(RNG_AEXCBINC); //make sure it's not 0
 *     rr = log(rr);
 * @endcode
 *
 * Uniform random double values are quite hard to generate. The code example at Ref. [TC14] provides a method to do it,
 * which is to interpret a random stream of bits as fractional part of the binary expansion of a number in [0, 1].
 * The file also goes into details about why other methods are troublesome if we rely on quality, e. g. sensitive random searches
 * on a plane due to the gap size of `1/2^(bits)`.
 *
 * In short, the doub() function returns a double within `[0, 1]` that is rounded to
 * the nearest `1/2^64`th.
 * To get around this, you can set `true_double` to use an algorithm that calculates doubles in `[0, 1]` without the usual limitation of
 * floating point randoms being on a grid (but then you must use `RNG_AINCBINC` and guard that in your own code
 * if it must not happen even once).
 * @code
 *     double rr;
 *     do {
 *         rr = RNG.doub(RNG_AINCBINC, true); //get a random float on continuous axis
 *     } while (rr == 0.); //make sure it's not 0
 *     rr = log(rr);
 * @endcode
 *
 *
 * @section rng_distribution Distributions / Random deviates
 * For doubles, you can select from a number of distribution functions.
 *
 * For example, you can draw a Student-t variate like this:
 * @code
 * RNG.setDistribution(mio::RandomNumberGenerator::RNG_STUDENTT);
 * double rr = RNG.doub();
 * @endcode
 * If you don't set any distribution parameters, they will be defaulted (cf. \ref
 * rng_distributionparams).
 *
 * So far, the following deviates are available, defined by their probability density:
 *
 * - <a href="https://en.wikipedia.org/wiki/Uniform_distribution">Uniform</a>, `RNG_UNIFORM`:
 * \f[
 * f(x) = \frac 1{b-a} \quad \mathrm{for} \quad a \le x \le b,
 * \f]
 *
 * - <a href="https://en.wikipedia.org/wiki/Normal_distribution">Gauss</a> (= Normal), `RNG_GAUSS`:
 * \f[
 * f(x \mid\mu,\sigma^2)=\frac{1}{\sqrt{2\pi\sigma^2}}\exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
 * \quad \mathrm{for} \quad -\infty < x < +\infty
 * \f]
 * - <a href="https://en.wikipedia.org/wiki/Gamma_distribution">Gamma</a>, `RNG_GAMMA`:
 * \f[
 * f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x} \quad \mathrm{for} \quad x > 0
 * \f]
 * \f[
 * \alpha > 0, \quad \beta > 0
 * \f]
 * - <a href="https://en.wikipedia.org/wiki/Student%27s_t-distribution">Student-t</a>, `RNG_STUDENTT`:
 * \f[
 * f_\nu(x)=\frac{\Gamma\left(\frac{\nu+1}{2}\right)}
 * {\sqrt{\nu\pi}~\Gamma\left(\frac{\nu}{2}\right)}
 * \left(1+\frac{x^{2}}{\nu}\right)^{-\frac{\nu+1}{2}}
 * \quad \mathrm{for} \quad -\infty < x < +\infty
 * \f]
 * \f[
 * \nu > 0
 * \f]
 * - <a href="https://en.wikipedia.org/wiki/Chi-squared_distribution">Chi-squared</a>, `RNG_CHISQUARED`:
 * \f[
 * Z\sim\mathcal{N}(0,1) \rightarrow Z^2\sim\chi^2(1) \rightarrow
 * Y=\chi^2(r_1)+\ldots+\chi^2(r_\nu)\sim\chi^2(r_1+\ldots+r_\nu)
 * \f]
 * \f[
 * \nu > 0
 * \f]
 * - <a href="https://en.wikipedia.org/wiki/Beta_distribution">Beta</a>, `RNG_BETA`:
 * \f[
 * f(x)=\frac{1}{\mathrm{B}(\alpha, \beta)} x^{\alpha-1}(1-x)^{\beta-1}
 * \quad \mathrm{for} \quad 0 < x < 1
 * \f]
 * \f[
 * \alpha > 0, \quad \beta > 0
 * \f]
 * - <a href="https://en.wikipedia.org/wiki/F-distribution">F</a>, `RNG_F`:
 * \f[
 * f(x|\nu _1, \nu _2) = m^{\frac{\nu _1}{2}} \nu _2^{\frac{\nu _2}{2}} \cdot
 * \frac{\Gamma (\frac{\nu _1}{2}+\frac{\nu _2}{2})}{\Gamma (\frac{\nu _1}{2})
 * \Gamma(\frac{\nu _2}{2})} \cdot \frac{x^{\frac{\nu _1}{2}-1}}
 * {(\nu _1 x+\nu _2)^\frac{\nu _1+\nu _2}{2}}
 * \quad \mathrm{for} \quad x \geq 0
 * \f]
 *
 * @note Integer values are always drawn with Uniform distribution, so if you really need a
 * different one, for now you'll have to draw a double value and multiply accordingly.
 *
 * @section rng_distributionparams Distribution parameters
 *
 * When you change the distribution, you switch to a completely new one.
 * The distribution parameters <i>have to be</i> provided each time, or they will be defaulted.
 *
 * You can set distribution parameters in two ways:
 *  -# By setting (getting) them one by one after a distribution has been set:
 *  @code
 *        RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAUSS);
 *        RNG.setDistributionParameter("mean", 5.);
 *        RNG.setDistributionParameter("sigma", 2.);
 *        double mean_out = RNG.getDistributionParameter("mean");
 *  @endcode
 *  -# Via accessing the `DistributionParameters` vector directly.
 *    A `std::vector<double>` is provided with input parameters, or it has the output stored to it.
 *    This vector is given to the setDistribution() or getDistribution() call (the latter also returns the distribution type).
 *
 *    Set distribution and parameters:
 *        @code
 *        std::vector<double> distribution_params;
 *        const double alpha = 1.2, beta = 1.;
 *        distribution_params.push_back(alpha);
 *        distribution_params.push_back(beta);
 *        RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAMMA, distribution_params);
 *        @endcode
 *
 *    Get distribution and parameters:
 *        @code
 *        distribution_params.clear();
 *        const mio::RandomNumberGenerator::RNG_DISTR dist_t = RNG.getDistribution(distribution_params);
 *        const double alpha_out = distribution_params.at(0); //check doc for indices
 *        @endcode
 *
 * Here are the names within MeteoIO, arguments and default values of the distributions described above:
 *
 * <table>
 *  <tr><th>Distribution</th><th>Index</th><th>Parameter</th><th>Description</th><th>Default value</th></tr>
 *  <tr><td>RNG_UNIFORM</td><td>-</td><td>-</td><td>no parameteres</td><td>-</td></tr>
 *  <tr><td rowspan="2">RNG_GAUSS = RNG_NORMAL</td><td>1</td><td>mean</td><td>center of curve</td><td>0</td></tr>
 *  <tr><td>2</td><td>sigma</td><td>standard deviation</td><td>1</td></tr>
 *  <tr><td rowspan="2">RNG_GAMMA</td><td>1</td><td>alpha</td><td>shape parameter 1</td><td>1</td></tr>
 *  <tr><td>2</td><td>beta</td><td>shape parameter 2</td><td>1</td></tr>
 *  <tr><td>RNG_CHISQUARED</td><td>1</td><td>nu</td><td>number of degrees of freedom</td><td>1</td></tr>
 *  <tr><td rowspan="3">RNG_STUDENTT</td><td>1</td><td>nu</td><td>number of degrees of freedom</td><td>1</td></tr>
 *  <tr><td>2</td><td>mean</td><td>center of curve</td><td>0</td></tr>
 *  <tr><td>3</td><td>sigma</td><td>standard deviation</td><td>1</td></tr>
 *  <tr><td rowspan="2">RNG_BETA</td><td>1</td><td>alpha</td><td>shape parameter 1</td><td>1</td></tr>
 *  <tr><td>2</td><td>beta</td><td>shape parameter 2</td><td>1</td></tr>
 *  <tr><td rowspan="2">RNG_F</td><td>1</td><td>nu1</td><td>degrees of freedom in numerator</td><td>1</td></tr>
 *  <tr><td>2</td><td>nu2</td><td>degrees of freedom in denominator</td><td>1</td></tr>
 * </table>
 *
 * @section rng_pdf PDFs and CDFs
 * So far, the probability density function and cumulative distribution function are
 * available for the Gauss distribution like this:
 * @code
 * mio::RandomNumberGenerator RNG;
 * RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAUSS);
 * const double rg = RNG.doub();
 * std::cout << "Drew: " << rg << std::endl;
 * std::cout << std::setprecision(4) << "Probability to hit a number close to this one: "
 *           << RNG.pdf(rg)*100 << " %" << std::endl;
 * std::cout << "Probability to hit below this number: " << RNG.cdf(rg)*100 << " %" << std::endl;
 * @endcode
 * @section rng_range Range calculations
 * You can draw 32 and 64 bit integers in a given range like this:
 * @code
 * uint32_t rn = RNG.range32(10, 20);
 * uint64_t rn = RNG.range64(100, 2000);
 * @endcode
 *
 * Note that whatever you do, for an arbitrary count of random numbers you cannot downscale them and keep the distribution
 * completely intact (although "non-trivial" methods are under investigation) due to the
 * <a href="https://en.wikipedia.org/wiki/Pigeonhole_principle">Pigeonhole principle</a>.
 * The only way not to distort the (uniform) distribution is to generate lots of numbers and reject
 * out of boundary values. This is done by the trueRange32() function with a default `1e6` tries before
 * resorting to downscaling (indicated by the return boolean). You can crank this up, but to state the obvious
 * if the range gets small this gets costly quickly. The bitshift-methods above only avoid the slow modulo and
 * its possible inherent bias.
 * @code
 * uint32_t rt;
 * const bool true_range_success = RNG.trueRange32(100, 3000, rt);
 * @endcode
 *
 * @section rng_algorithms Random number generator algorithms
 * Three algorithms are available, namely the Mersenne Twister, a "classical" combined generator,
 * and a rather new algorithm with promising statistical qualities.
 * You can set them when initializing the RNG like this:
 * @code
 * mio::RandomNumberGenerator RNG(mio::RandomNumberGenerator::RNG_XOR); //default
 * mio::RandomNumberGenerator MTW(mio::RandomNumberGenerator::RNG_MTW);
 * mio::RandomNumberGenerator PCG(mio::RandomNumberGenerator::RNG_PCG);
 * @endcode
 *
 * @subsection rng_mtw Mersenne Twister
 * - Implementation of the wide-spread Mersenne Twister algorithm by M. Matsumoto and T. Nishimura (Ref. [MN98]).
 * - By using 624 internal states, the period is extremely long.
 * - This does not make it crypto-secure (the state can be derived from 624 random numbers), but it
 *   passes many statistical tests and is the standard RNG in numerous well-known software packages.
 * - It needs a few kB buffer size, which is relatively large compared to the other generators.
 * - Facts:
 *   size: 32 bit, period: 2^19937-1 (Mersenne prime) ~ 4.3e6001
 *
 * @subsection rng_xor Combined generator
 * - Generator with xor, shift and multiplication
 * - This is a fast combined generator that should be suitable for all but very special Monte Carlo applications.
 *   Since more than one internal states are being propagated and combined to the output, this makes it
 *   somewhat less predictable than similar generators.
 * - Seed with any value except vv.
 * - Facts:
 *   size: 64 bit, period: ~3.138e57
 *
 * @subsection rng_pcg PCG
 * - Permuted linear congruential generator by Prof. Melissa O'Neill
 * - Range is overestimated, and this generator performs very well in statistical tests, i. e. it is less
 *   predictable than related generators. Even smaller versions with only 32 bit entropy pass SmallCrunch,
 *   which is only barely theoretically possible.
 *   The key element is the hashing function from the internal states to the random number.
 *   The algorithm author describes this RNG family in her paper (Ref. [MO14]) and offers a huge
 *   sophisticated <a href="https://github.com/imneme/pcg-c">C-library</a> for free download with from tiny to 128 bit generators.
 * - You should seed true 64 bit values or discard the first numbers.
 * - If drawing 64 bit naturally is slow on your machine, try this one.
 * - Facts:
 *   size: 32 bit, period: ~2^64 ~ 1.8e19
 *
 * @note The numbers were subjected to the dieharder random number test suite (Ref. [RB03]), passing
 * most tests (while `rand()` fails horribly). Here is a quick benchmark (cf. \ref rng_appendix_B
 * for dieharder results):
 *
 * Generator | Bit  |  #  | time (s) | # per second
 * ----------|------|-----|----------|-----------------
 * XOR       |  64  | 1e8 |    7.4   |   13.49e6
 * PCG       |  32  | 1e8 |    4.3   |   23.29e6
 * MTW       |  32  | 1e8 |    6.2   |   15.97e6
 * XOR-doub  |      | 1e8 |   18.3   |    5.46e6
 * " Gauss   |      | 1e8 |   59.0   |    1.69e6 (310% slower)
 * hardware seed  |  64 | 1e6 |  199.8   |    5005.2
 * time seed      |  64 | 1e7 |   17.6   |   567.0e3
 * rand()     |  32 | 1e8 |    8.0   |   12.53e6
 * std::MTW   |  64 | 1e8 |   10.6   |    9.45e6
 * std::Gauss |     | 1e8 |   94.2   |    1.06e6 (900% slower)
 *
 * @section rng_seeding Seeding
 * Each time a RNG is constructed, it auto-seeds from hardware noise, or if that fails by hashing the
 * system time. If the system time is the same when you initialize more generators at once, it will still
 * seed differently. Successful hardware noise can be checked with getHardwareSeedSuccess(), and it's
 * also noted in the toString() info.
 * Manually seeding the generator is done after the fact with setState(), for example, to resume
 * experiments after the state was saved via getState().
 * Finally, we offer the getUniqueSeed() function, so if you have set up your
 * calculations with a grandfathered in, better, faster, ... RNG we can at least help with the seeding.
 *
 * Example: By default, the Mersenne Twister initializes its 624 states with a linear congruential
 * generator and then mixes that together with 64 hardware noise (resp. time hash) values. If you wanted
 * to seed all 624 internal states with hardware noise (or time hashes) you could do it like this:
 * @code
 * mio::RandomNumberGenerator MTW(mio::RandomNumberGenerator::RNG_MTW);
 * std::vector<uint64_t> seed_array;
 * for (size_t i = 0; i < 624; ++i)
 *     seed_array.push_back(MTW.getUniqueSeed());
 * MTW.setState(seed_array);
 * @endcode
 *
 * You can retrieve the generator's state to later resume experiments at exactly this point:
 * @code
 * std::vector<uint64_t> out_seed;
 * RNG.getState(out_seed);
 * mio::RandomNumberGenerator RN2;
 * RN2.setState(out_seed);
 * @endcode
 *
 * @section rng_developer Developer's guide
 * For developers of statistical filters it may be important to be able to implement custom probability distributions,
 * for example for an empirical nonlinear sensor response. This class tries to be easy to expand in that regard.
 * There are comment markers in the header and source files leading with "`CUSTOM_DIST step #`: ..." in the 7 places
 * you need to register your custom distribution functions at. These 7 steps are:
 *  -# Give your distribution a name within MeteoIO
 *  -# Put your functions' prototypes in the header
 *  -# Point to your distribution function in the generic setDistribution() function,
 *     and use the interface to the caller to set your distribution parameters
 *  -# Give a small output info string
 *  -# Write your distribution function, its pdf and cdf (if only to throw a not-implemented error)
 *  -# If you want, you can map your parameters to names in the get- and setDistributionParameter() functions.
 *  -# Map a string shorthand to the name of your distribution.
 *
 * @section rng_bibliography Bibliography
 * - [AS73] Abramowitz, Stegun.
 *        <i>Handbook of Mathematical Functions.</i>
 *        Applied Mathematics Series 55, 10th edition, 1973.
 * - [DK81] Donald E. Knuth.
 *        <i>The art of computer programming 2.</i>
 *        Addison-Wesley series in computer science and information processing, 2nd edition, 1981.
 * - [GM03] George Marsaglia.
 *        <i>Xorshift RNGs.</i>
 *        Journal of Statistical Software, Articles, 8/14, 2003.
 * - [MN98] Makoto Matsumoto and Takuji Nishimura.
 *        <i>Mersenne Twister: A 623-dimensionally equidistributed uniform pseudo-random number generator.</i>
 *        ACM Transactions on Modeling and Computer Simulation, 8/1, 1998.
 *        http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 * - [MO14] Melissa O'Neill.
 *        <i>PCG: A family of simple fast space-efficient statistically good algorithms
 *        for random number generation.</i>
 *        Harvey Mudd College, 2014.
 *        http://www.pcg-random.org
 * - [MT00] G. Marsaglia and W. Tsang.
 *        <i>A simple method for generating gamma variables.</i>
 *        ACM Transactions on Mathematical Software, 26/3, 2000.
 * - [NR3] Press, Teukolsky, Vetterling, Flannery.
 *        <i>Numerical Recipes. The Art of Scientific Computing.</i>
 *        Cambridge University Press, 3rd edition, 2007.
 * - [PE97] Pierre L'Ecuyer.
 *        <i>Distribution properties of multiply-with-carry random number generators.</i>
 *        Mathematics of Computation, 66/218i, 1997.
 * - [PE99] Pierr L'Ecuyer.
 *        <i>Tables of linear congruential generators of different sizes and good
 *        lattice structure.</i>
 *        Mathematics of Computation, 68/225, 1999.
 *        <a href="https://www.iro.umontreal.ca/~lecuyer/myftp/papers/latrules99Errata.pdf">Errata</a>
 *        for the paper, read on 18-10-23)
 * - [RB03] Robert G. Brown.
 *        <i>Dieharder: A random number test suite.</i>
 *        http://webhome.phy.duke.edu/~rgb/General/dieharder.php
 * - [TC14] Taylor R. Campbell.
 *        http://mumble.net/~campbell/tmp/random_real.c (read on 18-10-23)
 *
 * @section rng_appendix_A Appendix A
 * Why is
 * @code
 * srand( time(NULL) );
 * return rand() % range;
 * return rand() / double(RAND_MAX + 1);
 * @endcode
 * bad?
 * - A purely linear congruential RNG has purely bad statistical qualities
 * - Quiet type collision between time() and srand()
 * - `(% range)` distorts the distribution at the borders and `(range + 1)` should be used
 * - Careful not to hit RAND_MAX = INT_MAX, maybe `((double)RAND_MAX) + 1`.
 *
 * \image html rng_planes.png "Triplets of random numbers generated by a poor linear congruential generator."
 *
 * Why is
 * @code
 * std::random_device RNG;
 * std::seed_seq seed{RNG()};
 * std::mt19937 RNG_MT(seed);
 * @endcode
 * not good?
 * - A 624 state Mersenne Twister is seeded with a single 32 bit value, not a sequence
 * - Leads to statistical flaws; some numbers are never drawn
 * - `std::seed_seq` isn't a bijection like it's supposed to be (as of C++17)
 * - Can produce zero-state
 *
 * @section rng_appendix_B Appendix B
 * Random number quality summary:
 * The RNG performs as expected and passes statistical tests within reason.
 *
 * The generators were subjected to the test suite dieharder, alongside with the
 * hardware device, rand(), and a state of the art crypto-generator. Please refer
 * to the <a href="https://linux.die.net/man/1/dieharder">man page for this project</a>
 * for an interpretation of the results. In short:
 * - The p-number denotes how likely it is that a perfect generator would produce
 *   this sequence.
 * - Values below 0.05 and above 0.95 are usually considered bad. However, if a
 *   generator does not produce p-values below 0.05 in 5% of the tests, this is
 *   equally bad, since the p-value itself is a uniform test statistic.
 * - This means that in a full test run, a handful of "weak" results are expected!
 *
 * Even with default values, dieharder uses up massive amounts of random numbers,
 * and it is designed to be able to push all generators to failure and make a
 * stronger assessment than the ambiguous "weak". However, this also pushes the
 * runtime from hours to many hours.
 * So far, only a small number of numbers (about 100 MB per test) was used to check
 * the integrity of the code itself, and not so much the algorithms behind it. For
 * the latter, you will find material in the bibliography. Some "weak" results are
 * definitely due to this.
 * The hardware seed was piped to dieharder for a continuous flow of random words.
 *
 *
 * @code
 * #=============================================================================#
 * #            dieharder version 3.31.1 Copyright 2003 Robert G. Brown          #
 * #=============================================================================#
 * # ---------- C's rand() function ----------
 *    rng_name    |           filename             |rands/second|
 *      file_input|          dieharder_rand_1e7.dat|  7.38e+05  |
 * #=============================================================================#
 *         test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 * # The file file_input was rewound 32 times
 *    diehard_birthdays|   0|       100|     100|0.16952200|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.00196919|   WEAK
 *   diehard_rank_32x32|   0|     40000|     100|0.00000000|  FAILED
 *     diehard_rank_6x8|   0|    100000|     100|0.01267242|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.00000000|  FAILED
 * #=============================================================================#
 * # ---------- MeteoIO's XOR generator ----------
 *   rng_name    |           filename             |rands/second|
 *     file_input|           dieharder_xor_1e7.dat|  7.35e+05  |
 * #=============================================================================#
 *         test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 * # The file file_input was rewound 88 times
 *    diehard_birthdays|   0|       100|     100|0.92189928|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.00000567|   WEAK
 *   diehard_rank_32x32|   0|     40000|     100|0.33514633|  PASSED
 *     diehard_rank_6x8|   0|    100000|     100|0.11124602|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.12018532|  PASSED
 *         diehard_opso|   0|   2097152|     100|0.26877388|  PASSED
 *          diehard_dna|   0|   2097152|     100|0.72380963|  PASSED
 * diehard_count_1s_str|   0|    256000|     100|0.31986270|  PASSED
 * diehard_count_1s_byt|   0|    256000|     100|0.00619749|  PASSED
 *  diehard_parking_lot|   0|     12000|     100|0.45123033|  PASSED
 *     diehard_2dsphere|   2|      8000|     100|0.52288444|  PASSED
 *     diehard_3dsphere|   3|      4000|     100|0.06412948|  PASSED
 * #=============================================================================#
 * # ---------- MeteoIO's PCG generator ----------
 *    rng_name    |           filename             |rands/second|
 *      file_input|           dieharder_pcg_1e7.dat|  7.36e+05  |
 * #=============================================================================#
 *         test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 * # The file file_input was rewound 112 times
 *    diehard_birthdays|   0|       100|     100|0.92445127|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.00479393|  PASSED
 *   diehard_rank_32x32|   0|     40000|     100|0.18435280|  PASSED
 *     diehard_rank_6x8|   0|    100000|     100|0.20921510|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.45293843|  PASSED
 *         diehard_opso|   0|   2097152|     100|0.08307580|  PASSED
 *          diehard_dna|   0|   2097152|     100|0.58915277|  PASSED
 * diehard_count_1s_str|   0|    256000|     100|0.57353220|  PASSED
 * diehard_count_1s_byt|   0|    256000|     100|0.01000397|  PASSED
 *  diehard_parking_lot|   0|     12000|     100|0.54459758|  PASSED
 *     diehard_2dsphere|   2|      8000|     100|0.80960785|  PASSED
 *     diehard_3dsphere|   3|      4000|     100|0.79071682|  PASSED
 * #=============================================================================#
 * # ---------- MeteoIO's Mersenne Twister ----------
 *   rng_name    |           filename             |rands/second|
 *      file_input|            dieharder_mt_1e7.dat|  7.26e+05  |
 * #=============================================================================#
 *         test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 * # The file file_input was rewound 3 times
 *    diehard_birthdays|   0|       100|     100|0.96357283|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.01695490|  PASSED
 *   diehard_rank_32x32|   0|     40000|     100|0.88275675|  PASSED
 *     diehard_rank_6x8|   0|    100000|     100|0.55027055|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.57085617|  PASSED
 * #=============================================================================#
 * # ---------- MeteoIO's hardware seed through /dev/urandom ----------
 *    rng_name    |rands/second|   Seed   |
 * stdin_input_raw|  3.60e+06  |2936014326|
 * #=============================================================================#
 *        test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 *    diehard_birthdays|   0|       100|     100|0.32694395|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.42749194|  PASSED
 *   diehard_rank_32x32|   0|     40000|     100|0.69443998|  PASSED
 *     diehard_rank_6x8|   0|    100000|     100|0.97392304|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.20992531|  PASSED
 *         diehard_opso|   0|   2097152|     100|0.72376741|  PASSED
 *          diehard_dna|   0|   2097152|     100|0.00581678|  PASSED
 * diehard_count_1s_str|   0|    256000|     100|0.44478974|  PASSED
 * diehard_count_1s_byt|   0|    256000|     100|0.99778492|   WEAK
 *  diehard_parking_lot|   0|     12000|     100|0.27195005|  PASSED
 *     diehard_2dsphere|   2|      8000|     100|0.61223463|  PASSED
 *     diehard_3dsphere|   3|      4000|     100|0.42724315|  PASSED
 *      diehard_squeeze|   0|    100000|     100|0.12032198|  PASSED
 *         diehard_sums|   0|       100|     100|0.01395721|  PASSED
 *         diehard_runs|   0|    100000|     100|0.58328875|  PASSED
 *        diehard_craps|   0|    200000|     100|0.53868076|  PASSED
 *  marsaglia_tsang_gcd|   0|  10000000|     100|0.46751336|  PASSED
 *          sts_monobit|   1|    100000|     100|0.43425440|  PASSED
 *             sts_runs|   2|    100000|     100|0.93593627|  PASSED
 *           sts_serial|   1|    100000|     100|0.06718051|  PASSED
 * #=============================================================================#
 * # ---------- Modern AES crypto generator ----------
 *    rng_name    |rands/second|   Seed   |
 *         AES_OFB|  6.37e+06  |3232005462|
 * #=============================================================================#
 *         test_name   |ntup| tsamples |psamples|  p-value |Assessment
 * #=============================================================================#
 *    diehard_birthdays|   0|       100|     100|0.37189795|  PASSED
 *       diehard_operm5|   0|   1000000|     100|0.73717774|  PASSED
 *   diehard_rank_32x32|   0|     40000|     100|0.82597218|  PASSED
 *     diehard_rank_6x8|   0|    100000|     100|0.46301177|  PASSED
 *    diehard_bitstream|   0|   2097152|     100|0.56940920|  PASSED
 *         diehard_opso|   0|   2097152|     100|0.49114028|  PASSED
 *         diehard_oqso|   0|   2097152|     100|0.48034789|  PASSED
 *          diehard_dna|   0|   2097152|     100|0.94132992|  PASSED
 * diehard_count_1s_str|   0|    256000|     100|0.26244149|  PASSED
 * diehard_count_1s_byt|   0|    256000|     100|0.63672623|  PASSED
 *  diehard_parking_lot|   0|     12000|     100|0.25282458|  PASSED
 *     diehard_2dsphere|   2|      8000|     100|0.28551551|  PASSED
 *     diehard_3dsphere|   3|      4000|     100|0.62322167|  PASSED
 *      diehard_squeeze|   0|    100000|     100|0.51563931|  PASSED
 *         diehard_sums|   0|       100|     100|0.15600758|  PASSED
 *         diehard_runs|   0|    100000|     100|0.16471462|  PASSED
 *         diehard_runs|   0|    100000|     100|0.51138070|  PASSED
 *        diehard_craps|   0|    200000|     100|0.79803071|  PASSED
 *        diehard_craps|   0|    200000|     100|0.87598421|  PASSED
 * @endcode
 * (This table is up for revision on a better machine once the numbers are actually
 * being used within MeteoIO.)
 *
 * @ingroup stats
 * @author Michael Reisecker
 * @date 2018-10
*/

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
		unsigned int countLeadingZeros(const uint64_t& nn) const;
};

class RandomNumberGenerator : private RngCore {
	public:
		enum RNG_TYPE //computation method
		{
			RNG_XOR, //!< Combined generator
			RNG_PCG, //!< Permuted linear congruential generator
			RNG_MTW //!< Mersenne Twister generator
		};
		
//CUSTOM_DIST step 1/7: Give your distribution a name in this enum
		enum RNG_DISTR //desired distribution, only used for doubles!
		{
			RNG_UNIFORM, //!< Uniform deviates
			RNG_GAUSS, //!< Gaussian deviates
			RNG_NORMAL, //!< = RNG_GAUSS
			RNG_GAMMA, //!< Gamma deviates
			RNG_CHISQUARED, //!< Chi-Squared deviates
			RNG_STUDENTT, //!< Student-t deviates
			RNG_BETA, //!< Beta deviates
			RNG_F //!< Fisher deviates
		};
		enum RNG_BOUND //return uniform double respecting these boundaries
		{
			RNG_AINCBINC, //!< [0, 1]
			RNG_AINCBEXC, //!< [0, 1)
			RNG_AEXCBINC, //!< (0, 1]
			RNG_AEXCBEXC  //!< (0, 1)
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
		double cdf(const double& xx); //cumulative distribution function

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

		static RNG_TYPE strToRngtype(const std::string& str); //get an RNG_TYPE from a string shorthand
		static RNG_DISTR strToRngdistr(const std::string& str); //get an RNG_DISTR from a string shorthand

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
		
//CUSTOM_DIST step 2/7: Add your distribution function, its pdf and cdf here, matching exactly this type:
		double doubUniform();
		double pdfUniform(const double& xx) const;
		double cdfUniform(const double& xx) const;
		double doubGauss(); //=normal
		double pdfGauss(const double& xx) const;
		double cdfGauss(const double& xx) const;
		double doubGamma();
		double doubChiSquare();
		double doubStudentT();
		double doubBeta();
		double doubF();

		double pdfNotImplemented(const double& xx) const;
		double cdfNotImplemented(const double& xx) const;

		double doubGaussKernel(const double& mean, const double& sigma); //internal calls with specific params
		double doubGammaKernel(const double& alpha, const double& beta);
		double doubBetaKernel(const double& alpha, const double& beta);
};

class RngXor : public RngCore { //combined generator with xor, shift and multiply
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

class RngPcg : public RngCore { //Permuted linear congruential generator
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

class RngMtw : public RngCore { //Mersenne Twister
	public:
		RngMtw();
		uint64_t int64();
		uint32_t int32( );
		void getState(std::vector<uint64_t>& ovec_seed) const;
		void setState(const std::vector<uint64_t>& ivec_seed);

	private:
		const unsigned int MT_NN; //number of states
		const unsigned int MT_MM; //middle word / offset
		unsigned int current_mt_index;
		
		std::vector<uint32_t> vec_states;

		bool initAllStates();
};

class RngFactory { //factory for the generator algorithm
	public: //(so that memory dedicated to the states lives only as long as the RNG)
		static RngCore* getCore(const RandomNumberGenerator::RNG_TYPE& algorithm);
};

} //namespace

#endif

