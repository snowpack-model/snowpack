/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file SeaIce.cc
 * @version 23.09
 * @bug     -
 * @brief This module contains the sea ice specific routines for the 1d snowpack model
 */

#include <snowpack/Constants.h>
#include <snowpack/Utils.h>
#include <snowpack/snowpackCore/Metamorphism.h>
#include <snowpack/snowpackCore/SeaIce.h>
#include <snowpack/snowpackCore/Snowpack.h>
#include <snowpack/snowpackCore/ReSolver1d.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

using namespace mio;
using namespace std;

/**
 * @page sea_ice Sea Ice
 *
 */

/************************************************************
 * static section                                           *
 ************************************************************/

//Threshold that defines ice
const double SeaIce::SeaWaterFreezingTemp = IOUtils::C_TO_K(-1.95);
const double SeaIce::SeaIceDensity = ReSolver1d::max_theta_ice * Constants::density_ice;
const double SeaIce::ice_threshold = 800.;
const double SeaIce::mu = 0.054;		// Freezing point coefficient
const double SeaIce::betaS = 0.824;		// Density coefficient (see: Appendix A Sea Water Density According to UNESCO Formula, https://link.springer.com/content/pdf/bbm%3A978-3-319-18908-6%2F1.pdf)
const double SeaIce::ThicknessFirstIceLayer = 0.01;
const double SeaIce::InitRg = 5.;
const double SeaIce::InitRb = 2.5;
const double SeaIce::OceanSalinity = 35.;
const double SeaIce::InitSeaIceSalinity = 5.;
const double SeaIce::InitSnowSalinity = 0.;


/************************************************************
 * non-static section                                       *
 ************************************************************/

SeaIce::SeaIce():
	SeaLevel(0.), ForcedSeaLevel(IOUtils::nodata), FreeBoard (0.), IceSurface(0.), IceSurfaceNode(0), OceanHeatFlux(0.), BottomSalFlux(0.), TopSalFlux(0.), check_initial_conditions(false), salinityprofile(SINUSSAL), thermalmodel(ASSUR1958) {}

SeaIce& SeaIce::operator=(const SeaIce& source) {
	if(this != &source) {
		SeaLevel = source.SeaLevel;
		ForcedSeaLevel = source.ForcedSeaLevel;
		FreeBoard = source.FreeBoard;
		IceSurface = source.IceSurface;
		IceSurfaceNode = source.IceSurfaceNode;
		OceanHeatFlux = source.OceanHeatFlux;
	}
	return *this;
}

SeaIce::~SeaIce() {}

void SeaIce::ConfigSeaIce(const SnowpackConfig& i_cfg) {
	// Read salinity profile
	std::string tmp_salinityprofile;
	i_cfg.getValue("SALINITYPROFILE", "SnowpackSeaice", tmp_salinityprofile, mio::IOUtils::nothrow);
	if (tmp_salinityprofile=="NONE") {
		salinityprofile=NONE;
	} else if (tmp_salinityprofile=="CONSTANT") {
		salinityprofile=CONSTANT;
	} else if (tmp_salinityprofile=="COXANDWEEKS") {
		salinityprofile=COXANDWEEKS;
	} else if (tmp_salinityprofile=="LINEARSAL") {
		salinityprofile=LINEARSAL;
	} else if (tmp_salinityprofile=="LINEARSAL2") {
		salinityprofile=LINEARSAL2;
	} else if (tmp_salinityprofile=="SINUSSAL") {
		salinityprofile=SINUSSAL;
	} else {
		prn_msg( __FILE__, __LINE__, "err", Date(), "Unknown salinity profile (key: SALINITYPROFILE).");
		throw;
	}

	// Read thermal model for sea ice to use
	std::string tmp_thermalmodel;
	i_cfg.getValue("THERMALMODEL", "SnowpackSeaice", tmp_thermalmodel, mio::IOUtils::nothrow);
	if(!tmp_thermalmodel.empty()) {
		if (tmp_thermalmodel=="IGNORE") {
			thermalmodel=IGNORE;
		} else if (tmp_thermalmodel=="ASSUR1958") {
			thermalmodel=ASSUR1958;
		} else if (tmp_thermalmodel=="VANCOPPENOLLE2019") {
			thermalmodel=VANCOPPENOLLE2019;
		} else if (tmp_thermalmodel=="VANCOPPENOLLE2019_M") {
			thermalmodel=VANCOPPENOLLE2019_M;
		} else {
			prn_msg( __FILE__, __LINE__, "err", Date(), "Unknown thermal model (key: THERMALMODEL).");
			throw;
		}
	}

	// Read whether or not to check the initial conditions
	i_cfg.getValue("CHECK_INITIAL_CONDITIONS", "SnowpackSeaice", check_initial_conditions, mio::IOUtils::nothrow);
	return;
}

std::iostream& operator<<(std::iostream& os, const SeaIce& data)
{
	os.write(reinterpret_cast<const char*>(&data.FreeBoard), sizeof(data.FreeBoard));
	os.write(reinterpret_cast<const char*>(&data.IceSurface), sizeof(data.IceSurface));
	os.write(reinterpret_cast<const char*>(&data.IceSurfaceNode), sizeof(data.IceSurfaceNode));
	os.write(reinterpret_cast<const char*>(&data.OceanHeatFlux), sizeof(data.OceanHeatFlux));
	return os;
}

std::iostream& operator>>(std::iostream& is, SeaIce& data)
{
	is.read(reinterpret_cast<char*>(&data.FreeBoard), sizeof(data.FreeBoard));
	is.read(reinterpret_cast<char*>(&data.IceSurface), sizeof(data.IceSurface));
	is.read(reinterpret_cast<char*>(&data.IceSurfaceNode), sizeof(data.IceSurfaceNode));
	is.read(reinterpret_cast<char*>(&data.OceanHeatFlux), sizeof(data.OceanHeatFlux));
	return is;
}

/**
 * @brief Determines the salinity and associated melting temperature
 * @param Xdata SnowStation to apply the salinity profile to
 */
void SeaIce::compSalinityProfile(SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();
	findIceSurface(Xdata);

	switch ( salinityprofile ) {

	case NONE:
		{
			break;
		}

	case CONSTANT:
		{
			for (size_t e = Xdata.SoilNode; e < nE; e++) {
				Xdata.Edata[e].salinity = 35.;			// Default: 35 g/kg
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			/*size_t e = Xdata.SoilNode;
			for (; e < IceSurfaceNode; e++) {
				Xdata.Edata[e].salinity = 35.;			// Default: 35 g/kg
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			for (; e < nE; e++) {
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}*/
			break;
		}

	case COXANDWEEKS:
		{
			for (size_t e = Xdata.SoilNode; e < nE; e++) {
				if(Xdata.Ndata[e].z >= findIceSurface(Xdata)) {
					// For snow
					Xdata.Edata[e].salinity = 1.;
				} else {
					// For ice
					if(Xdata.Ndata[e].z < 0.4) {
						Xdata.Edata[e].salinity = 14.24 - 19.39 * Xdata.Ndata[e].z;
					} else {
						Xdata.Edata[e].salinity = 7.88 - 1.59 * Xdata.Ndata[e].z;
					}
				}
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			break;
		}

	case LINEARSAL:
		{
			const double topSal = 1.;
			const double botSal = 5.;
			for (size_t e = Xdata.SoilNode; e < nE; e++) {
				Xdata.Edata[e].salinity = ((topSal - botSal) / (Xdata.Ndata[IceSurfaceNode].z - Xdata.Ndata[0].z)) * 0.5 * (Xdata.Ndata[e].z + Xdata.Ndata[e+1].z);		// linear gradient between 1 psu (top) to 4 psu (bottom)
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			break;
		}

	case LINEARSAL2:
		{
			const double topSal = 1.;
			const double botSal = 5.;
			// define salinity in ice
			size_t e = Xdata.SoilNode;
			for (; e <  IceSurfaceNode ; e++) {
				Xdata.Edata[e].salinity = ((topSal - botSal) / (Xdata.Ndata[IceSurfaceNode].z - Xdata.Ndata[0].z)) * 0.5 * (Xdata.Ndata[e].z + Xdata.Ndata[e+1].z);		// linear gradient between 1 psu (top) to 4 psu (bottom)
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			// define salinity in snow
			for (; e <  nE ; e++) {
				Xdata.Edata[e].salinity = 1;
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			break;
		}

	// C shaped salinity profile
	case SINUSSAL:
		{
			const double topSal = 12.;
			const double ampSal = 8.;
			const double PI = 3.141592653589793;
			// define salinity in ice
			size_t e = Xdata.SoilNode;
			for (; e <  IceSurfaceNode ; e++) {
				Xdata.Edata[e].salinity = ampSal* sin((Xdata.Ndata[e].z / (Xdata.Ndata[IceSurfaceNode].z - Xdata.Ndata[0].z))*PI+PI)+topSal;		// c shaped salinity profile in sea ice
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			// define salinity in snow
			for (; e <  nE ; e++) {
				Xdata.Edata[e].salinity = 1; // 8 after Massom et al. 1997
				Xdata.Edata[e].updDensity();
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			break;
		}
		default:
			InvalidArgumentException("Unknown salinity profile provided", AT);

	}
}

/**
 * @brief Updates the freeboard variable (i.e., sea level with respect to ice surface)\n
 *        positive: sea level below ice surface\n
 *        negative: sea level above ice surface (flooding)\n
 * @version 16.08
 * @param Xdata SnowStation object to use in calculation
 */
void SeaIce::updateFreeboard(SnowStation& Xdata)
{
	Xdata.compSnowpackMasses();
	SeaLevel = (ForcedSeaLevel!=IOUtils::nodata) ? (ForcedSeaLevel) : (Xdata.swe / (Constants::density_water + SeaIce::betaS * SeaIce::OceanSalinity));
	const double FreeBoard_snow = Xdata.cH - SeaLevel;	// This is the freeboard relative to snow surface
	FreeBoard = (findIceSurface(Xdata) - (Xdata.cH - FreeBoard_snow));
	return;
}

/**
 * @brief Find snow/ice transition for sea ice simulations\n
 * @version 16.08
 * @param Xdata SnowStation object to use in calculation
 */
double SeaIce::findIceSurface(SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();

	// Now find ice/snow transition
	if(nE == 0) {
		IceSurface = 0.;
		IceSurfaceNode = 0;
		return IceSurface;
	}
	// Deal with the case that the top element is ice
	if (Xdata.Edata[nE-1].theta[ICE] * Constants::density_ice > ice_threshold) {
		IceSurface = Xdata.Ndata[nE].z;
		IceSurfaceNode = nE;
		return IceSurface;
	}
	// Go from top to bottom. Note that ice layers inside the snowpack may fool this simple search.
	for (size_t e = nE-1; e-- > 0;) {
		if (Xdata.Edata[e].theta[ICE] * Constants::density_ice > ice_threshold && Xdata.Edata[e+1].theta[ICE] * Constants::density_ice < ice_threshold) {
			IceSurface = Xdata.Ndata[e+1].z;
			IceSurfaceNode = e+1;
			return IceSurface;
		}
	}
	IceSurfaceNode = 0;
	IceSurface = 0.;
	return IceSurface;
}

/**
 * @brief Apply flooding\n
 * @version 16.08
 * @param Xdata SnowStation object to use in calculation
 * @param Sdata
 */
void SeaIce::compFlooding(SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	size_t iN = 0;
	while (iN < Xdata.getNumberOfElements() && Xdata.Ndata[iN].z + 0.5 * Xdata.Edata[iN].L < SeaLevel) {
		const double dth_w = std::max(0., Xdata.Edata[iN].theta[AIR] * (Constants::density_ice / Constants::density_water) - Xdata.Edata[iN].theta[WATER] * (Constants::density_water / Constants::density_ice - 1.));
		Xdata.Edata[iN].theta[WATER] += dth_w;
		Xdata.Edata[iN].theta[AIR] -= dth_w;
		Xdata.Edata[iN].salinity += SeaIce::OceanSalinity * dth_w;
		Xdata.Edata[iN].salinity = std::min(SeaIce::OceanSalinity, Xdata.Edata[iN].salinity);
		Sdata.mass[SurfaceFluxes::MS_FLOODING]-=Xdata.Edata[iN].Rho * Xdata.Edata[iN].L;
		Xdata.Edata[iN].updDensity();
		Sdata.mass[SurfaceFluxes::MS_FLOODING]+=Xdata.Edata[iN].Rho * Xdata.Edata[iN].L;
		Xdata.Edata[iN].M = Xdata.Edata[iN].Rho * Xdata.Edata[iN].L;
		calculateMeltingTemperature(Xdata.Edata[iN]);
		iN++;
	}
	return;
}


/**
 * @brief Calculate melting temperature as function of brine salinity
 * @version 16.08
 * @param Edata Layer element to use in calculation
 */
void SeaIce::calculateMeltingTemperature(ElementData& Edata)
{
	// See: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	//      who is citing: Assur, A., Composition of sea ice and its tensile strength, in Arctic Sea Ice, N.  A.  S. N.  R.  C. Publ., 598, 106-138, 1958.
	Edata.meltfreeze_tk = (Edata.theta[WATER] + Edata.theta[WATER_PREF] > Constants::eps2 && Edata.salinity > Constants::eps2) ? (SeaIce::calculateMeltingTemperature(Edata.salinity / (Edata.theta[WATER] + Edata.theta[WATER_PREF]))) : (Constants::meltfreeze_tk);
	return;
}


/**
 * @brief Calculate brine salinity as a function of temperature
 * @version 17.12: initial version
 * @param T: Temperature (K)
 */
double SeaIce::calculateBrineSalinity(const double& T)
{
	if (thermalmodel == IGNORE) {
		return IOUtils::nodata;
	} else if (thermalmodel == ASSUR1958) {
		// See: Assur, A., Composition of sea ice and its tensile strength, in Arctic Sea Ice, N.  A.  S. N.  R.  C. Publ., 598, 106-138, 1958.
		const double tc = IOUtils::K_TO_C(T);
		return tc/-SeaIce::mu;
	} else if (thermalmodel == VANCOPPENOLLE2019) {
		// See Eq. 10 in: Vancoppenolle, M., Madec, G., Thomas, M., & McDougall, T. J. (2019). Thermodynamics of sea ice phase composition revisited. Journal of Geophysical Research: Oceans, 124, 615–634. doi: 10.1029/2018JC014611 
		const double a1 = -0.00535;
		const double a2 = -0.519;
		const double a3 = -18.7;
		const double tc = IOUtils::K_TO_C(T);
		return std::min(300., a1*tc*tc*tc + a2*tc*tc + a3*tc);
	} else if (thermalmodel == VANCOPPENOLLE2019_M) {
		// A quadratic fit to Eq. 10 in Vancoppenolle et al. (2019)
		const double a1 = -0.16055612425953938;
		const double a2 = -13.296596377964793;
		const double tc = IOUtils::K_TO_C(T);
		return std::min(270., a1 * tc * tc + a2 * tc);
	} else {
		throw; return IOUtils::nodata;
	}
}



/**
 * @brief Calculate melting temperature as function of brine salinity
 * @version 17.12: initial version
 * @param Sal: Brine salinity (PSU, which is g/kg)
 */
double SeaIce::calculateMeltingTemperature(const double& Sal)
{
	if (thermalmodel == IGNORE) {
		return IOUtils::C_TO_K(0.);
	} else if (thermalmodel == ASSUR1958) {
		return IOUtils::C_TO_K(-SeaIce::mu * Sal);
	} else if (thermalmodel == VANCOPPENOLLE2019) {
		const double a1 = -0.00535;
		const double a2 = -0.519;
		const double a3 = -18.7;

		// Normalize the cubic equation:
		double p = a2 / a1;
		double q = a3 / a1;
		double r = -((Sal > 300.)?(300.):(Sal)) / a1;

		// Calculate the depressed cubic:
		double a = q - (p * p / 3.0);
		double b = 2.0 * p * p * p / 27.0 - (p * q / 3.) + r;

		// Calculate discriminant
		double disc = (b / 2.) * (b / 2.) + (a / 3.) * (a / 3.) * (a / 3.);

		double rt1 = 0.;
		if (disc > 0.) {
			// One real root. Use Cardano's formula:
			double C = cbrt(-b / 2. + sqrt(disc));
			double D = cbrt(-b / 2. - sqrt(disc));

			rt1 = C + D - (p / 3.);
			return IOUtils::C_TO_K(rt1);
		} else {
			// We should not end up here...
			throw;
		}
	} else if (thermalmodel == VANCOPPENOLLE2019_M) {
		const double a1 = -0.16055612425953938;
		const double a2 = -13.296596377964793;
		const double t = -((sqrt(4.*a1*std::min(270.,Sal)+a2*a2)+a2)/a1)/2.;
		return IOUtils::C_TO_K(t);
	} else {
		throw; return IOUtils::nodata;
	}
}


/**
 * @brief Returns the tangent line to the salinity-melting point curve
 * General tangent line equation:
 *     T = mu0 + mu1 * Sal,
 * with mu1 = dTm/dS
 */
std::pair<double, double> SeaIce::getMu(const double& Sal)
{
	double mu0 = IOUtils::nodata;
	double mu1 = IOUtils::nodata;
	if (thermalmodel == IGNORE) {
		mu1 = 0.;
	} else if (thermalmodel == ASSUR1958) {
		mu1 = -SeaIce::mu;
	} else if (thermalmodel == VANCOPPENOLLE2019) {
		// Numerical differentiation:
		if(Sal < Constants::eps) {
			// Left window
			mu1 = (this->calculateMeltingTemperature(Sal - Constants::eps) - this->calculateMeltingTemperature(Sal - 2.*Constants::eps)) / Constants::eps;
		} else {
			// Centered window
			mu1 = (this->calculateMeltingTemperature(Sal + Constants::eps) - this->calculateMeltingTemperature(Sal - Constants::eps)) / (2. * Constants::eps);
		}
	} else if (thermalmodel == VANCOPPENOLLE2019_M) {
		const double a1 = -0.16055612425953938;
		const double a2 = -13.296596377964793;
		mu1 = (-1./sqrt(4.*a1*std::min(270.,Sal)+a2*a2));
	}
	mu0 = this->calculateMeltingTemperature(Sal) - mu1 * Sal;
	return std::make_pair(mu0, mu1);
}


/**
 * @brief Heat capacity of sea ice, for the combined system ICE + WATER (brine).
 * @version 16.08: initial version
 * @param T: Temperature (K)
 * @param Sal: Salinity (PSU, which is g/kg)
 * @return Heat capacity for sea ice (J / kg / K)
 */
double SeaIce::compSeaIceHeatCapacity(const double& T, const double& Sal)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 1 and 2
	const double L0 = Constants::lh_fusion;
	const double c0 = Constants::specific_heat_ice;
	return c0 + (SeaIce::mu * L0 * Sal) / (T * T);
}


/**
 * @brief Heat conduction in sea ice, for the combined system ICE + WATER (brine)
 * @version 16.08: initial version
 * @param Edata
 * @return Thermal conductivity for sea ice (W K-1 m-1)
 */
double SeaIce::compSeaIceThermalConductivity(const ElementData& Edata)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 9
	const double beta =  0.1172;		// W/m^2/permille
	const double k0 = 2.034;		// W/m/K, note that this is the thermal conductivity of fresh ice, and it may be coupled to the value in Constants.h
	// Note the conversion from kg/kg to permille for salinity
	return (k0 + ((beta * Edata.salinity) / Edata.Te));
}


/**
 * @brief Latent heat of melting for sea ice, for the combined system ICE + WATER (brine)
 * @version 16.08: initial version
 * @param T: Temperatur (K)
 * @param Sal: Salinity (PSU, which is g/kg)
 * @return Latent heat of fusion for sea ice (J / kg)
 */
double SeaIce::compSeaIceLatentHeatFusion(const double& T, const double& Sal)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 5
	const double L0 = Constants::lh_fusion;
	return L0 * (1. + (SeaIce::mu * Sal) / T);
}


/**
 * @brief Latent heat of melting for sea ice, for the combined system ICE + WATER (brine)
 * @version 16.08: initial version
 * @param Edata
 * @return Latent heat of fusion for sea ice (J / kg)
 */
double SeaIce::compSeaIceLatentHeatFusion(const ElementData& Edata)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 5
	const double L0 = Constants::lh_fusion;
	const double c0 = Constants::specific_heat_ice;
	return c0 * (Edata.meltfreeze_tk - Edata.Te) + L0 * (1. + (SeaIce::mu * Edata.salinity) / Edata.Te);
}


/**
 * @brief Calculate ice formation and decay at the bottom
 * @version 16.08: initial version
 * @param Xdata SnowStation object to use in calculation
 * @param Mdata Meteo data
 * @param sn_dt Time step (s)
 * @param Sdata
 */
void SeaIce::bottomIceFormation(SnowStation& Xdata, const CurrentMeteo& Mdata, const double& sn_dt, SurfaceFluxes& Sdata)
{
  	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	size_t nE = Xdata.getNumberOfElements();
	if (Xdata.getNumberOfElements() == 0 && Xdata.Ndata[Xdata.SoilNode].T >= SeaIce::SeaWaterFreezingTemp) {
		// Just open ocean
		return;
	}

	double netBottomEnergy = 0.;
	double dM = 0.;
	if (nE > 0 ) {
		// With at least one element, calculate net energy flux
		// Here: netBottomEnergy has units W/m^2
		netBottomEnergy = OceanHeatFlux
			+ compSeaIceThermalConductivity(EMS[Xdata.SoilNode]) * ( (NDS[Xdata.SoilNode+1].T - NDS[Xdata.SoilNode].T) / (NDS[Xdata.SoilNode+1].z - NDS[Xdata.SoilNode].z));
		dM = (-netBottomEnergy * sn_dt) / compSeaIceLatentHeatFusion(EMS[Xdata.SoilNode]);
	} else {
		// First time freezing, create first ice layer
		// Here: netBottomEnergy has units W / kg
		//netBottomEnergy = (Xdata.Ndata[Xdata.SoilNode].T - SeaIce::SeaWaterFreezingTemp) * compSeaIceHeatCapacity(Xdata.Ndata[Xdata.SoilNode].T, SeaIce::OceanSalinity) / sn_dt;
		// Convert netBottomEnergy to W / m^2, assuming a 1 cm freezing layer
		// TODO: insert density ocean water
		// TODO: we don't know the ocean temperature profile, so we cannot accurately know how thick the first ice layer is
		//const double freezing_ocean_depth = 0.1;
		//netBottomEnergy *= Constants::density_water * freezing_ocean_depth;
		//dM = (-netBottomEnergy * sn_dt) / compSeaIceLatentHeatFusion(Xdata.Ndata[Xdata.SoilNode].T, SeaIce::OceanSalinity);
		dM = ThicknessFirstIceLayer * SeaIceDensity;
	}
	ApplyBottomIceMassBalance(Xdata, Mdata, dM, Sdata);
}


/**
 * @brief Apply mass gain/loss at the bottom (dM)
 * @version 16.08: initial version
 * @param Xdata
 * @param Mdata
 * @param dM: mass change from phase changes at the bottom of the sea ice (kg/m^2), positive=gain, negative=loss. Note, dM only concerns change in ice mass.
 * @param Sdata
 */
void SeaIce::ApplyBottomIceMassBalance(SnowStation& Xdata, const CurrentMeteo& Mdata, double dM, SurfaceFluxes& Sdata)
{
	// Dereference pointers
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	size_t nE = Xdata.getNumberOfElements();

	// Apply mass change:
	double dz = 0.;
	if ( dM > 0 ) {
		// dM > 0: mass gain
		if ( nE == 0 || EMS[Xdata.SoilNode].Rho < ice_threshold ) {

			const double dH = dM / SeaIceDensity;								// Total height to be added. Note that dM represents the amount of ice resulting from the phase changes from the energy balance at the bottom of the ice.
			const size_t nAddE = 1;										// Number of elements
			const double dL = (dH / double(nAddE));								// Height of each individual layer
			for ( size_t j = 0; j < nAddE; j++ ) {
				dz += dL;
				nE++;
				Xdata.resize(nE);
				if(nE > 1) {
					// Shift all existing elements up in the domain
					for(size_t ee = nE-1; ee > Xdata.SoilNode; ee--) {
						EMS[ee]=EMS[ee-1];
						NDS[ee+1]=NDS[ee];
						NDS[ee]=NDS[ee-1];
					}
				} else {
					// Set upper node for very first element in the domain that will be newly created
					NDS[nE].T = SeaIce::calculateMeltingTemperature(OceanSalinity);
				}
				// Set the new ice element
				EMS[Xdata.SoilNode].depositionDate = Mdata.date;
				EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = dL;
				EMS[Xdata.SoilNode].theta[SOIL] = 0.;
				EMS[Xdata.SoilNode].theta[ICE] = (SeaIceDensity/Constants::density_ice);
				EMS[Xdata.SoilNode].theta[WATER] = (1. - EMS[Xdata.SoilNode].theta[ICE]) * (Constants::density_ice/Constants::density_water);
				EMS[Xdata.SoilNode].theta[WATER_PREF] = 0.;
				EMS[Xdata.SoilNode].theta[AIR] = 1.0 - EMS[Xdata.SoilNode].theta[WATER] - EMS[Xdata.SoilNode].theta[WATER_PREF] - EMS[Xdata.SoilNode].theta[ICE] - EMS[Xdata.SoilNode].theta[SOIL];
				EMS[Xdata.SoilNode].salinity = OceanSalinity * EMS[Xdata.SoilNode].theta[WATER];
				EMS[Xdata.SoilNode].updDensity();
				Sdata.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING]+=dM;


				for (unsigned short ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[Xdata.SoilNode].conc[ICE][ii]   = Mdata.conc[ii]*Constants::density_ice/Constants::density_water;
					EMS[Xdata.SoilNode].conc[WATER][ii] = Mdata.conc[ii];
					EMS[Xdata.SoilNode].conc[AIR][ii]   = 0.;
					EMS[Xdata.SoilNode].conc[SOIL][ii]  = 0.;
				}

				// Constitutive Parameters
				EMS[Xdata.SoilNode].k[TEMPERATURE] = EMS[Xdata.SoilNode].k[SEEPAGE] = EMS[Xdata.SoilNode].k[SETTLEMENT]= 0.;
				EMS[Xdata.SoilNode].heatCapacity();
				EMS[Xdata.SoilNode].c[SEEPAGE] = EMS[Xdata.SoilNode].c[SETTLEMENT]= 0.;
				EMS[Xdata.SoilNode].soil[SOIL_RHO] = EMS[Xdata.SoilNode].soil[SOIL_K] = EMS[Xdata.SoilNode].soil[SOIL_C] = 0.;
				EMS[Xdata.SoilNode].snowResidualWaterContent();

				//new snow micro-structure
				EMS[Xdata.SoilNode].sw_abs = 0.;
				EMS[Xdata.SoilNode].rg = InitRg;
				EMS[Xdata.SoilNode].dd = 0.;
				EMS[Xdata.SoilNode].sp = 1.;
				EMS[Xdata.SoilNode].rb = InitRb;
				EMS[Xdata.SoilNode].N3 = Metamorphism::getCoordinationNumberN3(EMS[Xdata.SoilNode].Rho);
				EMS[Xdata.SoilNode].opticalEquivalentGrainSize();
				EMS[Xdata.SoilNode].mk = 7;
				EMS[Xdata.SoilNode].metamo = 0.;
				EMS[Xdata.SoilNode].snowType(); // Snow classification
				EMS[Xdata.SoilNode].dth_w = 0.;
				EMS[Xdata.SoilNode].Qmf = 0.;
				EMS[Xdata.SoilNode].QIntmf = 0.;
				EMS[Xdata.SoilNode].dEps = 0.;
				EMS[Xdata.SoilNode].Eps = EMS[Xdata.SoilNode].Eps_e = EMS[Xdata.SoilNode].Eps_v = EMS[Xdata.SoilNode].Eps_Dot = EMS[Xdata.SoilNode].Eps_vDot = EMS[Xdata.SoilNode].E = 0.;
				EMS[Xdata.SoilNode].S = 0.;
				EMS[Xdata.SoilNode].C = EMS[Xdata.SoilNode].CDot = 0.;
				EMS[Xdata.SoilNode].ps2rb = 0.;
#ifndef SNOWPACK_CORE
				EMS[Xdata.SoilNode].s_strength = 0.;
#endif
				EMS[Xdata.SoilNode].hard = 0.;
#ifndef SNOWPACK_CORE
				EMS[Xdata.SoilNode].S_dr = IOUtils::nodata;
				EMS[Xdata.SoilNode].crit_cut_length = Constants::undefined;
#endif
				EMS[Xdata.SoilNode].VG.theta_r = 0.;
				EMS[Xdata.SoilNode].lwc_source = 0.;
				EMS[Xdata.SoilNode].PrefFlowArea = 0.;
#ifndef SNOWPACK_CORE
				EMS[Xdata.SoilNode].dsm = 0.;
#endif

				EMS[Xdata.SoilNode].h = EMS[Xdata.SoilNode+1].h + .5 * dL;

				// Initial nodal properties
				NDS[Xdata.SoilNode].u = 0.;                     // Initial displacement is 0
				NDS[Xdata.SoilNode].hoar = 0.;                  // The new snow surface hoar is set to zero
				NDS[Xdata.SoilNode].udot = 0.;                  // Settlement rate is also 0
				NDS[Xdata.SoilNode].f = 0.;                     // Unbalanced forces are 0
#ifndef SNOWPACK_CORE
				NDS[Xdata.SoilNode].S_n = IOUtils::nodata;
				NDS[Xdata.SoilNode].S_s = IOUtils::nodata;
#endif
				NDS[Xdata.SoilNode].z = 0.;

				BottomSalFlux += EMS[Xdata.SoilNode].salinity * dL;
			}
		} else {
			// In this case, increase existing element
			const double dL = dM / (EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice);
			dz += dL;
			const double L0 = EMS[Xdata.SoilNode].L;
			EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = (L0 + dL);
			EMS[Xdata.SoilNode].updDensity();
			EMS[Xdata.SoilNode].h += .5 * dL;
			BottomSalFlux += EMS[Xdata.SoilNode].salinity * dL;
			Sdata.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING]+=dM;

		}
	} else {
		// dM < 0: Mass loss
		while (dM < 0. && nE > 0) {
			if(EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice * EMS[Xdata.SoilNode].L + dM > Constants::eps2) {
				// Reduce element length
				const double dL = dM / (EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice);
				EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = EMS[Xdata.SoilNode].L + dL;
				EMS[Xdata.SoilNode].updDensity();
				BottomSalFlux += EMS[Xdata.SoilNode].salinity * dL;
				dz += dL;
				Sdata.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING]+=dM;
				dM = 0.;
			} else {
				// Remove element
				dM += EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice * EMS[Xdata.SoilNode].L;
				Sdata.mass[SurfaceFluxes::MS_ICEBASE_MELTING_FREEZING]-=EMS[Xdata.SoilNode].Rho*EMS[Xdata.SoilNode].L;
				dz += -EMS[Xdata.SoilNode].L;
				// TODO: put mass in SNOWPACK runoff!
				// Add salinity to BottomSalFlux
				BottomSalFlux += EMS[Xdata.SoilNode].salinity * -EMS[Xdata.SoilNode].L;
				if(nE > Xdata.SoilNode) {
					if(EMS[Xdata.SoilNode+1].VG.defined) {
						if(EMS[Xdata.SoilNode+1].h > EMS[Xdata.SoilNode+1].VG.h_e) {
							EMS[Xdata.SoilNode+1].h = EMS[Xdata.SoilNode].h;
						}
					}
					// Shift all existing elements down in the domain
					for(size_t ee = Xdata.SoilNode; ee < nE-1; ee++) {
						EMS[ee]=EMS[ee+1];
						NDS[ee]=NDS[ee+1];
						NDS[ee+1]=NDS[ee+2];
					}
				}
				nE--;
				Xdata.resize(nE);
			}
		}
	}

	// Adjust domain
	NDS[Xdata.SoilNode].z = 0.;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		NDS[e + 1].z = NDS[e].z + EMS[e].L;
	}

	// Adjust externally forced sea level:
	if (ForcedSeaLevel != IOUtils::nodata) {
		ForcedSeaLevel += dz;
	}

	// Ocean water is infinite, so as much ice will be created as energy available, i.e., the bottom node is at meltfreeze_tk!
	calculateMeltingTemperature(EMS[Xdata.SoilNode]);
	if (nE > 0) NDS[Xdata.SoilNode].T = SeaIce::calculateMeltingTemperature(SeaIce::OceanSalinity);
	EMS[Xdata.SoilNode].Te = 0.5 * (NDS[Xdata.SoilNode].T + NDS[Xdata.SoilNode+1].T);
	EMS[Xdata.SoilNode].gradT = (NDS[Xdata.SoilNode+1].T - NDS[Xdata.SoilNode].T) / EMS[Xdata.SoilNode].L;
	return;
}


/**
 * @brief Returns the average bulk salinity (g / kg)
 * @param Xdata Snow cover data
 */
double SeaIce::getAvgBulkSalinity(const SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();
	double ret = 0.;
	double dH = 0.;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		ret += Xdata.Edata[e].salinity * Xdata.Edata[e].Rho * Xdata.Edata[e].L;
		dH += Xdata.Edata[e].Rho * Xdata.Edata[e].L;
	}
	return (dH>0.) ? (ret/dH) : (IOUtils::nodata);
}


/**
 * @brief Returns the average brine salinity (g / kg)
 * @param Xdata Snow cover data
 */
double SeaIce::getAvgBrineSalinity(const SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();
	double ret = 0.;
	double dH = 0.;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		ret += Xdata.Edata[e].theta[WATER] * Xdata.Edata[e].L * (((Xdata.Edata[e].theta[WATER] + Xdata.Edata[e].theta[WATER_PREF]) != 0.) ? (Xdata.Edata[e].salinity / (Xdata.Edata[e].theta[WATER] + Xdata.Edata[e].theta[WATER_PREF])) : (0.));
		dH += Xdata.Edata[e].theta[WATER] * Xdata.Edata[e].L;
	}
	return (dH>0.) ? (ret/dH) : (IOUtils::nodata);
}


/**
 * @brief Returns the total salinity (g / m^2)
 * @param Xdata Snow cover data
 */
double SeaIce::getTotSalinity(const SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();
	double ret = 0.;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		ret += (Xdata.Edata[e].theta[WATER] + Xdata.Edata[e].theta[WATER_PREF]) * Constants::density_water * Xdata.Edata[e].L * Xdata.Edata[e].salinity;
	}
	return ret;
}


/**
 * @brief Initializes a SnowStation object for appropriate sea ice conditions \n
 * First, water and ice content is calculated, while maintaining initial bulk salinity and temperature
 * After that, initialize pressure head consistently with the displaced ocean water
 * @version 21.06: initial version
 * @author Nander Wever
 * @param Xdata The SnowStation object to initialize
 */
void SeaIce::InitSeaIce(SnowStation& Xdata)
{
	const size_t nE = Xdata.getNumberOfElements();
	if (nE==0 || !check_initial_conditions) return;	// Nothing to do...

	double totM = 0.;	// Tracks total mass

	// Set thermodynamical properties consistently (temperature, salinity, etc):
	bool inconsistent_layer = false;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		std::cout << "[i] Initializing sea ice layer " << e << ". Original values: " << " h=" << Xdata.Edata[e].h << " theta[ICE]=" << Xdata.Edata[e].theta[ICE] << " theta[WATER]=" << Xdata.Edata[e].theta[WATER] << " theta[AIR]=" << Xdata.Edata[e].theta[AIR] << " Bulk salinity: " << Xdata.Edata[e].salinity << std::endl;
		// If a layer is reported as dry, no salinity can be present:
		if (Xdata.Edata[e].theta[WATER]<Constants::eps) {
			Xdata.Edata[e].salinity = 0.;
		} else {
			// A given temperature corresponds to a specific brine salinity
			const double BrineSal = calculateBrineSalinity(Xdata.Edata[e].Te);
			if (BrineSal <= 0. && Xdata.Edata[e].salinity > 0.) {
				prn_msg( __FILE__, __LINE__, "err", Date(), "Inconsistent initial condition at layer %d / %d.\n    For the initial layer temperature (%lf), brine salinity is %lf, while prescribed bulk salinity is %lf. This cannot be made consistent.", e, nE, Xdata.Edata[e].Te, BrineSal, Xdata.Edata[e].salinity);
				inconsistent_layer=true;
			}
			// Given brine salinity and provided bulk salinity, we can derive theta[WATER]
			const double theta_water_new = Xdata.Edata[e].salinity / BrineSal;
			// We can now estimate theta[ICE]
			const double theta_ice_new = Xdata.Edata[e].theta[ICE] - ((theta_water_new - Xdata.Edata[e].theta[WATER]) * (Constants::density_water / Constants::density_ice));
			Xdata.Edata[e].theta[ICE] = theta_ice_new;
			Xdata.Edata[e].theta[WATER] = theta_water_new;
			Xdata.Edata[e].theta[AIR] = 1. - Xdata.Edata[e].theta[WATER] - Xdata.Edata[e].theta[ICE];
			if (Xdata.Edata[e].theta[ICE] < Snowpack::min_ice_content) {
				prn_msg( __FILE__, __LINE__, "err", Date(), "Inconsistent initial condition at layer %d / %d.\n    For the initial layer temperature (%lf), brine salinity is %lf, while prescribed bulk salinity is %lf. This would require theta[WATER] = %lf, which for this layer results in a too low ice content of %lf.", e, nE, Xdata.Edata[e].Te, BrineSal, Xdata.Edata[e].salinity, Xdata.Edata[e].theta[WATER], Xdata.Edata[e].theta[ICE]);
				inconsistent_layer=true;
			}
		}
		Xdata.Edata[e].updDensity();
		totM += Xdata.Edata[e].M;
	}

	if (inconsistent_layer) {
		prn_msg( __FILE__, __LINE__, "err", Date(), "Inconsistent initial conditions. See error message(s) above for details.");
		throw;
	}

	// Set pressure head consistently:
	const double hbottom = totM / (Constants::density_water + (SeaIce::betaS * SeaIce::OceanSalinity));	// pressure head at bottom
	double z = .5 * Xdata.Edata[0].L;
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		Xdata.Edata[e].h = (hbottom - Xdata.Edata[e].VG.h_e) - z;
		z += Xdata.Edata[e].L;
	}

	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		std::cout << "[i] Initializing sea ice layer " << e << ". Initialized values: " << " h=" << Xdata.Edata[e].h << " theta[ICE]=" << Xdata.Edata[e].theta[ICE] << " theta[WATER]=" << Xdata.Edata[e].theta[WATER] << " theta[AIR]=" << Xdata.Edata[e].theta[AIR] << " Bulk salinity: " << Xdata.Edata[e].salinity << std::endl;
	}
}

/**
 * @brief The sea ice module\n
 * This function runs the sea ice module of SNOWPACK. \n
 * @version 16.08: initial version
 * @author Nander Wever
 * @param Xdata
 * @param Mdata
 * @param Bdata
 * @param sn_dt SNOWPACK time step (s)
 * @param Sdata
 */
void SeaIce::runSeaIceModule(SnowStation& Xdata, const CurrentMeteo& Mdata, BoundCond& Bdata, const double& sn_dt, SurfaceFluxes& Sdata)
{
	Xdata.Seaice->compSalinityProfile(Xdata);
	Xdata.Seaice->OceanHeatFlux=(Bdata.qg == Constants::undefined)?(0.):(Bdata.qg);
	Xdata.Seaice->bottomIceFormation(Xdata, Mdata, sn_dt, Sdata);
	Xdata.Seaice->compSalinityProfile(Xdata);
	Xdata.Seaice->updateFreeboard(Xdata);
}
