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
 * @file Snowpack.cc
 * @version 11.06
 * @bug     -
 * @brief This module contains the driving routines for the 1d snowpack model
 */

#include <snowpack/Constants.h>
#include <snowpack/snowpackCore/Metamorphism.h>
#include <snowpack/snowpackCore/SeaIce.h>
#include <snowpack/snowpackCore/ReSolver1d.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

//Threshold that defines ice
const double SeaIce::SeaWaterFreezingTemp = IOUtils::C_TO_K(-1.95);
const double SeaIce::SeaIceDensity = ReSolver1d::max_theta_ice * Constants::density_ice;
const double SeaIce::ice_threshold = 800.;
const double SeaIce::mu = 0.054;
const double SeaIce::ThicknessFirstIceLayer = 0.01;
const double SeaIce::InitRg = 5.;
const double SeaIce::InitRb = 2.5;
const double SeaIce::OceanSalinity = 35.;
const double SeaIce::InitSeaIceSalinity = 5.;
const double SeaIce::InitSnowSalinity = 1.;


/************************************************************
 * non-static section                                       *
 ************************************************************/

SeaIce::SeaIce():
	SeaLevel(0.), FreeBoard (0.), IceSurface(0.), IceSurfaceNode(0), OceanHeatFlux(0.), salinityprofile(SINUSSAL) {}

SeaIce& SeaIce::operator=(const SeaIce& source) {
	if(this != &source) {
		SeaLevel = source.SeaLevel;
		FreeBoard = source.FreeBoard;
		IceSurface = source.IceSurface;
		IceSurfaceNode = source.IceSurfaceNode;
		OceanHeatFlux = source.OceanHeatFlux;
	}
	return *this;
}

SeaIce::~SeaIce() {}

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
 * @param Xdata Snow cover data
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
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			/*size_t e = Xdata.SoilNode;
			for (; e < IceSurfaceNode; e++) {
				Xdata.Edata[e].salinity = 35.;			// Default: 35 g/kg
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			for (; e < nE; e++) {
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
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			// define salinity in snow
			for (; e <  nE ; e++) {
				Xdata.Edata[e].salinity = 1;
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
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			// define salinity in snow
			for (; e <  nE ; e++) {
				Xdata.Edata[e].salinity = 1; // 8 after Massom et al. 1997
				calculateMeltingTemperature(Xdata.Edata[e]);
			}
			break;
		}
	}
}

/**
 * @brief Updates the freeboard variable (i.e., sea level with respect to ice surface)\n
 *        positive: sea level below ice surface\n
 *        negative: sea level above ice surface (flooding)\n
 * @version 16.08
 */
void SeaIce::updateFreeboard(SnowStation& Xdata)
{
	Xdata.compSnowpackMasses();
	SeaLevel = (Xdata.swe / Constants::density_water);
	const double FreeBoard_snow = Xdata.cH - SeaLevel;	// This is the freeboard relative to snow surface
	FreeBoard = (findIceSurface(Xdata) - (Xdata.cH - FreeBoard_snow));
	return;
}

/**
 * @brief Find marked layer level\n
 * @version 16.08
 */

/**
 * @brief Find snow/ice transition for sea ice simulations\n
 * @version 16.08
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
	if (Xdata.Edata[Xdata.SoilNode].theta[ICE] * Constants::density_ice < ice_threshold) {
		IceSurface = 0.;
		IceSurfaceNode = 0;
		return IceSurface;
	}
	for (size_t e = Xdata.SoilNode; e < nE-1; e++) {
		if (Xdata.Edata[e].theta[ICE] * Constants::density_ice > ice_threshold && Xdata.Edata[e+1].theta[ICE] * Constants::density_ice < ice_threshold) {
			IceSurface = Xdata.Ndata[e+1].z;
			IceSurfaceNode = e+1;
			return IceSurface;
		}
	}
	IceSurfaceNode = nE;
	IceSurface = Xdata.Ndata[nE].z;
	return IceSurface;
}

/**
 * @brief Apply flooding\n
 * @version 16.08
 */
void SeaIce::compFlooding(SnowStation& Xdata)
{
	size_t iN = 0;
	while (Xdata.Ndata[iN].z < SeaLevel && iN < Xdata.getNumberOfElements()) {
		const double dth_w = std::max(0., Xdata.Edata[iN].theta[AIR] * (Constants::density_ice / Constants::density_water) - Xdata.Edata[iN].theta[WATER] * (Constants::density_water / Constants::density_ice - 1.));
		Xdata.Edata[iN].theta[WATER] += dth_w;
		Xdata.Edata[iN].theta[AIR] -= dth_w;
		Xdata.Edata[iN].M += dth_w * Constants::density_water * Xdata.Edata[iN].L;
		Xdata.Edata[iN].Rho = Xdata.Edata[iN].M / Xdata.Edata[iN].L;
		Xdata.Edata[iN].salinity += SeaIce::OceanSalinity * dth_w;
		Xdata.Edata[iN].salinity = std::min(SeaIce::OceanSalinity, Xdata.Edata[iN].salinity);
		calculateMeltingTemperature(Xdata.Edata[iN]);
		iN++;
	}
	return;
}


/**
 * @brief Calculate melting temperature as function of salinity
 * @version 16.08
 * @param Edata
 */
void SeaIce::calculateMeltingTemperature(ElementData& Edata)
{
	// See: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	//      who is citing: Assur, A., Composition of sea ice and its tensile strength, in Arctic Sea Ice, N.  A.  S. N.  R.  C. Publ., 598, 106-138, 1958.
	Edata.melting_tk = Edata.freezing_tk = SeaIce::calculateMeltingTemperature(Edata.salinity);
	return;
}


/**
 * @brief Calculate melting temperature as function of salinity
 * @version 17.12: initial version
 * @param Sal: Salinity (PSU, which is g/kg)
 */
const double SeaIce::calculateMeltingTemperature(const double& Sal)
{
	// See: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	//      who is citing: Assur, A., Composition of sea ice and its tensile strength, in Arctic Sea Ice, N.  A.  S. N.  R.  C. Publ., 598, 106-138, 1958.
	return IOUtils::C_TO_K(-SeaIce::mu * Sal);
}


/**
 * @brief Heat capacity of sea ice.
 * @version 16.08: initial version
 * @param T: Temperature (K)
 * @param Sal: Salinity (PSU, which is g/kg)
 * @return Heat capacity for sea ice (J / kg / K)
 */
const double SeaIce::compSeaIceHeatCapacity(const double& T, const double& Sal)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 1 and 2
	const double L0 = Constants::lh_fusion;
	const double c0 = Constants::specific_heat_ice;
	return c0 + (SeaIce::mu * L0 * Sal) / (T * T);
}


/**
 * @brief Heat conduction in sea ice.
 * @version 16.08: initial version
 * @param Edata
 * @return Thermal conductivity for sea ice (W K-1 m-1)
 */
const double SeaIce::compSeaIceThermalConductivity(const ElementData& Edata)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 9
	const double beta =  0.1172;		// W/m^2/permille
	const double k0 = 2.034;		// W/m/K, note that this is the thermal conductivity of fresh ice, and it may be coupled to the value in Constants.h
	// Note the conversion from kg/kg to permille for salinity
	return (k0 + ((beta * Edata.salinity) / Edata.Te));
}


/**
 * @brief Latent heat of melting for sea ice.
 * @version 16.08: initial version
 * @param T: Temperatur (K)
 * @param Sal: Salinity (PSU, which is g/kg)
 * @return Latent heat of fusion for sea ice (J / kg)
 */
const double SeaIce::compSeaIceLatentHeatFusion(const double& T, const double& Sal)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 5
	const double L0 = Constants::lh_fusion;
	return L0 * (1. + (SeaIce::mu * Sal) / T);
}


/**
 * @brief Latent heat of melting for sea ice.
 * @version 16.08: initial version
 * @param Edata
 * @return Latent heat of fusion for sea ice (J / kg)
 */
const double SeaIce::compSeaIceLatentHeatFusion(const ElementData& Edata)
{
	// From: Bitz, C. M., and W. H. Lipscomb (1999), An energy-conserving thermodynamic model of sea ice, J. Geophys. Res., 104(C7), 15669–15677, doi:10.1029/1999JC900100.
	// See Eq. 5
	const double L0 = Constants::lh_fusion;
	const double c0 = Constants::specific_heat_ice;
	return c0 * (Edata.melting_tk - Edata.Te) + L0 * (1. + (SeaIce::mu * Edata.salinity) / Edata.Te);
}


/**
 * @brief Calculate ice formation and decay at the bottom
 * @version 16.08: initial version
 * @param Edata
 */
void SeaIce::bottomIceFormation(SnowStation& Xdata, const CurrentMeteo& Mdata, const double& sn_dt)
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

	// Apply mass change:
	if (dM > 0) {
		// dM > 0: mass gain
		const double dL = dM / SeaIceDensity;
		if ( nE == 0 || EMS[Xdata.SoilNode].Rho < ice_threshold ) {
			// In these case, add new element with ice density
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
				NDS[nE].T = SeaIce::SeaWaterFreezingTemp;
			}
			// Set the new ice element
			EMS[Xdata.SoilNode].depositionDate = Mdata.date;
			EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = dL;
			EMS[Xdata.SoilNode].theta[SOIL] = 0.;
			EMS[Xdata.SoilNode].theta[WATER] = 0.;
			EMS[Xdata.SoilNode].theta[WATER_PREF] = 0.;
			EMS[Xdata.SoilNode].theta[ICE] = (SeaIceDensity/Constants::density_ice);
			EMS[Xdata.SoilNode].theta[AIR] = 1.0 - EMS[Xdata.SoilNode].theta[WATER] - EMS[Xdata.SoilNode].theta[WATER_PREF] - EMS[Xdata.SoilNode].theta[ICE] - EMS[Xdata.SoilNode].theta[SOIL];
			EMS[Xdata.SoilNode].Rho = (EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice) + ((EMS[Xdata.SoilNode].theta[WATER]+EMS[Xdata.SoilNode].theta[WATER_PREF])
					*Constants::density_water) + (EMS[Xdata.SoilNode].theta[SOIL]
					  * EMS[Xdata.SoilNode].soil[SOIL_RHO]);
			EMS[Xdata.SoilNode].M = dM;

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
			EMS[Xdata.SoilNode].mk = 7.;
			EMS[Xdata.SoilNode].metamo = 0.;
			EMS[Xdata.SoilNode].snowType(); // Snow classification
			EMS[Xdata.SoilNode].salinity = SeaIce::InitSeaIceSalinity;
			EMS[Xdata.SoilNode].dth_w = 0.;
			EMS[Xdata.SoilNode].Qmf = 0.;
			EMS[Xdata.SoilNode].QIntmf = 0.;
			EMS[Xdata.SoilNode].dEps = 0.;
			EMS[Xdata.SoilNode].Eps = EMS[Xdata.SoilNode].Eps_e = EMS[Xdata.SoilNode].Eps_v = EMS[Xdata.SoilNode].Eps_Dot = EMS[Xdata.SoilNode].Eps_vDot = EMS[Xdata.SoilNode].E = 0.;
			EMS[Xdata.SoilNode].S = 0.;
			EMS[Xdata.SoilNode].C = EMS[Xdata.SoilNode].CDot = 0.;
			EMS[Xdata.SoilNode].ps2rb = 0.;
			EMS[Xdata.SoilNode].s_strength = 0.;
			EMS[Xdata.SoilNode].hard = 0.;
			EMS[Xdata.SoilNode].S_dr = INIT_STABILITY;
			EMS[Xdata.SoilNode].crit_cut_length = Constants::undefined;
			EMS[Xdata.SoilNode].VG.theta_r = 0.;
			EMS[Xdata.SoilNode].lwc_source = 0.;
			EMS[Xdata.SoilNode].PrefFlowArea = 0.;
			EMS[Xdata.SoilNode].dsm = 0.;
			
			// Initial nodal properties
			NDS[Xdata.SoilNode].u = 0.;                     // Initial displacement is 0
			NDS[Xdata.SoilNode].hoar = 0.;                  // The new snow surface hoar is set to zero
			NDS[Xdata.SoilNode].udot = 0.;                  // Settlement rate is also 0
			NDS[Xdata.SoilNode].f = 0.;                     // Unbalanced forces are 0
			NDS[Xdata.SoilNode].S_n = INIT_STABILITY;
			NDS[Xdata.SoilNode].S_s = INIT_STABILITY;
			NDS[Xdata.SoilNode].z = 0.;
		} else {
			// In this case, increase existing element
			const double L0 = EMS[Xdata.SoilNode].L;
			EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = (L0 + dL);

			EMS[Xdata.SoilNode].theta[WATER] *= L0 / (L0 + dL);
			EMS[Xdata.SoilNode].theta[WATER_PREF] *= L0 / (L0 + dL);
			EMS[Xdata.SoilNode].theta[ICE] = (EMS[Xdata.SoilNode].theta[ICE] * L0 + dL * (SeaIceDensity/Constants::density_ice)) / (L0 + dL);
			EMS[Xdata.SoilNode].theta[AIR] = 1.0 - EMS[Xdata.SoilNode].theta[WATER] - EMS[Xdata.SoilNode].theta[WATER_PREF] - EMS[Xdata.SoilNode].theta[ICE] - EMS[Xdata.SoilNode].theta[SOIL];
			EMS[Xdata.SoilNode].Rho = (EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice) + ((EMS[Xdata.SoilNode].theta[WATER]+EMS[Xdata.SoilNode].theta[WATER_PREF])
					*Constants::density_water) + (EMS[Xdata.SoilNode].theta[SOIL] * EMS[Xdata.SoilNode].soil[SOIL_RHO]);
			EMS[Xdata.SoilNode].M += dM;
		}
	} else {
		// dM < 0: Mass loss
		while (dM < 0. && nE > 0) {
			const double dL = dM / SeaIceDensity;
			if(EMS[Xdata.SoilNode].M + dM > Constants::eps2 && EMS[Xdata.SoilNode].L + dL > Constants::eps2) {
				// Reduce element length
				EMS[Xdata.SoilNode].Rho = (EMS[Xdata.SoilNode].theta[ICE] * Constants::density_ice) + ((EMS[Xdata.SoilNode].theta[WATER]+EMS[Xdata.SoilNode].theta[WATER_PREF])
						*Constants::density_water) + (EMS[Xdata.SoilNode].theta[SOIL]
						  * EMS[Xdata.SoilNode].soil[SOIL_RHO]);
				EMS[Xdata.SoilNode].M += dM;
				EMS[Xdata.SoilNode].L0 = EMS[Xdata.SoilNode].L = EMS[Xdata.SoilNode].M / EMS[Xdata.SoilNode].Rho;
				dM = 0.;
			} else {
				// Remove element
				dM += EMS[Xdata.SoilNode].M;
				// TODO: put mass in SNOWPACK runoff!
				if(nE > 1) {
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

	// Ocean water is infinite, so as much ice will be created as energy available, i.e., the bottom node is at melting_tk!
	calculateMeltingTemperature(EMS[Xdata.SoilNode]);
	if (nE > 0) NDS[Xdata.SoilNode].T = SeaIce::calculateMeltingTemperature(SeaIce::OceanSalinity);
	EMS[Xdata.SoilNode].Te = 0.5 * (NDS[Xdata.SoilNode].T + NDS[Xdata.SoilNode+1].T);
	EMS[Xdata.SoilNode].gradT = (NDS[Xdata.SoilNode+1].T - NDS[Xdata.SoilNode].T) / EMS[Xdata.SoilNode].L;
	return;
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
 */
void SeaIce::runSeaIceModule(SnowStation& Xdata, const CurrentMeteo& Mdata, BoundCond& Bdata, const double& sn_dt)
{
	Xdata.Seaice->compSalinityProfile(Xdata);
	Xdata.Seaice->OceanHeatFlux=(Bdata.qg == Constants::undefined)?(0.):(Bdata.qg);
	Xdata.Seaice->bottomIceFormation(Xdata, Mdata, sn_dt);
	Xdata.Seaice->compSalinityProfile(Xdata);
	Xdata.Seaice->updateFreeboard(Xdata);
}
