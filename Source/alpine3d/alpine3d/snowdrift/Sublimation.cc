/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
/*************************************************************************/
/* 	calculate sublimation drifting snow    					             */
/*************************************************************************/
/*------------------------------------------------------------------------+
 |  													    	          |
 +-----------------------------------------------------------------------*/
/*************************************************************************/
#include <alpine3d/snowdrift/SnowDrift.h>

using namespace mio;

/**
 * @brief Sublimation
 * Calculates drifting snow sublimation at every node of the concentration field
 */
void SnowDriftA3D::Sublimation()
{
	nodes_Subl.grid3D = 0.;

	for (unsigned int ix=0; ix<nx;ix++) {
		for (unsigned int iy=0; iy<ny; iy++) {
			for (unsigned int iz=1; iz<nz; iz++)	{ //start at second level since below > BC saturation
				//calculate magnitude windspeed instead of u,v,w
				const double WindVelocity = sqrt(mio::Optim::pow2(nodes_u.grid3D(ix,iy,iz))+mio::Optim::pow2(nodes_v.grid3D(ix,iy,iz))+mio::Optim::pow2(nodes_w.grid3D(ix,iy,iz)));
				if ( nodes_c.grid3D(ix,iy,iz)>1e-9 && nodes_RH.grid3D(ix,iy,iz)<1.0){
					const double repRadius = 6.25e-5; //representative radius for sublimation of a particle size distribution
					const double dmdt_ini = calcSubldM(repRadius, nodes_Tair.grid3D(ix,iy,iz), nodes_RH.grid3D(ix,iy,iz), WindVelocity, nodes_z.grid3D(ix,iy,iz));
					nodes_Subl.grid3D(ix,iy,iz) = calcS(nodes_c.grid3D(ix,iy,iz),repRadius,dmdt_ini);
				}
			}
		}
	}
}

/**
 * @brief Calculate sublimation mass change
 * This function calculates the mass change of a single ice particle for the given conditions
 * @param Radius radius of the ice particle
 * @param AirTemperature
 * @param RH relative humidity
 * @param WindSpeed magnitude of the wind
 */
double SnowDriftA3D::calcSubldM(const double Radius,const double AirTemperature,const double RH, const double WindSpeed, const double altitude)
{
	const double sigma = RH - 1.; //Undersaturation of the surrounding air

	//Nusselt number Nu
	double Nu;
	const double Re = reynoldsNumberforFallingParticles(Radius, WindSpeed, AirTemperature, RH, altitude); //Reynolds number
	if (Re<=10.) {
		Nu = 1.79+0.606*sqrt(Re);
	} else {
		Nu = 1.88+0.580*sqrt(Re);
	}

	//Sherwood number
	const double Sh = Nu;

	//diffusivity of water vapor in the atmosphere
	const double diffusivity = 2.11e-5 * pow((AirTemperature/273.15), 1.75);

	//water vapor density in the atmosphere
	double const SatVaporDensity = Atmosphere::waterVaporDensity(AirTemperature,Atmosphere::vaporSaturationPressure(AirTemperature));

	//Calculate dm/dt according to Thorpe and Mason (we split it in 3 parts)
	const double temp0 = (Constants::lh_sublimation*molecularWeightofWater) / (Constants::gas_constant_mol*AirTemperature) - 1.;
	const double temp1 = 2.*Constants::pi*Radius*sigma;
	const double temp2 = Constants::lh_sublimation / (thermalConductivityofAtm*AirTemperature*Nu) * temp0 +
	                     1. / (diffusivity * SatVaporDensity * Sh) ;

	//And return the value for dm/dt for 1 particle
	const double dmdt = temp1/temp2;
	return dmdt;
}

/**
 * @brief Calculate sublimation
 * The mass change of a single particle is now extended to all particles
 * @param concentration snow concentration
 * @param sublradius radius of a single particle
 * @param dmdt mass change of a single particle due to sublimation
 */
double SnowDriftA3D::calcS(const double concentration,const double sublradius, const double dmdt)
{
	const double subllossrate = (concentration / (4./3.*Constants::pi*pow(sublradius,3)*Constants::density_ice)) * dmdt;
	return subllossrate;
}

/**
 * @brief Reynolds number
 * Calculate the Reynolds number for a falling particle
 * @param Radius radius of the particle
 * @param Windspeed magnitude of the wind
 * @param AirTemperature
 * @param RH relative humidity of the air
 */
double SnowDriftA3D::reynoldsNumberforFallingParticles(const double Radius, const double Windspeed,const double AirTemperature,const double RH, const double altitude)
{
	const double tmp_venVelocity = ventilationVelocity(Radius, Windspeed, AirTemperature, RH, altitude);
	const double ReynoldsNumber = (2.*Radius*tmp_venVelocity) / kinematicViscosityAir;
	return ReynoldsNumber;
}

/**
 * @brief Ventilation velocity
 * @param Radius radius of the particle
 * @param Windspeed magnitude of the wind
 * @param AirTemperature
 * @param RH relative humidity of the air
 * @param altitude Altitude above sea level
 */
double SnowDriftA3D::ventilationVelocity(const double Radius, const double Windspeed, const double AirTemperature, const double RH, const double altitude)
{
	const double MeanTerminalFallVelocity = terminalFallVelocity(Radius, AirTemperature, RH, altitude);
	const double FluctuatingComponent = 7.5e-3*sqrt(2)*pow(Windspeed, 1.36);
	const double VentilationVelocity = MeanTerminalFallVelocity + FluctuatingComponent;

	return VentilationVelocity;
}

/**
 * @brief Terminal fall velocity
 * NB assumed to be 0.5
 * @param Radius radius of the particle
 * @param Temperature air temperature (K)
 * @param RH relative humidity of the air (between 0 and 1)
 * @param altitude Altitude above sea level
 */
double SnowDriftA3D::terminalFallVelocity(const double /*Radius*/, const double /*Temperature*/, const double /*RH*/ , const double /*altitude*/)
{
	//problems for radius larger than 0.000054>> try constant velocity 0.5
	return 0.5;//settling velocity assumed 0.5 m/s

	/*const double tmp1 = (6.203*kinematicViscosityAir) / Radius;
	const double Psat = Atmosphere::vaporSaturationPressure(Temperature);
	const double sigma = Constants::density_ice / ( Atmosphere::stdDryAirDensity(altitude, Temperature) +
	                     Atmosphere::waterVaporDensity(Temperature, Psat*RH) );
	const double vr = .5 * (-tmp1 + sqrt((tmp1*tmp1)-4.*(-1.379)*mio::Cst::gravity*sigma*Radius) );
	return vr;*/
}

/**
 * @brief Relative humidity
 * Gives the relative humidity for a given air temperature and specific humidity
 * @param AirTemp Air temperature
 * @param qi specific humidity
 * @param altitude Altitude above sea level
 */
double SnowDriftA3D::RH_from_q(const double AirTemp, const double qi, const double altitude)
{
	const double SatVaporDensity = Atmosphere::waterVaporDensity(AirTemp,Atmosphere::vaporSaturationPressure(AirTemp));
	const double tmp = (qi/(1.-qi)) * Atmosphere::stdDryAirDensity(altitude, AirTemp);

	const double RH = std::min(tmp/SatVaporDensity , 1.5);
	return RH;
}

