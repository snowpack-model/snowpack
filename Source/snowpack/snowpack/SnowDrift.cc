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

#include <snowpack/SnowDrift.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>

#include <vector>
#include <assert.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

///Deviation from geometrical factors defined by Schmidt
const double SnowDrift::schmidt_drift_fudge = 1.0;

///Enables erosion notification
const bool SnowDrift::msg_erosion = false;



/************************************************************
 * non-static section                                       *
 ************************************************************/

static bool get_bool(const SnowpackConfig& cfg, const std::string& key, const std::string& section)
{
	bool value;
	cfg.getValue(key, section, value);
	return value;
}

static std::string get_erosion(const SnowpackConfig& cfg)
{
	std::string erosion = "NONE";
	cfg.getValue("SNOW_EROSION", "SnowpackAdvanced", erosion);
	std::transform(erosion.begin(), erosion.end(), erosion.begin(), ::toupper);	// Force upper case
	if (erosion != "NONE" && erosion != "VIRTUAL" && erosion != "HS_DRIVEN" && erosion != "FREE" && erosion != "REDEPOSIT") {
		if (erosion == "TRUE") {
			// SNOW_EROSION==TRUE is deprecated and now interpreted as HS_DRIVEN.
			// TODO/HACK: remove this hack in the future.
			erosion="HS_DRIVEN";
		} else if (erosion == "FALSE") {
			// SNOW_EROSION==FALSE is deprecated and now interpreted as NONE.
			// TODO/HACK: remove this hack in the future.
			erosion="NONE";
		} else {
			std::stringstream msg;
			msg << "Value provided for SNOW_EROSION (" << erosion << ") is not valid. Choose either NONE, VIRTUAL, HS_DRIVEN, FREE or REDEPOSIT.";
			throw UnknownValueException(msg.str(), AT);
		}
	}
	return erosion;
}

static bool get_redistribution(const SnowpackConfig& cfg)
{
	bool redistribution = false;
	const int nSlopes = cfg.get("NUMBER_SLOPES", "SnowpackAdvanced");

	// Defines whether real snow erosion at main station or/and redistribution on virtual slopes (default in operational mode)
	// should happen under blowing snow conditions.
	//cfg.getValue("SNOW_EROSION", "SnowpackAdvanced", snow_erosion);
	if (nSlopes>1)
		cfg.getValue("SNOW_REDISTRIBUTION", "SnowpackAdvanced", redistribution);
	
	return redistribution;
}

static double get_sn_dt(const SnowpackConfig& cfg) 
{
	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	return M_TO_S(calculation_step_length);
}

static double get_fetch_length(const SnowpackConfig& cfg)
{
	double fetch_length = Constants::undefined;
	cfg.getValue("SNOW_EROSION_FETCH_LENGTH", "SnowpackAdvanced", fetch_length, IOUtils::nothrow);
	return fetch_length;
}

static double get_erosion_limit(const SnowpackConfig& cfg)
{
	double tmp_erosion_limit = Constants::undefined;
	cfg.getValue("SNOW_EROSION_LIMIT", "SnowpackAdvanced", tmp_erosion_limit, IOUtils::nothrow);
	return tmp_erosion_limit;
}

static double get_redeposit_avg_depth(const SnowpackConfig& cfg)
{
	double tmp_redeposit_avg_depth = Constants::undefined;
	cfg.getValue("SNOW_EROSION_REDEPOSIT_AVG_DEPTH", "SnowpackAdvanced", tmp_redeposit_avg_depth, IOUtils::nothrow);
	return tmp_redeposit_avg_depth;
}

double SnowDrift::get_tau_thresh(const ElementData& Edata)
{
	// Compute basic quantities that are needed: friction velocity, z0, threshold vw
	// For now assume logarithmic wind profile; TODO change this later
	const double weight = 0.02 * Constants::density_ice * (Edata.sp + 1.) * Constants::g * MM_TO_M(Edata.rg);
	// weight = Edata.Rho*(Edata.sp + 1.)*Constants::g*MM_TO_M(Edata.rg);
	const double sig = 300.;
	const double binding = 0.0015 * sig * Edata.N3 * Optim::pow2(Edata.rb/Edata.rg);
	const double tau_thresh = SnowDrift::schmidt_drift_fudge * (weight + binding);  // Original value for fudge: 1. (Schmidt)
	//const double ustar_thresh = sqrt(tau_thresh / Constants::density_air);
	unsigned int wet_marker = (int(Edata.mk/10)%10);
	const double fudge = (wet_marker==1 || wet_marker==2) ? (10.) : (1.);            // For wet/crusted snow, a pragmatic fudge factor makes the snow not easily erodible.
	return fudge * tau_thresh;
}

double SnowDrift::get_ustar_thresh(const ElementData& Edata)
{
	const double ustar_thresh = sqrt(SnowDrift::get_tau_thresh(Edata) / Constants::density_air);
	return ustar_thresh;
}


SnowDrift::SnowDrift(const SnowpackConfig& cfg) : saltation(cfg),
                     enforce_measured_snow_heights( get_bool(cfg, "ENFORCE_MEASURED_SNOW_HEIGHTS", "Snowpack") ), snow_redistribution( get_redistribution(cfg) ), snow_erosion( get_erosion(cfg) ), alpine3d( get_bool(cfg, "ALPINE3D", "SnowpackAdvanced") ),
                     sn_dt( get_sn_dt(cfg) ), fetch_length( get_fetch_length(cfg) ), erosion_limit( get_erosion_limit(cfg) ), redeposit_avg_depth( get_redeposit_avg_depth(cfg) ), forcing("ATMOS")
{
	cfg.getValue("FORCING", "Snowpack", forcing);
}

/**
 * @brief Computes the local mass flux of snow
 * @bug Contribution from suspension not considered yet!
 * @param *Edata
 * @param ustar Shear wind velocity (m s-1)
 * @param angle Slope angle (deg)
 * @return Saltation mass flux (kg m-1 s-1)
 */
double SnowDrift::compMassFlux(const std::vector<ElementData>& EMS, const double& ustar, const double& slope_angle) const
{
	size_t nE = EMS.size();
	const double tau = Constants::density_air * Optim::pow2(ustar);
	double tau_thresh = 0.;
	double rg = 0.;
	if (snow_erosion == "REDEPOSIT") {
		double dH = 0.;
		double tausum = 0.;
		double rgsum = 0.;
		for(size_t e = nE; e --> 0; ) {
			if (EMS[e].theta[SOIL] > 0.) break;
			dH += EMS[e].L;
			tausum += SnowDrift::get_tau_thresh(EMS[e]) * EMS[e].L;
			rgsum += EMS[e].rg * EMS[e].L;
			if (redeposit_avg_depth == Constants::undefined || dH > redeposit_avg_depth) break;
		}
		tau_thresh = (dH > 0.) ? (tausum / dH) : (0.);
		rg = (dH > 0.) ? (rgsum / dH) : (0.);
	} else {
		tau_thresh = SnowDrift::get_tau_thresh(EMS[nE-1]);
		rg = EMS[nE-1].rg;
	}

	// First, look whether there is any transport at all: use formulation of Schmidt
	if ( tau_thresh > tau ) {
		return (0.0);
	}

	// Compute the saltation mass flux (after Pomeroy and Gray)
	double Qsalt = 0., Qsusp = 0., c_salt; // The mass fluxes in saltation and suspension (kg m-1 s-1)
	if (!saltation.compSaltation(tau, tau_thresh, slope_angle, MM_TO_M(2.*rg), Qsalt, c_salt)) {
		prn_msg(__FILE__, __LINE__, "err", Date(), "Saltation computation failed");
		throw IOException("Saltation computation failed", AT);
	}

	if (Qsalt > 0.) {
		Qsusp = 0.; // TODO What about sd_IntegrateFlux(ustar, ustar_thresh, Qsalt, z0); ???
	} else {
		Qsalt = 0.;
		Qsusp = 0.;
	}

	return (Qsalt + Qsusp);
}

/**
 * @brief Erodes Elements from the top and computes the associated mass flux
 * Even so the code is quite obscure, it should cover all of the following cases:
 * -# Externally provided eroded mass (for example, by Alpine3D); see parameter forced_massErode
 * -# SNOW_REDISTRIBUTION is true: using vw_drift to erode the snow surface on the windward virtual slope
 * -# SNOW_EROSION is true: using vw to erode the snow surface at the main station (flat field or slope).
 * 	@note However, erosion will also be controlled by mH and thus a measured snow depth (HS1) is required
 * 	@note If either measured snow depth is missing or the conditions for a real erosion are not fulfilled,
 *          the possibility of a virtual erosion will be considered using the ErosionLevel marker (virtual erodible layer).
 * @param Mdata
 * @param Xdata
 * @param Sdata
 * @param forced_massErode if greater than 0, force the eroded mass to the given value (instead of computing it)
*/
void SnowDrift::compSnowDrift(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& forced_massErode) const
{
	size_t nE = Xdata.getNumberOfElements();
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	const bool no_snow = ((nE < Xdata.SoilNode+1) || (EMS[nE-1].theta[SOIL] > 0.));
	const bool no_wind_data = (Mdata.vw_drift == mio::IOUtils::nodata);
	Xdata.ErosionMass = 0.;
	Xdata.ErosionLength = 0.;
	Xdata.ErosionAge = 0.;
	if (no_snow || no_wind_data) {
		if (no_snow) {
			Xdata.ErosionLevel = Xdata.SoilNode;
			Sdata.drift = 0.;
			return;
		} else
			Sdata.drift = Constants::undefined;
	}

	// Real erosion either on windward virtual slope, from Alpine3D, or at main station.
	// At main station, measured snow depth controls whether erosion is possible or not
	const bool windward = !alpine3d && snow_redistribution && Xdata.windward; // check for windward virtual slope
	const bool erosion = ( 
		( snow_erosion == "FREE" || snow_erosion == "REDEPOSIT"	) 
		|| (
			snow_erosion == "HS_DRIVEN" && (Xdata.mH > (Xdata.Ground + Constants::eps)) && ((Xdata.mH + 0.02) < Xdata.cH)
		)  );
	double ustar = 0.;
	if (windward || erosion || (forced_massErode < -Constants::eps2) || (forcing == "MASSBAL")) {
		double massErode=0.; // Mass loss due to erosion
		if (forcing == "MASSBAL") {
			if (forced_massErode != IOUtils::nodata) {
				massErode = std::max(0., -forced_massErode); //negative mass is erosion
			}
		} else if (forced_massErode < -Constants::eps2) {
			massErode = std::max(0., -forced_massErode); //negative mass is erosion
		} else {
			const double ustar_max = (Mdata.vw>0.1) ? Mdata.ustar * Mdata.vw_drift / Mdata.vw : 0.; // Scale Mdata.ustar
			try {
				if (enforce_measured_snow_heights && !windward) {
					Sdata.drift = compMassFlux(EMS, Mdata.ustar, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, main station, local vw && nE-1
					ustar = Mdata.ustar;
				} else if (windward || (Xdata.meta.getSlopeAngle() < Constants::min_slope_angle)) { // only erode windward or flat (if erosion=1) (maybe rewrite this pile of if-statements entirely?)
					Sdata.drift = compMassFlux(EMS, ustar_max, Xdata.meta.getSlopeAngle()); // kg m-1 s-1, windward slope && vw_drift && nE-1
					ustar = ustar_max;
				} else{	
					Sdata.drift = 0.;
				}
			} catch(const exception&) {
					prn_msg(__FILE__, __LINE__, "err", Mdata.date, "Cannot compute mass flux of drifting snow!");
					throw;
			}
			massErode = Sdata.drift * sn_dt / ((fetch_length != Constants::undefined) ? (fetch_length) : (Hazard::typical_slope_length)); // Convert to eroded snow mass in kg m-2
		}

		unsigned int nErode=0; // number of eroded elements
		for(size_t e = nE; e --> Xdata.SoilNode; ) {
			// Check for limits to halt erosion
			if (erosion_limit != Constants::undefined && erosion_limit != 0.) {
				// A positive value of erosion_limit denotes a snow density above which layer erosion is halted.
				// A negative value of erosion_limit denotes that erosion is halted if a layer threshold friction velocity exceeds the actual friction velocity
				if (erosion_limit > 0.) {
					if (EMS[e].Rho > erosion_limit) break;
				} else {
					if (SnowDrift::get_ustar_thresh(EMS[e]) > ustar) break;
				}
			}

			// Continue with mass erosion:
			if (massErode >= 0.95 * EMS[e].M) {
				// Erode at most one element with a maximal error of +- 5 % on mass ...
				if (windward)
					Xdata.rho_hn = EMS[e].Rho;
				nE--;
				Xdata.cH -= EMS[e].L;
				NDS[e].hoar = 0.;
				Xdata.ErosionMass += EMS[e].M;
				Xdata.ErosionLength -= EMS[e].L;
				Xdata.ErosionLevel = std::min(e, Xdata.ErosionLevel);
				Xdata.ErosionAge += EMS[e].depositionDate.getJulian() * EMS[e].L;
				nErode++;
				massErode -= EMS[e].M;
				forced_massErode = -massErode;
			} else if (massErode > Constants::eps) { // ... or take away massErode from top element - partial real erosion
				if (fabs(EMS[e].L * EMS[e].Rho - EMS[e].M) > 0.001) {
					prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Inconsistent Mass:%lf   L*Rho:%lf", EMS[e].M,EMS[e].L*EMS[e].Rho);
					EMS[e].M = EMS[e].L * EMS[e].Rho;
					assert(EMS[e].M>=0.); //mass must be positive
				}
				if (windward)
					Xdata.rho_hn = EMS[e].Rho; // Density of drifting snow on virtual luv slope
				const double dL = -massErode / (EMS[e].Rho);
				NDS[e+1].z += dL;
				EMS[e].L0 = EMS[e].L = EMS[e].L + dL;
				Xdata.cH += dL;
				NDS[e+1].z += NDS[e+1].u;
				NDS[e+1].u = 0.0;
				NDS[e+1].hoar = 0.;
				EMS[e].M -= massErode;
				assert(EMS[e].M>=0.); //mass must be positive
				Xdata.ErosionMass += massErode;
				Xdata.ErosionLength += dL;
				Xdata.ErosionAge += EMS[e].depositionDate.getJulian() * -dL;
				massErode = 0.;
				forced_massErode = 0.;
				break;
			} else {
				Xdata.ErosionMass += 0.;
				Xdata.ErosionLength += 0.;
				break;
			}
			if (snow_erosion == "HS_DRIVEN") break;	// To be consistent with legacy SNOWPACK where only one element at a time can erode.
		}
		if (nErode > 0)
			Xdata.resize(nE);

		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionMass > 0.) {
				if (windward)
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (windward: azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
				else
					prn_msg(__FILE__, __LINE__, "msg+", Mdata.date, "Eroding %d layer(s) w/ total mass %.3lf kg/m2 (azi=%.1lf, slope=%.1lf)", nErode, Xdata.ErosionMass, Xdata.meta.getAzimuth(), Xdata.meta.getSlopeAngle());
			}
		}
	// ... or, in case of no real erosion, check whether you can potentially erode at the main station.
	// This will never contribute to the drift index VI24, though!
	} else if (snow_erosion == "VIRTUAL" && (Xdata.ErosionLevel > Xdata.SoilNode)) {
		std::vector<ElementData> subEMS; subEMS.push_back(EMS[Xdata.ErosionLevel]);  // Subset EMS vector with only the ErosionLevel ElementData
		double virtuallyErodedMass = compMassFlux(subEMS, Mdata.ustar, Xdata.meta.getSlopeAngle());  // kg m-1 s-1, main station, local vw && erosion level
		virtuallyErodedMass *= sn_dt / Hazard::typical_slope_length; // Convert to eroded snow mass in kg m-2
		if (virtuallyErodedMass > Constants::eps) {
			// Add (negative) value stored in Xdata.ErosionMass
			if (Xdata.ErosionMass < -Constants::eps)
				virtuallyErodedMass -= Xdata.ErosionMass;
			Sdata.mass[SurfaceFluxes::MS_WIND] = std::min(virtuallyErodedMass, EMS[Xdata.ErosionLevel].M); // use MS_WIND to carry virtually eroded mass
			// Now keep track of mass that either did or did not lead to erosion of full layer
			if (virtuallyErodedMass > EMS[Xdata.ErosionLevel].M) {
				virtuallyErodedMass -= EMS[Xdata.ErosionLevel].M;
				Xdata.ErosionLevel--;
			}
			Xdata.ErosionMass = -virtuallyErodedMass;
			Xdata.ErosionLevel = std::max(Xdata.SoilNode, std::min(Xdata.ErosionLevel, nE-1));
		} else {
			Xdata.ErosionMass = 0.;
		}
		if (!alpine3d && SnowDrift::msg_erosion) { //messages on demand but not in Alpine3D
			if (Xdata.ErosionLevel > nE-1)
				prn_msg(__FILE__, __LINE__, "wrn", Mdata.date, "Virtual erosion: ErosionLevel=%d did get messed up (nE-1=%d)", Xdata.ErosionLevel, nE-1);
		}
	} else {
		Xdata.ErosionMass = 0.;
	}
	Sdata.mass[SurfaceFluxes::MS_EROSION_DHS] += Xdata.ErosionLength;
	if (Xdata.ErosionLength < 0.) {
		Xdata.ErosionAge /= -Xdata.ErosionLength;
	} else {
		Xdata.ErosionAge = Constants::undefined;
	}
	return;
}
