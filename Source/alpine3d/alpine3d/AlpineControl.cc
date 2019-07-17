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
#include <alpine3d/AlpineControl.h>
#include <alpine3d/AlpineMain.h>
#include <meteoio/MeteoIO.h>

using namespace mio;
using namespace std;

/**
 * @brief Constructs and initialise the time steps loop.
 * This module calls all the necessary other modules (ie thoese that have been enabled) and runs through the time steps.
 * @param mysnowpack pointer to the initialized SNOWPACK Manager
 * @param mysnowdrift pointer to the initialized Snowdrift Manager
 * @param myeb pointer to the initialized radiation manager
 * @param myda pointer to the initialized data assimilation Manager
 * @param myrunoff pointer to the initialized runoff Manager
 * @param cfg User configuration keys
 * @param dem DEM defining the simulation
 */
AlpineControl::AlpineControl(SnowpackInterface *mysnowpack, SnowDriftA3D *mysnowdrift, EnergyBalance *myeb, DataAssimilation *myda, Runoff *myrunoff, const Config& cfg, const DEMObject& dem)
              : meteo(cfg, dem), snowpack(mysnowpack), snowdrift(mysnowdrift), eb(myeb), da(myda), runoff(myrunoff),
                snow_days_between(0.), max_run_time(-1.), enable_simple_snow_drift(false), nocompute(false), out_snow(true)
{
	cfg.getValue("SNOW_WRITE", "Output", out_snow);
	if (out_snow) {
		cfg.getValue("SNOW_DAYS_BETWEEN", "Output", snow_days_between);
	}

	//check if simple snow drift is enabled
	enable_simple_snow_drift = false;
	cfg.getValue("SIMPLE_SNOW_DRIFT", "Alpine3D", enable_simple_snow_drift, IOUtils::nothrow);

	//check if maximum run time is specified
	cfg.getValue("MAX_RUN_TIME", "Alpine3D", max_run_time, IOUtils::nothrow);
}

void AlpineControl::Run(Date i_startdate, const unsigned int max_steps)
{// This function organizes the whole simulation:
// initializations, time loop, getting input data, computing each time step and finally writing outputs
	Date calcDate(i_startdate); //Date object initialized to julian 0.0 in local time zone with no DST
	const double timeStep = dt_main/86400.;
	Timer elapsed;
	std::vector<MeteoData> vecMeteo; // to transfer meteo information
	mio::Grid2DObject p, psum, psum_ph, vw, vw_drift, dw, rh, ta, tsg, ilwr;
	const bool isMaster = MPIControl::instance().master();

	if (isMaster) {
		cout << "\n**** Done initializing\n";
		cout << "**** Starting Calculation on date: " << calcDate.toString(Date::ISO) << " using Alpine3D version " << A3D_VERSION << "\n";
		if (nocompute) 
			cout << "**** Performing dry run (--no-compute option)\n";
		cout << "\n";
		
		if (nocompute) {
			const Grid2DObject maskGlacier( snowpack->getGrid(SnGrids::GLACIER) );
			meteo.setGlacierMask(maskGlacier);
		}
	}

	//if the meteo data would need to be resampled, we try to fill the buffer with a date a little bit before
	if (snowdrift) meteo.setSkipWind(true); //do not fill grids if met3D
	meteo.prepare(i_startdate);
	
	elapsed.start();
	for (unsigned int t_ind=0; t_ind<max_steps; t_ind++) { //main computational loop
		const double elapsed_start = elapsed.getElapsed();
		const double est_completion = elapsed_start * ((double)max_steps/(double)(t_ind+1) - 1.);

		if (isMaster) {
			cout << "\nSimulation step " << t_ind+1 << "/" << max_steps << " at time step " << calcDate.toString(mio::Date::ISO_TZ) << "\n";
			cout << std::fixed << "Elapsed time: " << setprecision(1) << elapsed_start << " seconds\nEstimated completion in " << est_completion/3600. << " hours\n";
		}
		
		//for --no-compute, simply check the data and move on
		if (nocompute) {
			meteo.prepare(calcDate); //prepare the current timestep (because it could not be prepared before)
			meteo.checkMeteoForcing(calcDate);
			calcDate += timeStep; //move to next time step
			continue;
		}

		//get 1D and 2D meteo for the current time step
		try {
			meteo.get(calcDate, vecMeteo);
			meteo.get(calcDate, ta, tsg, rh, psum, psum_ph, vw, vw_drift, dw, p, ilwr);
		} catch (IOException&) {
			//saving state files before bailing out
			if (isMaster) {
				if (out_snow && t_ind>0 && snowpack) 
					snowpack->writeOutputSNO(calcDate-1./24.); //output for last hour
				throw;
			}
		}
		
		if (t_ind < (max_steps-1)) {
			meteo.prepare(calcDate+timeStep); //prepare next timestep
		}

		if (eb) {
			eb->setStations(vecMeteo);
		}

		if (!snowdrift) { //otherwise snowdrift calls snowpack.setMeteo()
			if (snowpack) snowpack->setMeteo(psum, psum_ph, vw, dw, rh, ta, tsg, calcDate);
		}
		if (snowpack && enable_simple_snow_drift) snowpack->setVwDrift(vw_drift, calcDate);

		try { //Snowdrift
			if (snowdrift) {
				//meteo1d(MeteoData::TA) : big TODO, see with Christine if we could get rid of it
				snowdrift->setMeteo(t_ind, psum, psum_ph, p, vw, rh, ta, tsg, ilwr, calcDate, vecMeteo);
				snowdrift->Compute(calcDate);
			}
		} catch (std::exception& e) {
			cout << "[E] Exception: Snowdrift compute\n";
			cout << e.what() << endl;
			throw;
		}

		try {
			if (eb && !snowdrift) { //otherwise snowdrift calls eb.setMeteo()
				eb->setMeteo(ilwr, ta, rh, p, calcDate);
			}
		} catch (std::bad_alloc&) {
			cout << "[E] AlpineControl : Virtual memory exceeded\n";
		} catch (std::exception& e) {
			cout << "[E] Exception: Ebalance compute\n";
			cout << e.what() << endl;
			throw;
		}

		try { //Data Assimilation
			if (da) da->Compute(calcDate);
		} catch (std::exception& e) {
			cout << "[E] Exception: Data Assimilation compute\n";
			cout << e.what() << endl;
			throw;
		}

		// Check if elapsed time exceeds specified maximum run time
		const double tmp_elapsed = elapsed.getElapsed();
		const bool ForceStop = (max_run_time > 0. && tmp_elapsed > max_run_time);
		if (!(max_run_time < 0.)) {
			if (ForceStop) {
				cout << std::fixed << "[W] !!! Elapsed time (" << setprecision(1) << tmp_elapsed << " seconds) exceeds specified maximum run time (" << setprecision(1) << max_run_time << " seconds) !!!\n        ---> Force writing restart files (if requested) and exiting...\n";
			} else {
				cout << std::fixed << "[i] Maximum run time set to: " << setprecision(1) << max_run_time << " seconds ---> time remaining: " << (max_run_time - tmp_elapsed)/3600. << " hours\n";
			}
		}

		try { //do some outputs (note, if ForceStop == true, output will be done outside of the loop)
			if ( snowpack && out_snow && (t_ind > 0) && !ForceStop ) {
				const unsigned int output_step = static_cast<unsigned int>( Optim::round(snow_days_between*24.) );
				if (snow_days_between>0 && (t_ind%output_step)==0)
					snowpack->writeOutputSNO(calcDate);
			}
		} catch (std::exception& e) {
			printf("[e] Exception: time dependent outputs\n");
			cout << e.what() << endl;
			throw;
		}

		if (isMaster) {
			cout << "[i] timing (seconds spent): " << std::setprecision(3) << std::fixed<< "\n";
			cout << "\tmeteo=" << meteo.getTiming() << "  ";

			if (eb) cout << "ebalance=" << eb->getTiming() << "  ";
			if (snowdrift) cout << "snowdrift=" << snowdrift->getTiming() << "  ";
			if (snowpack) cout << "snowpack=" << snowpack->getTiming() << "  ";
			if (runoff) cout << "runoff=" << runoff->getTiming() << " ";

			cout << "\n\ttotal=" << elapsed.getElapsed()-elapsed_start << endl;
		}

		if (ForceStop) break;
		calcDate += timeStep; //move to next time step

	} /* For all times max_steps */

	//Finish the program: Write SNO Files and put final output on the screen
	if (snowpack && out_snow && !nocompute){
		if (isMaster) cout << "[i] Simulation finished, writing output files...\n";
		snowpack->writeOutputSNO(calcDate);
	 }
}

