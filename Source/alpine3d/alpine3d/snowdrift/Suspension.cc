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
/*------------------------------------------------------------------------------------------+
 |  This module is the FINITE ELEMENT solution to Snow Diffusion                            |
 +------------------------------------------------------------------------------------------*/
/********************************************************************************************/
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/AlpineMain.h>
#include <iostream>

using namespace mio;
using namespace std;

#define DRIFT_MARC 1

void SnowDriftA3D::Diffusion(double deltaT,double &diff_max, double t1)
{//OLD VERSION
	(void)deltaT;
	(void)diff_max;
	(void)t1;
	std::cerr<<"[E] This method is obsolete, it should not be used!!"<<std::endl;
}
 /* End of Diffusion */

void SnowDriftA3D::Suspension()
{
  const double eps=0.00005;
  double diff_max[2];
  diff_max[1] = 0.0;
  double t_diff=0.0;
  static int tmpindex=0;

  if (DRIFT_MARC)
    {
	int timeStep;
	const double maxTime = 0.8; // the time and upper limit of the time interval of integration
	const int maxTimeStep = (int) ceil( maxTime/dt_diff );

	// Inititalize system
	initializeSystem( colA, rowA, sA,sB,rhs,f,Psi,c,c00, CON);
	prepareSolve();
	timeStep=0;
	SolveEquation(timeStep, maxTimeStep,CON);
	values_elements_to_nodes(nodes_c,c);
	double nodes_min=nodes_c.grid3D.getMin();
	if (nodes_min<0){
		cout<<"!!!WARNING!!! at least one point with negative snow concentration!  (minimum c="<<nodes_min<<")"<<endl;
	}

	// compute deposition flux
	computeDepositionFlux( c, theta );

	if (SUBLIMATION && (nodes_c.grid3D.getMean()>0.00000)){
		std::cout<<"Start sublimation routine"<< std::endl;

		//now find the steady state sublimation
		iterativeSublimationCalculation(timeStep, maxTimeStep);

		//steady state sublimation is found, if the concentration feedback is switched off, the new concentration field will have to be calculated here to be able to consider the effect of sublimation on concentration and deposition
		if (C_FB== 0){
		    //recalculate concentration field
		    initializeSystem( colA, rowA,sA,sB, rhs, f, Psi, c,c00, CON);
		    //change source f to steady state sublimation
		    Grid3DObject tmpSink_c;
		    tmpSink_c=nodes_Subl_ini;
		    tmpSink_c.grid3D=nodes_Subl_ini.grid3D*(-1);
		    values_nodes_to_elements(tmpSink_c,f);

		    prepareSolve();
		    SolveEquation(0, maxTimeStep, CON);
		    values_elements_to_nodes(nodes_c,c);
		    nodes_min=nodes_c.grid3D.getMin();
		    if (nodes_min<0){
			cout<<"!!WARNING!! at least one point with negative concentration! (minimum c="<<nodes_min<<")"<<endl;
		    }
		}
		computeDepositionFluxSublimation( c, theta );
	}
    }//driftmarc
     else {
		do {
			t_diff += dt_diff;
			diff_max[0] = diff_max[1];
			/* Calculate the diffusion */
			Diffusion(dt_diff, diff_max[1], t_diff);
			printf("\n Local c: %f Diff Time: %f s diff_max:%f\n", nodes_c.grid3D(2,2,2),t_diff,diff_max[1]);

			tmpindex++;
			//DEBUG("SEQ FLUX (ix=31): x=%f, y=%f, z=%f",
			//	checksum(flux_x,31, nx-1),checksum(flux_y,31, nx-1),checksum(flux_z,31, nx-1));
			//DEBUG("SEQ FLUX (ix=33): x=%f, y=%f, z=%f",
			//	checksum(flux_x,33, nx-1),checksum(flux_y,33, nx-1),checksum(flux_z,33, nx-1));

		} while ( (diff_max[1] > eps) && (t_diff/dt_diff < 100.0) );

		printf("\n Converged after time %f\n", t_diff);
	}
	//DEBUG("DRIFT FLUX #%d: x=%f, y=%f, z=%f", tmpindex, checksum(flux_x), checksum(flux_y), checksum(flux_z));
}


/**
 * @brief Solve the advection diffusion equation
 * @param timeStep (parameter only for non stationary)
 * @param maxTimeStep (parameter only for non stationary)
 * @param param variable that should be solved
 */
void SnowDriftA3D::SolveEquation(int timeStep, int maxTimeStep,const param_type param)
{
	// Solve equation
	if ( !(STATIONARY) ) {

	    while ( timeStep < maxTimeStep ) {
		if ( timeStep == 0 ) {
		    //extract the initial conditions for time=0!
		    for (unsigned int kk = 0; kk<= nz-3; kk++) {
			// kick all the internal node layers!
			for (unsigned int jj = 0; jj <= ny-3; jj++ ) {
			    for (unsigned int ii = 1; ii <= nx-2; ii++ ) {
				c[kk*(nx-2)*(ny-2)+jj*(nx-2)+ii] = c00[(kk+1)*nx*ny+(jj+1)*nx+ii];
			    }
			}
		    }
		}

		if (param==HUM){
			matmult(rhs,q,sB,colA,rowA);
		} else if (param==CON){
			matmult(rhs,c,sB,colA,rowA);
		} else if (param==TEM){
			matmult(rhs,T,sB,colA,rowA);
		}

		//construct the rhs of the system to solve: Ac=rhs
		for (unsigned int i = 0; i < nDOF; i++) {
			rhs[i] += ( -Psi[i] );
		}

		//solve system
		double test;
		if (param==HUM){
			bicgStab(q,rhs,sA,colA,rowA,900,1e-5, test);
		} else if (param==CON){
			bicgStab(c,rhs,sA,colA,rowA,900,1e-5, test);
		} else if (param==TEM){
			bicgStab(T,rhs,sA,colA,rowA,900,1e-5, test);
		}
		//update time
		timeStep++;
	    }
	} else {
		//construct the right hand side rhs of the system A*c = rhs
		for (unsigned int i = 0; i < nDOF; i++ ) {
			rhs[i] += ( -Psi[i] );

			//LH_BC
			//  	      rhs[i] = 0.0005;
			//rhs[i] = res1[i] - Psi[i] + res2[i];//Here real line
		}

		double testres=10.;
		CDoubleArray saveRhs=rhs;
		CIntArray saveRowA=rowA;
		if (param==HUM){
		    std::cout<<"Solving equation for q with rhs=" <<rhs[18]<<std::endl;
		    bicgStab(q,rhs,sA,colA,rowA,3000,1e-11,testres);//call the biconjgrad method,ori 1e-10
		} else if (param==CON){
		    std::cout<<"Solving equation for c with rhs=" <<rhs[18]<<std::endl;
		    /*std::string fname = std::string("../results/") + calcDate.toString(Date::NUM)+ "0" + std::string(".c00");
			writefield(fname,HUM);*/
		    bicgStab(c,rhs,sA,colA,rowA,3000,1e-10,testres); //ori 1e-10, 1000
		} else if (param==TEM){
		    std::cout<<"Solving equation for T with rhs=" <<rhs[18]<<std::endl;
		    bicgStab(T,rhs,sA,colA,rowA,3000,1e-11,testres);
		}
	}
}

/**
@brief
Calculate the steady state sublimation in several steps.
Single feedbacks can be switched on or off in SnowDrift.h
*/
void SnowDriftA3D::iterativeSublimationCalculation(int timeStep, int maxTimeStep)
{
	double factorf=(-1); //fraction of initial sublimation that should be used to solve the equation. -1 by definition

	//Initial field of concentration has been calculated, now calculate the initial sublimation
	Sublimation();

	//initialize some parameters
	int count=0;             //count iterations
	double maxdif_s=0.;      //the maximum difference in sublimation between 2 iterations
	double maxdif_c=0.;      //the maximum difference in concentration between 2 iterations
	double test=0;           //controls whether iteration stops
	double mean_startsubl=0.;//mean sublimation at start iteration
	double mean_endsubl=0;   //mean final sublimation

	//initialize  init sublimation
	nodes_Subl_ini.grid3D=0.;

	//Save initial sublimation and adjust to start iteration with a value probably closer to steady-state sublimation
	if (C_FB==1 || Q_FB==1 || T_FB==1){
	    nodes_Subl_ini.grid3D=nodes_Subl.grid3D*0.8;
	}else{
	    nodes_Subl_ini=nodes_Subl;
	}
	mean_startsubl=nodes_Subl.grid3D.getMean();

	//now start the iteration
	while (count<8 && test==0 && mean_startsubl!=0.){ //Usually not more than 7 iterations are needed but can be adjusted. Test is defined at the end of the loop and will test changes in sublimation and concentration
		count=count+1;
		std::cout <<"Sublimation iteration number: "<< count << std::endl;

		//adapt initial sublimation to sublimation in next step (take the average)
 		if (count>1){
 			nodes_Subl_ini.grid3D = nodes_Subl_ini.grid3D+nodes_Subl.grid3D;
			nodes_Subl_ini.grid3D = nodes_Subl_ini.grid3D/2.;
 		}

		//save ini concentration to compare to next step
		values_elements_to_nodes(nodes_tmp_c,c);

		//Combine initial sublimation and initial humidity field and solve new steady state for humidity
		if (Q_FB==1){ //switch humidity feedback
			initializeSystem( colA, rowA,sA,sB,rhs,f,Psi,q,q00, HUM);
			Grid3DObject tmpSource_q;
			tmpSource_q.set(nx,ny,nz,ta.cellsize,ta.llcorner);
			 for (size_t ii=0; ii<(nx*ny*nz); ii++) {
				tmpSource_q(ii) = -1.*nodes_Subl_ini(ii) / (Atmosphere::waterVaporDensity(nodes_Tair(ii), Atmosphere::vaporSaturationPressure(nodes_Tair(ii)))+Atmosphere::stdDryAirDensity(nodes_z(ii), nodes_Tair(ii))*factorf);
			}
			values_nodes_to_elements(tmpSource_q,f);
			prepareSolve();
			SolveEquation(0, maxTimeStep, HUM);

			//recalculate relative humidity field
			values_elements_to_nodes(nodes_q, q);
			for (size_t ii=0; ii<(nx*ny*nz); ii++) {
				nodes_RH(ii) = RH_from_q(nodes_Tair(ii),nodes_q(ii), nodes_z(ii));
				if (nodes_RH(ii)>=1.5){
					nodes_q(ii) = Atmosphere::relToSpecHumidity(nodes_z(ii), nodes_Tair(ii),nodes_RH(ii));
				}
			}
			double tmpmin=nodes_q.grid3D.getMin();
			if (tmpmin<0){
				cout<<"WARNING, at least one point with negative specific humidity! (minimum q="<<tmpmin<<")"<<endl;
			}
		}

		if (C_FB== 1){
			//recalculate concentration field
			initializeSystem( colA, rowA,sA,sB, rhs, f, Psi, c,c00, CON);
			Grid3DObject tmpSink_c;
			tmpSink_c.set(nx,ny,nz,ta.cellsize,ta.llcorner);
			tmpSink_c.grid3D = nodes_Subl_ini.grid3D*factorf;
			values_nodes_to_elements(tmpSink_c,f);
			prepareSolve();
			timeStep=0;
			SolveEquation(timeStep, maxTimeStep, CON);
			values_elements_to_nodes(nodes_c,c);
			const double tmpmin = nodes_c.grid3D.getMin();
			if (tmpmin<0){
				std::cout<<"!!!WARNING!!! at least one point with negative snow concentration!  (minimum c="<<tmpmin*1e3<<" g/m3)"<<std::endl;
			}
		}

		//temperature feedback
		if (T_FB==1){
			initializeSystem( colA, rowA,sA,sB,rhs,f,Psi,T,T00, TEM);
			Grid3DObject tmpSink_T;
			tmpSink_T.set(nx,ny,nz,ta.cellsize,ta.llcorner);
			for (size_t ii=0; ii<(nx*ny*nz); ii++) {
				tmpSink_T(ii) = (nodes_Subl_ini(ii)*factorf) * Constants::lh_sublimation * 1/(Constants::specific_heat_air*(1+0.84*nodes_q(ii))*Atmosphere::stdDryAirDensity(nodes_z(ii), nodes_Tair(ii)));
			}
			values_nodes_to_elements(tmpSink_T,f);
			prepareSolve();
			timeStep=0;
			SolveEquation(timeStep, maxTimeStep, TEM);

			//get temperature from potential temperature and calculate relative humidity
			Grid3DObject tmp_potT;
			tmp_potT.set(nx,ny,nz,ta.cellsize,ta.llcorner);
			values_elements_to_nodes(tmp_potT,T);
			for (size_t ii=0; ii<(nx*ny*nz); ii++) {
				nodes_Tair(ii) = tmp_potT(ii)-(Cst::gravity/Constants::specific_heat_air)*nodes_z(ii);
				nodes_RH(ii) = RH_from_q(nodes_Tair(ii),nodes_q(ii), nodes_z(ii));
				if (nodes_RH(ii)>=1.5){
					nodes_q(ii)=Atmosphere::relToSpecHumidity(nodes_z(ii), nodes_Tair(ii),nodes_RH(ii));
				}
			}
			const double tmpmin = nodes_Tair.grid3D.getMin();
			if (tmpmin<0){
				std::cout<<"WARNING, at least one point with negative temperature! (minimum T="<<tmpmin<<")"<<std::endl;
			}
		}

		Sublimation();

		//find maximum difference and calculate mean final sublimation
		Grid3DObject temp_dif_subl;
		Grid3DObject temp_dif_c;
		temp_dif_subl.set(nx,ny,nz,ta.cellsize,ta.llcorner);
		temp_dif_c.set(nx,ny,nz,ta.cellsize,ta.llcorner);

		temp_dif_subl.grid3D=(nodes_Subl.grid3D-nodes_Subl_ini.grid3D);
		temp_dif_subl.grid3D.abs();

		temp_dif_c.grid3D=nodes_c.grid3D-nodes_tmp_c.grid3D;
		temp_dif_c.grid3D.abs();

		maxdif_s=temp_dif_subl.grid3D.getMax();
		maxdif_c=temp_dif_c.grid3D.getMax();
		mean_endsubl=nodes_Subl.grid3D.getMean();

		//test whether changes in c and sublimation are too large
		if (C_FB==1){
		    if (maxdif_s>1.2e-6 && maxdif_c>1.1e-6){
			test=0; //go on with iteration
		    } else {
			test=1; //found steady state, stop iteration
		    }
		} else {//test without concentration as this will not change when feedback is switched off
		    if (maxdif_s>1.e-6){
			test=0;
		    }else{
			test=1;
		    }
		}
		std::cout <<"Maximum difference in sublimation= " << maxdif_s*1e6 << "e-6, maximum difference in c= " << maxdif_c*1e3 << "e-3"<<std::endl;
	}//close while

	std::cout<<"Mean final sublimation=" <<mean_endsubl<<", fraction of initial:"<<mean_endsubl/(mean_startsubl+1e-20)<<std::endl;

	if (test==1){
	  std::cout <<"Found steady state sublimation after "<<count<<" iterations."<< std::endl;
	}else{
	  std::cout<<"!!!WARNING!!! No steady state found yet (maximum difference in sublimation = "<<maxdif_s<<"), max number of iterations reached, using sublimation field anyway"<<std::endl;
	}

}

/**
 * @brief prepareSolve
 * Some preparations to solve a 3D-field.
 */
void SnowDriftA3D::prepareSolve(){

    prepareSparseMatrix( colA,rowA,adjA);
    resetArray( sA );
    resetArray( sB );
    resetArray( rhs );
    assembleSystem( colA, rowA,sA,sB,Psi,f,dt_diff);
    //-------------------------------------------------------------------------------
    // Apply BC-----------NOT implemented!!!!!!!!!!!!!!!!!!!!!!!
    //-------------------------------------------------------------------------------
    //applyBoundaryValues(c00,Psi);
}


