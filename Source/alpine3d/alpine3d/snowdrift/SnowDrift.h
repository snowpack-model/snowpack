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
#ifndef SNOWDRIFTA3D_H
#define SNOWDRIFTA3D_H

#define WS0 0.5
#define TKE 0
#define SALTATION 1		// switch for saltation simulation
#define SUBLIMATION 0 		// switch for drifting snow sublimation
#define FIELD3D_OUTPUT 0	// output with all three-dimensional fields (only possible for sublimation)
#define SUBLIMATION_OUTPUT 0	// debug output of drifting snow sublimation
#define T_FB 1 //switch for feedback between sublimation and air temperature
#define Q_FB 1 //switch for feedback between sublimation and humidity
#define C_FB 1 //switch for feedback between sublimation and snow concentration
#define READK 0  //define as 1 if you have K from ARPS wind fields INCLUDING turbulence
#define WRITE_DRIFT_FLUXES 0 //set to 1 in order to write snow drift fluxes
#define dt_diff 0.5   /* Small calculation step length for snow diffusion */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <meteoio/MeteoIO.h>
#include <meteoio/plugins/ARPSIO.h>

typedef mio::Array2D<int> CElementArray;
typedef mio::Array1D<double> CDoubleArray;
typedef mio::Array1D<int> CIntArray;

#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/snowdrift/checksum.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"

class SnowpackInterface;
class EnergyBalance;

typedef enum DRIFT_OUTPUT_ENUM {OUT_CONC, OUT_SUBL} DRIFT_OUTPUT;
typedef enum PARAM_TYPES {CON,HUM,SUB,TEM,SUB2} param_type;
typedef enum ASPECT_TYPES {OTHER,BOTTOM} aspect_type;

struct WIND_FIELD {unsigned int start_step;std::string wind;};

/**
 * @page snowdrift Snowdrift
 * This module computes the preferential deposition and redistribution of snow by the wind (see \ref principles_snowdrift). It is
 * enabled using the "--enable-drift" command line option.
 * The 3D wind fields are read using a GRID3D <A HREF="https://models.slf.ch/p/meteoio">MeteoIO</A> plugin such as ARPS:
 * @code
 * GRID3D		= ARPS
 * GRID3DPATH	= ../input/wind_fields/
 * @endcode
 *
 * The WINDFIELDS key must be defined in the [Input] section and contains a space delimited list of wind fields files (here within GRID3DPATH) and
 * associated number of time steps:
 * @code
 * WINDFIELDS = sw3.asc 1 nw3.asc 3 ww0.asc 2 nw9.asc 5 nw6.asc 10 ww0.asc 5 sw3.asc 6 nw3.asc 1
 * @endcode
 *
 */
class SnowDriftA3D {
	public:
		SnowDriftA3D(const mio::DEMObject& dem, const mio::Config& cfg);

		virtual ~SnowDriftA3D();

		virtual void setSnowSurfaceData(const mio::Grid2DObject& cH_in, const mio::Grid2DObject& sp_in, const mio::Grid2DObject& rg_in,
		                                const mio::Grid2DObject& N3_in, const mio::Grid2DObject& rb_in);

		void setSnowPack(SnowpackInterface &mysnowpack);
		void setEnergyBalance(EnergyBalance &myeb);

		virtual void Compute(const mio::Date& calcDate);
		bool isNewWindField(const unsigned int current_step) /*const*/;
		void setMeteo (const unsigned int& steps, const mio::Grid2DObject& new_psum, const mio::Grid2DObject& new_psum_ph, const mio::Grid2DObject& new_p, const mio::Grid2DObject& new_vw,
		                     const mio::Grid2DObject& new_rh, const mio::Grid2DObject& new_ta, const mio::Grid2DObject& new_tsg, const mio::Grid2DObject& new_ilwr, const mio::Date& calcDate,
		                     const std::vector<mio::MeteoData>& vecMeteo);

		void GetTResults(double  outtime_v[15], double outtime_tau[15], double outtime_salt[15], double outtime_diff[15]);

		double getTiming() const;

		void Destroy();
		
		std::string getGridsRequirements() const;

	protected:
		void Initialize();
		void ConstructElements();
		void InitializeNodes(const mio::Grid3DObject& z_readMatr);
		void CompleteNodes();

		virtual void compSaltation(bool setbound);
		virtual void SnowMassChange(bool setbound, const mio::Date& calcDate);
		virtual void Suspension();
		virtual void Diffusion(double deltaT, double &diff_max, double t);


		//sublimation functions begin
		virtual void Sublimation();
		double terminalFallVelocity(const double Radius, const double Temperature, const double RH, const double altitude);
		double calcS(const double concentration, const double sublradius, const double dmdt);
		double calcSubldM(const double Radius, const double AirTemperature, const double RH,const double WindSpeed, const double altitude);
		double reynoldsNumberforFallingParticles(const double Radius, const double Windspeed, const double AirTemperature, const double RH, const double altitude);
		double ventilationVelocity(const double Radius, const double Windspeed,const double AirTemperature, const double RH, const double altitude);
		double waterVaporDensity(const double Temperature,const double VaporPressure);
		double RH_from_q(const double AirTemp, const double q, const double altitude);
		void initializeTRH();
		//sublimation functions end

		//---------------------------------------------------------------------
		//functions which are required for the fe numerics
		//defined in SnowDriftFENumerics.cc
		//---------------------------------------------------------------------
		//bare numerics for element computations, i.e. integrations, jacobian and auxiliary stuff
		virtual void Jacobian(double *DETERMINANTJ,double J[][3],const int element, const double *P, int k, const int ix, const int iy, const int iz);
		virtual void J0fun(double J0[3][3], const double J[3][3]);
		virtual double GQIntB(double *DETERMINANTJ,const int i, const int j);
		virtual double GQIntC(double * DETERMINANTJ, const double J0M[3][3][8], const int i, const int j, const double b[3],const double K[3][3]);
		virtual double GQIntApdx(double DETERMINANTJ[],const double J0M[3][3][8],const int i, const int j, double b[],const double deltak);
		virtual double GQIntAdxdx(double *DETERMINANTJ, double J0M[3][3][8],const int i, const int j, double *b,const double deltak);
		virtual void TTfun(double TT[3][8],const double P[]);	// P is a pointer pointing at an adress containing a double
		virtual void phi(double*PHI, double*P);
		virtual void setQuadraturePoints();

		//functions required for solving the linear system
		virtual void matmult(CDoubleArray& res, const CDoubleArray& x, double* sm, int* ijm);
		virtual void matmult(CDoubleArray& res, const CDoubleArray& x, const CDoubleArray& sA, const CIntArray& colA, CIntArray& rowA);
		virtual void transmult(CDoubleArray& res, const CDoubleArray& x,double* sm, int* ijm);
		virtual void SolveEquation(int timeStep, int maxTimeStep, const param_type param );
		virtual void bicgStab(CDoubleArray& result, CDoubleArray& rhs, const CDoubleArray& sA, const CIntArray& colA, CIntArray& rowA, const int nmax, const double tol, double& testres);


		//---------------------------------------------------------------------
		// finite Element functions, basically wrappers for the bare
		// numerics.  Basically two different "classes" of functions, the
		// first (assembleSystem, applyBoundaryValues, computeDepositionFlux)
		// contain loops over all or a particular subset of elements and then
		// invoke the other class of functions which wrap the numerics for
		// each element (computeDiffusionTensor, computeDriftVector,
		// computeElementParameter,
		// computeElementSystem,computeDirichletBoundaryValues,
		// addElementMatrix
		//
		// all these functions are defined in SnowDriftFEControl.cc
		//---------------------------------------------------------------------
		virtual void assembleSystem( CIntArray& colA, CIntArray& rowA, CDoubleArray& sA, CDoubleArray& sB, CDoubleArray& Psi, CDoubleArray& f, const double dt);
		virtual void applyBoundaryValues(CDoubleArray& c00, CDoubleArray& Psi);
		virtual void prepareSolve();


		//functions which are required for building the element matrices,
		virtual void computeDiffusionTensor(double K[3][3], const unsigned int ix, const unsigned int iy, const unsigned int iz);
		virtual void computeDriftVector(double b[3], const unsigned int ix, const unsigned int iy, const unsigned int iz );
		virtual void computeElementParameters(const int& element, double DETERMINANTJ[8], double J0M[3][3][8], double J0[3][3], double J[3][3], double b[3], double K[3][3], double& deltak, double& qualla, const int ix, const int iy, const int iz);

		virtual void computeElementSystem(int &element, int &nDofNodes, int* dofNode, double Ael[9][9], double Del[9][9], bool stationary, double DETERMINANTJ[8], double J0M[3][3][8], double b[3], double K[3][3], double &deltak, const double &dt, CDoubleArray& f, CDoubleArray& Psi);

		virtual void computeDirichletBoundaryValues(int element,double DETERMINANTJ[8],double J0M[3][3][8], double J0[3][3], double b[3], double K[3][3], double deltak, int spec[8], int length_spec, int length_complSpec, CDoubleArray& c00, CDoubleArray& Psi);

		virtual void addElementMatrix( CDoubleArray& sA, const CIntArray& colInd, const CIntArray& rowPtr,const double Bel[9][9], const int element, const int* spec,const int length_spec);

		virtual void computeDepositionFlux(const CDoubleArray& c, const double theta);
		virtual void computeDepositionFluxSublimation(const CDoubleArray& c, const double theta);


		//---------------------------------------------------------------------
		// Functions for initializing the finite element procedure
		// defined in SnowDriftFEInit.cc
		//---------------------------------------------------------------------
		void setBC_BottomLayer(CDoubleArray& var00,  const param_type param);
		void values_nodes_to_elements(const mio::Grid3DObject& nodesGrid, CDoubleArray& elementsArray );
		void values_elements_to_nodes(mio::Grid3DObject& nodesGrid, const CDoubleArray& elementsArray );
		void setRobinBoundaryCondition(const aspect_type aspect, const double gamma_val, const int ix, const int iy, const int iz, CDoubleArray& var00, const param_type param);
		virtual void InitializeFEData();
		virtual void prepareSparseMatrix( CIntArray& colA, CIntArray& rowA, CDoubleArray& adjA);
		virtual void initializeSystem( CIntArray& colA,CIntArray& rowA,CDoubleArray& sA,CDoubleArray& sB, CDoubleArray& rhs,CDoubleArray& f, CDoubleArray& Psi,CDoubleArray& var,CDoubleArray& var00, const param_type param);
		virtual void resetArray(CDoubleArray& sA);
		virtual void resetArray(CIntArray& sA);
		virtual void classifySubdomain();
		virtual int numberOfNonzeros();
		void iterativeSublimationCalculation(int timeStep, int maxTimeStep);

	protected:
		Saltation saltation_obj;
		
		//Time dependent data output after each computation step (SnowDrift::Compute)
		double  time_v[15], time_tau[15], time_salt[15], time_diff[15];
		int DOITERATION;
		double auxLayerHeight;
		double station_altitude;

		mio::IOManager io;
		SnowpackInterface *snowpack;
		EnergyBalance *eb;
		mio::Timer timer;

		mio::Grid2DObject cH;
		mio::Grid2DObject sp;
		mio::Grid2DObject rg;
		mio::Grid2DObject N3;
		mio::Grid2DObject rb;

		//the dimensions of the rectangular domain
		unsigned int nx, ny, nz;

		// Variables for the finite element solution
		unsigned int nDOF; //the total number of degrees of freedom
		unsigned int nNodes; //the total number of nodes
		unsigned int nElements; //the total number of elements
		unsigned int nNZ; //the total number of nonzero elements of the system matrix

		//create system matrix in compressed sparse row (CSR) storage format
		CDoubleArray sA;
		CDoubleArray sB;
		CIntArray colA;
		CIntArray rowA;

		//auxiliary vector for incorporating inhomogeneous Dirichlet
		//boundary conditions
		CDoubleArray Psi;

		//vector of nodal concentrations, solution of the linear system
		CDoubleArray c; //snow concentration in suspension
		CDoubleArray q; //specific humidity
		CDoubleArray T; //potential temperature

		//vector which contains the right hand side of the linear system
		CDoubleArray rhs;

		//vector of sink/source terms
		CDoubleArray f;

		//vector which contains boundary and initial conditions
		CDoubleArray c00;
		CDoubleArray q00;
		CDoubleArray T00;

		//vector which contains boundary and initial conditions
		CDoubleArray precond;
		//mio::Grid3DObject newElements_precond;
		//LH_BC
		CDoubleArray gNeumann;
		CDoubleArray gDirichlet;
		CDoubleArray gamma;
		void zeroRow(int node);

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//LH_DEBUG: Algorithm test BEGIN
		CIntArray nnzA;
		CDoubleArray adjA;
		//LH_DEBUG: Algorithm test END
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		double qualla;// paramater to vary the deltak coefficient of the SUPG approach

		double theta; // here you can vary the heigth above the point used in

		//the face,bar,corner interior arrays are used for the loops over
		//all elements and are set in classifySubdomain, see there for details
		mio::Array2D<int> nz_face;
		mio::Array2D<int> ny_face;
		mio::Array2D<int> nx_face;

		mio::Array2D<int> nz_bar;
		mio::Array2D<int> ny_bar;
		mio::Array2D<int> nx_bar;

		mio::Array2D<int> nz_interior;
		mio::Array2D<int> ny_interior;
		mio::Array2D<int> nx_interior;

		CIntArray n_corner;

		//for the quadrature points of the numerical integration scheme
		mio::Array2D<double> qPoint;

		//array for mapping the local node indices onto the global nodes
		//indices (actually this is the analog of the array elems of the old code)
		CElementArray nodeMap;

		//array for mapping the local node indices onto the global
		//degree-of-freedom (dof) indices
		CElementArray dofMap;

		CDoubleArray flux_x, flux_y, flux_z, flux_x_subl, flux_y_subl, flux_z_subl;

		bool new_wind_status;

		CElementArray elems;
		mio::Array2D<double> saltation, c_salt, mns_subl, mns_nosubl,dif_mns_subl;

		//Meteo 2D data
		mio::Grid2DObject mns, vw, rh, ta, tsg, p, psum, psum_ph/*, iswr, ea*/; // TODO ISWR activate, TODO EA activate

		//Meteo 1D data
		double ta_1D;
		
		mio::Date skip_date; //time step to skip because there would be no drift

		//3D
		mio::Grid3DObject nodes_x, nodes_y, nodes_z, nodes_sy, nodes_sx, nodes_slope, nodes_wstar;
		mio::Grid3DObject nodes_u, nodes_v, nodes_w, nodes_K;
		mio::Grid3DObject nodes_e, nodes_c;
		mio::Grid3DObject nodes_Tair,nodes_Tair_ini, nodes_q, nodes_q_ini, nodes_RH, nodes_Subl, nodes_Subl_ini, nodes_WindVel, nodes_tmp_c;
		bool STATIONARY;

	protected:
		void buildWindFieldsTable(const std::string& wind_field_string);
		std::vector<struct WIND_FIELD> wind_fields;
		int wind_field_index;

		void debugOutputs(const mio::Date& calcDate, const std::string& fname, const DRIFT_OUTPUT& filetype);
		void writeOutput(const std::string& fname); //HACK: this should be done by MeteoIO

		//sublimation constants
		static const double kinematicViscosityAir, USTAR, molecularWeightofWater, thermalConductivityofAtm;

		// constants originally from Snowpack
		static const double c_red, grain_size, tau_thresh, z0;
		static const bool thresh_snow;
};

#pragma GCC diagnostic pop

#endif





