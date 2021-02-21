
#include <alpine3d/ebalance/SolarPanel.h>
#include <alpine3d/ebalance/ViewFactorsHelbig.h> //for testing and comparison


//file operations
#include <iostream>
#include <fstream>
#include <stdlib.h>

//for string conversion
#include <string>
#include <stdio.h>
#include <unistd.h>

using namespace mio;


SolarPanel::SolarPanel(const mio::Config& cfg, const mio::DEMObject &dem_in, const std::vector<std::vector<double> > &pv_pts):
						dimx(dem_in.getNx()), dimy(dem_in.getNy()), dem(dem_in), pv_points(pv_pts), BRDFobject(cfg), timer(), Shadelist()
{
	// PRECISION PARAMETERS
	// ########################
	M_epsilon_panel=30;				// Number of small circles in Basic Set of Panels [MT fig. 2.1]
	M_phi_panel=30;					// Number of vectors per small circle of Basic Set [MT fig. 2.1]
	// ########################

	S_panel=M_epsilon_panel*M_phi_panel;	// Number of vectors per Basic Set [MT fig. 2.1]

	// ENABLE HERE ALL KIND OF OUTPUT
	// ##############################################################################################################################
	generate_PVP_sum=false; 					// Write Averaged Radiation to files. ( Make sure you have folder ../output/PVP )
	sun_tracker=false;							// Writes File with irradiation of Tracker Following the Sun [MT 3.4.3 Tracker]
	optimal_tracker=false;						// Writes File with irradiation of Tracker searching the optimum all the time (very expensive!) [MT 3.4.3 Tracker]

	// readSumPVP();							// Reads the 3 sum-files (direct.sum, diffuse.sum, terrain.sum; must be placed in ../input/PVP_SUM
	//
	// ->   this enables other output (the following are examples)
	// PrintVector(optimize(get_ii(x),get_jj(x),pv_points[x][2],100,&SolarPanel::minfun_MonoStatic));	// Prints, for PVP nr. x+1, optimal angles and corresponding radiation (for sum-period)
	// RadiationMap(get_ii(x),get_jj(x),pv_points[x][2]);		// Write file with Irradiation for all directions of the basic set, for PVP nr. x+1,
	// GridRadiationMap(3)										// Write File with optimal angles and corresponding radiation for whole DEM (panels @ a height of 3m)
	// ##############################################################################################################################

	// Abort if TerrainRadiationComplex is required but not activated
	std::string method=cfg.get("Terrain_Radiation_Method", "EBalance");
	if (method!="COMPLEX" && (generate_PVP_sum || sun_tracker || optimal_tracker) ){
		std::cout<<"[E] In SolarPanel: TerrainRadiationComplex is required for <generate_PVP_sum>, <sun_tracker> or <optimal_tracker>.\n";
		abort();
	}

	// Check for required Key
	if (cfg.keyExists("PV_shadowing", "Ebalance")) cfg.getValue("PV_shadowing", "Ebalance", if_shadowing);
		else std::cout<<"[i] In SolarPanel: No flag <<PV_shadowing>> set in [Ebalance]. Use default PV_shadowing = "<<if_shadowing<<" .\n";

	// Initialize the embedding of the panels in the terrain
	initBasicSetHorizontal();
	initBasicSetRotated();
	initViewListPanel();
	if (if_shadowing) initShadelist();
	writeHeader();

	std::cout<<"[i] SolarPanel object created\n";
}



/**
* @brief Called by TerrainRadiationComplex, updates Terrain Lists. A Terrain List (TList) stores Radiance in all S directions for whole DEM.
*        TList_ms is full terrain radiance while TList_direct is for shading, TList_sky's are for identifying anisotropy and multiple scattering effect
* @param[in] TList1-TList4
* @param[out] -
*
*/
void SolarPanel::setTLists(mio::Array4D<double> TList1, mio::Array4D<double> TList2, mio::Array4D<double> TList3, mio::Array4D<double> TList4){

	TList_direct=TList1;
	TList_sky_iso=TList2;
	TList_sky_aniso=TList3;
	TList_ms=TList4;
}



/**
* @brief Initializing interface between TerrainRadiationComplex and SolarPanel.
* @param[in] TList1-TList4
* @param[out] -
*
*/
void SolarPanel::initTerrain(size_t M_epsilon_terrain_in, size_t M_phi_terrain_in){

	M_epsilon_terrain=M_epsilon_terrain_in;
	M_phi_terrain=M_phi_terrain_in;
	S_terrain=M_epsilon_terrain*M_phi_terrain;

	Terrain_complex_mode=true;

	initViewListTerrain();
	if (generate_PVP_sum) initSumPVP();

}


/**
* @brief Updates grid Radiation (incoming SWR for all grid points): preparing TerrainRadiation if TerrainRadiationAlg!=COMPLEX
* @param[in] grid-albedo
* @param[in] grid-direct SWR
* @param[in] grid-diffuse SWR
* @param[in] grid-direct-horizontal-unshaded (for projection on PVP/triangles)
* @param[out] -
*
*/
void SolarPanel::setGridRadiation(const mio::Array2D<double>& in_albedo, const mio::Array2D<double>& in_direct, const mio::Array2D<double>& in_diffuse, const mio::Array2D<double>& in_direct_unshaded_horizontal, const double solarAzimuth, const double solarElevation){

	albedo = in_albedo;

	d_diffuse = in_diffuse;
	d_direct_unshaded_horizontal=in_direct_unshaded_horizontal;
	d_direct_A = in_direct;
	d_direct_B = in_direct;

	v_globalsun=getVectorSun(solarAzimuth, solarElevation);

	// Calculate Direct Radiation for all triangular pixels [MT Figure 2.2]. If TerrainRadiationAlg==COMPLEX, this is already done
	if (Terrain_complex_mode==false)
	{
		#pragma omp parallel for

		for (size_t ii=1; ii<dem.getNx()-1; ii++) {

			for (size_t jj=1; jj<dem.getNy()-1; jj++) {

				// Triangle Type A; Geometric projection: horizontal radiation -> beam radiation -> triangle radiation
				if (v_globalsun[2]>0 && VectorScalarProduct(TriangleNormal(ii,jj,1),v_globalsun)>0 && in_direct(ii,jj)!=0) // Criteria: 1) Day, 2) Self-shading, 3) horizon-shading
				{
					double proj_to_ray, proj_to_triangle;
					std::vector<double> triangle_normal=TriangleNormal(ii,jj,1);

					proj_to_ray=1./VectorScalarProduct(v_globalsun, {0,0,1});
					proj_to_triangle=VectorScalarProduct(v_globalsun, triangle_normal);

					d_direct_A(ii,jj)=in_direct_unshaded_horizontal(ii,jj)*proj_to_ray*proj_to_triangle;
				}
				else d_direct_A(ii,jj)=0;

				// Triangle Type B; Geometric projection: horizontal radiation -> beam radiation -> triangle radiation
				if (v_globalsun[2]>0 && VectorScalarProduct(TriangleNormal(ii,jj,0),v_globalsun)>0 && in_direct(ii,jj)!=0)
				{
					double proj_to_ray, proj_to_triangle;
					std::vector<double> triangle_normal=TriangleNormal(ii,jj,0);

					proj_to_ray=1./VectorScalarProduct(v_globalsun, {0,0,1});
					proj_to_triangle=VectorScalarProduct(v_globalsun, triangle_normal);

					d_direct_B(ii,jj)=in_direct_unshaded_horizontal(ii,jj)*proj_to_ray*proj_to_triangle;
				}
				else d_direct_B(ii,jj)=0;
			}
		}
	}
}

/**
* @brief Writes output for SolarPanels (PVP files), and updates the sum if generate_PVP_sum==true
* @param[in] Timeobject for timestamp
* @param[in] grid-albedo
* @param[out] -
*
*/
void SolarPanel::setSP(const mio::Date timestamp, const double solarAzimuth, const double solarElevation){
	std::cout << "Key1";
	getRadfield(); // Trigger Calculation of Iradiation on Solar Panels

	std::vector<std::ofstream> PVP_files(pv_points.size());
	std::string path = "../output/PVP/";

	// Open Static Panel Files
 	for (size_t i = 0; i < pv_points.size(); ++i)
	{
		std::string filename=path;
		filename.append(std::to_string(i+1)).append(".pvp");

		PVP_files[i].open(filename, std::ios_base::app | std::ios_base::out);

		if (!PVP_files[i]) {
		    std::cout << "Exception opening PVP file";
		}
	}

	// Write Static Panel Files
	for (size_t i = 0; i < pv_points.size(); ++i)
	{
		PVP_files[i]<<std::fixed <<std::setprecision(6) <<timestamp.toString(Date::ISO)<<","<<timestamp.getJulian()<<","<<solarAzimuth<<","<<solarElevation<<","<<direct[i]<<","<<diffuse[i];
		PVP_files[i]<<std::fixed <<std::setprecision(6) <<","<<terrain_iso[i]<<","<<terrain_aniso[i]<<","<<terrain_ms[i]<<","<<terrain_ms_noshadow[i]<<","<<albedo(get_ii(i),get_jj(i))<<","<<direct_beam[i]<<"\n";
		PVP_files[i].close();
	}

	// Write Tracker Files
	if(sun_tracker) WriteOptimumTrackerRadiation(0, "OptimumTracker.pvp", timestamp);
	if(sun_tracker) WriteSunTrackerRadiation(0, "SunTracker.pvp", timestamp);

	// Update sum [MT 3.4.4 Radiation Maps], [MT 2.2.4 Solar Panels, fig. 2.16]
	if (generate_PVP_sum && v_globalsun[2]>=0) // no accumulation during the night
	{
		TList_sum+=TList_ms;

		#pragma omp parallel for
		for (size_t ii = 0; ii < dimx; ++ii)
		{
			for (size_t jj = 0; jj < dimy; ++jj)
			{
				Diffuse_sum(ii, jj) += d_diffuse(ii,jj);

				double proj_to_ray;
				size_t solidangle_sun_flat;
				proj_to_ray=1./VectorScalarProduct(v_globalsun, {0,0,1});
				solidangle_sun_flat=vectorToSPixel(v_globalsun, M_epsilon_panel, M_phi_panel);

				Direct_sum(ii, jj, solidangle_sun_flat) += d_direct_unshaded_horizontal(ii,jj)*proj_to_ray;
			}
		}
	}

	alpine_steps++;
	std::cout<<"[i] SolarPanel updated - step "<<alpine_steps<<" \n";

}


/**
* @brief Updates Radiation (diffuse, direct, Terrain) for Static Solarpanels [MT 2.2.4 Solar Panels]
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::getRadfield(){

	direct.clear();					// Direct Incident Radiation
	diffuse.clear();				// Diffuse Radiation
	direct_beam.clear();			// Beam Radiation
	terrain_iso.clear();			// Single-Scattered terrain Radiation (Isotropic)
	terrain_aniso.clear();			// Single-Scattered terrain Radiation (Anisotropic)
	terrain_ms.clear();				// Multiple-Scattered terrain Radiation
	terrain_ms_noshadow.clear();	// Multiple-Scattered terrain Radiation (if there were no panel induced shadows)

	bool day=true;

	if(v_globalsun[2]<=0) day=false; // <0 : not visible;  ==0 : impossible back-projection

	for (size_t number_pvp = 0; number_pvp < ViewList_panel.size(); ++number_pvp)
	{
		double dir_horinzontal;
		std::vector<double> PVP_normal;
		int ii_PVP =get_ii(number_pvp);
		int jj_PVP =get_jj(number_pvp);

		PVP_normal=RotationAnglesToNormalVector(pv_points[number_pvp][4], pv_points[number_pvp][3]);

		/////////////////////////////
		///// DIRECT [MT eq. 2.103]

		dir_horinzontal=d_direct_unshaded_horizontal(ii_PVP,jj_PVP);

		if (day && VectorScalarProduct(PVP_normal,v_globalsun) > 0)
		{

			double proj_to_ray, proj_to_PVP;
			size_t solidangle_sun;

			proj_to_ray=1./VectorScalarProduct(v_globalsun, {0,0,1});
			proj_to_PVP=VectorScalarProduct(v_globalsun, PVP_normal);
			solidangle_sun=vectorToSPixel(v_globalsun, number_pvp, M_epsilon_panel, M_phi_panel);
			direct_beam.push_back(dir_horinzontal*proj_to_ray);

			if (ViewList_panel[number_pvp][solidangle_sun][0]==-999) direct.push_back(dir_horinzontal*proj_to_ray*proj_to_PVP);
			else direct.push_back(0);
		}
		else{
			direct.push_back(0);
			direct_beam.push_back(0);
		}


		/////////////////////////////
		///// DIRECT [MT eq. 2.105]

		double total_terrain_iso=0, total_terrain_aniso=0, total_terrain_ms=0, total_terrain_ms_noshadow=0;
		size_t solidangle_sun_flat=0;
		if (day) solidangle_sun_flat=vectorToSPixel(v_globalsun, M_epsilon_panel, M_phi_panel); // No day, no shading

		if (Terrain_complex_mode==false)
		{

			int ii_dem, jj_dem, which_triangle;
			std::vector<double> normal_pixel,ray_out;
			double diffuse_solidangle, direct_solidangle, skyview_solidangle;
			double cth_i, cth_v, cphi,R;

			for (size_t number_solidangle = 0; number_solidangle < ViewList_panel[number_pvp].size(); ++number_solidangle)
			{
				if (ViewList_panel[number_pvp][number_solidangle][0]==-999) continue;
				ii_dem=ViewList_panel[number_pvp][number_solidangle][1];
				jj_dem=ViewList_panel[number_pvp][number_solidangle][2];
				which_triangle=ViewList_panel[number_pvp][number_solidangle][3];

				skyview_solidangle=mio::DEMAlgorithms::getCellSkyViewFactor(dem, ii_dem, jj_dem);

				if (which_triangle==1) direct_solidangle=d_direct_A(ii_dem, jj_dem);
				else direct_solidangle=d_direct_B(ii_dem, jj_dem);

				if (if_shadowing && Shadelist(number_pvp, solidangle_sun_flat, number_solidangle)) direct_solidangle=0; //Self-shadowing
				diffuse_solidangle=d_diffuse(ii_dem,jj_dem)*skyview_solidangle;

				ray_out=VectorStretch(BasicSet_rotated[number_pvp][number_solidangle],-1);
				normal_pixel=TriangleNormal(ii_dem, jj_dem, which_triangle);

				cth_i=VectorScalarProduct(v_globalsun,normal_pixel);
				cth_v=VectorScalarProduct(normal_pixel,ray_out);
				cphi=VectorScalarProduct(normalizeVector(ProjectVectorToPlane(v_globalsun,normal_pixel)),normalizeVector(ProjectVectorToPlane(ray_out,normal_pixel)));
				R=1;

				if (albedo(ii_dem, jj_dem)<0) albedo(ii_dem, jj_dem)=0.21; //For Non-Data
				total_terrain_iso=total_terrain_iso+(direct_solidangle*R+diffuse_solidangle)*albedo(ii_dem, jj_dem)/(S_panel);
				if (albedo(ii_dem, jj_dem)>0.5) R=BRDFobject.get_RF(cth_i, cphi, cth_v);
				total_terrain_aniso=total_terrain_aniso+(direct_solidangle*R+diffuse_solidangle)*albedo(ii_dem, jj_dem)/(S_panel);

			}
		}

		if (Terrain_complex_mode==true)
		{

			for (size_t number_solidangle = 0; number_solidangle < ViewList_panel[number_pvp].size(); ++number_solidangle)
			{
				if (ViewList_panel[number_pvp][number_solidangle][0]==-999) continue;

				size_t ii_source=ViewList_panel[number_pvp][number_solidangle][1];
				size_t jj_source=ViewList_panel[number_pvp][number_solidangle][2];
				size_t which_triangle_source=ViewList_panel[number_pvp][number_solidangle][3];
				size_t solidangle_source=ViewList_panel[number_pvp][number_solidangle][4];

				total_terrain_iso+=TList_sky_iso(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
				total_terrain_aniso+=TList_sky_aniso(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
				total_terrain_ms+=TList_ms(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
				total_terrain_ms_noshadow+=TList_ms(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;

				if (if_shadowing && Shadelist(number_pvp, solidangle_sun_flat, number_solidangle))
				{
					total_terrain_iso-=TList_direct(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
					total_terrain_aniso-=TList_direct(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
					total_terrain_ms-=TList_direct(ii_source, jj_source, which_triangle_source, solidangle_source)/ViewList_panel[number_pvp].size()*Cst::PI;
				}
			}
		}

		terrain_iso.push_back(total_terrain_iso);
		terrain_aniso.push_back(total_terrain_aniso);
		terrain_ms.push_back(total_terrain_ms);
		terrain_ms_noshadow.push_back(total_terrain_ms_noshadow);


		/////////////////////////////
		///// DIRECT [MT eq. 2.104]
		diffuse.push_back(d_diffuse(ii_PVP, jj_PVP)*getSkyViewFactor(number_pvp));


	}
}


/**
* @brief Writes unprojected average radiation to files
* @param[in] max_steps timesteps of simulation for mean
* @param[out] -
*
*/
void SolarPanel::writeSP(const unsigned int max_steps)
{
	if (generate_PVP_sum){
		std::ofstream F_DIFF, F_DIR, F_TER;

		F_DIFF.open("../output/PVP/diffuse.sum");
		F_DIFF<<"[HEADER]\n";
		F_DIFF<<"Hemispherical Resolution (NxM) : \n"<< M_epsilon_panel<<"\t"<<M_phi_panel<<"\n[DATA]\n";
		F_DIFF<<"ii\tjj\tRadiation\n";

		for (size_t ii = 0; ii < dimx; ++ii)
		{
			for (size_t jj = 0; jj < dimy; ++jj)
			{
				F_DIFF<<ii<<"\t"<<jj<<"\t"<<Diffuse_sum(ii, jj)/(double)max_steps<<"\n";
			}
		}
		F_DIFF.close();

		F_DIR.open("../output/PVP/direct.sum");
		F_DIR<<"[HEADER]\n";
		F_DIR<<"Hemispherical Resolution (NxM) : \n"<< M_epsilon_panel<<"\t"<<M_phi_panel<<"\n[DATA]\n";
		F_DIR<<"ii\tjj\tsolidangle\tRadiation\n";

		for (size_t ii = 0; ii < dimx; ++ii)
		{
			for (size_t jj = 0; jj < dimy; ++jj)
			{
				for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
				{
					F_DIR<<ii<<"\t"<<jj<<"\t"<<number_solidangle<<"\t"<<Direct_sum(ii, jj, number_solidangle)/(double)max_steps<<"\n";
				}
			}
		}
		F_DIR.close();

		F_TER.open("../output/PVP/terrain.sum");
		F_TER<<"[HEADER]\n";
		F_TER<<"Hemispherical Resolution (NxM) : \n"<< M_epsilon_terrain<<"\t"<<M_phi_terrain<<"\n[DATA]\n";
		F_TER<<"ii\tjj\twhichtriangle\tsolidangle\tRadiation\n";

		for (size_t ii = 0; ii < dimx; ++ii)
		{
			for (size_t jj = 0; jj < dimy; ++jj)
			{
				for (size_t which_triangle = 0; which_triangle < 2; ++which_triangle)
				{
					for (size_t number_solidangle = 0; number_solidangle < S_terrain; ++number_solidangle)
					{
						F_TER<<ii<<"\t"<<jj<<"\t"<<which_triangle<<"\t"<<number_solidangle<<"\t"<<TList_sum(ii, jj,which_triangle, number_solidangle)/(double)max_steps<<"\n";
					}
				}
			}
		}
		F_TER.close();
	}
}


//########################################################################################################################
//                                             INITIALISATION FUNCTIONS
//########################################################################################################################


/**
* @brief Write Header for PVP-output
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::writeHeader(){


	std::vector<std::ofstream> PVP_files(pv_points.size());
	std::string path = "../output/PVP/";


	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp)
	{
		std::string filename=path;
		filename.append(std::to_string(number_pvp+1)).append(".pvp");

		PVP_files[number_pvp].open(filename);

		if (!PVP_files[number_pvp]) {
		    std::cout << "Exception opening PVP file";
		}
	}

	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp)
	{
		PVP_files[number_pvp]<<"[HEADER]\n";
		PVP_files[number_pvp]<<"Station_id = PVP_"<< number_pvp+1 << "\n";
		PVP_files[number_pvp]<<"Easting: "<< dem.llcorner.getEasting()+dem.cellsize*get_ii(number_pvp) << "\n";
		PVP_files[number_pvp]<<"Northing: "<< dem.llcorner.getNorthing()+dem.cellsize*get_jj(number_pvp) << "\n";
		PVP_files[number_pvp]<<"Height over ground: "<< pv_points[number_pvp][2] << "m\n";
		PVP_files[number_pvp]<<"Inclination: "<< pv_points[number_pvp][3] << "\n";
		PVP_files[number_pvp]<<"Azimuth: "<< pv_points[number_pvp][4] << "\n";
		PVP_files[number_pvp]<<"Height: "<< pv_points[number_pvp][5] << "m\n";
		PVP_files[number_pvp]<<"Width: "<< pv_points[number_pvp][6] << "m\n";
		PVP_files[number_pvp]<<"LandviewFactor = "<<getLandViewFactor(number_pvp)<<"\n";
		PVP_files[number_pvp]<<"SkyviewFactor = "<<getSkyViewFactor(number_pvp)<<"\n";
		PVP_files[number_pvp]<<"Hemispherical Resolution :  "<< M_epsilon_panel<<"x"<<M_phi_panel<<" = "<< S_panel<<"\n";
		PVP_files[number_pvp]<<"[DATA]\n";
		PVP_files[number_pvp]<<"timestamp, Julian Date, solarAzimuth, solarElevation, ISR direct, ISR sky-diffuse, ISWR terrain (ISO), ISWR terrain (ANISO), ISWR terrain (MS), ISWR terrain (MS-NOSHADOW),  albedo, direct beam\n";
		PVP_files[number_pvp].close();
	}

	if(sun_tracker){
		std::ofstream File;
		std::string filename=path;
		filename.append("SunTracker.pvp");
		File.open(filename);

		File<<"[HEADER]\n";
		File<<"Station_id = SunTracker\n";
		File<<"Hemispherical Resolution :  "<< M_epsilon_panel<<"x"<<M_phi_panel<<" = "<< S_panel<<"\n";
		File<<"[DATA]\n";
		File<<"timestamp, Julian Date, OtimalAzimuth, OptimalElevation, ISR direct, ISR sky-diffuse, ISWR terrain, ISWR global,  albedo, direct beam\n";
		File.close();
	}
	if(optimal_tracker){
		std::ofstream File;
		std::string filename=path;
		filename.append("OptimumTracker.pvp");
		File.open(filename);

		File<<"[HEADER]\n";
		File<<"Station_id = OptimumTracker\n";
		File<<"Hemispherical Resolution :  "<< M_epsilon_panel<<"x"<<M_phi_panel<<" = "<< S_panel<<"\n";
		File<<"[DATA]\n";
		File<<"timestamp, Julian Date, solarAzimuth, solarElevation, ISR direct, ISR sky-diffuse, ISWR terrain, ISWR global,  albedo, direct beam\n";
		File.close();
	}
}

/**
* @brief  Initializes set of Vectors that point to equal solid angels for horizontal hemisphere (BasicSet) [MT 2.1.1 Basic Set]
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initBasicSetHorizontal(){

	double psi_0=0,psi_1=0,epsilon=0;
	double delta=0;
	double phi=0;
	int i=0;

	BasicSet_horizontal.clear();

	for (size_t n = 0; n < M_epsilon_panel; ++n)
	{
		psi_0=psi_1;
		if (n==(M_epsilon_panel-1))
		{
			psi_1=3.1415926/2;
		}else{
		delta=acos(-2/float(M_epsilon_panel)+cos(2*psi_0))/2-psi_0;											// [~ MT eq. 2.5]
		psi_1=psi_0+delta;
		}
		epsilon=(psi_1+psi_0)/2;																			// [MT eq. 2.4]

		for (size_t m = 0; m < M_phi_panel; ++m)
		{
			phi=2*3.1415926*(m+0.5)/M_phi_panel;															// [MT eq. 2.3]
			i+=1;
			BasicSet_horizontal.push_back({cos(epsilon)*sin(phi), cos(epsilon)*cos(phi), sin(epsilon)}); 	// [MT eq. 2.1]
		}
	}
}



/**
* @brief Rotates set of Vectors for horizontal hemisphere in Plane of PVP [MT 2.1.3 View-List]
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initBasicSetRotated(){

	double inclination_panel,azimuth_panel;

	std::vector<std::vector<double> > SVector_temp;

	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		SVector_temp.push_back({0,0,0});
	}

	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp)
	{

		inclination_panel=pv_points[number_pvp][3]*Cst::to_rad;
		azimuth_panel=pv_points[number_pvp][4]*Cst::to_rad;

		#pragma omp parallel for
		for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
		{
			std::vector<double> v_rotated=BasicSet_horizontal[number_solidangle];
			std::vector<double> normal_panel=RotationAnglesToNormalVector(azimuth_panel*Cst::to_deg, inclination_panel*Cst::to_deg);
			std::vector<double> axis=VectorCrossProduct({0,0,1}, normal_panel); // [MT eq. 2.39]
			if (NormOfVector(axis)==0) axis={1,0,0};							// [MT see text after eq. 2.41]
			if ( ( ( (int)(inclination_panel*Cst::to_deg) % 360)+360) % 360 > 180) axis=VectorStretch(axis, -1); // prevent flip of axis

			v_rotated=RotN(axis, v_rotated, inclination_panel);					// [MT eq. 2.38]

			SVector_temp[number_solidangle]=v_rotated;
		}

		BasicSet_rotated.push_back(SVector_temp);
	}
}


/**
* @brief Assigns a pixel (or sky) to each space vector of each PVP very similar to TerrainRadiationComplex::initViewList
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initViewListPanel(){

	std::vector<std::vector<double> > SList_temp;
	for (size_t  number_solidangle = 0; number_solidangle < S_panel; ++ number_solidangle) // loop over solid angles
	{
		SList_temp.push_back({0,0,0,0});
	}
	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp) // loop over PVP
	{

		size_t ii_PVP= get_ii(number_pvp);
		size_t jj_PVP= get_jj(number_pvp);
		double offset_PVP=pv_points[number_pvp][2];

		#pragma omp parallel for
		for (size_t  number_solidangle = 0; number_solidangle < BasicSet_rotated[number_pvp].size(); ++ number_solidangle) // loop over solid angles
		{
			size_t ii_temp=999999;
			size_t jj_temp=999999;
			size_t which_triangle_temp=9;
			double distance, minimal_distance=dem.cellsize*(dem.getNx()+dem.getNy());

			size_t ii_dem=ii_PVP, jj_dem=jj_PVP, nb_cells=0;
			std::vector<double> projected_ray=ProjectVectorToPlane(BasicSet_rotated[number_pvp][number_solidangle],{0,0,1});
			if (NormOfVector(projected_ray) != 0 ) projected_ray=normalizeVector(projected_ray);
			else projected_ray={1,0,0}; // if it goes straight up there will be sky (no caves for 2D DEM)
			while ( !(ii_dem<1 || ii_dem>dem.getNx()-2 || jj_dem<1 || jj_dem>dem.getNy()-2) ) {

				for (int k = -1; k < 2; ++k)
				{
					for (int kk = -1; kk < 2; ++kk)
					{

						ii_dem = ii_PVP + (int)round( ((double)nb_cells)*projected_ray[0] ) + k;
						jj_dem = jj_PVP + (int)round( ((double)nb_cells)*projected_ray[1] ) + kk;
						if((ii_dem<1 || ii_dem>dem.getNx()-2 || jj_dem<1 || jj_dem>dem.getNy()-2)) continue;

						distance=IntersectionRayTriangle(BasicSet_rotated[number_pvp][number_solidangle],ii_PVP,jj_PVP,offset_PVP,ii_dem,jj_dem,0);
						if (distance!=-999 && distance<minimal_distance)
						{
							minimal_distance=distance;
							ii_temp=ii_dem;
							jj_temp=jj_dem;
							which_triangle_temp=0;
						}
						distance=IntersectionRayTriangle(BasicSet_rotated[number_pvp][number_solidangle],ii_PVP,jj_PVP,offset_PVP,ii_dem,jj_dem,1);
						if (distance!=-999 && distance<minimal_distance)
						{
							minimal_distance=distance;
							ii_temp=ii_dem;
							jj_temp=jj_dem;
							which_triangle_temp=1;
						}

					}
				}
				ii_dem = ii_PVP + (int)round( ((double)nb_cells)*projected_ray[0] );
				jj_dem = jj_PVP + (int)round( ((double)nb_cells)*projected_ray[1] );
				nb_cells+=2.7;
			}

			if (minimal_distance==dem.cellsize*(dem.getNx()+dem.getNy())) minimal_distance=-999;
			SList_temp[number_solidangle]={minimal_distance, (double)ii_temp, (double)jj_temp, (double)which_triangle_temp};
		}

		ViewList_panel.push_back(SList_temp);
	}
}

/**
* @brief For terraincomplex-Mode: Assigns index of back-staring vector
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initViewListTerrain(){

	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp) // loop over PVP
	{
		for (size_t  number_solidangle = 0; number_solidangle < BasicSet_rotated[number_pvp].size(); ++ number_solidangle) // loop over solid angles
		{
			if (ViewList_panel[number_pvp][number_solidangle][0]==-999){
				ViewList_panel[number_pvp][number_solidangle].push_back(99999);
			}
			else{
				std::vector<double> v_out=VectorStretch(BasicSet_rotated[number_pvp][number_solidangle], -1);
				size_t ii_source=ViewList_panel[number_pvp][number_solidangle][1];
				size_t jj_source=ViewList_panel[number_pvp][number_solidangle][2];
				size_t which_triangle_source=ViewList_panel[number_pvp][number_solidangle][3];

				size_t solidangle_source=vectorToSPixel(v_out, ii_source, jj_source, which_triangle_source, M_epsilon_terrain, M_phi_terrain);

				ViewList_panel[number_pvp][number_solidangle].push_back(solidangle_source);
			}
		}
	}
}


/**
* @brief initialize shadows based on BasicSet of sun positions
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initShadelist()
{
	Shadelist.resize(pv_points.size(), S_panel, S_panel, 0);

	for (size_t number_pvp = 0; number_pvp < pv_points.size(); ++number_pvp) // loop over PVP
	{
		#pragma omp parallel for
		for (size_t  solidangle_sun = 0; solidangle_sun < S_panel; ++ solidangle_sun) // loop over solid angles
		{
			std::vector<double> v_sun=BasicSet_horizontal[solidangle_sun];

			for (size_t  solidangle_out = 0; solidangle_out < S_panel; ++ solidangle_out) // loop over solid angles
			{
				Shadelist(number_pvp, solidangle_sun, solidangle_out)=doesPanelShadowPixel(v_sun, number_pvp, solidangle_out);
			}
		}
	}
}


//########################################################################################################################
//                                                 SUM & TRACKER FUNCTIONS
//########################################################################################################################


/**
* @brief Initialisation of Sum-objects
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::initSumPVP(){

	TList_sum.resize(dimx, dimy, 2, S_terrain, 0);
	Direct_sum.resize(dimx, dimy, S_panel, 0);
	Diffuse_sum.resize(dimx, dimy, 0);

}


/**
* @brief Read Sum from files
* @param[in] -
* @param[out] -
*
*/
void SolarPanel::readSumPVP(){

	generate_PVP_sum=false; // Or else it would add to these numbers

	// DIFFUSE
	std::ifstream F_DIFF ("../input/PVP_SUM/diffuse.sum");
    std::string garbage;

    for (int i = 0; i < 5; ++i) getline(F_DIFF, garbage);

    Diffuse_sum.resize(dimx, dimy, 0);

    size_t ii,jj, number_solidangle;
	double diffuse_temp;

    while (F_DIFF >> ii >> jj >> diffuse_temp)
	{
		Diffuse_sum(ii, jj) = diffuse_temp;
	}

	// DIRECT

	std::ifstream F_DIR ("../input/PVP_SUM/direct.sum");

    for (int i = 0; i < 2; ++i) getline(F_DIR, garbage);
    F_DIR >> M_epsilon_panel >> M_phi_panel;
	for (int i = 0; i < 3; ++i) getline(F_DIR, garbage);
	S_panel=M_epsilon_panel*M_phi_panel;

    Direct_sum.resize(dimx, dimy, S_panel, 0);

	double direct_temp;

    while (F_DIR >> ii >> jj >> number_solidangle >> direct_temp)
	{
		Direct_sum(ii, jj, number_solidangle) = direct_temp;
	}

	// TERRAIN

	std::ifstream F_TER("../input/PVP_SUM/terrain.sum");

    for (int i = 0; i < 2; ++i) getline(F_TER, garbage);
    F_TER >> M_epsilon_terrain >> M_phi_terrain;
	S_terrain=M_epsilon_terrain*M_phi_terrain;

	for (int i = 0; i < 3; ++i) getline(F_TER, garbage);

    TList_sum.resize(dimx, dimy, 2, S_terrain, 0);

	double terrain_temp;
	int which_triangle;

    while (F_TER >> ii >> jj >> which_triangle >> number_solidangle >> terrain_temp)
	{
		TList_sum(ii, jj, which_triangle, number_solidangle) = terrain_temp;
	}

	if(M_epsilon_terrain!=M_epsilon_panel || M_phi_terrain!=M_phi_panel) throw std::invalid_argument( "readSumPVP: S_terrain and S_panel dont agree.\n" );

	std::cout<<"Done readSumPVP \n";

}


/**
* @brief prepare geometry for 1-panel sum
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain szrface [m]
* @param[in] SVector_temp Rotated BasicSet
* @param[out] SList_SumPVP ViewList for this one panel
*
*/
mio::Array2D<double> SolarPanel::initSListSumPVP(size_t ii, size_t jj, double height, std::vector<std::vector<double> > SVector_temp)
{
	mio::Array2D<double> SList_SumPVP;
	SList_SumPVP.resize(S_panel, 5, -999);

	size_t ii_PVP= ii;
	size_t jj_PVP= jj;
	double offset_PVP=height;

	#pragma omp parallel for
	for (size_t  number_solidangle = 0; number_solidangle < S_panel; ++ number_solidangle) // loop over solid angles
	{
		size_t solidangle_source;
		size_t ii_temp=999999;
		size_t jj_temp=999999;
		size_t which_triangle_temp=9;
		double distance, minimal_distance=dem.cellsize*(dem.getNx()+dem.getNy());

		size_t ii_dem=ii_PVP, jj_dem=jj_PVP, nb_cells=0;
		std::vector<double> projected_ray=ProjectVectorToPlane(SVector_temp[number_solidangle],{0,0,1});
		if (NormOfVector(projected_ray) != 0 ) projected_ray=normalizeVector(projected_ray);
		else projected_ray={1,0,0}; // if it goes straight up there will be sky (no caves for 2D DEM)
		while ( !(ii_dem<1 || ii_dem>dem.getNx()-2 || jj_dem<1 || jj_dem>dem.getNy()-2) ) {

			for (int k = -1; k < 2; ++k)
			{
				for (int kk = -1; kk < 2; ++kk)
				{
					ii_dem = ii_PVP + (int)round( ((double)nb_cells)*projected_ray[0] ) + k;
					jj_dem = jj_PVP + (int)round( ((double)nb_cells)*projected_ray[1] ) + kk;
					if((ii_dem<1 || ii_dem>dem.getNx()-2 || jj_dem<1 || jj_dem>dem.getNy()-2)) continue;

					distance=IntersectionRayTriangle(SVector_temp[number_solidangle],ii_PVP,jj_PVP,offset_PVP,ii_dem,jj_dem,0);
					if (distance!=-999 && distance<minimal_distance)
					{
						minimal_distance=distance;
						ii_temp=ii_dem;
						jj_temp=jj_dem;
						which_triangle_temp=0;
					}
					distance=IntersectionRayTriangle(SVector_temp[number_solidangle],ii_PVP,jj_PVP,offset_PVP,ii_dem,jj_dem,1);
					if (distance!=-999 && distance<minimal_distance)
					{
						minimal_distance=distance;
						ii_temp=ii_dem;
						jj_temp=jj_dem;
						which_triangle_temp=1;
					}
				}
			}
			ii_dem = ii_PVP + (int)round( ((double)nb_cells)*projected_ray[0] );
			jj_dem = jj_PVP + (int)round( ((double)nb_cells)*projected_ray[1] );
			nb_cells++;
		}

		if (minimal_distance==dem.cellsize*(dem.getNx()+dem.getNy())) minimal_distance=-999;

		if (minimal_distance==-999) solidangle_source=999;
		else{
			std::vector<double> v_out=VectorStretch(SVector_temp[number_solidangle], -1);
			solidangle_source=vectorToSPixel(v_out, ii_temp, jj_temp, which_triangle_temp, M_epsilon_terrain, M_phi_terrain);
		}

		SList_SumPVP(number_solidangle,0)=minimal_distance;
		SList_SumPVP(number_solidangle,1)=(double)ii_temp;
		SList_SumPVP(number_solidangle,2)=(double)jj_temp;
		SList_SumPVP(number_solidangle,3)=(double)which_triangle_temp;
		SList_SumPVP(number_solidangle,4)=(double)solidangle_source;

	}

	return SList_SumPVP;
}


/**
* @brief calculate Average-Radiation on specific panel; very similar to SolarPanel::getRadfield but for sum
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain szrface [m]
* @param[in] azimuth of panel [deg]
* @param[in] inclination of panel [deg]
* @param[out] {direct_p,diffuse_p,terrain_p} Radiation components
*/
std::vector<double> SolarPanel::projectSum(size_t ii, size_t jj, double height, double azimuth, double inclination)
{
	// generate Hemispherical Vectors
	std::vector<std::vector<double> > SVector_temp;
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		SVector_temp.push_back({0,0,0});
	}
	#pragma omp parallel for
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		std::vector<double> v_rotated=BasicSet_horizontal[number_solidangle];
		std::vector<double> normal_panel=RotationAnglesToNormalVector(azimuth, inclination);
		std::vector<double> axis=VectorCrossProduct({0,0,1}, normal_panel);
		if (NormOfVector(axis)==0) axis={1,0,0};
		if ( ( ( (int)(inclination) % 360)+360) % 360 > 180) axis=VectorStretch(axis, -1); // prevent flip of axis

		v_rotated=RotN(axis, v_rotated, inclination*Cst::to_rad);

		SVector_temp[number_solidangle]=v_rotated;
	}

	mio::Array2D<double> SList_SumPVP=initSListSumPVP(ii, jj, height, SVector_temp);
	double direct_p=0, diffuse_p=0, terrain_p=0;

	// Diffuse
	int count=0;
	double skyview_factor;
	#pragma omp parallel for reduction(+:count)
	for (size_t i = 0; i < S_panel; ++i)
	{
		if(SList_SumPVP(i,0)!=-999) count++;
	}
	skyview_factor=1-(double)count/(double)S_panel;

	diffuse_p=Diffuse_sum(ii,jj)*skyview_factor;

	// Terrain
	#pragma omp parallel for reduction(+:terrain_p)
	for (size_t number_solidangle = 0; number_solidangle < S_terrain; ++number_solidangle)
	{
		if (SList_SumPVP(number_solidangle,0)==-999) continue;

		size_t ii_source=SList_SumPVP(number_solidangle,1);
		size_t jj_source=SList_SumPVP(number_solidangle,2);
		size_t which_triangle_source=SList_SumPVP(number_solidangle,3);
		size_t solidangle_source=SList_SumPVP(number_solidangle,4);
		terrain_p+=TList_sum(ii_source, jj_source, which_triangle_source, solidangle_source)/S_panel*Cst::PI;
	}
	// Direct
	std::vector<double> PVP_normal=RotationAnglesToNormalVector(azimuth, inclination);

	#pragma omp parallel for reduction(+:direct_p)
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		if (Direct_sum(ii,jj, number_solidangle)==0) continue;
		std::vector<double> v_sun=BasicSet_horizontal[number_solidangle];

		double cos_sun_panel= VectorScalarProduct(v_sun, PVP_normal );

		if (cos_sun_panel<0) continue;
		size_t solidangle_sun=vectorToSPixel(v_sun, inclination, azimuth, M_epsilon_panel, M_phi_panel);

		if (SList_SumPVP(solidangle_sun,0)!=-999) continue; // Horizon Check

		direct_p+=Direct_sum(ii,jj, number_solidangle)*cos_sun_panel;
	}

	return {direct_p,diffuse_p,terrain_p};
}


/**
* @brief  calculate hourly radiation on specific panel; very similar to SolarPanel::getRadfield but just for 1 specific panel
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain szrface [m]
* @param[in] azimuth of panel [deg]
* @param[in] inclination of panel [deg]
* @param[out] {direct_p,diffuse_p,terrain_p} Radiation components
*/
std::vector<double> SolarPanel::projectTracker(size_t ii, size_t jj, double height, double azimuth, double inclination){

	// generate Hemispherical Vectors
	std::vector<std::vector<double> > SVector_temp;
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		SVector_temp.push_back({0,0,0});
	}

	#pragma omp parallel for
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		std::vector<double> v_rotated=BasicSet_horizontal[number_solidangle];
		std::vector<double> normal_panel=RotationAnglesToNormalVector(azimuth, inclination);
		std::vector<double> axis=VectorCrossProduct({0,0,1}, normal_panel);
		if (NormOfVector(axis)==0) axis={1,0,0};
		if ( ( ( (int)(inclination) % 360)+360) % 360 > 180) axis=VectorStretch(axis, -1); // prevent flip of axis

		v_rotated=RotN(axis, v_rotated, inclination*Cst::to_rad);

		SVector_temp[number_solidangle]=v_rotated;
	}


	mio::Array2D<double> SList_SumPVP=initSListSumPVP(ii, jj, height, SVector_temp);
	double direct_p=0, diffuse_p=0, terrain_p=0;

	// Diffuse
	int count=0;
	double skyview_factor;
	#pragma omp parallel for reduction(+:count)
	for (size_t i = 0; i < S_panel; ++i)
	{
		if(SList_SumPVP(i,0)!=-999) count+=1;
	}
	skyview_factor=1-(double)count/(double)S_panel;

	diffuse_p=d_diffuse(ii,jj)*skyview_factor;

	// Terrain
	#pragma omp parallel for reduction(+:terrain_p)
	for (size_t number_solidangle = 0; number_solidangle < S_panel; ++number_solidangle)
	{
		if (SList_SumPVP(number_solidangle,0)==-999) continue;

		size_t ii_source=SList_SumPVP(number_solidangle,1);
		size_t jj_source=SList_SumPVP(number_solidangle,2);
		size_t which_triangle_source=SList_SumPVP(number_solidangle,3);
		size_t solidangle_source=SList_SumPVP(number_solidangle,4);

		terrain_p+=TList_ms(ii_source, jj_source, which_triangle_source, solidangle_source)/S_panel*Cst::PI;
	}

	// Direct
	size_t solidangle_sun;
	double cos_sun_panel, cos_sun_horizontal;
	std::vector<double> PVP_normal=RotationAnglesToNormalVector(azimuth, inclination);

	cos_sun_panel= VectorScalarProduct(v_globalsun, PVP_normal);
	cos_sun_horizontal=VectorScalarProduct(v_globalsun, {0,0,1});

	if (cos_sun_horizontal>0) cos_sun_horizontal=1/cos_sun_horizontal;
	else cos_sun_horizontal=0;
	if (cos_sun_panel<0) {
		cos_sun_panel=0;
		solidangle_sun=0; // Need some finite number, direct radiation is anyway zero.
	}
	else solidangle_sun=vectorToSPixel(v_globalsun, inclination, azimuth, M_epsilon_panel, M_phi_panel);
	if (solidangle_sun > S_panel) solidangle_sun=0; // Hard fix to prevent crash in case of numerical rotation error
	if (SList_SumPVP(solidangle_sun,0)!=-999) cos_sun_panel=0; // shading

	direct_p=d_direct_unshaded_horizontal(ii,jj)*cos_sun_panel*cos_sun_horizontal;

	return {direct_p,diffuse_p,terrain_p};
}


/**
* @brief  Nelder-Mead blackbox optimisation algorithm, can be given different functions to minimize
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain szrface [m]
* @param[in] rounds upper limit of iterations
* @param[in] f_min function to minimize
* @param[out] {f[k_min], x[k_min][0], x[k_min][1]} min of f_min [W/m^2], optimal inclination [deg], optimal azimuth [deg]
*/
std::vector<double> SolarPanel::optimize(size_t ii, size_t jj, double height, size_t rounds, minfun f_min){
	const double a=1.0, b=1.0, g=0.5, h=0.5, tolerance=1;
	bool shrink=false;

	std::vector<double> f={0,0,0};
	std::vector<double> center_old, center_new, x_r, x_e, x_c;
	center_new={99,99}; center_old=center_new;
	double f_r, f_e, f_c;

	std::vector<std::vector<double>> x; // [inclination, azimuth]
	size_t k_min, k_max, k_n;

	x.push_back({10,0});
	x.push_back({60,-90});
	x.push_back({60,20});


	for (size_t count = 0; count < rounds; ++count)
	{
		if (count==0)
		{
			for (int i = 0; i < 3; ++i)
			{
				f[i]=(this->*f_min)(ii, jj, height, x[i]);
			}
		}
		if (count>0 && !shrink) f[k_max]=(this->*f_min)(ii, jj, height, x[k_max]);
		if (shrink){
			for (size_t i = 0; i < 3; ++i)
			{
				if(i==k_min) continue;
				f[i]=(this->*f_min)(ii, jj, height, x[i]);
			}
		}
		shrink=false;


		k_min=0; k_max=0; k_n=0;

		//Order
		for(size_t i=0;i<x.size();++i){
			if(f[i]<f[k_min]){
			  k_min=i;
			}
			if(f[i]>f[k_max]){
			  k_max=i;
			}
		}
		k_n=k_min;
		for (size_t i = 0; i < x.size(); ++i)
		{
			if (f[i]<f[k_max] && f[i]>f[k_n]) k_n=i;
		}

		// centroid & termination
		center_new=VectorStretch(VectorSum(x[k_n],x[k_min]),0.5);


		double simplex_area=fabs(0.5*((x[1][0]-x[0][0])*(x[2][1]-x[0][1])-(x[2][0]-x[0][0])*(x[1][1]-x[0][1])));
		if (simplex_area<tolerance) break;

		if(count%5==1){
			for (size_t i = 0; i < x.size(); ++i)
			{
				if(i==k_min) continue;
				x[i]=VectorSum(VectorStretch(x[i], 1.2), {10*(double)rand()/RAND_MAX-5,10*(double)rand()/RAND_MAX-5});
			}
			continue;
		}


		// Reflection
		x_r=VectorSum(center_new, VectorStretch(VectorDifference(center_new,x[k_max]),a));
		f_r=(this->*f_min)(ii, jj, height, x_r);
		if (f_r>f[k_min] && f_r<f[k_n])
		{
			x[k_max]=x_r;
			continue;
		}

		// Expansion
		if (f_r<f[k_min])
		{
			x_e=VectorSum(center_new, VectorStretch(VectorDifference(x_r,center_new),b));
			f_e=(this->*f_min)(ii, jj, height, x_e);
			if(f_e < f_r){
				x[k_max]=x_e;
				continue;
			} else{
				x[k_max]=x_r;
				continue;
			}
		}

		// Contraction
		x_c=VectorSum(center_new,VectorStretch(VectorDifference(x[k_max],center_new),g));
		f_c=(this->*f_min)(ii, jj, height, x_c);

		if (f_c<f[k_max])
		{
			x[k_max]=x_c;
			continue;
		}

		// Shrink
		for (size_t i = 0; i < x.size(); ++i)
		{
			if(i==k_min) continue;
			x[i]=VectorSum(x[k_min],VectorStretch(VectorDifference(x[i],x[k_min]),h));
		}
		shrink=1;
	}

	return {f[k_min], x[k_min][0], x[k_min][1]};
}


/**
* @brief minimisation function of Monofacial tracker
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain surface [m]
* @param[in] x Vector defining point on hemnisphere
* @param[out] total radiation [W/m^2]
*/
double SolarPanel::minfun_MonoTracker(size_t ii, size_t jj, double height, std::vector<double> x){

	std::vector<double> rad=projectTracker(ii, jj, height, x[1], x[0]);
	return -(rad[0]+rad[1]+rad[2]);
}


/**
* @brief minimisation function of Monofacial static Panel
* @param[in] ii grid coordinate (East)
* @param[in] jj grid coordinate (North)
* @param[in] height offset over terrain surface [m]
* @param[in] x Vector defining point on hemnisphere
* @param[out] total radiation [W/m^2]
*/
double SolarPanel::minfun_MonoStatic(size_t ii, size_t jj, double height, std::vector<double> x){

	std::vector<double> rad=projectSum(ii, jj, height, x[1], x[0]);
	return -(rad[0]+rad[1]+rad[2]);
}


//########################################################################################################################
//                                               ESSENTIAL FUNCTIONS
//########################################################################################################################


/**
* @brief intersection finder of direction vector and triangular surface object
*   See TerrainRadiationComplex::IntersectionRayTriangle for More information
* @param[in] ray normalized direction-vector
* @param[in] ii_PVP grid coordinate (East) of ray
* @param[in] jj_PVP grid coordinate (North) of ray
* @param[in] offset_PVP height over terrain surface [m]
* @param[in] ii_dem grid coordinate (East) of object
* @param[in] jj_dem grid coordinate (North) of object
* @param[in] which_triangle triangle A (1) or B (0)
* @param[out] distance to intersection [m] (-999 if no intersection)
*/
double SolarPanel::IntersectionRayTriangle(std::vector<double> ray, size_t ii_PVP, size_t jj_PVP, double offset_PVP, size_t ii_dem, size_t jj_dem, int which_triangle)
{
	std::vector<double> aufpunkt_ray, balance_point_par, intersection;
	std::vector<double> e_x,e_y,e_xT,e_yT,n;


	double cellsize=dem.cellsize;
	double distance, P_3, r_3;
	double intersection_x,intersection_y;

	if (which_triangle==1)
	{
		e_x={-cellsize,0, dem(ii_dem-1,jj_dem)-dem(ii_dem,jj_dem)};
		e_y={0,cellsize, dem(ii_dem,jj_dem+1)-dem(ii_dem,jj_dem)};
	}
	else
	{
		e_x={cellsize,0, dem(ii_dem+1,jj_dem)-dem(ii_dem,jj_dem)};
		e_y={0,-cellsize, dem(ii_dem,jj_dem-1)-dem(ii_dem,jj_dem)};
	}
	n=TriangleNormal(ii_dem, jj_dem, which_triangle);

	if (VectorScalarProduct(n,ray)>0) return -999; // Must hit frontal

	aufpunkt_ray={ii_PVP*cellsize ,jj_PVP*cellsize , dem(ii_PVP,jj_PVP)+offset_PVP};
	balance_point_par={ii_dem*cellsize ,jj_dem*cellsize , dem(ii_dem,jj_dem)};


	P_3=VectorScalarProduct(VectorDifference(aufpunkt_ray,balance_point_par), n);
	r_3=VectorScalarProduct(ray, n);

	if (fabs(r_3)<0.0000001) return -999;

	distance=-P_3/r_3;

	if (distance<0) return -999;

	intersection=VectorDifference(VectorSum(aufpunkt_ray, VectorStretch(ray, distance)), balance_point_par);
	e_xT=VectorCrossProduct(n,e_x);
	e_yT=VectorCrossProduct(n,e_y);

	intersection_x=VectorScalarProduct(intersection, e_yT)/VectorScalarProduct(e_x, e_yT);
	intersection_y=VectorScalarProduct(intersection, e_xT)/VectorScalarProduct(e_y, e_xT);

	if (intersection_x>=0 && intersection_y>=0 && (intersection_x+intersection_y)<=1) return distance;
	else return -999;

}


/**
* @brief Assigns index s of Basic-Set to given vector
* @param[in] vec_in normalized vector
* @param[in] inclination [deg]
* @param[in] azimuth [deg]
* @param[out] list_index
*/
size_t SolarPanel::vectorToSPixel(std::vector<double> vec_in, double inclination, double azimuth, size_t N, size_t M){

	double azimuth_flat;
	double delta, phi_temp=0, phi;
	int m,n=-1;
	size_t list_index;

	std::vector<double> vec_horizontal, vec_projected_xy;
	std::vector<double> normal_panel=RotationAnglesToNormalVector(azimuth, inclination);
	std::vector<double> axis=VectorCrossProduct({0,0,1}, normal_panel);

	if (NormOfVector(axis)==0) axis={1,0,0};
	if ( ( ( (int)(inclination) % 360)+360) % 360 > 180) axis=VectorStretch(axis, -1); // prevent flip of axis
	vec_horizontal=RotN(axis, vec_in, -inclination*Cst::to_rad);

	if (vec_horizontal[2]<0){
		std::cout<<"Vector lower than horizon SolarPanel::vectorToSPixel\n";
		return -999;
	}

	vec_projected_xy={vec_horizontal[0],vec_horizontal[1],0};

	if (NormOfVector(vec_projected_xy)!=0) azimuth_flat=AngleBetween2Vectors(vec_projected_xy,{0,1,0});
	else{ // Deal with zenith-Vector
		azimuth_flat=0;
		vec_projected_xy={0, -1, 0};
	}


	if(vec_projected_xy[0]<0) azimuth_flat=2*Cst::PI-azimuth_flat; // AngleBetween2Vectors is not uniquely defined;

	phi=AngleBetween2Vectors(vec_horizontal,vec_projected_xy);

	m=(int) M*(azimuth_flat)/2/Cst::PI; // Inverse operation Of azi=2*3.1415926*(m+0.5)/M; the 0.5 is off for a proper rounding
	size_t i=0;
	while(n==-1)
	{

		if (i!=(N-1)){
			delta=acos(-2/float(N)+cos(2*phi_temp))/2-phi_temp;
			phi_temp+=delta;
		}
		else phi_temp=Cst::PI/2;
		if(phi<=phi_temp) n=i;
		i+=1;
		if (i>N) std::cout<<"[E] in SolarPanel::vectorToSPixel\n";

	}

	list_index=n*M+m;
	return list_index;
}


/**
* @brief Assigns index s of Basic-Set to given vector of surface object
* @param[in] vec_in normalized vector
* @param[in] ii_dem grid coordinate (East) of object
* @param[in] jj_dem grid coordinate (North) of object
* @param[in] which_triangle triangle A (1) or B (0)
* @param[out] list_index
*/
size_t SolarPanel::vectorToSPixel(std::vector<double> vec_in, size_t ii_dem, size_t jj_dem, size_t which_triangle, size_t N, size_t M){

	std::vector<double> triangle_normal=TriangleNormal(ii_dem, jj_dem, which_triangle);
	std::vector<double> angles=NormalVectorToRotationAngles(triangle_normal);

	double azimuth=angles[0];
	double inclination=angles[1];

	return vectorToSPixel(vec_in, inclination, azimuth, N, M);
}


/**
* @brief Assigns index s of Basic-Set to given vector of panel p
* @param[in] vec_in normalized vector
* @param[in] number_pvp
* @param[out] list_index
*/
size_t SolarPanel::vectorToSPixel(std::vector<double> vec_in, size_t number_pvp, size_t N, size_t M){

	double inclination_panel, azimuth_panel;

	inclination_panel=pv_points[number_pvp][3];
	azimuth_panel=pv_points[number_pvp][4];

	return vectorToSPixel(vec_in,inclination_panel,azimuth_panel, N, M);
}


/**
* @brief Assigns index s of flat Basic-Set to given vector
* @param[in] vec_in normalized vector
* @param[out] list_index
*/
size_t SolarPanel::vectorToSPixel(std::vector<double> vec_in,size_t N, size_t M){

	return vectorToSPixel(vec_in,0,0, N, M);
}


/**
* @brief yes/no intersection of given vector of flat basic set with shadow on terrain
* @param[in] v_sun normalized vector pointing to the sun
* @param[in] number_pvp number of pv-panel as initialized
* @param[in] number_solidangle index of BasicSet
* @param[out] yes/no intersection
*/
bool SolarPanel::doesPanelShadowPixel(std::vector<double> v_sun, size_t number_pvp, size_t number_solidangle)
{

	std::vector<double> ray=BasicSet_rotated[number_pvp][number_solidangle];
	double distance=ViewList_panel[number_pvp][number_solidangle][0];
	std::vector<double> r_0={pv_points[number_pvp][0],pv_points[number_pvp][1], dem(get_ii(number_pvp), get_jj(number_pvp))+pv_points[number_pvp][2]};

	for (size_t pp = 0; pp < pv_points.size(); ++pp)
	{

		//if (pp==number_pvp) continue;
		double d_intersect, intersection_x,intersection_y;
		double panel_height=pv_points[pp][5];
		double panel_width=pv_points[pp][6];
		double inclination_panel=pv_points[pp][3]*Cst::to_rad;
		double azimuth_panel=pv_points[pp][4]*Cst::to_rad;

		std::vector<double> normal_pp=RotationAnglesToNormalVector(pv_points[pp][4], pv_points[pp][3]);
		std::vector<double> r_pp={pv_points[pp][0],pv_points[pp][1], dem(get_ii(pp), get_jj(pp))+pv_points[pp][2]};
		std::vector<double> delta_r=VectorDifference(r_pp, r_0);
		std::vector<double> e_x,e_y,n, intersection;

		d_intersect=(VectorScalarProduct(delta_r, normal_pp) -distance*VectorScalarProduct(ray, normal_pp))/VectorScalarProduct(v_sun, normal_pp);
		if (d_intersect<0) return false;

		intersection=VectorDifference(VectorSum(VectorStretch(ray,distance),VectorStretch(v_sun,d_intersect)), delta_r);

		e_x=RotX({0,-1,0}, inclination_panel);
		e_x=RotZ(e_x, -azimuth_panel);
		e_y=RotX({1,0,0}, inclination_panel);
		e_y=RotZ(e_y, -azimuth_panel);

		intersection_x=VectorScalarProduct(intersection, e_x);
		intersection_y=VectorScalarProduct(intersection, e_y);

		if (fabs(intersection_x)<=panel_height/2 && fabs(intersection_y)<=panel_width/2) return true;
	}

	return false;
}


/**
* @brief land-view factor of panel
* @param[in] name_pvp number of pv-panel as initialized
* @param[out] land-view factor €[0,1]
*/
double SolarPanel::getLandViewFactor(size_t name_pvp)
{
	std::vector<std::vector<double> > SList=ViewList_panel[name_pvp];
	double sum=0;
	double factor;

	for (size_t l = 0; l < SList.size(); ++l)
	{
		if (SList[l][0]!=-999)
		{
			sum++;
		}
	}

	factor=sum/M_epsilon_panel/M_phi_panel;
	return factor;
}


/**
* @brief sky-view factor of panel
* @param[in] name_pvp number of pv-panel as initialized
* @param[out] sky-view factor €[0,1]
*/
double SolarPanel::getSkyViewFactor(size_t name_pvp)
{
	return 1-getLandViewFactor(name_pvp);
}


/**
* @brief Surface-normal of surface object of triangulation
* @param[in] ii_dem grid coordinate (East) of object
* @param[in] jj_dem grid coordinate (North) of object
* @param[in] which_triangle triangle A (1) or B (0)
* @param[out] Surface-normal
*/
std::vector<double> SolarPanel::TriangleNormal(size_t ii_dem, size_t jj_dem, int which_triangle)
{
	std::vector<double> e_x,e_y,n;
	double cellsize=dem.cellsize;

	if (which_triangle==1)
	{
		e_x={-cellsize,0, dem(ii_dem-1,jj_dem)-dem(ii_dem,jj_dem)};
		e_y={0,cellsize, dem(ii_dem,jj_dem+1)-dem(ii_dem,jj_dem)};
	}
	else
	{
		e_x={cellsize,0, dem(ii_dem+1,jj_dem)-dem(ii_dem,jj_dem)};
		e_y={0,-cellsize, dem(ii_dem,jj_dem-1)-dem(ii_dem,jj_dem)};
	}
	n=VectorCrossProduct(e_y, e_x);

	return normalizeVector(n);

}


//########################################################################################################################
//                                               AUXILARY FUNCTIONS
//########################################################################################################################


/**
* @brief gives angles of vector of Basic set
* @param[in] index of vector of BasicSet
* @param[out] {azimuth, phi_eff} in [rad]
*/
std::vector<double> SolarPanel::listindexToAngles(size_t index)
{
	size_t m, n;
	double azimuth, phi_eff, phi1,phi0;
	double delta=0;

	m=index%M_phi_panel;
	n=(index-m)/M_phi_panel;


	azimuth=2*Cst::PI*(m+0.5)/M_phi_panel;
	phi1=0;
	phi0=0;


	for (size_t k = 0; k < n+1; ++k)
	{
		phi0=phi1;
		if (k==M_epsilon_panel) std::cout<<"problem in listindexToAngles\n";
		if (k==(M_epsilon_panel-1))
		{
			phi1=Cst::PI/2;
		}else{
		delta=acos(-2/float(M_epsilon_panel)+cos(2*phi0))/2-phi0;
		phi1=phi0+delta;
		}
		phi_eff=(phi1+phi0)/2;
	}

	return {azimuth, phi_eff};

}


/**
* @brief gives rotation angles of normal vector
* @param[in] normal normalized normal vector
* @param[out] {azimuth, phi} in [deg]
*/
std::vector<double> SolarPanel::NormalVectorToRotationAngles(std::vector<double> normal)
{
	double azimuth, phi;
	std::vector<double> n_projected={normal[0],normal[1],0};

	phi=acos( VectorScalarProduct(normal,{0,0,1})/NormOfVector(normal) );

	if (NormOfVector(n_projected)==0) return {0, phi*Cst::to_deg};

	azimuth=acos( VectorScalarProduct(n_projected, {0,-1,0})/NormOfVector(n_projected) );
	if (normal[0]>0) azimuth=-azimuth;


	return {azimuth*Cst::to_deg, phi*Cst::to_deg};
}


/**
* @brief calculates Normal vector of panel given the rotation angles
* @param[in] azimuth angle in [deg]
* @param[in] phi (inclination) angle in[deg]
* @param[out] v_out normalized vector
*/
std::vector<double> SolarPanel::RotationAnglesToNormalVector(double azimuth, double phi)
{
	std::vector<double> v_out={0,0,0};
	v_out[0]=-sin(phi*Cst::to_rad)*sin(azimuth*Cst::to_rad);
	v_out[1]=-sin(phi*Cst::to_rad)*cos(azimuth*Cst::to_rad);
	v_out[2]=cos(phi*Cst::to_rad);

	return v_out;
}

int SolarPanel::get_ii(int number_pvp){

			Coords coord_temp;

			coord_temp.setEPSG(21781);
			//add 1/2*cellsize because gridify always rounds to next lower pixel
			coord_temp.setXY(pv_points[number_pvp][0]+dem.cellsize/2, pv_points[number_pvp][1]+dem.cellsize/2, pv_points[number_pvp][2]);
			dem.gridify(coord_temp);

			return coord_temp.getGridI();
}

int SolarPanel::get_jj(int number_pvp){

			Coords coord_temp;

			coord_temp.setEPSG(21781);
			//add 1/2*cellsize because gridify always rounds to next lower pixel
			coord_temp.setXY(pv_points[number_pvp][0]+dem.cellsize/2, pv_points[number_pvp][1]+dem.cellsize/2, pv_points[number_pvp][2]);
			dem.gridify(coord_temp);

			return coord_temp.getGridJ();
}

std::vector<double> SolarPanel::getVectorSun(const double solarAzimuth, const double solarElevation)
{
	std::vector<double> v_out={0,0,0};
	v_out[0]=cos(solarElevation*Cst::to_rad)*sin(solarAzimuth*Cst::to_rad);
	v_out[1]=cos(solarElevation*Cst::to_rad)*cos(solarAzimuth*Cst::to_rad);
	v_out[2]=sin(solarElevation*Cst::to_rad);

	return v_out;

}



//########################################################################################################################
//                                               ELEMENTARY FUNCTIONS
//########################################################################################################################

double SolarPanel::AngleBetween2Vectors(std::vector<double> vec1, std::vector<double> vec2)
{
	double sum=0;
	double angle=0;
	double norm1,norm2;
	if(vec1.size()!=vec2.size()) return -999;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		sum+=vec1[i]*vec2[i];
	}
	norm1=NormOfVector(vec1);
	norm2=NormOfVector(vec2);

	double val=sum/norm1/norm2;
	if(val>1) val=1;
	if(val<-1) val=-1;

	angle=acos(val);
	return angle;
}

double SolarPanel::NormOfVector(std::vector<double> vec1)
{
	return sqrt(VectorScalarProduct(vec1,vec1));
}

std::vector<double> SolarPanel::normalizeVector(std::vector<double> vec1)
{
	double norm=NormOfVector(vec1);
	if (norm==0){
		std::cout<<"Exception in SolarPanel::normalizeVector : norm is zero\n";
		return vec1;
	}

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		vec1[i]=vec1[i]/norm;
	}

	return vec1;

}


std::vector<double> SolarPanel::ProjectVectorToPlane(std::vector<double> vec1, std::vector<double> plane_normal)
{
	std::vector<double> normal=plane_normal, v_out=vec1;
	normal=normalizeVector(normal);
	double proj=VectorScalarProduct(normal, vec1);

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out[i]=vec1[i]-proj*normal[i];
	}

	return v_out;
}

double SolarPanel::VectorScalarProduct(std::vector<double> vec1, std::vector<double> vec2)
{
	double sum=0;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		sum+=vec1[i]*vec2[i];
	}

	return sum;
}

std::vector<double> SolarPanel::VectorCrossProduct(std::vector<double> vec1, std::vector<double> vec2)
{
	std::vector<double> v_out;
	v_out={0,0,0};

	v_out[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
	v_out[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
	v_out[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];

	return v_out;
}

std::vector<double> SolarPanel::VectorSum(std::vector<double> vec1, std::vector<double> vec2)
{
	std::vector<double> v_out;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out.push_back(vec1[i]+vec2[i]);
	}

	return v_out;
}

std::vector<double> SolarPanel::VectorDifference(std::vector<double> vec1, std::vector<double> vec2)
{
	std::vector<double> v_out;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out.push_back(vec1[i]-vec2[i]);
	}

	return v_out;
}

std::vector<double> SolarPanel::VectorStretch(std::vector<double> vec1, double factor)
{
	std::vector<double> v_out;

	for (size_t i = 0; i < vec1.size(); ++i)
	{
		v_out.push_back(vec1[i]*factor);
	}

	return v_out;
}

std::vector<double> SolarPanel::RotX(std::vector<double> vec_in, double rad)
{
	std::vector<double> vec_out(3);
	vec_out[0]=vec_in[0];
	vec_out[1]=cos(rad)*vec_in[1]-sin(rad)*vec_in[2];
	vec_out[2]= sin(rad)*vec_in[1]+cos(rad)*vec_in[2];

	return vec_out;
}

std::vector<double> SolarPanel::RotY(std::vector<double> vec_in, double rad)
{
	std::vector<double> vec_out(3);
	vec_out[0]=cos(rad)*vec_in[0]+sin(rad)*vec_in[2];
	vec_out[1]=vec_in[1];
	vec_out[2]= -sin(rad)*vec_in[0]+cos(rad)*vec_in[2];

	return vec_out;
}

std::vector<double> SolarPanel::RotZ(std::vector<double> vec_in, double rad)
{
	std::vector<double> vec_out(3);
	vec_out[0]=cos(rad)*vec_in[0]-sin(rad)*vec_in[1];
	vec_out[1]=sin(rad)*vec_in[0]+cos(rad)*vec_in[1];
	vec_out[2]= vec_in[2];

	return vec_out;
}

std::vector<double> SolarPanel::RotN(std::vector<double> axis, std::vector<double> vec_in, double rad)
{
	std::vector<double> vec_out(3);
	double c=cos(rad);
	double s=sin(rad);

	double n1=axis[0]/NormOfVector(axis);
	double n2=axis[1]/NormOfVector(axis);
	double n3=axis[2]/NormOfVector(axis);

	vec_out[0]=(n1*n1*(1-c)+c)*vec_in[0]+(n1*n2*(1-c)-n3*s)*vec_in[1]+(n1*n3*(1-c)+n2*s)*vec_in[2];
	vec_out[1]=(n2*n1*(1-c)+n3*s)*vec_in[0]+(n2*n2*(1-c)+c)*vec_in[1]+(n2*n3*(1-c)-n1*s)*vec_in[2];
	vec_out[2]=(n3*n1*(1-c)-n2*s)*vec_in[0]+(n3*n2*(1-c)+n1*s)*vec_in[1]+(n3*n3*(1-c)+c)*vec_in[2];

	return vec_out;
}


//########################################################################################################################
//                                                     TEST OUTPUT
//########################################################################################################################



void SolarPanel::PrintVector(std::vector<double> vec1)
{
	for (size_t i = 0; i < vec1.size()-1; ++i)
	{
		std::cout<<std::fixed <<std::setprecision(3)<<vec1[i]<<"  ";
	}
	std::cout<<std::fixed <<std::setprecision(3)<<vec1[vec1.size()-1]<<"\n";
}

// Hemispherical radiation Map for 1 gridpoint and Summed radiation
void SolarPanel::RadiationMap(size_t ii, size_t jj, double elevation){

	mio::Array2D<double> Rad_List;
	Rad_List.resize(S_panel, 3);

	std::ofstream VF_file;
	VF_file.open("Radiation.map");


	VF_file<<"index, azimuth, phi, radiation\n";

	#pragma omp parallel for
	for (size_t solidangle = 0; solidangle < S_panel; ++solidangle)
	{
		std::vector<double> angles=listindexToAngles(solidangle);
		angles[0]=angles[0]*Cst::to_deg-180; // to rotation default (azimuth)
		angles[1]=90-angles[1]*Cst::to_deg; // to rotation default (inclination))


		std::vector<double> radiation=projectSum(ii, jj, elevation, angles[0],angles[1]);
		double radiation_tot=radiation[0]+radiation[1]+radiation[2];

		Rad_List(solidangle,0)=angles[0];
		Rad_List(solidangle,1)=angles[1];
		Rad_List(solidangle,2)=radiation_tot;

		std::cout<<"RadiationMap: done Solidangle"<<solidangle<<"\n";
	}


	for (size_t solidangle = 0; solidangle < S_panel; ++solidangle)
	{

		VF_file<<solidangle<<", "<<Rad_List(solidangle,0)<<", "<<Rad_List(solidangle,1)<<", "<<Rad_List(solidangle,2)<<"\n";

	}

	VF_file.close();
}


// Optimize radiation and angles for whole DEM
void SolarPanel::GridRadiationMap(double offset){

	mio::Array3D<double> Rad_List;
	Rad_List.resize(dimx, dimy,3, -999);

	std::ofstream VF_file;
	VF_file.open("Grid_Radiation.map");

	VF_file<<"ii, jj, altitude, radiation, inclination, azimuth\n";

	for (size_t ii = 200; ii < 225; ii+=1) // for (size_t ii = 1; ii < dimx-1; ii+=1)
	{
		#pragma omp parallel for
		for (size_t jj = 1; jj < dimy-1; jj+=1)
		{
			std::vector<double> radiation=optimize(ii, jj, offset, 40, &SolarPanel::minfun_MonoStatic);

			Rad_List(ii,jj,0)=radiation[0];
			Rad_List(ii,jj,1)=radiation[1];
			Rad_List(ii,jj,2)=radiation[2];

		}
		std::cout<<"Done ii="<<ii<<"\n";
	}

	for (size_t ii = 1; ii < dimx-1; ii+=1)
	{
		for (size_t jj = 1; jj < dimy-1; jj+=1)
		{
			VF_file<<ii<<", "<<jj<<", "<<dem(ii,jj)<<", "<<Rad_List(ii,jj,0)<<", "<<Rad_List(ii,jj,1)<<", "<<Rad_List(ii,jj,2)<<"\n";
		}
	}

	VF_file.close();
	std::cout<<"Done GridRadiationMap \n";
}



// Write OptimumTracker radiation for 1 hour
void SolarPanel::WriteOptimumTrackerRadiation(size_t number_pvp, std::string filename, const mio::Date timestamp){

	if (Terrain_complex_mode!=true) throw std::invalid_argument( "[E] SolarPanel::WriteOptimumTrackerRadiation only runs in TerrainRadiationComplex-mode \n");

	double height=3;

	std::vector<double> optimum_tracker, radiation_tracker;
	if (v_globalsun[2]>0)
	{
		optimum_tracker=optimize(get_ii(number_pvp),get_jj(number_pvp),height, 20, &SolarPanel::minfun_MonoTracker);
		radiation_tracker=projectTracker(get_ii(number_pvp),get_jj(number_pvp),height, optimum_tracker[2], optimum_tracker[1]);
	}
	else{
		optimum_tracker={0,0,0};
		radiation_tracker={0,0,0};
	}

	std::ofstream TF;
	TF.open(filename, std::ios_base::app | std::ios_base::out);
	TF<<std::fixed <<std::setprecision(6) <<timestamp.toString(Date::ISO)<<","<<timestamp.getJulian()<<","<<optimum_tracker[2]<<","<<optimum_tracker[1]<<","<<radiation_tracker[0]<<","<<radiation_tracker[1]<<","<<radiation_tracker[2]<<","<<-optimum_tracker[0]<<","<<albedo(get_ii(number_pvp),get_jj(number_pvp))<<","<<direct_beam[number_pvp]<<"\n";
	TF.close();
}


// Write SunTracker radiation for 1 hour
void SolarPanel::WriteSunTrackerRadiation(size_t number_pvp, std::string filename, const mio::Date timestamp){

	if (Terrain_complex_mode!=true) throw std::invalid_argument( "[E] SolarPanel::WriteSunTrackerRadiation only runs in TerrainRadiationComplex-mode \n");

	double height=3;
	std::vector<double> angles=NormalVectorToRotationAngles(v_globalsun);

	std::vector<double> optimum_tracker, radiation_tracker;
	if (v_globalsun[2]>0)
	{
		radiation_tracker=projectTracker(get_ii(number_pvp),get_jj(number_pvp),height, angles[0], angles[1]);
	}
	else{
		radiation_tracker={0,0,0};
	}
	double radiation_tot=radiation_tracker[0]+radiation_tracker[1]+radiation_tracker[2];

	std::ofstream TF;
	TF.open(filename, std::ios_base::app | std::ios_base::out);
	TF<<std::fixed <<std::setprecision(6) <<timestamp.toString(Date::ISO)<<","<<timestamp.getJulian()<<","<<angles[0]<<","<<angles[1]<<","<<radiation_tracker[0]<<","<<radiation_tracker[1]<<","<<radiation_tracker[2]<<","<<radiation_tot<<","<<albedo(get_ii(number_pvp),get_jj(number_pvp))<<","<<direct_beam[number_pvp]<<"\n";
	TF.close();
}
