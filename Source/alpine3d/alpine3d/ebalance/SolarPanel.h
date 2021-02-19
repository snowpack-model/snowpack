#ifndef SOLARPANEL_H
#define SOLARPANEL_H


#include <meteoio/MeteoIO.h>
#include <alpine3d/ebalance/RadiationField.h>
#include <meteoio/dataClasses/Grid2DObject.h>

#include <alpine3d/ebalance/SnowBRDF.h>

class SolarPanel{

	public:
		SolarPanel(){};

		SolarPanel(const mio::Config& cfg, const mio::DEMObject &dem_in, const std::vector<std::vector<double> > &pv_pts);

		void setSP(const mio::Date timestamp, const double solarAzimuth, const double solarElevation);
		void setGridRadiation(const mio::Array2D<double>& in_albedo, const mio::Array2D<double>& in_direct,
                          const mio::Array2D<double>& in_diffuse, const mio::Array2D<double>&
                          in_direct_unshaded_horizontal, const double solarAzimuth, const double solarElevation);

		void initTerrain(size_t N_terrain_in, size_t M_terrain_in);
		void setTLists(mio::Array4D<double> TList1, mio::Array4D<double> TList2, mio::Array4D<double> TList3, mio::Array4D<double> TList4);
		void writeSP(const unsigned int max_steps);

	private:
		typedef double (SolarPanel::*minfun)(size_t, size_t, double, std::vector<double>);


		// Initialisation Functions Functions
		void getRadfield();
		void writeHeader();

		void initBasicSetHorizontal();
		void initBasicSetRotated();
		void initViewListPanel();
		void initViewListTerrain();
		void initShadelist();


		// Essential Functions
		size_t vectorToSPixel(std::vector<double> vec_in, size_t N, size_t M);
		size_t vectorToSPixel(std::vector<double> vec_in, size_t number_pvp, size_t N, size_t M);
		size_t vectorToSPixel(std::vector<double> vec_in, double inclination_panel, double azimuth_panel, size_t N, size_t M);
		size_t vectorToSPixel(std::vector<double> vec_in, size_t ii_dem, size_t jj_dem, size_t which_triangle, size_t N, size_t M);

		double IntersectionRayTriangle(std::vector<double> ray, size_t ii_PVP, size_t jj_PVP, double offset_PVP, size_t ii_dem, size_t jj_dem, int which_triangle);
		bool doesPanelShadowPixel(std::vector<double> v_sun, size_t number_pvp, size_t number_solidangle);
		double getLandViewFactor(size_t name_pvp);
		double getSkyViewFactor(size_t name_pvp);


		// Sum & Tracker Functions
		void initSumPVP();
		void readSumPVP();
		mio::Array2D<double> initSListSumPVP(size_t ii, size_t jj, double height, std::vector<std::vector<double> > SVector_temp);

		std::vector<double> projectSum(size_t ii, size_t jj, double height, double azimuth, double inclination);

		std::vector<double> projectTracker(size_t ii, size_t jj, double height, double azimuth, double inclination);

		std::vector<double> optimize(size_t ii, size_t jj, double height, size_t rounds, minfun f_min);
		double minfun_MonoTracker(size_t ii, size_t jj, double height, std::vector<double> x);
		double minfun_MonoStatic(size_t ii, size_t jj, double height, std::vector<double> x);


		// auxiliary functions
		std::vector<double> listindexToAngles(size_t index);
		std::vector<double> NormalVectorToRotationAngles(std::vector<double> normal);
		std::vector<double> RotationAnglesToNormalVector(double azimuth, double phi);
		int get_ii(int number_pvp);
		int get_jj(int number_pvp);
		std::vector<double> TriangleNormal(size_t ii_dem, size_t jj_dem, int which_triangle);
		std::vector<double> getVectorSun(const double solarAzimuth, const double solarElevation);


		// Elementary functions
		double AngleBetween2Vectors(std::vector<double> vec1, std::vector<double> vec2);
		double NormOfVector(std::vector<double> vec1);
		std::vector<double> normalizeVector(std::vector<double> vec1);
		std::vector<double> ProjectVectorToPlane(std::vector<double> vec1, std::vector<double> plane_normal);
		double VectorScalarProduct(std::vector<double> vec1, std::vector<double> vec2);
		std::vector<double> VectorCrossProduct(std::vector<double> vec1, std::vector<double> vec2);

		std::vector<double> VectorSum(std::vector<double> vec1, std::vector<double> vec2);
		std::vector<double> VectorDifference(std::vector<double> vec1, std::vector<double> vec2);
		std::vector<double> VectorStretch(std::vector<double> vec1, double factor);

		std::vector<double> RotX(std::vector<double> vec_in, double rad);
		std::vector<double> RotY(std::vector<double> vec_in, double rad);
		std::vector<double> RotZ(std::vector<double> vec_in, double rad);
		std::vector<double> RotN(std::vector<double> axis, std::vector<double> vec_in, double rad);


		// Output functions for testing
		void PrintVector(std::vector<double> vec1);
		void RadiationMap(size_t ii, size_t jj, double elevation);
		void GridRadiationMap(double offset);
		void WriteOptimumTrackerRadiation(size_t number_pvp, std::string filename, const mio::Date timestamp);
		void WriteSunTrackerRadiation(size_t number_pvp, std::string filename, const mio::Date timestamp);


		size_t dimx, dimy;
		size_t M_epsilon_panel, M_phi_panel; // sperical pixelcoordinates
		size_t S_panel;
		size_t M_epsilon_terrain, M_phi_terrain;
		size_t S_terrain;
		bool Terrain_complex_mode=false;
		std::vector<double> v_globalsun;

		mio::DEMObject dem;
		std::vector<std::vector<double> > pv_points;
		SnowBRDF BRDFobject;
		mio::Timer timer;
		mio::Array3D<size_t> Shadelist;
		bool generate_PVP_sum;
		bool sun_tracker;
		bool optimal_tracker;

		std::vector<std::vector<double> > BasicSet_horizontal; //Spherical-pixel vectors horizontal case
		std::vector<std::vector<std::vector<double> > > BasicSet_rotated; //For each PVP a set of Spherical-pixel vectors according to inclination and azimuth of Pannel
		std::vector<std::vector<std::vector<double> > > ViewList_panel; //For each Spherical-pixel (and for all PVP's), the ii.jj of the DEM-grid are assigned as well as the distance


		std::vector<double> direct, diffuse, terrain_iso, terrain_aniso, terrain_ms, terrain_ms_noshadow, direct_beam;
		mio::Array4D<double> TList_direct, TList_sky_iso, TList_sky_aniso,  TList_ms;
		mio::Array4D<double> TList_sum;
		mio::Array3D<double> Direct_sum;
		mio::Array2D<double> Diffuse_sum;

		mio::Array2D<double> albedo;
		mio::Array2D<double> d_direct_A, d_direct_B, d_diffuse, d_direct_unshaded_horizontal;
		const mio::Date* dateobject;
		std::vector<std::string> filenames;

		bool if_shadowing=true;
		size_t alpine_steps=0;

};

#endif
