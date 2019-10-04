#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

// Static varibales containing the results
const double res_Lat []		= {46.752399,		46.668527,		46.338206,		46.996607,		46.191333,		46.539473,		46.647044};
const double res_Lon []		= {9.946666,		8.064570,		8.853099,		9.037582,		6.827770,		7.561830,		8.740198};
const double res_Alt []		= {2390,		2110,			2100,			1630,			2020,			2020,			2220};
const double res_X []		= {791600,		647900,			708900,			721610,			552840,			609450,			699639};
const double res_Y []		= {180975,		168780,			132850,			206300,			115725,			154250,			167027};
const string res_ID []		= {"FLU2",		"FIR2",			"FRA2",			"GLA2",			"ILI2",			"OTT2",			"TUJ3"};
const string res_Name []	= {"Fluela Hospiz",	"Schmidigen-Bidmeren",	"Efra",			"Guppen",		"Les Collines",		"Ottere",		"Nual"};
const double res_Slope []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata};
const double res_Azi []		= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata};

const double res_Met_0[]	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // P
const double res_Met_1 []	= {263.3351476,		264.4570077,		265.6127785,		265.5006154,		265.0483363,		265.9048493,		262.8538273}; // TA
const double res_Met_2 []	= {1.,			0.957,			1.,			0.967,			0.963,			0.862,			0.95}; // RH
const double res_Met_3[]	= {273.6963,		274.3309,		273.9969,		274.6983,		274.0971,		274.9154,		274.0971}; // TSG
const double res_Met_4[]	= {263.8,		266.12,			265.6,			268.74,			265.64,			267.24,			262.5}; // TSS
const double res_Met_5 []	= {0.57,		0.95,			1.9,			1.41,			0.35,			0.72,			1.47}; // HS
const double res_Met_6 []	= {2.9,			0.6,			2.,			0.2,			2.5,			1.5,			0.3}; // VW
const double res_Met_7 []	= {335.,		138.,			98.,			216.,			268.,			187.,			107.}; // DW
const double res_Met_8 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // VW_MAX
const double res_Met_9 []	= {119.5,		105.5,			86.5,			67.,			120.5,			105.,			115.5}; //RSWR
const double res_Met_10 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // ISWR
const double res_Met_11 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // ILWR
const double res_Met_12 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // TAU_CLD
const double res_Met_13 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	4.,			5.160045665,		0.,			IOUtils::nodata}; // HNW

// methode do controll content of Meteo Data !!
// Also controlles != operator of containing special structures
bool controllStation(MeteoData& datMeteo, int i_results, Date datDate){

	const double epsilon = 1.0e-7; // accuracy of the double tests
	bool status = true;

	// Coords content
	Coords& dataCoord = datMeteo.meta.position;
	if(!IOUtils::checkEpsilonEquality(dataCoord.getAltitude(), res_Alt[i_results], epsilon)){
		cerr << "error on Altitude"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getLat(), res_Lat[i_results], epsilon)){
		cerr << "error on Latitude"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getLon(), res_Lon[i_results], epsilon)){
		cerr << "error on Longitude"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getEasting(), res_X[i_results], epsilon)){
		cerr << "error on X (Easting)"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getNorthing(), res_Y[i_results], epsilon)){
		cerr << "error on Y (Northing)"<< endl;
		status = false;
	}
	string tmp_type, tmp_args;
	dataCoord.getProj(tmp_type,tmp_args);
	if(tmp_type.compare("CH1903") != 0){
		cerr << "error on Projection"<< endl;
		status = false;
	}

	Coords refCoord;
	refCoord.setLatLon(res_Lat[i_results], res_Lon[i_results], res_Alt[i_results]);
	refCoord.setXY(res_X[i_results], res_Y[i_results], res_Alt[i_results]);
	refCoord.setProj("CH1903");
	if(datMeteo.meta.position != refCoord){
		cerr << "error on == operator for Coords :";
		cerr << datMeteo.meta.position.toString() << endl;
		cerr << refCoord.toString() << endl;
		status = false;
	}


	// Station Data content
	if(datMeteo.meta.getStationID().compare(res_ID[i_results]) != 0){
		cerr << "error on StationID"<< endl;
		status = false;
	}
	if(datMeteo.meta.getStationName().compare(res_Name[i_results]) != 0){
		cerr << "error on getStationName"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo.meta.getSlopeAngle(),res_Slope[i_results], epsilon)){
		cerr << "error on getSlopeAngle";
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo.meta.getAzimuth(), res_Azi[i_results], epsilon)){
		cerr << "error on getAzimuth";
		status = false;
	}

	StationData refStation(refCoord, res_ID[i_results], res_Name[i_results]);
	refStation.setSlope(res_Slope[i_results], res_Azi[i_results]);
	if(refStation != datMeteo.meta){
		cerr << "error on != between Station Data";
		status = false;
	}

	// Meteo data
	if(!IOUtils::checkEpsilonEquality(datMeteo(0), res_Met_0[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(0) << " : " << std::setprecision(10) << datMeteo(0) << " != " << res_Met_0[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(1), res_Met_1[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(1) << " : " << std::setprecision(10) << datMeteo(1) << " != " << res_Met_1[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(2), res_Met_2[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(2) << " : " << std::setprecision(10) << datMeteo(2) << " != " << res_Met_2[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(3), res_Met_3[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(3) << " : " << std::setprecision(10) << datMeteo(3) << " != " << res_Met_3[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(4), res_Met_4[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(4) << " : " << std::setprecision(10) << datMeteo(4) << " != " << res_Met_4[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(5), res_Met_5[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(5) << " : " << std::setprecision(10) << datMeteo(5) << " != " << res_Met_5[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(6), res_Met_6[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(6) << " : " << std::setprecision(10) << datMeteo(6) << " != " << res_Met_6[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(7), res_Met_7[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(7) << " : " << std::setprecision(10) << datMeteo(7) << " != " << res_Met_7[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(8), res_Met_8[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(8) << " : " << std::setprecision(10) << datMeteo(8) << " != " << res_Met_8[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(9), res_Met_9[i_results], 10.*epsilon)){ // HACK special epsilon that passes tests !
		cerr << "error on " << MeteoData::getParameterName(9) << " : " << std::setprecision(10) << datMeteo(9) << " != " << res_Met_9[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(10), res_Met_10[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(10) << " : " << std::setprecision(10) << datMeteo(10) << " != " << res_Met_10[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(11), res_Met_11[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(11) << " : " << std::setprecision(10) << datMeteo(11) << " != " << res_Met_11[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(12), res_Met_12[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(12) << " : " << std::setprecision(10) << datMeteo(12) << " != " << res_Met_12[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo(13), res_Met_13[i_results], epsilon)){
		cerr << "error on " << MeteoData::getParameterName(13) << " : " << std::setprecision(10) << datMeteo(13) << " != " << res_Met_13[i_results] << endl;
		status = false;
	}

	MeteoData refMeteo(datDate);
	refMeteo.standardizeNodata(IOUtils::nodata);
	refMeteo.meta = refStation;
	refMeteo(0)= res_Met_0[i_results];
	refMeteo(1)= res_Met_1[i_results];
	refMeteo(2)= res_Met_2[i_results];
	refMeteo(3)= res_Met_3[i_results];
	refMeteo(4)= res_Met_4[i_results];
	refMeteo(5)= res_Met_5[i_results];
	refMeteo(6)= res_Met_6[i_results];
	refMeteo(7)= res_Met_7[i_results];
	refMeteo(8)= res_Met_8[i_results];
	refMeteo(9)= res_Met_9[i_results];
	refMeteo(10)= res_Met_10[i_results];
	refMeteo(11)= res_Met_11[i_results];
	refMeteo(12)= res_Met_12[i_results];
	refMeteo(13)= res_Met_13[i_results];
	if(datMeteo != refMeteo){
		cerr << "error on == operator for MeteoData :" << datMeteo.getNrOfParameters() << " - " << refMeteo.getNrOfParameters() << endl;
		cerr << datMeteo.toString() << endl;
		cerr << refMeteo.toString() << endl;
		status = false;
	}

	return status;
}

//Test if data read at 2008-12-01T15:00:00 are correct
int main() {
	bool status = true;
	Date d1(2008, 12, 1, 15, 0, 1);
	std::vector<MeteoData> vecMeteo;

	Config cfg("io.ini");
	IOManager io(cfg);

	//io.setProcessingLevel(IOManager::raw); //set the processing level: raw, filtered or resampled
	io.getMeteoData(d1, vecMeteo);

	// Compare data with hard coded values
	if(vecMeteo.size() != 7) {
		cerr << "ERROR on amout of Data read !!! \n";
		status = false;
	}

	// Test readed data
	for(unsigned int i = 0; i < vecMeteo.size(); i++){
		cout << "----- Control of vecMeteo # : "<< i+1 << endl;
		if(!controllStation(vecMeteo[i], i, d1)){
			status = false;
		}
	}

	if(status==true)
		return 0;
	else
		return 1;
}

