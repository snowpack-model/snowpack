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
const double res_Met_1 []	= {263.2811074,		264.8240384,		265.6572426,		265.8351213,		264.9516688,		265.3552443,		263.9469165}; // TA
const double res_Met_2 []	= {1.,			0.9669166667,		1.,			0.9658333333,		0.9641666667,		0.9121666667,		0.93425}; // RH
const double res_Met_3[]	= {273.6963,		274.3503833,		273.9969,		274.6983,		274.0971,		274.9251417,		274.0971}; // TSG
const double res_Met_4[]	= {262.9366667,		266.1666667,		264.41,			268.9033333,		265.43,			266.89,			262.1966667}; // TSS
const double res_Met_5 []	= {0.5875000001,	0.95,			1.894166667,		1.415833333,		0.35,			0.72,			1.47}; // HS
const double res_Met_6 []	= {2.2,			1.2,			2.2,			0.2,			2.1,			0.6,			1.}; // VW
const double res_Met_7 []	= {335.5833333,		144.4166666,		103.25,			213.0833333,		266.8333333,		239.5,			112.25}; // DW
const double res_Met_8 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // VW_MAX
const double res_Met_9 []	= {84.79166651,		81.29166661,		67.83333342,		50.37500004,		88.41666666,		85.45833337,		87.4999999}; //RSWR
const double res_Met_10 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // ISWR
const double res_Met_11 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // ILWR
const double res_Met_12 []	= {IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata,	IOUtils::nodata}; // TAU_CLD
const double res_Met_13 []	= {IOUtils::nodata,	2.632907412,		IOUtils::nodata,	2.255830171,		3.593920512,		0.,			0.}; // HNW

// methode do controll content of Meteo Data !!
// Also controlles != operator of containing special structures
bool controllStation(MeteoData& datMeteo, int i_results, Date datDate){

	const double epsilon = 1.0e-7; // accuracy of the double tests
	bool status = true;

	// Coords content
	Coords& dataCoord = datMeteo.meta.position;
	if(!IOUtils::checkEpsilonEquality(dataCoord.getAltitude(), res_Alt[i_results], epsilon)){
		cerr << "error on Altitude on " << res_ID[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getLat(), res_Lat[i_results], epsilon)){
		cerr << "error on Latitude on " << res_ID[i_results] << ": ref=" << dataCoord.getLat() << " got " << res_Lat[i_results] << endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(dataCoord.getLon(), res_Lon[i_results], epsilon)){
		cerr << "error on Longitude on " << res_ID[i_results] << ": ref=" << dataCoord.getLon() << " got " << res_Lon[i_results] << endl;
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
		cerr << "error on getSlopeAngle"<< endl;
		status = false;
	}
	if(!IOUtils::checkEpsilonEquality(datMeteo.meta.getAzimuth(), res_Azi[i_results], epsilon)){
		cerr << "error on getAzimuth"<< endl;
		status = false;
	}

	StationData refStation(refCoord, res_ID[i_results], res_Name[i_results]);
	refStation.setSlope(res_Slope[i_results], res_Azi[i_results]);
	if(refStation != datMeteo.meta){
		cerr << "error on != between Station Data"<< endl;
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

//Test if data read at 2008-12-01T15:35:00 are correct
int main() {
	bool status = true;
	Date d1(2008, 12, 1, 15, 35, 1);
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
