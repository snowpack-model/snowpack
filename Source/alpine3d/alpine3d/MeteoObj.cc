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
#include <alpine3d/MeteoObj.h>

using namespace mio;
using namespace std;

MeteoObj::MeteoObj(const mio::Config& in_config, const mio::DEMObject& in_dem)
                   : timer(), config(in_config), io(in_config), dem(in_dem),
                     ta(in_dem, IOUtils::nodata), rh(in_dem, IOUtils::nodata), psum(in_dem, IOUtils::nodata), 
                     psum_ph(in_dem, IOUtils::nodata), vw(in_dem, IOUtils::nodata), p(in_dem, IOUtils::nodata), ilwr(in_dem, IOUtils::nodata),
                     sum_ta(), sum_rh(), sum_rh_psum(), sum_psum(), sum_psum_ph(), sum_vw(), sum_ilwr(),
                     vecMeteo(), date(), glaciers(NULL), count_sums(0), count_precip(0), skipWind(false) {}

MeteoObj::~MeteoObj() 
{
	if (glaciers!=NULL) delete glaciers;
}

void MeteoObj::setSkipWind(const bool& i_skipWind) {
	skipWind = i_skipWind;
}

void MeteoObj::prepare(const mio::Date& in_date)
{
	if (!MPIControl::instance().master())  // Only master reads data
		return;

	date = in_date;
	getMeteo(date);
}

void MeteoObj::get(const mio::Date& in_date, mio::Grid2DObject& out_ta, mio::Grid2DObject& out_rh, mio::Grid2DObject& out_psum,
                   mio::Grid2DObject& out_psum_ph, mio::Grid2DObject& out_vw, mio::Grid2DObject& out_p, mio::Grid2DObject& out_ilwr)
{
	timer.restart(); //this method is called first, so we initiate the timing here
	
	if (MPIControl::instance().master()) {
		if (date.isUndef()) {
			date = in_date;
			getMeteo(date); //it will throw an exception if something goes wrong
		}
		if (in_date != date) {
			cerr << "[w] Meteo data was prepared for " << date.toString(Date::ISO);
			cerr << ", requested for " << in_date.toString(Date::ISO) << ", this is not optimal...\n";
			date = in_date;
			getMeteo(date); //it will throw an exception if something goes wrong
		}
	}
	
	//this acts as a barrier and forces MPI synchronization
	MPIControl::instance().broadcast(ta);
	MPIControl::instance().broadcast(rh);
	MPIControl::instance().broadcast(psum);
	MPIControl::instance().broadcast(psum_ph);
	MPIControl::instance().broadcast(vw);
	MPIControl::instance().broadcast(p);
	MPIControl::instance().broadcast(ilwr);

	out_ta = ta;
	out_rh = rh;
	out_psum = psum;
	out_psum_ph = psum_ph;
	out_vw = vw;
	out_p = p;
	out_ilwr = ilwr;
	timer.stop();
}

void MeteoObj::get(const mio::Date& in_date, std::vector<mio::MeteoData>& o_vecMeteo)
{
	timer.start();
	
	if (MPIControl::instance().master()) {
		if (date.isUndef()) {
			date = in_date;
			getMeteo(date); //it will throw an exception if something goes wrong
		}
		if (in_date != date) {
			cerr << "[w] Meteo data was prepared for " << date.toString(Date::ISO);
			cerr << ", requested for " << in_date.toString(Date::ISO) << ", this is not optimal...\n";
			date = in_date;
			getMeteo(date); //it will throw an exception if something goes wrong
		}
	}

	//this acts as a barrier and forces MPI synchronization
	MPIControl::instance().broadcast(vecMeteo);

	o_vecMeteo = vecMeteo;
	timer.stop();
}

double MeteoObj::getTiming() const
{
	return timer.getElapsed();
}

void MeteoObj::checkInputsRequirements(std::vector<MeteoData>& vecData)
{
	//This function checks that the necessary input data are available for the current timestamp
	unsigned int nb_ta=0, nb_iswr=0, nb_rh=0, nb_ilwr=0;
	unsigned int nb_ilwr_ta=0;

	if (vecData.empty())
		throw IOException("Vector of input meteo data is empty!", AT);

	for (size_t ii=0; ii<vecData.size(); ii++) {
		if (vecData[ii](MeteoData::TA) != IOUtils::nodata) nb_ta++;
		
		if (vecData[ii](MeteoData::ISWR) != IOUtils::nodata) nb_iswr++;
		
		if (vecData[ii](MeteoData::RH) != IOUtils::nodata) nb_rh++;
		
		if (vecData[ii](MeteoData::ILWR) != IOUtils::nodata) {
			nb_ilwr++;
			//We need ILWR and TA at the same location (so that the emissivity can be computed)
			//But since we rely on having one meteo1D station in the code, we also need ISWR at this place
			//(since we use the location of the measurement for computations with ISWR)
			if (vecData[ii](MeteoData::TA) != IOUtils::nodata && vecData[ii](MeteoData::ISWR) != IOUtils::nodata) nb_ilwr_ta++;
		}
	}

	if ( nb_ta==0 || nb_iswr==0 || nb_rh==0 || nb_ilwr==0) {
		printf("nb(ta)=%d nb(iswr)=%d nb(rh)=%d nb(ilwr)=%d\n",nb_ta, nb_iswr, nb_rh, nb_ilwr);
                throw IOException("Not enough input meteo data on "+vecData[0].date.toString(Date::ISO), AT);
	}
	if ( nb_ilwr_ta==0 ) {
		cout << "[e] For emissivity calculation, at least one set of both TA and ILWR are needed at the same station!\n";
		throw IOException("Not enough input meteo data", AT);
	}
}

void MeteoObj::fillMeteoGrids(const Date& calcDate)
{
	//fill the meteo parameter grids (of the AlpineControl object) using the data from the stations
	try {
		io.getMeteoData(calcDate, dem, MeteoData::PSUM, psum);
		io.getMeteoData(calcDate, dem, MeteoData::PSUM_PH, psum_ph);
		io.getMeteoData(calcDate, dem, MeteoData::RH, rh);
		io.getMeteoData(calcDate, dem, MeteoData::TA, ta);
		if (!skipWind) io.getMeteoData(calcDate, dem, MeteoData::VW, vw);
		io.getMeteoData(calcDate, dem, MeteoData::P, p);
		io.getMeteoData(calcDate, dem, MeteoData::ILWR, ilwr);
		cout << "[i] 2D Interpolations done for " << calcDate.toString(Date::ISO) << "\n";
	} catch (long) {
		cout << "[e] at " << calcDate.toString(Date::ISO) << " Could not fill 2D meteo grids" << endl;
	} catch(std::exception& e) {
		cerr << e.what() << endl;
		throw;
	}
}

void MeteoObj::getMeteo(const Date& calcDate)
{
	//Note: in case of MPI simulation only master node is responsible for file I/O
	if (!MPIControl::instance().master()) return;

	// Collect the Meteo values at each stations
	io.getMeteoData(calcDate, vecMeteo);
	checkInputsRequirements(vecMeteo);

	// Now fill the 2D Meteo Fields. Keep in mind that snowdrift might edit these fields
	fillMeteoGrids(calcDate);
	
	cout << "[i] Success reading/preparing meteo data for date: " << calcDate.toString(Date::ISO) << endl;
}

void MeteoObj::checkLapseRate(const std::vector<mio::MeteoData>& i_vecMeteo, const mio::MeteoData::Parameters& param)
{
	std::vector<double> vecData, vecAltitudes;
	for (size_t ii=0; ii<i_vecMeteo.size(); ii++){
		const double& val = i_vecMeteo[ii](param);
		if (val != IOUtils::nodata){
			vecData.push_back( val );
			vecAltitudes.push_back( i_vecMeteo[ii].meta.position.getAltitude() );
		}
	}
	
	if (vecData.size()<2) return;
	if (param==MeteoData::PSUM) { //skip when there is no precip
		if (mio::Interpol2D::allZeroes(vecData)) return;
	}
	
	double A, B, R;
	std::string mesg;
	Interpol1D::NoisyLinRegression(vecAltitudes, vecData, A, B, R, mesg);
	const std::string date_str( i_vecMeteo[0].date.toString(Date::ISO) );
	const std::string param_str( MeteoData::getParameterName(param) );
	std::cout << "[check:Data_Lapse_Rate] " << date_str << " " << param_str << " " << std::fixed << std::setw(7) << std::setprecision(5) << A << " ";
	
	if (param==MeteoData::PSUM && A<-1e-3) { //when precip gradient is "too wrong", print values
		const std::streamsize prec = std::cout.precision();
		std::cout << "- ";
		for (size_t ii=0; ii<vecData.size(); ii++){
			std::cout << std::fixed << std::setw(4) << std::setprecision(2) << vecData[ii] << "@" << std::fixed << std::setw(6) << std::setprecision(1) << vecAltitudes[ii] << " ";
		}
		std::cout << std::setprecision(prec);
		std::cout .unsetf(ios_base::floatfield);
	}
	std::cout << "\n";
}

void MeteoObj::checkGridRange(const mio::Date& calcDate, const mio::Grid2DObject& grid, const mio::MeteoData::Parameters& param)
{
	if (param==MeteoData::RH) {
		const double min = grid.grid2D.getMin();
		if (min<=0.05)
			std::cout << "[check:Grids_Range_Check] " << calcDate.toString(Date::ISO) << " Rh_min=" << min << "\n";
	} else if (param==MeteoData::VW) {
		const double max = grid.grid2D.getMax();
		if (max>40.)
			std::cout << "[check:Grids_Range_Check] " << calcDate.toString(Date::ISO) << " VW_max=" << max << "\n";
	} else
		throw IOException("Parameter '"+MeteoData::getParameterName(param)+"' not supported here", AT); 
}

void MeteoObj::setGlacierMask(const Grid2DObject& glacierMask)
{
	if (glacierMask.grid2D.getCount()>0) { //at least one pixel is glaciated...
		glaciers = new Glaciers(config, dem);
		glaciers->setGlacierMap( glacierMask );
	}
}

//this should only be called when "--nocompute" was set. So we consider that 
//most of the other modules have NOT been called.
void MeteoObj::checkMeteoForcing(const mio::Date& calcDate)
{
	if (calcDate != date) {
		cerr << "[w] Meteo data was prepared for " << date.toString(Date::ISO);
		cerr << ", requested for " << calcDate.toString(Date::ISO) << ", this is not optimal...\n";
		date = calcDate;
		getMeteo(date); //it will throw an exception if something goes wrong
	}
	
	if (glaciers!=NULL) {
		glaciers->correctTemperatures( ta );
	}
	
	//produce monthly gridded sums
	int year, month, day, hour, minute, second;
	calcDate.getDate(year, month, day, hour, minute, second);
	const bool startOfMonth = (day==1 && hour==0 && minute==0 && second==0);
	if (startOfMonth || sum_ta.empty()) { //we need to (re-)initialize the sums
		if (startOfMonth && !sum_ta.empty()) {
			sum_ta /= static_cast<double>(count_sums);
			//sum_rh /= static_cast<double>(count_sums);
			if (count_precip>0) sum_rh_psum /= static_cast<double>(count_precip);
			//if (count_precip>0) sum_psum /= static_cast<double>(count_precip);
			if (count_precip>0) sum_psum_ph /= static_cast<double>(count_precip);
			sum_vw /= static_cast<double>(count_sums);
			sum_ilwr /= static_cast<double>(count_sums);
			
			io.write2DGrid(sum_ta, MeteoGrids::TA, calcDate);
			//io.write2DGrid(sum_rh, MeteoGrids::RH, calcDate+1./(3600.*24));
			io.write2DGrid(sum_vw, MeteoGrids::VW, calcDate);
			io.write2DGrid(sum_ilwr, MeteoGrids::ILWR, calcDate);
			if (count_precip>0) io.write2DGrid(sum_rh_psum, MeteoGrids::RH, calcDate);
			if (count_precip>0) io.write2DGrid(sum_psum, MeteoGrids::PSUM, calcDate);
			if (count_precip>0) io.write2DGrid(sum_psum_ph, MeteoGrids::PSUM_PH, calcDate);
		}
		
		count_sums=0;
		count_precip=0;
		sum_ta.set(dem, 0.);
		//sum_rh.set(dem, 0.);
		sum_rh_psum.set(dem, 0.);
		sum_psum.set(dem, 0.);
		sum_psum_ph.set(dem, 0.);
		sum_vw.set(dem, 0.);
		sum_ilwr.set(dem, 0.);
	}
	
	count_sums++;
	sum_ta += ta;
	//sum_rh += rh;
	sum_vw += vw;
	sum_ilwr += ilwr;
	if (psum.grid2D.getMin()>0.) { //there are some precip
		sum_psum += psum;
		sum_psum_ph += psum_ph;
		sum_rh_psum += rh;
		count_precip++;
	}
	
	//check range in the grids
	checkGridRange(calcDate, rh, MeteoData::RH);
	checkGridRange(calcDate, vw, MeteoData::VW);
	
	//now get and print the lapse rates
	vector<MeteoData> vecMd;
	io.getMeteoData(calcDate, vecMd);
	checkLapseRate(vecMd, MeteoData::TA);
	checkLapseRate(vecMd, MeteoData::RH);
	checkLapseRate(vecMd, MeteoData::VW);
	checkLapseRate(vecMd, MeteoData::PSUM);
}
