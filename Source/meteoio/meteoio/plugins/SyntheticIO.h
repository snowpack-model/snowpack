// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SYNTHETICIO_H
#define SYNTHETICIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {
/**
 * @class Synthesizer
 * @brief Generator to produce synthetic data for the SynthIO plugin
 * @author Mathias Bavay
 * @date   2022-09-19
 */
class Synthesizer {
	public:
		/*Synthesizer(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs) {(void)vecArgs;}*/
		virtual ~Synthesizer() {}
		
		virtual double generate(const Date& dt) const {(void)dt; return IOUtils::nodata;}
};

//object factory for the Synthesizer class
class SynthFactory {
	public:
		static Synthesizer* getSynth(std::string type, const std::string& station, const StationData& sd, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ);
};


/**
 * @class SynthIO
 * @brief This plugin generate synthetic data
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2022-09-19
 */
class SynthIO : public IOInterface {
	public:
		SynthIO(const std::string& configfile);
		SynthIO(const SynthIO&);
		SynthIO(const Config& cfgreader);
		
		~SynthIO();

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void init();
		std::map< std::string, Synthesizer* > getSynthesizer(const std::string& stationRoot, const StationData &sd) const;
		
		const Config cfg;
		std::map< std::string, std::map< std::string, Synthesizer* > > mapSynthesizers; //map[station][parmname] to contain the Synthesizer
		std::vector<StationData> vecStations;
		Date dt_start, dt_end;
		double dt_step, TZ;
};

///////////////////////////////////////////////////////
//below derived classes of Synthesizer to implement specific algorithms (or generators)
class CST_Synth : public Synthesizer {
	public:
		CST_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs);
		virtual double generate(const Date& dt) const;
	private:
		double value;
};

class STEP_Synth : public Synthesizer {
	public:
		STEP_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ);
		virtual double generate(const Date& dt) const;
	private:
		Date dt_step;
		double value_before, value_after;
};

class RECT_Synth : public Synthesizer {
	public:
		RECT_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ);
		virtual double generate(const Date& dt) const;
	private:
		Date step_start, step_stop;
		double value, value_step;
};

class STDPRESS_Synth : public Synthesizer {
	public:
		STDPRESS_Synth(const StationData& sd);
		virtual double generate(const Date& dt) const;
	private:
		double altitude;
};

} //namespace
#endif
