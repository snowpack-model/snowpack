#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the most basic example. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and
//it will retrieve the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** argv) {
	Date d1;
	std::vector<MeteoData> vecMeteo;

	Config cfg("io.ini");
	IOManager io(cfg);

	const double TZ = cfg.get("TIME_ZONE", "Input");
	IOUtils::convertString(d1,argv[1], TZ);
	//io.setProcessingLevel(IOManager::raw); //set the processing level: raw, filtered or resampled
	io.getMeteoData(d1, vecMeteo);

	std::cout << vecMeteo.size() << " stations with an average sampling rate of " << io.getAvgSamplingRate() << " or 1 point every " << 1./(io.getAvgSamplingRate()*60.+1e-12) << " minutes\n";
	//writing some data out in order to prove that it really worked!
	for (unsigned int ii=0; ii < vecMeteo.size(); ii++) {
		std::cout << "---------- Station: " << (ii+1) << " / " << vecMeteo.size() << std::endl;
		std::cout << vecMeteo[ii].toString() << std::endl;
	}
	
	return 0;
}
