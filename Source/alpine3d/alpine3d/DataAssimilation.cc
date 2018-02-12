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
#include <alpine3d/DataAssimilation.h>

using namespace std;

DataAssimilation::DataAssimilation(mio::IOManager& myio)
                : data(), io(&myio), snowpack(NULL), eb(NULL)
{
  Initialize();
  printf("DataAssimilation created\n");
}

DataAssimilation::~DataAssimilation()  {
	Destroy();
}

DataAssimilation& DataAssimilation::operator=(const DataAssimilation& source) {
	if (this != &source) {
		data = source.data;
		io = source.io;
		snowpack = source.snowpack;
		eb = source.eb;
	}
	return *this;
}

void DataAssimilation::Destroy() {}

std::string DataAssimilation::getGridsRequirements() const
{
	return "";
}

void DataAssimilation::Initialize()
{
	if (io==NULL){
		printf("No Input object. Ignore initialization\n");
		return;
	}
}

void DataAssimilation::setSnowPack(SnowpackInterface& mysnowpack)
{
  snowpack=&mysnowpack;
}

void DataAssimilation::SetEnergyBalance(EnergyBalance& myeb)
{
  eb=&myeb;
}


void DataAssimilation::Compute(mio::Date julian)
{

	printf("Try to read assimilation data...\n");
	try {
		io->readAssimilationData(julian,data);
		cout << "Size of DA-fields: nx=" << data.getNx() << ", ny=" << data.getNy() << "\n";

		if ( snowpack!=NULL ) {
			snowpack->assimilate(data, julian);
			cout <<"Data assimilation done at time " << julian.toString(mio::Date::ISO) << "\n";
		}
	}

	//in fact, GetAssimilationData throws a string exception, but catch anything here
	catch (...) {
		cout << "No assimilation data at time " << julian.toString(mio::Date::ISO) << "\n";
	}
}
