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
#ifndef DATAASSIMILATION_H
#define DATAASSIMILATION_H

#include <meteoio/MeteoIO.h>

class DataAssimilation;
class SnowpackInterface;
class EnergyBalance;

#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/ebalance/EnergyBalance.h>

class DataAssimilation{
	public:

		DataAssimilation(mio::IOManager &myio);
		DataAssimilation(const DataAssimilation&);
		~DataAssimilation();
		
		DataAssimilation& operator=(const DataAssimilation&); ///<Assignement operator, required because of pointer member
		void Destroy();

		void setSnowPack(SnowpackInterface &mysnowpack);
		void SetEnergyBalance(EnergyBalance &myeb);

		void Compute(mio::Date julian);
		void Initialize();
		std::string getGridsRequirements() const;

	private:
		mio::Grid2DObject data;
		mio::IOManager* io;
		SnowpackInterface* snowpack;
		EnergyBalance* eb;
};

#endif // DATAASSIMILATION_H
