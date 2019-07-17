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
#ifndef ALPINE3D_H
#define ALPINE3D_H

#include <alpine3d/DataAssimilation.h>
//#include <alpine3d/Alpine3D.h>
#include <alpine3d/AlpineControl.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/OMPControl.h>
//#include <alpine3d/MainPage.h>
#include <alpine3d/MeteoObj.h>
//#include <alpine3d/runoff/prevah_runoff/fortran_and_c.h>
//#include <alpine3d/runoff/prevah_runoff/Runoff.h>
//#include <alpine3d/runoff/prevah_runoff/f77-uscore.h>
//#include <alpine3d/runoff/prevah_runoff/wasim.h>
//#include <alpine3d/runoff/prevah_runoff/ascii.h>
#include <alpine3d/runoff/Runoff.h>
//#include <alpine3d/snowdrift/SnowDriftParallel.h>
#include <alpine3d/snowdrift/checksum.h>
#include <alpine3d/snowdrift/SnowDrift.h>
#include <alpine3d/snowdrift/Cell.h>
#include <alpine3d/ebalance/EBStruct.h>
#include <alpine3d/ebalance/TerrainRadiationHelbig.h>
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/ebalance/RadiationField.h>
#include <alpine3d/ebalance/VFSymetricMatrix.h>
#include <alpine3d/ebalance/TerrainRadiationSimple.h>
#include <alpine3d/ebalance/ViewFactorsHelbig.h>
#include <alpine3d/ebalance/ViewFactors.h>
#include <alpine3d/ebalance/ViewFactorsCluster.h>
#include <alpine3d/ebalance/TerrainRadiation.h>
#include <alpine3d/ebalance/ViewFactorsSectors.h>
#include <alpine3d/ebalance/ViewFactorsAlgorithm.h>
#include <alpine3d/ebalance/TerrainRadiationAlgorithm.h>
//#include <alpine3d/ebalance/TerrainRadiationPETSc.h>
//#include <alpine3d/snowdrift_par/clock.h>
//#include <alpine3d/snowdrift_par/SnowDriftParallel.h>
//#include <alpine3d/snowdrift_par/DriftData.h>
//#include <alpine3d/snowdrift_par/checksum.h>
//#include <alpine3d/snowdrift_par/SnowDrift.h>
//#include <alpine3d/snowdrift_par/AziSlope.h>
#include <alpine3d/Glaciers.h>
#include <alpine3d/TechSnowA3D.h>
#include <alpine3d/SnowpackInterface.h>
#include <alpine3d/SnowpackInterfaceWorker.h>
#include <alpine3d/AlpineMain.h>

#endif
