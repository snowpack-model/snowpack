# Description: Generates SNOWPACK input files from RACMO data (ice blocks)
#
# Details: temperature profile is either constant or with seasonal variation
#
# Author: Christian Steger, June 2016

# Load modules
import os
import sys
import platform
import numpy as np
import datetime as dt
from netCDF4 import Dataset

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"], # Desktop-Mac
             ["Christians-MacBook-Pro-2.local"],                                         # Macbook
             ["cca-login1"]]                                                             # ECMWF

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]):
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = root_IMAU + "Data/Greenland/RACMO2.3/"
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/sno/"
elif (platform.node() in platf_nam[1]):  
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = root_IMAU + "Data/Greenland/RACMO2.3/"
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/sno/"
elif (platform.node() in platf_nam[2]):
    path_in = "/scratch/ms/nl/rucs/raw_Data/"
    path_out = "/scratch/ms/nl/rucs/SP_Data/INPUT/sno/"
else:
    sys.exit("Unknown platform")

########################################################################################################################
# Settings
########################################################################################################################

# Path and name of file with indices (rlat, rlon; on RACMO2.3 grid)
path_file_ind = root_IMAU + "FDM_SP_runs/ZGRN_V5_ALL_ind.txt"

fn_part_out = "ZGRN_V5_all"

# Ice block settings
thick_min = 30. # minimal thickness of initial ice block [m] (110 + 50 m for first spin-up)
elem_thickn = 1.0 # initial element thickness [m]
rg = 13.254 # grain size [mm]
rb = 2.895 # bond size [mm]

# Miscellaneous
temp_prof = "seas" # either const or seas (winter profile)
date_sp = dt.datetime(1960, 1, 1) # snow profile date

# Constants
dens_ice = 917.0 # ice density [kg m-3]
spec_heat_ice = 2100.0 # specific heat ice at 0 degC (same value as in SNOWPACK) [J kg-1 K-1]
conduct_ice = 2.2 # thermal conductivity of ice (same value as in SNOWPACK) [W m-1 K-1]

########################################################################################################################
# Load data
########################################################################################################################

# Read text file with indices of grid boxes on RACMO2.3 domain
indices = np.loadtxt(path_file_ind, usecols = [0, 1], dtype = int) # (column 0: rlat, column 1: rlon)

# Load geographical coordinates and topography
ncfile = Dataset(path_in + "RACMO23_masks_ZGRN11.nc", "r")
lon = ncfile.variables["lon"][:]
lat = ncfile.variables["lat"][:]
topo = ncfile.variables["topography"][:]
ncfile.close()

# Load necessary variables and create .sno input files (1960 - 2014; 55 years)
var_names = ["tskin"] # [K]
data = np.empty((55 * 12, len(var_names)), dtype = np.float32)
#for k in range(indices.shape[0]):
for k in [662]:
    
    # Load RACMO data (monthly)
    for m in range(len(var_names)):
        ncfile = Dataset(path_in + "Monthly/" + var_names[m] + ".BN_1958_2014.MM.nc", "r")
        data[:,m] = ncfile.variables[var_names[m]][24:,indices[k,0],indices[k,1]]
        ncfile.close()

    # Climatological values for reference period (1960 - 1979)
    tskin_ref = data[:(20 * 12),var_names.index("tskin")].mean() # [K]
    tskin_seas_ref = data[:(20 * 12),var_names.index("tskin")].reshape(20, 12).mean(axis = 0)
    tskin_amp_ref = (tskin_seas_ref.max() - tskin_seas_ref.min()) / 2. # [K]  
      
    # Profile temperature
    numb_elem = int(np.ceil(thick_min / elem_thickn))
    elem_temp = np.empty(numb_elem, dtype = np.float32) # (surface -> bottom)
    if (temp_prof == "const"):
        elem_temp.fill(tskin_ref)
    elif (temp_prof == "seas"):
        depth = np.linspace(elem_thickn / 2., (numb_elem * elem_thickn) - elem_thickn / 2, numb_elem)
        alpha = conduct_ice / (dens_ice * spec_heat_ice) # thermal diffusivity [m2 s-1]
        omega = 1. / (3600. * 24. * 365.25) # [s-1]
        phas_shift = 1.5 * np.pi # assumption: minimal temperature at beginning of January
        elem_temp.fill(tskin_ref)
        for m in range(numb_elem):
            elem_temp[m] += tskin_amp_ref * np.exp(-depth[m] * np.sqrt((np.pi * omega) / alpha)) * \
                            np.sin(-depth[m] * np.sqrt((np.pi * omega) / alpha) + phas_shift)
            # Source: Paterson 2010: The Physics of Glaciers, p. 401
        elem_temp = elem_temp.clip(max = 273.15)
    else:
        sys.exit("Unknown value for temp_prof")

    # Check temperature profile
    #plt.figure()
    #plt.plot(elem_temp, -depth)

    # Array with element age (age difference layers: 1 day)
    time_axis = [(dt.datetime(1960, 1, 1) - dt.timedelta(days = (m + 1))) for m in range(numb_elem)]

    # Write text file in format .met
    fn_out = fn_part_out + "_" + str(k + 1) + ".sno"
    file = open(path_out + fn_out,"w")
    file.write("SMET 1.1 ASCII\n")
    file.write("[HEADER]\n")
    file.write("station_id       = temp\n")
    file.write("station_name     = temp\n")
    file.write("latitude         = " + "%.4f" % lat[indices[k,0],indices[k,1]] + "\n")
    file.write("longitude        = " + "%.4f" % lon[indices[k,0],indices[k,1]] + "\n")
    file.write("altitude         = " + "%.2f" % topo[indices[k,0],indices[k,1]] + "\n")
    file.write("nodata           = -999\n")
    file.write("source           = IMAU; CSteger, " + dt.datetime.now().strftime("%Y-%m-%d") + "\n")
    file.write("ProfileDate      = " + dt.datetime(1960, 1, 1, 1).isoformat()[:-3] + "\n")
    file.write("HS_Last          = " + "%.6f" % thick_min + "\n") # [m]
    file.write("SlopeAngle       = 0.00\n")
    file.write("SlopeAzi         = 0.00\n")
    file.write("nSoilLayerData   = 0\n")
    file.write("nSnowLayerData   = " + str(numb_elem) + "\n")
    file.write("SoilAlbedo       = 0.09\n")
    file.write("BareSoil_z0      = 0.020\n")
    file.write("CanopyHeight     = 0.00\n")
    file.write("CanopyLeafAreaIndex = 0.00\n")
    file.write("CanopyDirectThroughfall = 1.00\n")
    file.write("WindScalingFactor = 1.00\n")
    file.write("ErosionLevel     = "  + str(numb_elem - 1) + "\n")
    file.write("TimeCountDeltaHS = 0.000000\n")
    file.write("fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo\n")
    file.write("[DATA]\n")
    for m in reversed(range(numb_elem)):
        file.write(time_axis[m].isoformat()[:-3] + "   " +
                    "%6.3f" % elem_thickn + "   " +    # layer-thickness [m]
                    "%7.3f" % elem_temp[m] + "   " +   # temperature [K]
                    "%5.3f" % 1.0 + "   " +            # volume fraction (ice) [-]
                    "%5.3f" % 0.0 + "   " +            # volume fraction (water) [-]
                    "%5.3f" % 0.0 + "   " +            # volume fraction (air) [-]
                    "%5.3f" % 0.0 + "   " +            # volume fraction (soil) [-]
                    "%8.3f" % 2600.0 + "   " +         # soil dry dens. [kg m-3]
                    "%6.3f" % 2.6 + "   " +            # soil heat cond. [W m-1 K-1]
                    "%8.3f" % 1100.00 + "   " +        # soil heat capacity [J m-3 K-1]
                    "%9.3f" % rg + "   " +             # grain radius [mm]
                    "%9.3f" % rb + "   " +             # bond radius [mm]
                    "%5.3f" % 0.0 + "   " +            # dendricity (0: old snow) [-]
                    "%5.3f" % 1.0 + "   " +            # sphericity (1: rounded) [-]
                    str(0) + "   " +                   # grain marker [-] (set to 7 -> glacier ice ?) ##########################
                    str(0.0) + "   " +                 # Mass of surface hoar [kg m-2]
                    str(1) + "   " +                   # Numb. of fin. elem. in layer [-]
                    str(0.0) + "   " +                 # Stress rate [Pa s-1]
                    str(0.0) + "   " +  "\n")          # Keep track of metamorphism [-]
    file.close()
    
    print "Input file " + fn_out + " created"

# Comments
# - Grain marker: value same as for new snow -> Snowpack::setHydrometeorMicrostructure
# - Stress rate: value same as for new snow -> Snowpack::fillNewSnowElement
# - Keep track of metamorphism: value same as for new snow -> Snowpack::setHydrometeorMicrostructure
# - Influence of initial ice block in accumulation area on upper profile evolution -> only over heat flux (temperature)
#   -> heat conductivity of ice is constant (not dependent on micro-properties (e.g. grain and bond size))