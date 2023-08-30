# Description: Generates configuration files for SNOWPACK - ECMWF version
#
# Author: Christian Steger, July 2016

# Load modules
import os
import sys
import platform
import numpy as np

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"], # Desktop-Mac
             ["Christians-MacBook-Pro-2.local", "ChristiansMBP2.home"],                  # Macbook
             ["cca-login1", "cca-login3", "ccb-login3"]]                                 # ECMWF

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]):
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = root_IMAU + "FDM_SP_runs/ECMWF_DATA/"
    path_out = root_IMAU + "FDM_SP_runs/ECMWF_DATA/cfgfiles/"
elif (platform.node() in platf_nam[1]):  
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = root_IMAU + "FDM_SP_runs/ECMWF_DATA/"
    path_out = root_IMAU + "FDM_SP_runs/ECMWF_DATA/cfgfiles/"
elif (platform.node() in platf_nam[2]):
    path_in = "/home/ms/nl/rucs/Snowpack/DATA/"
    path_out = "/scratch/ms/nl/rucs/SP_Data/Run_temp/cfgfiles/"
else:
    sys.exit("Unknown platform")

# Read command line arguments
point_file = sys.argv[1] # text file with selection of points to run (e.g ZGRN_V5_all_sel_points_cca.txt)

########################################################################################################################
# Parameters with different values for spin-up and actual run
########################################################################################################################

# Spin-up, actual run
prof_days_between = ["10000.0", "28.0"]
ts_days_between = ["10000.0", "2.0"]
config_fn = ["config_spinup", "config"]

########################################################################################################################
# Write configuration files
########################################################################################################################

points = np.loadtxt(path_in + point_file, dtype = int)

for k in points:
    for m in range(len(config_fn)):
        file = open(path_out + config_fn[m] + "_" + str(k) + ".ini","w")
        file.write("[General]\n") ################################################################# [General] ##########
        file.write("BUFFER_SIZE = 370\n")
        file.write("BUFF_BEFORE = 1.5\n") 
        file.write("\n")    
        file.write("[Input]\n") ##################################################################### [Input] ##########
        file.write("COORDSYS = CH1903\n")
        file.write("ISWR_IS_NET = false\n")
        file.write("METEO = SMET\n")
        file.write("METEOPATH = /scratch/ms/nl/rucs/SP_Data/INPUT/XGRN11_3hourly\n") ##### set before run #####
        file.write("SNOW = SMET\n")
        file.write("SNOWFILE1 = temp\n")    
        file.write("SNOWPATH = /scratch/ms/nl/rucs/SP_Data/Run_temp/point_" + str(k) + "/input\n")
        file.write("STATION1 = ZGRN_V5_all_" + str(k) + "\n") ##### set before run #####
        file.write("TIME_ZONE = 1\n")
        file.write("\n")
        file.write("[Output]\n") ################################################################### [Output] ##########
        file.write("BACKUP_DAYS_BETWEEN = 100000.0\n")
        file.write("CLASSIFY_PROFILE = false\n")
        file.write("COORDSYS = CH1903\n")
        file.write("CUMSUM_MASS = false\n")
        file.write("FIRST_BACKUP = 100000.0\n")
        file.write("HARDNESS_IN_NEWTON = false\n")
        file.write("METEO = SMET\n")
        file.write("METEOPATH = /scratch/ms/nl/rucs/SP_Data/Run_temp/point_" + str(k) + "/output\n")
        file.write("OUT_CANOPY = false\n")
        file.write("OUT_HAZ = true\n")
        file.write("OUT_HEAT = true\n")
        file.write("OUT_LW = true\n")
        file.write("OUT_MASS = true\n")
        file.write("OUT_METEO = true\n")
        file.write("OUT_SOILEB = false\n")
        file.write("OUT_STAB = false\n")
        file.write("OUT_SW = true\n")
        file.write("OUT_T = true\n")
        file.write("PRECIP_RATES = true\n")
        file.write("PROFILE_FORMAT = PRO\n")
        file.write("PROF_DAYS_BETWEEN = " + prof_days_between[m] + "\n")
        file.write("PROF_START = 0.0\n")
        file.write("PROF_WRITE = true\n")
        file.write("SNOW = SMET\n")
    	file.write("TIME_ZONE = 1\n")
        file.write("TS_DAYS_BETWEEN = " + ts_days_between[m] + "\n")
        file.write("TS_START = 0.0\n")
        file.write("TS_WRITE = true\n")
        file.write("\n")
        file.write("[Snowpack]\n") ############################################################### [Snowpack] ##########
        file.write("ATMOSPHERIC_STABILITY = MONIN_OBUKHOV\n")
        file.write("CALCULATION_STEP_LENGTH = 30.0\n")
        file.write("CANOPY = false\n")
        file.write("CHANGE_BC = true\n")
        file.write("ENFORCE_MEASURED_SNOW_HEIGHTS = false\n")
        file.write("GEO_HEAT = 0.0\n")
        file.write("HEIGHT_OF_METEO_VALUES = 2.0\n")
        file.write("HEIGHT_OF_WIND_VALUE = 10.0\n")
        file.write("MEAS_TSS = true\n")
        file.write("ROUGHNESS_LENGTH = 0.002\n")
        file.write("SNP_SOIL = false\n")
        file.write("SOIL_FLUX = true\n")
        file.write("SW_MODE = BOTH\n")
        file.write("THRESH_CHANGE_BC = -1.0\n")
        file.write("\n")
        file.write("[SnowpackAdvanced]\n") ############################################### [SnowpackAdvanced] ##########
        file.write("ALBEDO_PARAMETERIZATION = LEHNING_2\n")
        file.write("ALPINE3D = false\n")
        file.write("COMBINE_ELEMENTS = true\n")
        file.write("DETECT_GRASS = false\n")
        file.write("FORCE_RH_WATER = true\n")
        file.write("FORCE_SW_MODE = false\n")
        file.write("FORCING = MASSBAL\n")
        file.write("HEIGHT_NEW_ELEM = 0.02\n")
        file.write("HN_DENSITY = MEASURED\n")
        file.write("HN_DENSITY_FIXEDVALUE = 250.0\n")
        file.write("HN_DENSITY_PARAMETERIZATION = LEHNING_NEW\n")
        file.write("HOAR_DENSITY_BURIED = 125.0\n")
        file.write("HOAR_DENSITY_SURF = 100.0\n")
        file.write("HOAR_MIN_SIZE_BURIED = 2.0\n")
        file.write("HOAR_MIN_SIZE_SURF = 0.5\n")
        file.write("HOAR_THRESH_RH = 0.97\n")
        file.write("HOAR_THRESH_TA = 1.2\n")
        file.write("HOAR_THRESH_VW = 3.5\n")
        file.write("JAM = false\n")
        file.write("LB_COND_WATERFLUX = FREEDRAINAGE\n")
        file.write("MEAS_INCOMING_LONGWAVE = false\n")
        file.write("METAMORPHISM_MODEL = DEFAULT\n")
        file.write("MINIMUM_L_ELEMENT = 0.01\n")
        file.write("MIN_DEPTH_SUBSURF = 0.07\n")
        file.write("NEW_SNOW_GRAIN_SIZE = 0.2\n")
        file.write("PREVAILING_WIND_DIR = 0.0\n")
        file.write("REDUCE_N_ELEMENTS = true\n")
        file.write("SNOW_ALBEDO = PARAMETERIZED\n")
        file.write("SNOW_EROSION = true\n")
        file.write("SNOW_REDISTRIBUTION = false\n")
        file.write("SOOT_PPMV = 0.0\n")
        file.write("SW_ABSORPTION_SCHEME = MULTI_BAND\n")
        file.write("THRESH_DTEMP_AIR_SNOW = 25.0\n")
        file.write("THRESH_RAIN = 1.2\n")
        file.write("THRESH_RH = 0.0\n")
        file.write("T_CRAZY_MAX = 340.0\n")
        file.write("T_CRAZY_MIN = 210.0\n")
        file.write("VARIANT = POLAR\n")
        file.write("VISCOSITY_MODEL = DEFAULT\n")
        file.write("WATERTRANSPORTMODEL_SNOW = BUCKET\n")
        file.write("WATERTRANSPORTMODEL_SOIL = BUCKET\n")
        file.write("WATER_LAYER = false\n")
        file.write("WIND_SCALING_FACTOR = 1.0\n")
        file.write("\n")
        file.write("[Filters]\n") ################################################################# [Filters] ##########
        file.write("\n")
        file.write("[Interpolations1D]\n") ############################################### [Interpolations1D] ##########
        file.write("WINDOW_SIZE = 86400\n")
        file.write("PSUM::ACCUMULATE = 1800\n")
        file.write("PSUM::RESAMPLE = accumulate\n")
        file.write("\n")
        file.write("[GENERATORS]\n") ########################################################### [GENERATORS] ##########
        file.close()
