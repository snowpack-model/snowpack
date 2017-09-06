# Description: Generates textfiles with meteorological input data for SNOWPACK from RACMO NetCDF data
#
# Source: RACMO2.3 NetCDF-files (vertically sliced, 1 to 78)
#
# Author: Christian Steger, June 2016

# Load modules
import os
import sys
import platform
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"], # Desktop-Mac
             ["Christians-MacBook-Pro-2.local"],                                         # Macbook
             ["cca-login1", "cca-login3"]]                                               # ECMWF

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]):
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = "/Users/steger/Data/Greenland/RACMO2.3/FDM_INPUT/"
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/"   
elif (platform.node() in platf_nam[1]):  
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in = "/Volumes/IMAU/Greenland/RACMO2.3/FDM_INPUT/"
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/"
elif (platform.node() in platf_nam[2]):
    path_in = "/scratch/ms/nl/rucs/FDM_Data/INPUT/"
    path_out = "/scratch/ms/nl/rucs/SP_Data/INPUT/"
else:
    sys.exit("Unknown platform")

########################################################################################################################
# Settings
########################################################################################################################

# Path and name of file with indices (rlat, rlon; on RACMO2.3 grid)
#path_file_ind = root_IMAU + "FDM_SP_runs/ZGRN_V5_all_ind.txt"
path_file_ind = "/home/ms/nl/rucs/FDM_SP_runs/ZGRN_V5_all_ind.txt"

fn_part_in = "_ZGRN11_1960_2014_3H_p" # part of input NetCDF-file (without variable, slice number and ending)
fn_part_out = "ZGRN_V5_all"
# number of first slice must start with 1

# General
ph_thresh = 0.01 # treshold below all precipitation is considered snowfall [-]
tskin_melt_thresh = 268.15 # surface temperature threshold below snow melt is set to 0. [K]

# Meteorological RACMO variables
var_names = ["evap", "ff10m", "lwsd", "precip", "psurf", "q2m", "sndiv", "snowfall", "snowmelt", "swsd", "swsu",
             "t2m", "tskin"]
# units: mass fluxes [kg m-2 s-1], wind [m s-1], radiations [W m-2], pressure [Pa], specific humidity [kg kg-1],
#        temperatures [K]

# Physical Constants
t_wtp = 273.16 # triple point water (temperature) [K]
p_wtp = 611.73 # triple point water (pressure) [Pa]

# Further constants
rho_hn_max = 470. # maximum density of surface snow [kg m-3]

########################################################################################################################
# Load general data, check if necessary file exist
########################################################################################################################

# Read text file with indices of grid boxes on RACMO2.3 domain
indices = np.loadtxt(path_file_ind, usecols = [0, 1], dtype = int) # (column 0: rlat, column 1: rlon)

# Load geographical coordinates and topography
ncfile = Dataset(path_in + "RACMO23_masks_ZGRN11.nc", "r")
lon = ncfile.variables["lon"][:]
lat = ncfile.variables["lat"][:]
topo = ncfile.variables["topography"][:]
ncfile.close()

# Loop over all indices to check if all files exist
file_numbers = (indices[:,0] / 4) + 1
for k in np.unique(file_numbers):
    for m in range(0, len(var_names)):
         if not os.path.isfile(path_in + "XGRN11_3hourly/" + var_names[m] + fn_part_in + str(k) + ".nc"):
             sys.exit("Error: File " + var_names[m] + fn_part_in + str(k) + ".nc does not exist")
 
########################################################################################################################
# Convert RACMO NetCDF to SNOWPACK text-files
########################################################################################################################    

# Get time axis of first slice and variable (equal for all locations and variables)
ncfile = Dataset(path_in + "XGRN11_3hourly/" + var_names[0] + fn_part_in + str(24) + ".nc", "r")
time_axis = num2date(ncfile.variables["time"][:], units = ncfile.variables["time"].units, 
                     calendar = ncfile.variables["time"].calendar)
ncfile.close()
                                            
# Loop over all locations
data = np.empty((len(time_axis), len(var_names)), dtype = np.float32)
for k in range(indices.shape[0]):

    # Calculate necessary NetCDF-slice and sub-index
    file_numb = (indices[k,0] / 4) + 1
    ind_sub = np.mod(indices[k,0] , 4)

    # Get data for different variables
    for m in range(len(var_names)):
        ncfile = Dataset(path_in + "XGRN11_3hourly/" + var_names[m] + fn_part_in + \
                         str(file_numb) + ".nc", "r")
        data[:,m] = ncfile.variables[var_names[m]][:,ind_sub,indices[k,1]]
        ncfile.close()
    
    # Convert units [kg m-2 s-1] to [kg m-2 h-1]
    data[:,var_names.index("evap")] *= 3600.
    data[:,var_names.index("precip")] *= 3600.
    data[:,var_names.index("sndiv")] *= -3600. # (RACMO drift: removal (+), depos. (-))
    data[:,var_names.index("snowfall")] *= 3600.  
    data[:,var_names.index("snowmelt")] *= 3600.      
                                    
    # Calculate ground temperature (set equal to mean surface temperature)
    TSG = np.mean(data[:,var_names.index("tskin")])

    # Calculate precipitation phase fraction (0 is fully solid; 1 is fully liquid)
    PSUM_PH = np.zeros(len(time_axis), dtype = np.float32)
    ind_rain = np.logical_and.reduce((data[:,var_names.index("precip")] > data[:,var_names.index("snowfall")],
                                      data[:,var_names.index("precip")] > 0.,
                                      data[:,var_names.index("t2m")] > 273.15)) # time steps with rain
    # temperature criteria necessary -> time step in winter with tiny amounts of precipitation occur ->
    # minor issue: if snow drift occurs simultaneously, it is partially added as rain if PSUM_PH > 0.                                 
    PSUM_PH[ind_rain] = 1. - data[ind_rain,var_names.index("snowfall")] / data[ind_rain,var_names.index("precip")]
    PSUM_PH[ph_thresh > PSUM_PH] = 0.
  
    # Calculate relative humidty (over water or ice; depending on air temperature)
    T = data[:,var_names.index("t2m")] # air temperature [K]
    e_s = np.empty(len(time_axis), dtype = np.float32) # saturation vapor pressure [Pa]
    e_s[T < t_wtp] = p_wtp * np.exp(21.8745584 *  (T[T < t_wtp] - t_wtp) / (T[T < t_wtp] - 7.66)) # over ice
    e_s[T >= t_wtp] = p_wtp * np.exp(17.2693882 *  (T[T >= t_wtp] - t_wtp) / (T[T >= t_wtp] - 35.86)) # over water
    rh = data[:,var_names.index("q2m")] / (0.622 * e_s / data[:,var_names.index("psurf")])
    # Used approximations: w = q, w_s = 0.622 * e_s / p_surf
    # Sources: python script Stef ((August-Roche-)Magnus formula)
  
    # Calculate new snow density (firn model parameterization) (averaged calculation over all years)
    rho_hn = np.empty(len(time_axis), dtype = np.float32)
    tskin_mean = np.mean(data[:,var_names.index("tskin")]) # [K]
    ff10m_mean = np.mean(data[:,var_names.index("ff10m")]) # [m s-1]
    accum_mean = (np.mean(data[:,var_names.index("precip")]) + \
                  np.mean(data[:,var_names.index("evap")]) + \
                  np.mean(data[:,var_names.index("sndiv")])) * (24. * 365.25) # [kg m-2 h-1] to [mm a-1]  
    #rho_hn[:] = -151.94 + 1.4266 * (73.6 + 1.06 * tskin_mean + 0.0669 * accum_mean + 4.77 * ff10m_mean) # old param.
    rho_hn[:] = 481.0 + 4.834 * (tskin_mean - 273.15) # new param. (Peter's paper)
    rho_hn[rho_hn > rho_hn_max] = rho_hn_max

    # Set reflected shortwave radiation to positive values
    data[:,var_names.index("swsu")] *= (-1)

    # Add snow drift deposition to precipitation
    ind_depos = (data[:,var_names.index("sndiv")] > 0.) # time steps with snow drift deposition
    data[ind_depos,var_names.index("precip")] += data[ind_depos,var_names.index("sndiv")]
    data[ind_depos,var_names.index("sndiv")] = 0. # sndiv is now only negative

    # Physical restrictions
    rh = np.clip(rh, 0., 1.)
    PSUM_PH = np.clip(PSUM_PH, 0., 1.) # precipitation phase  
    data[data[:,var_names.index("swsd")] < 0.,var_names.index("swsd")] = 0.
    data[data[:,var_names.index("swsu")] < 0.,var_names.index("swsu")] = 0.
    ind = (data[:,var_names.index("swsu")] > data[:,var_names.index("swsd")])
    data[ind,var_names.index("swsu")] = data[ind,var_names.index("swsd")] # shortwave out can't exceed shortwave in
    data[data[:,var_names.index("lwsd")] < 0.,var_names.index("lwsd")] = 0.
    data[data[:,var_names.index("precip")] < 0.,var_names.index("precip")] = 0.       
    data[data[:,var_names.index("snowmelt")] < 0.,var_names.index("snowmelt")] = 0.
    data[data[:,var_names.index("tskin")] < tskin_melt_thresh,var_names.index("snowmelt")] = 0.   
    # Enforce consistency further between skin temperature and melt?
    
    # Miscellaneous
    data[:,var_names.index("swsu")] = np.abs(data[:,var_names.index("swsu")]) # set possible -0.0 values to 0.0
    #data[:,var_names.index("sndiv")] = 0.0 # set snow drift to zero
                                                                                                            
    # Calculate mean of certain parameters
    TA_MEAN = np.mean(data[:,var_names.index("t2m")]) - 273.15 # [degC]
    TSS_MEAN = np.mean(data[:,var_names.index("tskin")]) - 273.15 # [degC]
    VW_MEAN = np.mean(data[:,var_names.index("ff10m")]) # [m s-1]
    SUBLI_MEAN = np.mean(data[:,var_names.index("evap")]) # [kg m-2 h-1]
    PINT_MEAN = np.mean(data[:,var_names.index("precip")]) # [kg m-2 h-1]
    SNOWD_MEAN = np.mean(data[:,var_names.index("sndiv")]) # [kg m-2 h-1]
    SMELT_MEAN = np.mean(data[:,var_names.index("snowmelt")]) # [kg m-2 h-1]
    
    # Write text file in format .met
    fn_out = fn_part_out + "_" + str(k + 1) + ".smet"
    current_date = dt.datetime.now().strftime("%Y-%m-%d")
    file = open(path_out + "XGRN11_3hourly/" + fn_out,"w")
    file.write("SMET 1.1 ASCII\n")
    file.write("[HEADER]\n")
    file.write("station_id       = temp\n") # determines the output file name
    file.write("latitude         = " + "%.4f" % lat[indices[k,0],indices[k,1]] + "\n")
    file.write("longitude        = " + "%.4f" % lon[indices[k,0],indices[k,1]] + "\n")
    file.write("altitude         = " + "%.2f" % topo[indices[k,0],indices[k,1]] + "\n")
    file.write("nodata           = -999\n")
    file.write("source           = IMAU; CSteger, " + current_date + "\n")
    file.write("comment          = AirT [degC]: " + "%.3f" % TA_MEAN + 
                                ", SurfT [degC]: " +  "%.3f" % TSS_MEAN + 
                                ", Wind [m s-1]: " + "%.3f" %  VW_MEAN +
                                ", Prec (+ pos. Snowdrift) [kg m-2 h-1]: " + "%.5f" % PINT_MEAN +
                                ", Melt [kg m-2 h-1]: " + "%.5f" % SMELT_MEAN +
                                ", Snowdrift (only neg.) [kg m-2 h-1]: " + "%.5f" %  SNOWD_MEAN +
                                ", Sublim [kg m-2 h-1]: " + "%.5f" %  SUBLI_MEAN + "\n")                
    file.write("fields           = timestamp TA RH TSG TSS VW OSWR ISWR ILWR PINT PSUM_PH SMELT SNOWD SUBLI RHO_HN\n")
    file.write("[DATA]\n")
    for m in range(0, len(time_axis)):
        file.write(time_axis[m].isoformat()[:-3] + "   " +                    
                   "%7.3f" % data[m,var_names.index("t2m")] + "   " +       # [K]
                   "%5.3f" % rh[m] + "   " +      # [-]
                   "%7.3f" % TSG + "   " +                                  # [K]
                   "%7.3f" % data[m,var_names.index("tskin")] + "   " +     # [K]
                   "%6.3f" % data[m,var_names.index("ff10m")] + "   " +     # [m s-1]
                   "%8.3f" % data[m,var_names.index("swsu")] + "   " +      # [W m-2]
                   "%8.3f" % data[m,var_names.index("swsd")] + "   " +      # [W m-2]
                   "%8.3f" % data[m,var_names.index("lwsd")] + "   " +      # [W m-2]
                   "%9.5f" % data[m,var_names.index("precip")] + "   " +    # [kg m-2 h-1]
                   "%5.3f" % PSUM_PH[m] + "   " +                           # [-]
                   "%9.5f" % data[m,var_names.index("snowmelt")] + "   " +  # [kg m-2 h-1]
                   "%9.5f" % data[m,var_names.index("sndiv")]  + "   " +    # [kg m-2 h-1] (always < 0.)
                   "%9.5f" % data[m,var_names.index("evap")] + "   " +      # [kg m-2 h-1]
                   "%7.3f" % rho_hn[m] + "\n")      # [kg m-3]            
    file.close()
    
    print "Input file " + fn_out + " created"