# Description: Divide RACMO2.3 domain (ZGRN11; lat: 312 and lon: 306) in lat. slices (78 files with 4 grid boxes) 
#              with continuous time series (1960 - 2014)
#
# Important: Re-chunk RACMO2.3 input data first to considerably speed-up reading of these NetCDF-files
#            (-> file RACMO_rechunk.sh)
#
# Author: Christian Steger, October 2015

# Load modules
import sys
import platform
import numpy as np
from netCDF4 import Dataset, num2date, date2num
import time

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"],        # Desktop-Mac
             ["Christians-MacBook-Pro-2.local", "ChristiansMBP2.home", "christiansmbp2.home"],  # Macbook
             ["cca-login1", "cca-login3"]]                                                      # ECMWF

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]):
    path_in = "/Users/steger/Data/Greenland/RACMO2.3/"
    path_out = "/Users/steger/Data/Greenland/RACMO2.3/FDM_INPUT/"
elif (platform.node() in platf_nam[1]):
    path_in = "/Volumes/IMAU/Greenland/RACMO2.3/"
    path_out = "/Volumes/IMAU/Greenland/RACMO2.3/FDM_INPUT/"
elif (platform.node() in platf_nam[2]):
    path_in = "/scratch/ms/nl/rucs/raw_Data/"
    path_out = "/scratch/ms/nl/rucs/FDM_Data/INPUT/"
else:
    sys.exit("Unknown platform")

# Constant values
part_fn = ["1958", "1961", "1971", "1981", "1991", "2001", "2011"]
var_names = ["evap", "ff10m", "lwsd", "precip", "psurf", "q2m", "sndiv", "snowfall", "snowmelt",
             "swsd", "swsu", "t2m", "tskin"]                                     
glob_att_names = ["Conventions", "source", "Domain", "Experiment", "institution", "CreationDate", "comment", "title"]
fold = {3: "3hourly", 6: "6hourly", 24: "daily"} # folder for RACMO data

########################################################################################################################
# Settings
########################################################################################################################

dt_in = 3 # temporal resolution of input data (integer) [h]

########################################################################################################################
# Process data
########################################################################################################################

# Get time axis for entire interval and time steps per block
time_axis = np.array([], dtype = object)
dt_per_block = np.array([], dtype = np.int)
for k in range(0, len(part_fn)):
    ncfile = Dataset(path_in + fold.get(dt_in) + "/" + var_names[0] + ".KNMI-" + part_fn[k] + \
                     ".ZGRN11.BN_1958_2013.3H_rc.nc", mode = "r")
    time_axis = np.append(time_axis, num2date(ncfile.variables["time"][:], units = "days since 1950-01-01 00:00:00", 
                                              calendar = "standard"))
    dt_per_block = np.append(dt_per_block, ncfile.variables[var_names[0]].shape[0])                              
    ncfile.close()
# Get indices of first time step in 1960   
ind = np.where(np.asarray([time_axis[:dt_per_block[0]][k].year for k in range(0, dt_per_block[0])]) >= 1960)[0][0]
dt_per_block[0] -= ind # remove years 1958 / 1959
time_axis = time_axis[ind:]

# Create data array
data = np.empty((np.sum(dt_per_block), 4, 306), dtype = np.float32)
ind_time = np.append(0, np.cumsum(dt_per_block))

# Load data (1960 - 2014)
for k in var_names:
    
    # Get certain global and variable attributes
    ncfile = Dataset(path_in + fold.get(dt_in) + "/" + k + ".KNMI-" + part_fn[0] + ".ZGRN11.BN_1958_2013.3H_rc.nc", \
                     mode = "r")
    glob_att = [str(ncfile.getncattr(glob_att_names[m])) for m in range(0, len(glob_att_names))]
    data_units = str(ncfile.variables[k].units)
    cell_methods = str(ncfile.variables[k].cell_methods)
    ncfile.close()
    
    for m in range(10, 70): # (Greenland binary land mass (10, 70), entire Greenland domain (1, 79))

        start = time.time()

        # Calculate indices for slice
        ind_beg = (m - 1) * 4
        ind_end = m * 4

        # Get data
        for n in range(0, len(part_fn)):
            ncfile = Dataset(path_in + fold.get(dt_in) + "/" + k + ".KNMI-" + part_fn[n] + \
                             ".ZGRN11.BN_1958_2013.3H_rc.nc", mode = "r")
            if (n == 0):
                data[ind_time[n]:ind_time[n + 1],:,:] = ncfile.variables[k][ind:,0,ind_beg:ind_end,:]
            else:
                data[ind_time[n]:ind_time[n + 1],:,:] = ncfile.variables[k][:,0,ind_beg:ind_end,:]
            ncfile.close()

        # Save data (chunksize optimized for accessing time series of a single location)
        fn_out = k + "_ZGRN11_1960_2014_3H_p" + str(m) + ".nc"
        ncfile = Dataset(path_out + "XGRN11_" + fold.get(dt_in) + "/" + fn_out, mode = "w", format = "NETCDF4")
        [str(ncfile.setncattr(glob_att_names[m], glob_att[m])) for m in range(0, len(glob_att_names))]
        ncfile.createDimension(dimname = "time", size = (np.sum(dt_per_block)))
        ncfile.createDimension(dimname = "rlat", size = 4)
        ncfile.createDimension(dimname = "rlon", size = 306)
        nc_time = ncfile.createVariable(varname = "time", datatype = "f", dimensions = ("time"))
        nc_time.units = "days since 1950-01-01 00:00:00"
        nc_time.calendar = "standard"
        nc_time[:] = date2num(time_axis, units = nc_time.units, calendar = nc_time.calendar)
        nc_data = ncfile.createVariable(varname = k, datatype = "f", dimensions = ("time", "rlat", "rlon"),
                                        chunksizes = (np.sum(dt_per_block), 1, 1)) 
        nc_data[:] = data
        nc_data.units = data_units
        nc_data.cell_methods = cell_methods
        ncfile.close()
        
        end = time.time()
        
        print "File " + fn_out + " created (" + "%.2f" % (end - start) + " seconds)"