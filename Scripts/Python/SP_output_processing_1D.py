# Description: Rearrangement (and interpolation) of 1D output data from SNOWPACK
#              (sublimation, evaporation, rainfall, melt)
#
# Author: Christian Steger, January 2017

# Load modules
import os
import sys
import platform
import glob
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import datetime as dt
from dateutil.relativedelta import relativedelta

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"],        # Desktop-Mac
             ["Christians-MacBook-Pro-2.local", "ChristiansMBP2.home", "christiansmbp2.home"],  # Macbook
             ["cca-login1", "cca-login3"]]                                                      # ECMWF

# Paths to folders (dependent on platform)
if ((platform.node() in platf_nam[0]) or (platform.node() in platf_nam[1])):
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_RACMO = root_IMAU + "Data/Greenland/RACMO2.3/"
    path_ind = root_IMAU + "FDM_SP_runs/"
    sys.path.append(root_IMAU + "Scripts/Python/Functions/")
elif (platform.node() in platf_nam[2]):
    path_RACMO = "/scratch/ms/nl/rucs/raw_Data/"
    path_ind = "/home/ms/nl/rucs/FDM_SP_runs/"
else:
    sys.exit("Unknown platform")

# Load required functions
from array_processing import consresamp

# Constants
density_ice = 917. # [kg m-3]

########################################################################################################################
# Settings
########################################################################################################################

# Process settings
interp_frame = 1 # frame for spatial interpolation (frame length: 2 * interp_frame + 1)
temp_interp = True # additional temporal interpolation (time step: quarter of a month)

# Output variables and units
var_names_out = ["sublim", "evap", "rainfall", "melt"]
units_out = ["kg m-2", "kg m-2", "kg m-2", "kg m-2"]

# Temporal interpolation: time axis (time step: quarter of a month)
temp = [(dt.datetime(1960, 1, 1) + relativedelta(months = k)) for k in range(55 * 12 + 1)]
temp = np.asarray([(k - dt.datetime(1960, 1, 1)).total_seconds() for k in temp], dtype = np.float32)
tb_out_num = np.interp(np.linspace(0., 660., (55 * 12 * 4 + 1)), np.arange(661.), temp)
tb_out = np.asarray([dt.datetime(1960, 1, 1) + dt.timedelta(seconds = k) for k in tb_out_num])
ta_out_num = tb_out_num[:-1] + np.diff(tb_out_num) / 2.
ta_out = np.asarray([dt.datetime(1960, 1, 1) + dt.timedelta(seconds = k) for k in ta_out_num])

# Calendar
ta_units = "days since 1960-01-01 00:00:00" 
ta_cal = "gregorian"

# Path, names and meta data for SNOWPACK
file_ind = ["ZGRN_V5_all_ind.txt"]
var_names_in = ["sublim", "evapor", "rainrate", "melt"] # [kg m-2], [kg m-2], [kg m-2 h-1], [kg m-2]
# Paths and file names
if ((platform.node() in platf_nam[0]) or (platform.node() in platf_nam[1])):
    path_in = [root_IMAU + "Data/Greenland/SP/ZGRN_V5_all/1D/"]
    path_out = root_IMAU + "Data/Greenland/SP/processed/"
else:
    path_in = ["/scratch/ms/nl/rucs/SP_Data/OUTPUT/ZGRN_V5_all/"]
    path_out = "/scratch/ms/nl/rucs/SP_Data/OUTPUT/processed/"
file_nam_in = ["ZGRN_V5_all"]
file_nam_out = "ZGRN_V5_SP"

# Missing value for NetCDF-files
miss_val = 9.96921e+36

########################################################################################################################
# Load general data
########################################################################################################################

# Indices with RACMO subdomain
ind_subdom = np.load(path_RACMO + "subdomain_indices.npz")
slice_rlat = slice(ind_subdom["ind_rlat"][0], ind_subdom["ind_rlat"][-1] + 1)
slice_rlon = slice(ind_subdom["ind_rlon"][0], ind_subdom["ind_rlon"][-1] + 1)

# Read text file with indices of grid boxes on RACMO2.3 domain
indices = []
for k in file_ind:
    ind_temp = np.loadtxt(path_ind + k, usecols = [0, 1], dtype = int) # (ind_rlat, ind_rlon)
    ind_temp[:,0] -= ind_subdom["ind_rlat"][0]
    ind_temp[:,1] -= ind_subdom["ind_rlon"][0]
    indices.append(ind_temp)

# Load geographical coordinates and ice mask
ncfile = Dataset(path_RACMO + "RACMO23_masks_ZGRN11.nc", "r")
ism = ncfile.variables["GrIS_caps_mask"][slice_rlat,slice_rlon].astype(bool)
lon = ncfile.variables["lon"][slice_rlat,slice_rlon]
lat = ncfile.variables["lat"][slice_rlat,slice_rlon]
ncfile.close()

########################################################################################################################
# Process model data (SNOWPACK; 1960 - 2014)
########################################################################################################################

# Get time axis
ncfile = Dataset(path_in[0] + file_nam_in[0] + "_103_met.nc", "r")               
ta_in = num2date(ncfile.variables["time"][:], units = ncfile.variables["time"].units, 
                 calendar = ncfile.variables["time"].calendar)
ncfile.close()
ta_in_num = np.asarray([(k - dt.datetime(1960, 1, 1)).total_seconds() for k in ta_in], dtype = np.float32)
tb_in_num = np.append(0., ta_in_num) # start at (1960, 1, 1, 0, 0)

# Compute seconds per time step
dt = np.array([k.total_seconds() for k in np.diff(ta_in)], dtype = np.float32)
dt = np.append(dt[0], dt) # [seconds]

# Check available data
data_avail = np.zeros_like(ism)
for k in range(len(path_in)):
    file_list = glob.glob(path_in[k] + "*met.nc")
    point_numbs = np.asarray([int(m.split("/")[-1].split("_")[-2]) for m in file_list], dtype = int)
    data_avail[indices[k][(point_numbs - 1),0],indices[k][(point_numbs - 1),1]] = True
    if (len(file_list) != len(indices[k])):
        print "Warning: not all files available from " + file_ind[k] + " (" + str(len(file_list)) + \
              " / " + str(len(indices[k])) + ")"

# Indicies for interpolation
data_interp_try = np.logical_and(ism, np.invert(data_avail)) # locations where data should be interpolated
ind_interp_0_try, ind_interp_1_try = np.where(data_interp_try)
data_interp = np.zeros_like(ism)

# Create array for data
data = np.empty((len(ta_in), len(ind_subdom["ind_rlat"]), len(ind_subdom["ind_rlon"])), dtype = np.float32)

# Loop over variables
interp_frame_mask = np.zeros_like(ism)
for k in range(len(var_names_out)):
    data.fill(np.nan)
    
    print "#################### Variable " + var_names_out[k] + " ####################"

    print "Importing files:"
    
    numb_files_imp = 0
    for m in range(len(file_ind)):

        file_list = glob.glob(path_in[m] + "*met.nc")
        for n in file_list:
            
            # Load model data
            ncfile = Dataset(n, "r")
            data_loc = ncfile.variables[var_names_in[k]][:len(ta_in)]
            ncfile.close()
            
            # Assign data to spatial array
            point_numb = int(n.split("/")[-1].split("_")[-2])
            data[:,indices[m][(point_numb - 1),0],indices[m][(point_numb - 1),1]] = data_loc         
            numb_files_imp += 1

            if (np.mod(numb_files_imp, 250) == 0):
                print "First " + str(numb_files_imp) + " files imported"
        
    # Convert data (if necessary)
    if (var_names_in[k] == "rainrate"):
        data *= dt[:,np.newaxis,np.newaxis] / 3600. # [kg m-2 h-1] to [kg m-2]
    
    # Spatial interpolation of data
    print "Data interpolation (spatial):"
    
    for m in range(len(ind_interp_0_try)):
          
        interp_frame_mask.fill(False)
        box_add = interp_frame
        while(True):

            slice_0 = slice(np.maximum((ind_interp_0_try[m] - box_add), 0), 
                            np.minimum((ind_interp_0_try[m] + box_add + 1), interp_frame_mask.shape[0]))
            slice_1 = slice(np.maximum((ind_interp_1_try[m] - box_add), 0), 
                            np.minimum((ind_interp_1_try[m] + box_add + 1), interp_frame_mask.shape[1]))
            interp_frame_mask[slice_0,slice_1] = True 
            ind_sel_0, ind_sel_1 = np.where(np.logical_and(interp_frame_mask, data_avail))
        
            if (len(ind_sel_0) > 0):   
                data[:,ind_interp_0_try[m],ind_interp_1_try[m]] = np.mean(data[:,ind_sel_0,ind_sel_1], axis = 1)
                data_interp[ind_interp_0_try[m],ind_interp_1_try[m]] = True
                break
            else:
                box_add += 1
                print "Warning (" + str(ind_interp_0_try[m]) + ", " + str(ind_interp_1_try[m]) + "): " + \
                      "no neighbouring grid boxes, increase frame size to " + str((2 * box_add + 1) ** 2)
                   
        if (np.mod(m, 1000) == 0):
            print "First " + str(m) + " of " + str(len(ind_interp_0_try)) + " grid boxes spatially interpolated"

    # Save data to NetCDF
    data[np.isnan(data)] = miss_val
    ncfile = Dataset(path_out + file_nam_out + "_" + var_names_out[k] + ".nc", "w")
    ncfile.createDimension("time", len(ta_in))
    ncfile.createDimension("rlat", len(ind_subdom["ind_rlat"]))
    ncfile.createDimension("rlon", len(ind_subdom["ind_rlon"]))    
    nc_lon = ncfile.createVariable("lon","f", ("rlat", "rlon"))
    nc_lon.units = "degrees_east"
    nc_lon[:] = lon
    nc_lat = ncfile.createVariable("lat","f", ("rlat", "rlon"))
    nc_lat.units = "degrees_north"
    nc_lat[:] = lat
    nc_time = ncfile.createVariable("time","f", ("time"))
    nc_time.units = ta_units
    nc_time.calendar = ta_cal
    nc_time[:] = date2num(ta_in, units = ta_units, calendar = ta_cal)
    nc_data = ncfile.createVariable(var_names_out[k].lower(),"f", ("time", "rlat", "rlon"), fill_value = miss_val)
    nc_data.unit = units_out[k]
    nc_data[:] = data
    nc_mask = ncfile.createVariable("data_avail","i", ("rlat", "rlon"))
    nc_mask[:] = data_avail.astype(int)
    nc_mask = ncfile.createVariable("data_interp","i", ("rlat", "rlon"))
    nc_mask[:] = data_interp.astype(int)
    ncfile.close()    
    
    print "File " + file_nam_out + "_" + var_names_out[k] + ".nc created"

########################################################################################################################
# Optional temporal interpolaton (quarter monthly values)
########################################################################################################################

    if (temp_interp): 

        # Temporal interpolation
        data[data == miss_val] = np.nan
        print "Data interpolation (temporal)"
        tb_out_num[-1] = tb_in_num[-1] # few days in last time step are missing
        data_temp_interp = consresamp(data, tb_in_num, tb_out_num, cons_quant = "sum", axis = 0)
        
        # Save data to NetCDF
        data[np.isnan(data)] = miss_val
        ncfile = Dataset(path_out + file_nam_out + "_" + var_names_out[k] + "_qm.nc", "w")
        ncfile.createDimension("time", len(ta_out))
        ncfile.createDimension("rlat", len(ind_subdom["ind_rlat"]))
        ncfile.createDimension("rlon", len(ind_subdom["ind_rlon"]))    
        nc_lon = ncfile.createVariable("lon","f", ("rlat", "rlon"))
        nc_lon.units = "degrees_east"
        nc_lon[:] = lon
        nc_lat = ncfile.createVariable("lat","f", ("rlat", "rlon"))
        nc_lat.units = "degrees_north"
        nc_lat[:] = lat
        nc_time = ncfile.createVariable("time","f", ("time"))
        nc_time.units = ta_units
        nc_time.calendar = ta_cal
        nc_time[:] = date2num(ta_out, units = ta_units, calendar = ta_cal)
        nc_data = ncfile.createVariable(var_names_out[k].lower(),"f", ("time", "rlat", "rlon"), fill_value = miss_val)
        nc_data.unit = units_out[k]
        nc_data[:] = data_temp_interp
        nc_mask = ncfile.createVariable("data_avail","i", ("rlat", "rlon"))
        nc_mask[:] = data_avail.astype(int)
        nc_mask = ncfile.createVariable("data_interp","i", ("rlat", "rlon"))
        nc_mask[:] = data_interp.astype(int)
        ncfile.close()    
 
        print "File " + file_nam_out + "_" + var_names_out[k] + "_qm.nc created"