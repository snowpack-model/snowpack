# Description: Computation of firn air content and refreezing capacity from 2D SNOWPACK output
#
# Author: Christian Steger, March 2017

# Load modules
import os
import sys
import platform
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import datetime as dt
from dateutil.relativedelta import relativedelta
import glob

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
import firn_utilities_cy

########################################################################################################################
# Settings
########################################################################################################################

# Depth threshold (lower limit for firn layer consideration) -> lower part may not be refreshed at all locations...
depth_tresh = -40. # [m]

# General settings
interp_frame = 1 # frame for spatial interpolation (frame length: 2 * interp_frame + 1)
file_ind = "ZGRN_V5_all_ind.txt"

# Output variables and units
vn_out = ["firn_air", "refrez_cap"] # refreezing capacity not really useful quantity... (not used)
units_out = ["m", "kg m-2"]
                                                                     
# Output time axis
ta_out = np.asarray([(dt.datetime(1960, 1, 16, 0, 0) + relativedelta(months = k)) for k in range(55 * 12)])
ta_units = "days since 1960-01-01 00:00:00" 
ta_cal = "gregorian"

# Path, names and meta data for SNOWPACK
vn_in = ["temper", "ice", "layerthickn"] # [deg C], [%], [%], [m]
if ((platform.node() in platf_nam[0]) or (platform.node() in platf_nam[1])):
    path_in = root_IMAU + "Data/Greenland/SP/ZGRN_V5_all/2D/"
    path_out = root_IMAU + "Data/Greenland/SP/processed/"
else:
    path_in = "/scratch/ms/nl/rucs/SP_Data/OUTPUT/ZGRN_V5_all/"
    path_out = "/scratch/ms/nl/rucs/SP_Data/OUTPUT/processed/"
file_nam_in = "ZGRN_V5_all"
file_nam_out = "ZGRN_V5_SP"

# Missing value for NetCDF-files
miss_val = 9.96921e+36

# Constants
dens_ice = 917. # density of ice [kg m-3]
dens_water = 1000. # density of water [kg m-3]
spec_heat_ice = 2100.0 # specific heat of ice at 273.15 K (equal to SNOWPACK) [J kg-1 K-1]
lat_heat_fus = 3.34e5 # latent heat of fusion (equal to SNOWPACK) [J kg-1]

########################################################################################################################
# Load general data
########################################################################################################################

# Indices with RACMO subdomain
ind_subdom = np.load(path_RACMO + "subdomain_indices.npz")
slice_rlat = slice(ind_subdom["ind_rlat"][0], ind_subdom["ind_rlat"][-1] + 1)
slice_rlon = slice(ind_subdom["ind_rlon"][0], ind_subdom["ind_rlon"][-1] + 1)

# Read text file with indices of grid boxes on RACMO2.3 domain
indices = np.loadtxt(path_ind + file_ind, usecols = [0, 1], dtype = int) # (ind_rlat, ind_rlon)
indices[:,0] -= ind_subdom["ind_rlat"][0]
indices[:,1] -= ind_subdom["ind_rlon"][0]

# Load geographical coordinates and ice mask
ncfile = Dataset(path_RACMO + "RACMO23_masks_ZGRN11.nc", "r")
ism = ncfile.variables["GrIS_caps_mask"][slice_rlat,slice_rlon].astype(bool)
ncfile.close()

########################################################################################################################
# Process model data (firn model or SNOWPACK; 1960 - 2014)
########################################################################################################################

# Get time axis
ncfile = Dataset(path_in + file_nam_in + "_103_pro.nc", "r")               
ta_in = num2date(ncfile.variables["time"][:], units = ncfile.variables["time"].units, 
                        calendar = ncfile.variables["time"].calendar)
ncfile.close()

# Check available data
data_avail = np.zeros_like(ism)
file_list = glob.glob(path_in + "*pro.nc")
point_numbs = np.asarray([int(m.split("/")[-1].split("_")[-2]) for m in file_list], dtype = int)
data_avail[indices[(point_numbs - 1),0],indices[(point_numbs - 1),1]] = True
if (len(file_list) != len(indices)):
    print "Warning: not all files available from " + file_ind + " (" + str(len(file_list)) + \
          " / " + str(len(indices)) + ")"

# Indicies for interpolation
data_interp_try = np.logical_and(ism, np.invert(data_avail)) # locations where data should be interpolated
ind_interp_0_try, ind_interp_1_try = np.where(data_interp_try)
data_interp = np.zeros_like(ism)

# Create array for data
data = np.empty((len(vn_out), len(ta_out), len(ind_subdom["ind_rlat"]), len(ind_subdom["ind_rlon"])), 
                dtype = np.float32)
data.fill(np.nan)

# Numerical time axis
rel_date = dt.datetime(1960, 1, 1)
ta_in_num = np.asarray([(k - rel_date).total_seconds() for k in ta_in], dtype = np.float32)
ta_out_num = np.asarray([(k - rel_date).total_seconds() for k in ta_out], dtype = np.float32)

# Loop over files
interp_frame_mask = np.zeros_like(ism)    
numb_files_imp = 0
file_list = glob.glob(path_in + "*pro.nc")
bound_out = np.array([depth_tresh, 0.], dtype = np.float32)
for k in file_list:
    
    # Load model data (set missing values to nan when array is masked)
    ncfile = Dataset(k, "r")
    data_loc = np.empty((len(vn_in),) + ncfile.variables[vn_in[0]].shape, dtype = np.float32) # (variables, layer, time)
    for m in range(len(vn_in)):
        temp = ncfile.variables[vn_in[m]][:]
        if (type(temp).__name__ == "MaskedArray"):
            temp = temp.filled(np.nan)
        data_loc[m,:,:] = temp
    ncfile.close()
               
    # Check SNOWPACK data for elements with thickness 0.
    if (np.any(data_loc[vn_in.index("layerthickn"),:,:] == 0.)):
        sum_nn = np.sum(np.invert(np.isnan(data_loc[vn_in.index("layerthickn"),:,:])), axis = 0)
        ind_zero_0, ind_zero_1 = np.where(data_loc[vn_in.index("layerthickn"),:,:] == 0.)
        for o in range(len(ind_zero_0)):
            if ((sum_nn[ind_zero_1[o]] - 1) == ind_zero_0[o]): # check if element is at the top
                data_loc[:,ind_zero_0,ind_zero_1] = np.nan
            else:
                sys.exit("Error: Found element with thickness 0. in SNOWPACK data is not at the top")

    # Calculate layer boundaries
    bound_in = np.vstack((np.repeat(0., len(ta_in)).astype(np.float32), 
                          np.cumsum(data_loc[vn_in.index("layerthickn"),:,:], axis = 0)))
    bound_in = bound_in - np.nanmax(bound_in, axis = 0)

    # Compute firn air content (liquid water content not considered)
    vol_ice_cont = np.asarray(firn_utilities_cy.firnproresamp(data_loc[vn_in.index("ice"),:,:], bound_in, bound_out, 
                                                              cons_quant = "weighted_mean"))[0,:]  
    firn_air = (1. - (vol_ice_cont / 100.)) * np.abs(depth_tresh) # [m]
    firn_air = firn_air.clip(min = 0., max = 40.)
    
    # Compute refreezing capacity (limited either by pore space or cold content) [kg m-2]
    # -> neglect heat capacity of air in pore space
    # -> liquid water not considered (refreezing capacity 0. anyway because no cold content available)
    vol_frac_av = (1.0 - (data_loc[vn_in.index("ice"),:,:] / 100.))
    space_lim = vol_frac_av * data_loc[vn_in.index("layerthickn"),:,:] * dens_ice # [kg m-2] 
    # -> consider expansion of water when refreezing
    delta_T = (0. - data_loc[vn_in.index("temper"),:,:]) # [K]
    delta_T = delta_T.clip(min = 0.)
    energ_abs = ((data_loc[vn_in.index("ice"),:,:] / 100.) * data_loc[vn_in.index("layerthickn"),:,:] * dens_ice) * \
                spec_heat_ice * delta_T # [J m-2]      
    energ_lim = energ_abs / lat_heat_fus # [kg m-2]
    lim = np.minimum(space_lim, energ_lim)
    refrez_cap = np.asarray(firn_utilities_cy.firnproresamp(lim, bound_in, bound_out, cons_quant = "sum"))[0,:]

    # Temporal interpolation (monthly values)
    point_numb = int(k.split("/")[-1].split("_")[-2])
    data[0,:,indices[(point_numb - 1),0],indices[(point_numb - 1),1]] = np.interp(ta_out_num, ta_in_num, firn_air)
    data[1,:,indices[(point_numb - 1),0],indices[(point_numb - 1),1]] = np.interp(ta_out_num, ta_in_num, refrez_cap)
      
    numb_files_imp += 1
    if (np.mod(numb_files_imp, 250) == 0):
        print "First " + str(numb_files_imp) + " files imported and temporally interpolated"

# Spatial interpolation of data
print "Data interpolation:"
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
            data[:,:,ind_interp_0_try[m],ind_interp_1_try[m]] = np.mean(data[:,:,ind_sel_0,ind_sel_1], axis = 2)
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
for k in range(len(vn_out)):
    file_nam = file_nam_out + "_" + vn_out[k] + "_top" + ("%.0f" % np.abs(depth_tresh)) + "m.nc"
    ncfile = Dataset(path_out + file_nam, "w")
    ncfile.createDimension("time", len(ta_out))
    nc_time = ncfile.createVariable("time","f", ("time"))
    nc_time.units = ta_units
    nc_time.calendar = ta_cal
    nc_time[:] = date2num(ta_out, units = ta_units, calendar = ta_cal)
    ncfile.createDimension("rlat", len(ind_subdom["ind_rlat"]))
    ncfile.createDimension("rlon", len(ind_subdom["ind_rlon"]))
    nc_data = ncfile.createVariable(vn_out[k].lower(),"f", ("time", "rlat", "rlon"), fill_value = miss_val)
    nc_data.unit = units_out[k]
    nc_data[:] = data[k,:,:,:]   
    nc_mask = ncfile.createVariable("data_avail","i", ("rlat", "rlon"))
    nc_mask[:] = data_avail.astype(int)
    nc_mask = ncfile.createVariable("data_interp","i", ("rlat", "rlon"))
    nc_mask[:] = data_interp.astype(int)
    ncfile.close()
    
    print "File " + file_nam_out + "_" + vn_out[k] + ".nc created"