# Description: Converts SNOWPACK output format .met to NetCDF
#
# Input arguments: file(s)
#
# Author: Christian Steger, October 2015

# Load modules
import sys
import numpy as np
import datetime as dt
from netCDF4 import Dataset, date2num
import time

# Missing value for NetCDF-files
miss_val = 9.96921e+36

# Settings
print_sum = True # print summary information

########################################################################################################################
# Table with variables and properties that should be saved in the NetCDF-file
########################################################################################################################

# Table with variables that should be saved in NetCDF (short name, long name)
variables = np.chararray((37, 2), itemsize = 65)
variables[:][:] = [["sensheat",         "Sensible heat"],
                   ["latheat",          "Latent heat"],
                   ["longwaveout",      "Outgoing longwave radiation"],
                   ["longwavein",       "Incoming longwave radiation"],
                   ["shortwaveout",     "Reflected shortwave radiation"],
                   ["shortwavein",      "Incoming shortwave radiation"],
                   ["modsurfalbedo",    "Modelled surface albedo"],
                   ["airtemp",          "Air temperature"],
                   ["modsurftemp",      "Modeled surface temperature"],
                   ["meassurftemp",     "Measured surface temperature"],
                   ["bottemp",          "Temperature at bottom of snow or soil pack"],
                   ["botheatflux",      "Heat flux at bottom of snow or soil pack"],
                   ["groundtemp",       "Ground surface temperature"],
                   ["groundheatflux",   "Heat flux at ground surface"],
                   ["heatadvec",        "Heat advected to the surface by liquid precipitation"],
                   ["meassurfalbedo",   "Measured surface albedo"],
                   ["relhumid",         "Relative humidity"],
                   ["windspeed",        "Wind speed"],
                   ["precradsolid",     "Precipitation rate at surface (solid only)"],
                   ["modsnowdepth",     "Modelled snow depth (vertical)"],
                   ["enfsnowdepth",     "Enforced snow depth (vertical)"],
                   ["surfhoarsize",     "Surface hoar size"],
                   ["SWE",              "SWE (of snowpack)"],
                   ["erodedmass",       "Eroded mass"],
                   ["rainrate    ",     "Rain rate"],                    
                   ["runoff",           "Snowpack runoff (virtual lysimeter)"],
                   ["sublim",           "Sublimation"],
                   ["evapor",           "Evaporation"],
                   ["liquidwatcont",    "Liquid Water Content (of snowpack)"],
                   ["soilrunoff",       "Soil runoff"],
                   ["energchange",      "Internal energy change"],                    
                   ["surfin",           "Surface input (sum fluxes)"],
                   ["measnewsnowdens",  "Measured new snow density"],
                   ["modnewsnowdens",   "Modeled new snow density"],                  
                   ["meltfreez",        "Melt freeze part of internal energy change"],
                   ["melt",             "Melted mass"],
                   ["refreeze",         "Refrozen mass"]]

########################################################################################################################
# Loop over file(s) provided as input argument(s)
########################################################################################################################

numb_arg = len(sys.argv) - 1 # number of files
for k in range(0, numb_arg):
    path_and_name = sys.argv[k + 1]     
                                                                                                
    # Check file ending
    if (path_and_name.split(".")[-1] != "met"):
        print "Error: File ending is not .met"
        continue
                                                 
    # Import data from text file
    file = open(path_and_name, "r")
    content = file.read()
    file.close()
    content = content.split("\n")

    # Get indices of selected variables, units of variables and model information
    var_avail = content[12].split(",") # line with long names of variables
    var_ind = np.asarray([var_avail.index(variables[n,1]) for n in range(0, variables.shape[0])], dtype = int)
    units = np.asarray(content[13].split(","))[var_ind] # line with units of variables
    model_info = content[10][1:]

    # Delete part with [STATION_PARAMETERS] and [HEADER]
    del content[:16] # line with start of data

    # Assign selected data to array
    data = np.empty((len(content), variables.shape[0]), dtype = np.float32)
    time_axis = []
    for n in range(0, len(content)):
        time_axis.append(content[n].split(",")[1])
        data_raw = np.asarray(content[n].split(","))[var_ind] # may contain empty elements
        ind_empty = (data_raw == "")
        data_raw[ind_empty] = "0"
        data[n,:] = data_raw.astype(float)
        data[n,ind_empty] = np.nan

    # Convert time axis to list of datetime-objects
    time_axis = [dt.datetime.strptime(time_axis[n], "%d.%m.%Y %H:%M:%S") for n in range(0, len(content))]

    # Convert variables supplied in [cm] to [m]
    ind_cm = (units == "cm")
    data[:,ind_cm] /= 100
    units[ind_cm] = "m"

    # Print summary (snow depth and SWE -> last time step, other variables: temporal sum or mean)
    if (print_sum):
        modsnowdepth = data[-1,19] # [m]
        enfsnowdepth = data[-1,20] # [m]
        SWE = data[-1,22] # [m]
        liquidwatcont = data[-1,28] # [kg m-2]
        ind_sum = [25, 26, 27, 29, 35, 36]
        data_sum = np.nansum(data[:,ind_sum], axis = 0)
        ind_mean = range(19) + [21, 23, 24, 28, 30, 31, 32, 33, 34]
        data_mean = np.nanmean(data[:,ind_mean], axis = 0)
        print "###############################################"
        print path_and_name.split("/")[-1]
        print "--------------- Last time step ----------------"
        print "modsnowdepth: " + "%.3f" % modsnowdepth + " [m]"
        print "enfsnowdepth: " + "%.3f" % enfsnowdepth + " [m]"
        print "SWE: " + "%.3f" % SWE + " [kg m-2]"
        print "liquidwatcont: "  + "%.3f" % liquidwatcont + " [kg m-2]"
        print "---------------- Temporal sum -----------------"
        for m in range(len(ind_sum)):
            print variables[ind_sum[m],0] + ": " + "%.3f" % data_sum[m] + " [" + units[ind_sum[m]] + "]"
        print "---------------- Temporal mean ----------------"
        for m in range(len(ind_mean)):
            print variables[ind_mean[m],0] + ": " + "%.5f" % data_mean[m] + " [" + units[ind_mean[m]] + "]"

    # Save data to NetCDF-file
    data[np.isnan(data)] = miss_val
    ncfile = Dataset(path_and_name[:-4] + "_met.nc", "w")
    ncfile.history = "Created " + time.ctime(time.time())
    ncfile.source = model_info.replace('"', "") # remove quotation marks from string
    ncfile.createDimension("time", data.shape[0])
    nc_time = ncfile.createVariable("time", "f", ("time"))
    nc_time.units = "days since 1950-01-01 00:00:00"
    nc_time.calendar = "gregorian"
    nc_time[:] = date2num(time_axis, units = nc_time.units, calendar = nc_time.calendar)
    for n in range(0, variables.shape[0]):
        nc_data = ncfile.createVariable(variables[n,0], "f", ("time"), fill_value = miss_val)
        nc_data.long_name = variables[n,1]
        nc_data.unit = units[n]
        nc_data[:] = data[:,n]
    ncfile.close()

if (print_sum):
    print "###############################################"