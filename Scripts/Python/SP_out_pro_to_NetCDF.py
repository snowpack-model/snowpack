# Description: Converts SNOWPACK output format .pro to NetCDF
#
# Input arguments: file(s)
#
# Comments: Viscosity-values are very large for ice (1.000e+90) and are Inf in an array of type np.float32
#
# Author: Christian Steger, February 2016

# Load modules
import sys
import numpy as np
import datetime as dt
from netCDF4 import Dataset, date2num
import time

# Settings
print_sum = True # print summary information

# Missing value for NetCDF-files
miss_val = 9.96921e+36

########################################################################################################################
# Table with variables and properties that should be saved in the NetCDF-file
########################################################################################################################

variables = np.chararray((19, 4), itemsize = 30)
variables[:][:] = [["0501",   "height",         "Height",                          "cm"],
                   ["0502",   "densit",         "Density",                         "kg m-3"],
                   ["0503",   "temper",         "Temperature",                     "degC"],
                   ["0506",   "water",          "Water volume fraction",           "%"],
                   ["0508",   "dendric",        "Dendricity",                      "1"],
                   ["0509",   "spheric",        "Sphericity",                      "1"],
                   ["0510",   "coordnumb",      "Coordination number",             "1"],
                   ["0511",   "bondsize",       "Bond size",                       "mm"],
                   ["0512",   "grainsize",      "Grain size",                      "mm"],
                   ["0515",   "ice",            "Ice volume fraction",             "%"],
                   ["0516",   "air",            "Air volume fraction",             "%"],
                   ["0517",   "stress",         "Stress",                          "kPa"],
                   ["0518",   "viscos",         "Viscosity",                       "GPa s"],
                   ["0519",   "soil",           "Soil volume fraction",            "%"],
                   ["0520",   "tempgrad",       "Temperature gradient",            "K m-1"],
                   ["0521",   "conduct",        "Thermal conductivity",            "W K-1 m-1"],
                   ["0522",   "absorbrad",      "Absorbed shortwave radiation",    "W m-2"],
                   ["0523",   "deformrate",     "Viscous deformation rate",        "1.e-6 s-1"],
                   ["0535",   "optgrainsize",   "Optical equivalent grain size",   "mm"]]

# Warning: variables 0513 (grain type) and 0530 (position (cm) and minimum stability indices) may have different
#          numbers of elements per timestamp than other variables!

########################################################################################################################
# Loop over file(s) provided as input argument(s)
########################################################################################################################

numb_arg = len(sys.argv) - 1 # number of files
for k in range(0, numb_arg):    
    path_and_name = sys.argv[k + 1]
    
    # Check file ending
    if (path_and_name.split(".")[-1] != "pro"):
        print "Error: File ending is not .pro"
        continue

    # Import data from text file
    file = open(path_and_name, "r")
    content = file.read()
    file.close()
    content = content.split("\n")
    
    # Get model information, delete part with [STATION_PARAMETERS] and [HEADER]
    model_info = content[9][1:]
    del content[:44] # line with start of data
        
    # Get maximal number of elements and number of time steps
    time_steps = np.sum(np.asarray([(content[n].split(",")[0] == "0500") for n in range(len(content))]))
    ind_vv = np.asarray([(content[n].split(",")[0] in variables[:,0]) for n in range(len(content))])
    max_num_elem = np.max((np.asarray([content[n].split(",")[1] for n in range(len(content))])[ind_vv]).astype(int))

    # Assign selected data to array
    data = np.empty((max_num_elem, time_steps, variables.shape[0]), dtype = np.float32)
    data.fill(np.nan)
    time_axis = []
    current_time = -1
    for n in range(len(content)):
        if (content[n].split(",")[0] == "0500"):
            time_axis.append(content[n].split(",")[1])
            current_time += 1
        elif (content[n].split(",")[0] in variables[:,0]):
            length = int(content[n].split(",")[1])
            index = list(variables[:,0]).index(content[n].split(",")[0])
            data[:length,current_time,index] = np.asarray(content[n].split(",")[2:]).astype(float)

    # Convert time axis to list of datetime-objects
    time_axis = [dt.datetime.strptime(time_axis[n], "%d.%m.%Y %H:%M:%S") for n in range(time_steps)]
    
    # Convert variables supplied in [cm] to [m]
    data[:,:,0] /= 100
    variables[0,3] = "m"

    # Conversion of height/depth to element thickness
    ne_soil = (data[:,0,0] < 0.).sum()
    if (not np.all((data[:,:,0] < 0.).sum(axis = 0) == ne_soil)):
        sys.exit("Error: Number of soil elements varies with time")
    if (ne_soil == 0):
        data[:,:,0] = np.diff(np.concatenate((np.zeros((1, time_steps), dtype = np.float32), data[:,:,0])), axis = 0)
    else: # (height [> 0: top, < 0: bottom of elem.]) and line with 0. in height/depth array
        data[:(max_num_elem - 1),:,0] = np.diff(data[:,:,0], axis = 0)
        data = data[:(max_num_elem - 1),:,:] # delete upmost line (only contains nan)
        max_num_elem -= 1
    variables[0,1] = "layerthickn"
    variables[0,2] = "Layer thickness"

    # Print summary of last time step (only consider snow elements)
    if (print_sum):
        snow_last = data[ne_soil:,-1,:]
        if ((len(snow_last) > 0) & (np.invert(np.isnan(snow_last[:,0])).sum() > 0)):
            ind_nn = np.invert(np.isnan(snow_last[:,0]))
            layer_thickn = snow_last[ind_nn,0] # [m]
            mass = np.nansum(layer_thickn * snow_last[ind_nn,1]) # vertically summed total mass [kg m-2]
            liq_mass = np.nansum(layer_thickn * (snow_last[ind_nn,3] / 100.)) * 1000. # vertic. summed water [kg m-2]
            ind = [1, 2, 3, 9, 10, 11, 12, 14, 15] # indices for mean
            data_mean = np.average(snow_last[ind_nn,:][:,ind], weights = layer_thickn, axis = 0) # weighted average
            print "###############################################"
            print path_and_name.split("/")[-1]
            print "---- Summary of last time step (only snow) ----"
            print "--------------------- Sum ---------------------"
            print "thick: " + "%.5f" % layer_thickn.sum() + " [m]"
            print "totalmass: " + "%.5f" % mass + " [kg m-2]" 
            print "watermass: " + "%.5f" % liq_mass + " [kg m-2]"
            print "--------------------- Mean --------------------"
            for m in range(len(ind)):
                print variables[ind[m],1] + ": " + "%.5f" % data_mean[m] + " [" + variables[ind[m],3] + "]"        
                        
    # Save data to NetCDF-file
    data[np.isnan(data)] = miss_val
    ncfile = Dataset(path_and_name[:-4] + "_pro.nc", "w")
    ncfile.history = "Created " + time.ctime(time.time())
    ncfile.source = model_info.replace('"', "") # remove quotation marks from string
    ncfile.createDimension("layer", max_num_elem)
    ncfile.createDimension("time", time_steps)
    nc_time = ncfile.createVariable("time", "f", ("time"))
    nc_time.units = "days since 1950-01-01 00:00:00"
    nc_time.calendar = "gregorian"
    nc_time[:] = date2num(time_axis, units = nc_time.units, calendar = nc_time.calendar)
    for k in range(variables.shape[0]):
        nc_data = ncfile.createVariable(variables[k,1], "f", ("layer", "time"), fill_value = miss_val)
        nc_data.long_name = variables[k,2]
        nc_data.unit = variables[k,3]
        nc_data[:] = data[:,:,k]
    ncfile.close()

if (print_sum):    
    print "###############################################"