# Description: Add ice at bottom
#
# Author: Christian Steger, July 2016

# Load modules
import os
import sys
import numpy as np

########################################################################################################################
# Settings
########################################################################################################################

# Path and name of file with mass changes
file_in = os.getenv("HOME") + "/Dropbox/IMAU/FDM_SP_runs/ECMWF_DATA/Charalampidis_NoR_and_mc.txt"
#file_in = "/home/ms/nl/rucs/Snowpack/DATA/ZGRN_V5_all_NoR_and_mc.txt"

# Ice block settings
firn_thick_min = 60. # minimal firn column thickness [m]
elem_thickn = 1.95 # element thickness for ice added at the bottom [m]

# Constants
dens_water = 1000. # density of water [kg m-3]

########################################################################################################################
# Compute mass loss during entire period from RACMO2.3 data
########################################################################################################################

# Reads file name and point number as input argument
path_and_name = sys.argv[1] # path and file name of snow file (.sno)
p_ind = int(sys.argv[2]) - 1 # indices point = (point number - 1)
period = sys.argv[3] # period over which SNOWPACK is run next (spinup or entire)

# Define necessary column
if (period == "reference"):
    usecol = 1 
elif (period == "entire"):
    usecol = 2
else:
    sys.exit("Error: Unknown input for period")

# Read text file with mass changes for reference and entire period
mass_change = np.loadtxt(file_in, usecols = (usecol,), dtype = np.float32)[p_ind]
# mass changes for reference [m weq (20 a-1)] and entire [m weq (55 a-1)] period

if (mass_change < 0.):
    numb_elem_add_mc = np.ceil(np.abs(mass_change) / elem_thickn).astype(int)
else:
    numb_elem_add_mc = 0

########################################################################################################################
# Load .sno-file, add ice at bottom
########################################################################################################################

# Import text file
file = open(path_and_name,"r")
content = file.read()
file.close()
content = content.split("\n")

# Find indices of relevant rows
ind_HS = [(k[:7] == "HS_Last") for k in content].index(True)
ind_DATA = [(k[:6] == "[DATA]") for k in content].index(True)
ind_numb_elem = [(k[:14] == "nSnowLayerData") for k in content].index(True)
ind_e_lev = [(k[:12] == "ErosionLevel") for k in content].index(True)

# Retrieve firn column thickness, compute numb_elem_add to bring the firn colum to the minimal thickness
firn_thick = float(content[ind_HS].split("=")[1])
if (firn_thick < firn_thick_min):
    numb_elem_add_min = np.ceil(np.abs(firn_thick_min - firn_thick) / elem_thickn).astype(int)
else:
    numb_elem_add_min = 0
numb_elem_add = numb_elem_add_mc + numb_elem_add_min # ensure that ablation always start at least from "base thickness"

# Add new elements
if (numb_elem_add > 0):

    # Retrieve number of elements
    numb_elem = int(content[ind_numb_elem].split("=")[1])

    # Add necessary number of elements
    lowest_elem = content[ind_DATA + 1].split()
    lowest_elem[1] = "%.3f" % elem_thickn
    lowest_elem[3] = "1.000" # volumetric ice content
    lowest_elem[4] = "0.000" # volumetric water content
    lowest_elem[5] = "0.000" # volumetric air content
    add_elem = "     ".join(lowest_elem)
    for k in range(numb_elem_add):
        content.insert((ind_DATA + 1), add_elem)

    # Update firn thickness, number of layers and erosion level
    content[ind_HS] = "HS_Last = " + ("%.6f" % (float(content[ind_HS].split("=")[1]) + numb_elem_add * elem_thickn))
    content[ind_numb_elem] = "nSnowLayerData = " + str(numb_elem + numb_elem_add)            
    content[ind_e_lev] = "ErosionLevel = " + str(numb_elem + numb_elem_add - 1)                              

# Output modified text file
file = open(path_and_name,"w")
for k in range(0, len(content) - 1):
    file.write(content[k] + "\n")  
file.close()