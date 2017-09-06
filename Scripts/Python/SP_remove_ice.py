# Description: Try to remove ice at the bottom if firn thickness gets to large
#
# Author: Christian Steger, July 2016

# Load modules
import sys

########################################################################################################################
# Settings
########################################################################################################################

firn_thick_max = 120. # thickness to check for excess ice [m]
theta_ice_thresh = 0.99 # threshold for volumetric ice content (907.83 kg m-3) [-]

########################################################################################################################
# Load .sno-file, remove ice at bottom if necessary and possible
########################################################################################################################

# Reads file name as input argument
path_and_name = sys.argv[1] # path and file name of snow file (.sno)

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

# Retrieve firn thickness and number of elements
firn_thick = float(content[ind_HS].split("=")[1]) # [m]
numb_elem = int(content[ind_numb_elem].split("=")[1])

# Try to remove ice at bottom
numb_elem_rem = 0
while (firn_thick > firn_thick_max):
    if (float(content[ind_DATA + 1].split()[3]) > theta_ice_thresh): # volumetric ice content [-]
        firn_thick -= float(content[ind_DATA + 1].split()[1]) # layer thickness [-]
        del content[ind_DATA + 1]
        numb_elem_rem += 1
    else:
        break # non-ice element

# Update firn thickness, number of layers and erosion level
content[ind_HS] = "HS_Last = " + ("%.6f" % firn_thick)
content[ind_numb_elem] = "nSnowLayerData = " + str(numb_elem - numb_elem_rem)           
content[ind_e_lev] = "ErosionLevel = " + str(numb_elem - numb_elem_rem - 1)                              

# Output modified text file
file = open(path_and_name,"w")
for k in range(0, len(content) - 1):
    file.write(content[k] + "\n")  
file.close()