# Description: Removes the second element of a CPU_work file and reduces the first number by one
#
# Author: Christian Steger, July 2016

# Load modules
import os
import sys

########################################################################################################################
# Settings
########################################################################################################################

# Root path to CPU-files
path_in = os.getenv("HOME") + "/Dropbox/IMAU/FDM_SP_runs/ECMWF_DATA/"
#path_in = "/home/ms/nl/rucs/Snowpack/DATA/"

########################################################################################################################
# Remove point and decrease number by one
########################################################################################################################

# Read command line arguments
CPU_numb = sys.argv[1] # CPU number
host_name = sys.argv[2] # either cca or ccb

# Modify point-CPU-file
file = open(path_in + host_name + "/points_CPU_" + CPU_numb + ".txt", "r")
content = file.read().split(" ")
file.close()
content = [str(int(content[0]) - 1)] + content[2:]
file = open(path_in + host_name + "/points_CPU_" + CPU_numb + ".txt", "w")
file.write(" ".join(content))
file.close()