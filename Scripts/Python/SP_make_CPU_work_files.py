# Description: Distribute workload approximately equal over available CPUs (for cca and ccb)
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
    path_in_out = root_IMAU + "FDM_SP_runs/ECMWF_DATA/"
elif (platform.node() in platf_nam[1]):  
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in_out = root_IMAU + "FDM_SP_runs/ECMWF_DATA/"
elif (platform.node() in platf_nam[2]):
    path_in_out = "/home/ms/nl/rucs/Snowpack/DATA/"
else:
    sys.exit("Unknown platform")

# Read command line arguments
noc = int(sys.argv[1]) # number of CPUs
host_name = sys.argv[2] # host name
point_file = sys.argv[3] # text file with selection of points to run (e.g ZGRN_V5_all_sel_points_cca.txt)
NoR_file = sys.argv[4] # text file with NoR (e.g. ZGRN_V5_all_NoR_and_mc.txt)

########################################################################################################################
# Generate files with points for different CPUs
########################################################################################################################

# Load points and computational weights
points = np.loadtxt(path_in_out + point_file, dtype = int)
comp_weight = np.loadtxt(path_in_out + NoR_file, dtype = int, usecols = (0,))[(points - 1)]
comp_weight += np.round(55. / 20.).astype(int)

# Distribute points approximately equal over available CPUs (ascending order -> start with fast points)
ind_sort = np.argsort(comp_weight)
CPUs = [[] for k in range(noc)]
work_load_per_CPU = np.empty(20, dtype = np.float32) # fraction [-]
for k in range(noc):
    CPUs[k] = points[ind_sort[slice(k, len(ind_sort), noc)]]
    work_load_per_CPU[k] = float(comp_weight[ind_sort[slice(k, len(ind_sort), noc)]].sum()) / float(comp_weight.sum())

# Workload per CPU
numb_points = 0
print "##### Workload per CPU [%] #####"
for k in range(noc):
    print "CPU " + str(k + 1) + ": " + "%.2f" % (work_load_per_CPU[k] * 100.)
    print "Number of points: " + str(len(CPUs[k]))
    numb_points += len(CPUs[k])
print "Total workload: " "%.2f" % (work_load_per_CPU.sum() * 100.)
print "Total number of points: " + str(numb_points)
print "################################"

# Create folder (if necessary)
if (not os.path.exists(path_in_out + host_name)):
    os.makedirs(path_in_out + host_name)

# Write to text-files
for k in range(len(CPUs)):
    file = open(path_in_out + host_name + "/points_CPU_" + str(k + 1) + ".txt", "w")
    file.write(str(len(CPUs[k])) + " " + " ".join(list(CPUs[k].astype(str))))
    file.close()

# Check if all points are distributed
points_check = np.empty(0, dtype = int)
for k in range(len(CPUs)):
    points_check = np.append(points_check, CPUs[k])
if ((len(points_check) != len(points)) or (not np.all(np.sort(points_check) == np.sort(points)))):
    sys.exit("Error: Distribution of points to CPUs erroneous")