# Description: Change the dates of a .sno-file
#
# Details: Datetime objects in python cannot be BC. Therefore, as soon as the date should be shifted below a certain
# threshold (default: year 100), the date is simply set to the previous date minus one hour
#
# Author: Christian Steger, July 2016

# Load modules
import sys
import datetime as dt

# Settings
profile_date_new = dt.datetime(1960, 01, 01, 01, 00) # new date of profile (YYYY, MM, DD, HH, MM)
date_thresh = dt.datetime(100, 1, 1) # threshold date (python datetime objects cannot be BC)

########################################################################################################################
# Modify date in .sno-file
########################################################################################################################

# Reads file name as input argument
path_and_name = sys.argv[1] # path and file name of snow file (.sno)

# Import text file
file = open(path_and_name,"r")
content = file.read()
file.close()
content = content.split("\n")

# Find indices of relevant rows
ind_pd = [(k[:11] == "ProfileDate") for k in content].index(True)
ind_DATA = [(k[:6] == "[DATA]") for k in content].index(True)
ind_numb_elem_snow = [(k[:14] == "nSnowLayerData") for k in content].index(True)
ind_numb_elem_soil = [(k[:14] == "nSoilLayerData") for k in content].index(True)

# Retrieve number of elements
numb_elem = int(content[ind_numb_elem_snow].split("=")[1]) + \
            int(content[ind_numb_elem_soil].split("=")[1])

# Calculate time difference between old and new date
line_pd = content[ind_pd]
profile_date_old = dt.datetime.strptime(line_pd.split()[-1], "%Y-%m-%dT%H:%M:%S")
delta = profile_date_new - profile_date_old

# Change dates in .sno-file according to time difference
content[ind_pd] = " ".join(line_pd.split()[:2]) + " " + profile_date_new.isoformat() # strftime only works year > 1900
for k in range((ind_DATA + numb_elem), ind_DATA, -1): # loop from bottom to top (from younger to older elements)
    date_old = dt.datetime.strptime(content[k].split()[0], "%Y-%m-%dT%H:%M:%S")
    try:
        date_new = date_old + delta # may cause OverflowError
        if (date_new > date_thresh):
            content[k] = date_new.isoformat() + content[k][19:]
        else:
            last_date = dt.datetime.strptime(content[(k + 1)].split()[0], "%Y-%m-%dT%H:%M:%S")
            content[k] = (last_date - dt.timedelta(hours = 1)).isoformat() + content[k][19:]
    except OverflowError:
        last_date = dt.datetime.strptime(content[(k + 1)].split()[0], "%Y-%m-%dT%H:%M:%S")
        content[k] = (last_date - dt.timedelta(hours = 1)).isoformat() + content[k][19:]
      
# Write modified text file to SNOWPACK-output-folder
file = open(path_and_name,"w")
for k in range(0, len(content) - 1):
    file.write(content[k] + "\n")  
file.close()