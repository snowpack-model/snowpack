# Description: Generates textfile with artifical meteorological input data for SNOWPACK
#
# No seasonal variation in meteorological parameters (only diurnal variation if days resolved)
#
# Author: Christian Steger, October 2015

# Load modules
import sys
import platform
import datetime as dt
import calendar as cal
import numpy as np

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"], # Desktop-Mac
             ["Christians-MacBook-Pro-2.local"]]  

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]): 
    path_out = "/Users/steger/Desktop/SNOWPACK_run_desktop/input/"
elif (platform.node() in platf_nam[1]): 
    path_out = "/Users/christiansteger/Desktop/SNOWPACK_run_desktop/input/"
else:
    sys.exit("Unknown platform")

########################################################################################################################
# Settings
########################################################################################################################

# File name
fn_part_out = "test"

# General
rain_thresh = 273.15 # treshold for rainfall [K]
dt_in = 3 # temporal resolution of input data (integer) [h] (allowed steps: 1, 2, 3, 4, 6, 8, 12, 24)
seas_cycle = True # diurnal cycle for surface and air temperature and shortwave radiation

# Meteorological values (name, mean and seasonal amplitude)
val_mean_amp = [["TA",     268.15,   0.],  # air temperature [K]
                ["RH",     0.9],            # relative  humidity [-]
                ["TSG",    268.15],         # ground temperature [K]
                ["TSS",    268.15,   0.],  # surface temperature [K]
                ["VW",     0.0],            # wind speed [m s-1]
                ["OSWR",   400.],           # outgoing shortwave rad. [W m-2]
                ["ISWR",   410.],           # incoming shortwave rad. [W m-2]
                ["ILWR",   300.],           # incoming longwave radiation [W m-2]
                ["PINT",   1.0],            # precipitation [kg m-2 h-1]
                ["SMELT",  0.0],            # snow melt [kg m-2 h-1]
                ["SNOWD",  0.0],            # snow drift [kg m-2 h-1]
                ["SUBLI",  -2.0],            # sublimation [kg m-2 h-1]
                ["RHO_HN", 369.0]]          # new snow density [kg m-3]
                
########################################################################################################################
# Constant settings
########################################################################################################################

# Meteorological variables
var_names = ["TA", "RH", "TSG", "TSS", "VW", "OSWR", "ISWR", "ILWR", "PINT", "PSUM_PH", "SMELT", "SNOWD", "SUBLI",
             "RHO_HN"]

########################################################################################################################
# Generate artifical meteorological data and save to SNOWPACK textfile
########################################################################################################################

# Create time axis
beg = dt.datetime(1960, 1, 1)
end = dt.datetime(1961, 1, 1)
numb_dt = (end - beg).days * (24 / dt_in)
time_axis = np.asarray([beg + dt.timedelta(hours = (dt_in * k)) for k in range(0, numb_dt)])

# Data array
data = np.zeros((len(time_axis), len(var_names)), dtype = np.float32)

# Generate data
if (not seas_cycle):
    for k in range(0, len(val_mean_amp)): 
        data[:,var_names.index(val_mean_amp[k][0])] = val_mean_amp[k][1]
else:
    year_frac = np.asarray([((time_axis[k] - dt.datetime(time_axis[k].year, 1, 1)).total_seconds() / \
                            ((365. + float(cal.isleap(time_axis[k].year))) * 24. * 3600.)) \
                            for k in range(len(time_axis))], dtype = np.float32)
    cos_func = - np.cos(year_frac * 2 * np.pi) # start: winter -> negative anomaly
    for k in range(0, len(val_mean_amp)):
        if (len(val_mean_amp[k]) == 2): # no seasonal variation
            data[:,var_names.index(val_mean_amp[k][0])] = val_mean_amp[k][1]
        else: # sinusoidal seasonal variations       
            data[:,var_names.index(val_mean_amp[k][0])] = val_mean_amp[k][1] + cos_func * val_mean_amp[k][2]

# Derive values for precipitation phase (0 is fully solid while 1 is fully liquid)
#data[(data[:,var_names.index("TA")] > 273.15),var_names.index("PSUM_PH")] = 1.

# Physical constrictions  
#data[(data[:,var_names.index("TSS")] < 273.15),var_names.index("SMELT")] = 0. # Set melt to zero if TSS < 273.15 K
 
# Changes values for testing ###########################################################################################
#data[480:,var_names.index("TA")] = 233.15
#data[-2:,var_names.index("RH")] -= 0.3
#data[-4:,var_names.index("TSG")] = 273.15
#data[480:,var_names.index("TSS")] = 233.15
#data[-2:,var_names.index("VW")] -= 5.
#data[-4:,var_names.index("OSWR")] += 100.
#data[-4:,var_names.index("ISWR")] += 200.
#data[-4:,var_names.index("ILWR")] += 100.
#data[-4:,var_names.index("PINT")] = 0.0
#data[:,var_names.index("PSUM_PH")] = 1.0
#data[-2:,var_names.index("SMELT")] += 1.0
#data[-4:,var_names.index("SNOWD")] += 1.0
#data[-8:,var_names.index("SUBLI")] -= 0.5
#data[-100:,var_names.index("RHO_HN")] -= 50.

#data[:,var_names.index("SUBLI")] = np.linspace(0., 2., len(time_axis))

#data[:(len(time_axis) / 2),var_names.index("PINT")] = 0.0
#data[:(len(time_axis) / 2),var_names.index("SNOWD")] = 0.0

           
# Calculate mean of certain paramters
TA_MEAN = np.mean(data[:,var_names.index("TA")]) - 273.15 # [degC]
TSS_MEAN = np.mean(data[:,var_names.index("TSS")]) - 273.15 # [degC]
VW_MEAN = np.mean(data[:,var_names.index("VW")]) # [m s-1]
SUBLI_MEAN = np.mean(data[:,var_names.index("SUBLI")]) # [kg m-2 h-1]
PINT_MEAN = np.mean(data[:,var_names.index("PINT")]) # [kg m-2 h-1]
SNOWD_MEAN = np.mean(data[:,var_names.index("SNOWD")]) # [kg m-2 h-1]
SMELT_MEAN = np.mean(data[:,var_names.index("SMELT")]) # [kg m-2 h-1]
    
# Write text file in format .met
station_id = "art_" + fn_part_out
current_date = dt.datetime.now().strftime("%Y-%m-%d")
file = open(path_out + station_id + ".smet","w")
file.write("SMET 1.1 ASCII\n")
file.write("[HEADER]\n")
file.write("station_id       = temp\n") # determines the output file name
file.write("latitude         = 0.0\n")
file.write("longitude        = 0.0\n")
file.write("altitude         = 0.0\n")
file.write("nodata           = -999\n")
file.write("source           = IMAU; CSteger, " + current_date + "\n")
file.write("comment          = AirT [degC]: " + "%.3f" % TA_MEAN + 
                            ", SurfT [degC]: " +  "%.2f" % TSS_MEAN + 
                            ", Wind [m s-1]: " + "%.2f" %  VW_MEAN +
                            ", Prec (+ pos. Snowdrift) [kg m-2 h-1]: " + "%.5f" % PINT_MEAN +
                            ", Melt [kg m-2 h-1]: " + "%.5f" % SMELT_MEAN +
                            ", Snowdrift (only neg.) [kg m-2 h-1]: " + "%.5f" %  SNOWD_MEAN +
                            ", Sublim [kg m-2 h-1]: " + "%.5f" %  SUBLI_MEAN + "\n")                
file.write("fields           = timestamp TA RH TSG TSS VW OSWR ISWR ILWR PINT PSUM_PH SMELT SNOWD SUBLI RHO_HN\n")
file.write("[DATA]\n")
for m in range(0, len(time_axis)):
    file.write(time_axis[m].isoformat()[:-3] + "   " +                    
               "%7.3f" % data[m,var_names.index("TA")] + "   " +
               "%5.3f" % data[m,var_names.index("RH")] + "   " +
               "%7.3f" % data[m,var_names.index("TSG")] + "   " +
               "%7.3f" % data[m,var_names.index("TSS")] + "   " +
               "%6.3f" % data[m,var_names.index("VW")] + "   " +
               "%8.3f" % data[m,var_names.index("OSWR")] + "   " +
               "%8.3f" % data[m,var_names.index("ISWR")] + "   " +
               "%8.3f" % data[m,var_names.index("ILWR")] + "   " +
               "%7.3f" % data[m,var_names.index("PINT")] + "   " + 
               "%5.3f" % data[m,var_names.index("PSUM_PH")] + "   " +
               "%7.3f" % data[m,var_names.index("SMELT")] + "   " +
               "%7.3f" % data[m,var_names.index("SNOWD")] + "   " + 
               "%7.3f" % data[m,var_names.index("SUBLI")] + "   " +
               "%7.3f" % data[m,var_names.index("RHO_HN")] + "\n")                 
file.close()