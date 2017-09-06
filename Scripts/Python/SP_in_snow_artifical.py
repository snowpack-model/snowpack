# Description: Generates textfile with artifical snow input data for SNOWPACK
#
# Generates a homogenous block of ice (where only temperature has a seasonal amplitude)
#
# To do:
# - Use better estimates for rg and rb
# - Use thermal diffusivity dependent on temperature to compute temperature profile
#
# Author: Christian Steger, November 2015

# Load modules
import sys
import platform
import datetime as dt
import numpy as np

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"], # Desktop-Mac
             ["Christians-MacBook-Pro-2.local"]]                                         # Macbook

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
fn_out = "art_test"

# Total layer thickness, temperature
total_thickn = 50. # total thickness of ice block [m]
tskin_mean = 273.15 # mean surface temperature [K]
temp_prof = "linear" # either constant, linear or seas_cycle
tskin_amp = 0.0 # surface temperature amplitude [K]
soil_elem = 0 # number of soil elements

# Range of element thicknesses and microproperties (surface -> bottom) (linear interpolated)
desity_range = [369., 917.] # [kg m-3]
temp_range = [273.15, 250.] # [K]
elem_thickn = 0.1 # element thickness [m]
rg = [16.25, 21.42] # grain size range [mm]
rb = [3.85, 4.15] # bond size range [mm]

########################################################################################################################
# Constant settings
########################################################################################################################

# Snow profile date
date_sp = dt.datetime(1960, 1, 1)

# Constants
rho_i = 917.0 # ice density [kg m-3]
c_i = 2100.0 # specific heat ice [J kg-1 K-1]
k_i = 2.2 # conductivity ice [W m-1 K-1]

########################################################################################################################
# Generate artifical snow data and save to SNOWPACK textfile
########################################################################################################################

# Element thickness array
numb_elem = int(np.ceil(total_thickn / elem_thickn))
elem_thickn = np.repeat(elem_thickn, numb_elem)

# Array with element age
time_axis = [(date_sp - dt.timedelta(days = k)) for k in range(numb_elem)] # age difference layers: 1 day

# Array with temperature
elem_temp = np.empty(numb_elem, dtype = np.float32)
if (temp_prof == "constant"):
    elem_temp.fill(tskin_mean)
if (temp_prof == "linear"):    
    elem_temp = np.linspace(temp_range[0], temp_range[1], numb_elem)
elif (temp_prof == "seas_cycle"):
    elem_temp.fill(tskin_mean)
    bound = np.append(0, np.cumsum(elem_thickn)) # element boundaries
    depth = bound[:-1] + elem_thickn / 2.
    alpha = k_i / (rho_i * c_i) # thermal diffusivity [m2 s-1]
    omega = 1. / (float((dt.datetime(date_sp.year, 12, 31) - dt.datetime(date_sp.year, 1, 1)).days + 1) * 3600. * 24.)
    t = (date_sp - dt.datetime(date_sp.year, 1, 1)).total_seconds() # seconds since start of year [s]
    phas_shift = -(4. / 12. * 2 * np.pi) # minimal temperature is at end of January [rad]
    elem_temp += tskin_amp * np.exp(-depth * np.sqrt((np.pi * omega) / alpha)) * \
                 np.sin(2 * np.pi * omega * t - depth * np.sqrt((np.pi * omega) / alpha) + phas_shift)
else:
    sys.exit("Error: Unknown value for temp_prof")

# Array with grain size (rg) and bond size (rb)
elem_rg = np.linspace(rg[0], rg[1], numb_elem)
elem_rb = np.linspace(rb[0], rb[1], numb_elem)

# Set volumetric contents
theta_ice = np.linspace(desity_range[0], desity_range[1], numb_elem) / rho_i
theta_water = np.empty(numb_elem, dtype = np.float32)
theta_water.fill(0.0)
theta_soil = np.empty(numb_elem, dtype = np.float32)
theta_soil.fill(0.0)
theta_air = np.empty(numb_elem, dtype = np.float32)
theta_air = 1. - (theta_ice + theta_water + theta_soil)

# Set properties of soil elements
if (soil_elem > 0):
    elem_rg[-soil_elem:] = 10000.
    elem_rb[-soil_elem:] = 0.
    theta_ice[-soil_elem:] = 0.050
    theta_soil[-soil_elem:] = 0.950

# Write text file in format .met
file = open(path_out + fn_out + ".sno","w")
file.write("SMET 1.1 ASCII\n")
file.write("[HEADER]\n")
file.write("station_id       = temp\n")
file.write("station_name     = temp\n")
file.write("latitude         = 0.0\n")
file.write("longitude        = 0.0\n")
file.write("altitude         = 0.0\n")
file.write("nodata           = -999\n")
file.write("source           = IMAU; CSteger, " + dt.datetime.now().strftime("%Y-%m-%d") + "\n")
file.write("ProfileDate      = " + (date_sp + dt.timedelta(hours = 1)).isoformat()[:-3] + "\n") # add 1 h to start date
file.write("HS_Last          = " + "%.6f" % np.sum(elem_thickn) + "\n") # [m]
file.write("SlopeAngle       = 0.00\n")
file.write("SlopeAzi         = 0.00\n")
file.write("nSoilLayerData   = " + str(soil_elem) + "\n")
file.write("nSnowLayerData   = " + str(numb_elem - soil_elem) + "\n")
file.write("SoilAlbedo       = 0.09\n")
file.write("BareSoil_z0      = 0.020\n")
file.write("CanopyHeight     = 0.00\n")
file.write("CanopyLeafAreaIndex = 0.00\n")
file.write("CanopyDirectThroughfall = 1.00\n")
file.write("WindScalingFactor = 1.00\n")
file.write("ErosionLevel     = "  + str(numb_elem - 1) + "\n")
file.write("TimeCountDeltaHS = 0.000000\n")                          
file.write("fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo\n")
file.write("[DATA]\n")
for k in reversed(range(numb_elem)):
    file.write(time_axis[k].isoformat()[:-3] + "   " +                    
               "%6.3f" % elem_thickn[k] + "   " +  # element thickness [m]
               "%7.3f" % elem_temp[k] + "   " +    # element temperature [K]
               "%5.3f" % theta_ice[k] + "   " +    # volume fraction ice [-]
               "%5.3f" % theta_water[k] + "   " +  # volume fraction water [-]
               "%5.3f" % theta_air[k] + "   " +    # volume fraction void [-]
               "%5.3f" % theta_soil[k] + "   " +   # volume fraction soil [-]
               "%8.3f" % 2600.0 + "   " +          # soil dry density [kg m-3] (Rho_S)
               "%6.3f" % 2.6 + "   " +             # soil heat conductivity [W m-1 K-1] (Conduc_S)
               "%8.3f" % 1100.00 + "   " +         # soil heat capacity [J m-3 K-1] (HeatCapac_S)
               "%9.3f" % elem_rg[k] + "   " +      # grain size [mm] (rg)
               "%9.3f" % elem_rb[k] + "   " +      # bond size [mm] (rb)
               "%5.3f" % 0.0 + "   " +             # dendricity [-] (dd) (0 = none, 1 = newsnow)
               "%5.3f" % 1.0 + "   " +             # sphericity [-] (sp) (1 = round, 0 = angular)
               str(0) + "   " +                    # grain marker (value same as for new snow -> Snowpack::setHydrometeorMicrostructure) [-] (mk) 
               str(0.0) + "   " +                  # Mass of surface hoar [kg m-2] (mass_hoar)
               str(1) + "   " +                    # Number of finite elements in the layer [-] (DataClasses.h -> LayerData) (ne)
               str(0.0) + "   " +                  # Stress rate (value same as for new snow -> Snowpack::fillNewSnowElement) [Pa s-1] (CDot)
               str(0.0) + "   " +  "\n")           # Keep track of metamorphism (value same as for new snow -> Snowpack::setHydrometeorMicrostructure) [-] (metamo)                
file.close()