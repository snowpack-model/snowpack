# Description: Generates SNOWPACK input files from previous SNOWPACK run
#
# Author: Christian Steger, July 2016

# Load modules
import os
import sys
import platform
import numpy as np
import datetime as dt
from netCDF4 import Dataset, num2date

# List with platform names
platf_nam = [["iMac-Steger.local", "imac-steger.soliscom.uu.nl", "hst44188.phys.uu.nl"],        # Desktop-Mac
             ["Christians-MacBook-Pro-2.local", "ChristiansMBP2.home", "christiansmbp2.home"],  # Macbook
             ["cca-login1", "cca-login3"]]                                                      # ECMWF

# Paths to folders (dependent on platform)
if (platform.node() in platf_nam[0]):
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in_RACMO = root_IMAU + "Data/Greenland/RACMO2.3/"
    path_in_SP = [root_IMAU + "Data/Greenland/SP/ZGRN_V5/2D/",
                  root_IMAU + "Data/Greenland/SP/ZGRN_V5_BO/2D/"]
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/sno/"
elif (platform.node() in platf_nam[1]):  
    root_IMAU = os.getenv("HOME") + "/Dropbox/IMAU/"
    path_in_RACMO = root_IMAU + "Data/Greenland/RACMO2.3/"
    path_in_SP = [root_IMAU + "Data/Greenland/SP/ZGRN_V5/2D/",
                  root_IMAU + "Data/Greenland/SP/ZGRN_V5_BO/2D/"]
    path_out = root_IMAU + "Data/Greenland/RACMO2.3/SP_INPUT/sno/"
elif (platform.node() in platf_nam[2]):
    path_in_RACMO = "/scratch/ms/nl/rucs/raw_Data/"
    path_in_SP = ["/scratch/ms/nl/rucs/SP_Data/OUTPUT/ZGRN_V5/NetCDF/",
                  "/scratch/ms/nl/rucs/SP_Data/OUTPUT/ZGRN_V5_BO/NetCDF/"]
    path_out = "/scratch/ms/nl/rucs/SP_Data/INPUT/sno/"
else:
    sys.exit("Unknown platform")

########################################################################################################################
# Settings
########################################################################################################################

# Path and name of file with indices (rlat, rlon; on RACMO2.3 grid)
path_file_ind = [root_IMAU + "FDM_SP_runs/ZGRN_V5_checkerboard_ind.txt",
                 root_IMAU + "FDM_SP_runs/ZGRN_V5_border_ind.txt"] # new run: all grid boxes together
#path_file_ind = ["/home/ms/nl/rucs/FDM_SP_runs/ZGRN_V5_checkerboard_ind.txt",
#                 "/home/ms/nl/rucs/FDM_SP_runs/ZGRN_V5_border_ind.txt"]

fn_part_in = ["ZGRN_V5", "ZGRN_V5_BO"]
fn_part_out = "ZGRN_V5_all" # output name

# Ice block settings
elem_thickn = 1.95 # element thickness for ice added at the bottom [m]
firn_thick_range = [60., 120.] # minimal and maximal initial firn thickness [m]
theta_ice_thresh = 0.99 # threshold for volumetric ice content (907.83 kg m-3) [-]

# Constants
dens_water = 1000. # density of water [kg m-3]

########################################################################################################################
# Load general data, check if necessary file exist
########################################################################################################################

# Read text file with indices of grid boxes on RACMO2.3 domain
indices = []
for k in path_file_ind:
    indices.append(np.loadtxt(k, usecols = [0, 1], dtype = int)) # (column 0: rlat, column 1: rlon)

# Load geographical coordinates and topography
ncfile = Dataset(path_in_RACMO + "RACMO23_masks_ZGRN11.nc", "r")
lon = ncfile.variables["lon"][:]
lat = ncfile.variables["lat"][:]
topo = ncfile.variables["topography"][:]
ncfile.close()

# Loop over all indices to check if all files exist
#for k in range(len(indices)):
#    for m in range(len(indices[k])):
#        if not os.path.isfile(path_in_SP[k] + fn_part_in[k] + "_" + str(m + 1) + "_pro.nc"):
#            sys.exit("Error: File " + fn_part_in[k] + "_" + str(m + 1) + "_pro.nc does not exist")

########################################################################################################################
# Generate .sno files
########################################################################################################################

# Load time axis of previous SNOWPACK output, get indices of last time step in year 1979
ncfile = Dataset(path_in_SP[0] + fn_part_in[0] + "_" + str(434) + "_pro.nc", "r")
time_axis_SP = num2date(ncfile.variables["time"][:], units = ncfile.variables["time"].units, 
                        calendar = ncfile.variables["time"].calendar)
ncfile.close()
year_SP = np.asarray([k.year for k in time_axis_SP])
ind_dt = np.where(year_SP == 1979)[0][-1] # last time step in year 1979

# Load necessary variables and create .sno input files
var_names_SP = ["layerthickn", "temper", "water", "dendric", "spheric", "bondsize", "grainsize", "ice", "air"]
#for k in range(len(indices)):
for k in [0]:
#    for m in range(len(indices[k])):
    for m in [273, 2101, 2166]:
#    for m in [314, 441]:
                
        # Load SNOWPACK file from previous run
        ncfile = Dataset(path_in_SP[k] + fn_part_in[k] + "_" + str(m + 1) + "_pro.nc", "r")
        if (type(ncfile.variables["layerthickn"][:,ind_dt]).__name__ == "MaskedArray"):
            numb_elem = np.invert(ncfile.variables["layerthickn"][:,ind_dt].mask).sum()  
        else:
            numb_elem = ncfile.variables["layerthickn"].shape[0]
        data_SP = np.empty((len(var_names_SP), numb_elem), dtype = np.float32)
        for n in range(len(var_names_SP)):
            data_SP[n,:] = ncfile.variables[var_names_SP[n]][:numb_elem,ind_dt]
        ncfile.close()

        # Conversion of units
        data_SP[var_names_SP.index("temper"),:] += 273.15 # [degC -> K]
        data_SP[var_names_SP.index("water"),:] /= 100. # [% -> -]
        data_SP[var_names_SP.index("ice"),:] /= 100. # [% -> -]
        data_SP[var_names_SP.index("air"),:] /= 100. # [% -> -]

        # Add ice at bottom if necessary
        firn_thick = data_SP[var_names_SP.index("layerthickn"),:].sum()
        if (firn_thick < firn_thick_range[0]): # Add ice at the bottom
            numb_elem_add = np.ceil(np.abs(firn_thick_range[0] - firn_thick) / elem_thickn).astype(int)
            data_SP_add = np.copy(data_SP[:,0])
            data_SP_add[var_names_SP.index("layerthickn")] = elem_thickn
            data_SP_add[var_names_SP.index("water")] = 0.
            data_SP_add[var_names_SP.index("air")] = 0.
            data_SP_add[var_names_SP.index("ice")] = 1.
            data_SP_add = np.repeat(data_SP_add[:,np.newaxis], numb_elem_add, axis = 1)
            data_SP = np.hstack((data_SP_add, data_SP))
            numb_elem += numb_elem_add
        
        # Try to remove ice at bottom
        firn_thick = data_SP[var_names_SP.index("layerthickn"),:].sum()
        numb_elem_rem = 0
        while (firn_thick > firn_thick_range[1]): # Try to remove ice at the bottom
            if (data_SP[var_names_SP.index("ice"),0] > theta_ice_thresh):
                firn_thick -= data_SP[var_names_SP.index("layerthickn"),0]
                data_SP = np.delete(data_SP, 0, axis = 1)
                numb_elem_rem += 1
            else:
                break # non-ice element
        numb_elem -= numb_elem_rem

        # Flip data (start at surface)
        data_SP = np.fliplr(data_SP)

        # Check total volumetric content (threshold 1 %; equal to ElementData::checkVolContent())
        ind = [var_names_SP.index("ice"), var_names_SP.index("water"), var_names_SP.index("air")]
        theta_tot = data_SP[var_names_SP.index("ice"),:] + \
                    data_SP[var_names_SP.index("water"),:] + \
                    data_SP[var_names_SP.index("air"),:]
        mask = np.logical_or(theta_tot >= 1.01, theta_tot <= 0.99)         
        if (np.any(mask)):
                for n in np.where(mask)[0]:
                    if (theta_tot[n] < 1.):
                        data_SP[var_names_SP.index("air"),n] = (1.0 - data_SP[var_names_SP.index("ice"),n] - \
                                                                data_SP[var_names_SP.index("water"),n])
                    else:
                        data_SP[ind[np.argmax(data_SP[ind,n])],n] -= (theta_tot[n] - 1.)

        # Check for negative values (could be caused by the checking routine above)
        if (np.any(data_SP < 0.)):
            sys.exit("Error: At least 1 negative value found in file " + fn_part_in[k] + "_" + str(m + 1) + "_pro.nc")

        # Array with element age (age difference layers: 10 day)
        time_axis = [(dt.datetime(1960, 1, 1) - dt.timedelta(days = (n + 1) * 10)) for n in range(numb_elem)]

        # Write text file in format .met
        out_file_numb = (k * len(indices[0])) + m + 1
        fn_out = fn_part_out + "_" + str(out_file_numb) + ".sno"
        tot_thick = data_SP[var_names_SP.index("layerthickn"),:].sum()
        file = open(path_out + fn_out,"w")
        file.write("SMET 1.1 ASCII\n")
        file.write("[HEADER]\n")
        file.write("station_id       = temp\n")
        file.write("station_name     = temp\n")
        file.write("latitude         = " + "%.4f" % lat[indices[k][m,0],indices[k][m,1]] + "\n")
        file.write("longitude        = " + "%.4f" % lon[indices[k][m,0],indices[k][m,1]] + "\n")
        file.write("altitude         = " + "%.2f" % topo[indices[k][m,0],indices[k][m,1]] + "\n")
        file.write("nodata           = -999\n")
        file.write("source           = IMAU; CSteger, " + dt.datetime.now().strftime("%Y-%m-%d") + "\n")
        file.write("ProfileDate      = " + dt.datetime(1960, 1, 1, 1).isoformat()[:-3] + "\n")
        file.write("HS_Last          = " + "%.6f" % tot_thick + "\n") # [m]
        file.write("SlopeAngle       = 0.00\n")
        file.write("SlopeAzi         = 0.00\n")
        file.write("nSoilLayerData   = 0\n")
        file.write("nSnowLayerData   = " + str(numb_elem) + "\n")
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
        for n in reversed(range(numb_elem)):
            file.write(time_axis[n].isoformat()[:-3] + "   " +
                        "%6.3f" % data_SP[var_names_SP.index("layerthickn"),n] + "   " +  # layer-thickness [m]
                        "%7.3f" % data_SP[var_names_SP.index("temper"),n] + "   " +       # temperature [K]
                        "%5.3f" % data_SP[var_names_SP.index("ice"),n] + "   " +          # volume fraction (ice) [-]
                        "%5.3f" % data_SP[var_names_SP.index("water"),n] + "   " +        # volume fraction (water) [-]
                        "%5.3f" % data_SP[var_names_SP.index("air"),n] + "   " +          # volume fraction (air) [-]
                        "%5.3f" % 0.0 + "   " +                                           # volume fraction (soil) [-]
                        "%8.3f" % 2600.0 + "   " +                                        # soil dry dens. [kg m-3]
                        "%6.3f" % 2.6 + "   " +                                           # soil heat cond. [W m-1 K-1]
                        "%8.3f" % 1100.00 + "   " +                                       # soil heat capacity [J m-3 K-1]
                        "%9.3f" % data_SP[var_names_SP.index("grainsize"),n] + "   " +    # grain radius [mm]
                        "%9.3f" % data_SP[var_names_SP.index("bondsize"),n] + "   " +     # bond radius [mm]
                        "%5.3f" % data_SP[var_names_SP.index("dendric"),n] + "   " +      # dendricity (0: old snow) [-]
                        "%5.3f" % data_SP[var_names_SP.index("spheric"),n] + "   " +      # sphericity (1: rounded) [-]
                        str(0) + "   " +                                                  # grain marker [-] (set to 7 -> glacier ice ?)
                        str(0.0) + "   " +                                                # Mass of surface hoar [kg m-2]
                        str(1) + "   " +                                                  # Numb. of fin. elem. in layer [-]
                        str(0.0) + "   " +                                                # Stress rate [Pa s-1]
                        str(0.0) + "   " +  "\n")                                         # Keep track of metamorphism [-]
        file.close()
    
        print "Input file " + fn_out + " created"

# Comments
# - Grain marker: value same as for new snow -> Snowpack::setHydrometeorMicrostructure
# - Stress rate: value same as for new snow -> Snowpack::fillNewSnowElement
# - Keep track of metamorphism: value same as for new snow -> Snowpack::setHydrometeorMicrostructure