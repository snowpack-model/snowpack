# No Shebang, execute shell_file in current shell (where path and libraries are defined)
set -ex

##########################################################################################
# Settings
##########################################################################################

point_current=14
run_name=Charalampidis
# set correct meteo-file in config.ini and config_spinup.ini (-> STATION1)!!!

# Output file name
fn_out=${run_name}_${point_current}

# Initial snow profile
cp input/${run_name}_${point_current}.sno input/temp.sno

# Number of repeats
nor=1 # 8

##########################################################################################
# Perform spin-up (1960 - 1979; 20 years) and actual run (1960 - 2014; 55 years)
##########################################################################################

for i in $(seq 1 $nor);
do

	# Add expected mass loss during spin-up to snow file
	python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/ECMWF/SP_add_ice.py input/temp.sno $point_current reference
	
	echo "#########################################################################################"
	echo " Spin up number ${i}"
	echo "#########################################################################################"
	snowpack -c cfgfiles/config_spinup.ini -e 1979-12-31T21:00

	# Update snow file
	cp output/temp.sno input/temp.sno
	#cp output/temp.sno output/${fn_out}_spinup_${i}.sno # backup snow file
	python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/ECMWF/SP_change_date.py input/temp.sno

	# Try to remove ice at bottom if firn column gets to thick
	python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/ECMWF/SP_remove_ice.py input/temp.sno
	
done

# Add expected mass loss during simulation to snow file
python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/ECMWF/SP_add_ice.py input/temp.sno $point_current entire

# Run Model
snowpack -c cfgfiles/config.ini -e 2014-12-31T21:00

##########################################################################################
# Convert output
##########################################################################################

# Remove irrelevant output files (only if exist)
rm -f output/temp.haz
rm -f output/temp_NO_EXP.ini
rm -f output/temp.sno

# Rename relevant output files
mv output/temp.met output/$fn_out.met
mv output/temp.pro output/$fn_out.pro

# Convert .met and .pro to NetCDF-files
python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/SP_out_met_to_NetCDF.py output/$fn_out.met
rm output/$fn_out.met
python $HOME/Dropbox/IMAU/Scripts/Python/SNOWPACK/SP_out_pro_to_NetCDF.py output/$fn_out.pro
rm output/$fn_out.pro