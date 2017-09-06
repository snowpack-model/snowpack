#!/bin/ksh -x

# Host name
host_name="cca"

# Number of cpu's (limit to number of points if smaller than number of cpu's!!!)
noc=7

# File with points to run of project
#file_name="ZGRN_V5_all_sel_points_"${host_name}".txt"
file_name="test_points.txt"
#file_name="test_cca.txt"

# File with NoR from all points of project
#file_name_NoR="ZGRN_V5_all_NoR_and_mc.txt"
file_name_NoR="ZGRN_V5_all_NoR_and_mc_test_lim_2.txt"

# Output name (without _ at end)
#run_name="ZGRN_V5_all_SPr48"
run_name="test"

##########################################################################################

# Modifications also necessary in:
# - SP_multiple_load.sh
# 		- path and name of input.sno-files
# - SP_cfgfiles_generator.py
#		- path and name of input .smet-files
#		- general settings
# - SP_update_CPU_work_files.py
#		- root path to CPU-files
# - SP_add_ice.py
#		- path and name of file with mass changes

# Run SNOWPACK simultaneously on cca and ccb
# - ensure that there are no overlapping points (same point assigned to cca and ccb)

# Data which can be removed before run (do not remove when SNOWPACK is running on other host!)
# - /home/ms/nl/rucs/Snowpack/DATA -> folder cca and ccb
# - /scratch/ms/nl/rucs/SP_Data/Run_temp/cfgfiles -> configuration files
# - /scratch/ms/nl/rucs/SP_Data/Run_temp -> folders point_*
# - /scratch/ms/nl/rucs/SP_Data/loadscripts -> loadscripts
# - /scratch/ms/nl/rucs/SP_Data/logfiles -> logfiles
# - /scratch/ms/nl/rucs/SP_Data/OUTPUT -> old data (.met, .pro, .sno)

##########################################################################################

# Export path and input file names
export p2input=$HOME/Snowpack/DATA/$file_name
export p2input_NoR=$HOME/Snowpack/DATA/$file_name_NoR

# Path for loadscripts (location on scratch)
export p2loads=$SCRATCH/SP_Data/loadscripts/
if [ ! -r $p2loads ] ; then
  mkdir -p $p2loads
fi

# Path for logfiles (location on scratch)
export p2logs=$SCRATCH/SP_Data/logfiles/
if [ ! -r $p2logs ] ; then
  mkdir -p $p2logs
fi

# Path for points_CPU-files
export p2cpufiles=$HOME/Snowpack/DATA/${host_name}/

# Output path name
export p2output=$SCRATCH/SP_Data/OUTPUT/

# Output file name
export filename_part1=${run_name}

# Path to executable
export p2exe=$HOME/Snowpack/snowpackbin/bin/

# Make CPU files (in Snowpack/DATA/${host_name})
python /home/ms/nl/rucs/Scripts/Python/SP_make_CPU_work_files.py $noc $host_name $file_name $file_name_NoR

# Generate configuration files for SNOWPACK (one file per location)
python /home/ms/nl/rucs/Scripts/Python/SP_config_file_generator.py $file_name

# Generate loadscripts for the different CPU's
k=1
while [[ $k -le $noc ]] ; do
 cpu_numb=$k
 ./SP_multiple_load.sh $cpu_numb $host_name
  
 k=$(($k + 1))
done
