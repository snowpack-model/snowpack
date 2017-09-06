#!/bin/ksh -x

# Read in arguments
cpu_numb=$1
host_name=$2

# Read in current point number
p2cpuinput=${p2cpufiles}points_CPU_${cpu_numb}.txt
line1=`head -1 $p2cpuinput | tail -1`
rem_numb=`echo $line1 | cut -d' ' -f1`
point_current=`echo $line1 | cut -d' ' -f2`

# Read in spin-up number of point
line2=`head -$point_current $p2input_NoR | tail -1`
nor=`echo $line2 | cut -d' ' -f1`

# Generate loadscript
loadscript=${p2loads}loadscript_$point_current
cat <<EOS_L > $loadscript
#!/bin/ksh -f
#PBS -S /usr/bin/ksh
#PBS -N dens_${point_current}
#PBS -q ns
#PBS -j oe
#PBS -m a
#PBS -V
#PBS -l EC_billing_account=nlcko
#PBS -o /scratch/ms/nl/rucs/SP_Data/logfiles/log_SNOWPACK_${point_current}.out
#PBS -l EC_memory_per_task=1000mb
#PBS -l EC_threads_per_task=1
#PBS -l walltime=48:00:00

# Move to directory of SNOWPACK executable
cd $p2exe

echo "Start loadscript jobno $cpu_numb"
echo "Start simulation of point $point_current"

# Generate directory for storing temporary data
mkdir -p /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input
mkdir -p /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/output

# Copy snow file to temporary directory
cp /scratch/ms/nl/rucs/SP_Data/INPUT/sno/ZGRN_V5_all_${point_current}.sno /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno ##### set before run #####

# Perform spin-up's (1960 - 1979; 20 years)
for i in {1..$nor}
do

	# Add expected mass loss during spin-up to snow file
	python /home/ms/nl/rucs/Scripts/Python/SP_add_ice.py /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno $point_current reference
	
	echo "#########################################################################################"
	echo "Point ${point_current}: Spin up number ${i}"
	echo "#########################################################################################"
	snowpack -c /scratch/ms/nl/rucs/SP_Data/Run_temp/cfgfiles/config_spinup_${point_current}.ini -e 1979-12-31T21:00

	# Update snow file
	cp /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/output/temp.sno /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno
	python /home/ms/nl/rucs/Scripts/Python/SP_change_date.py /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno

	# Try to remove ice at bottom if firn column gets to thick
	python /home/ms/nl/rucs/Scripts/Python/SP_remove_ice.py /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno
	
done

# Add expected mass loss during simulation to snow file
python /home/ms/nl/rucs/Scripts/Python/SP_add_ice.py /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/input/temp.sno $point_current entire

# Perform simulation (1960 - 2014; 55 years)
echo "#########################################################################################"
echo "Point ${point_current}: Simulation over 1960 - 2014"
echo "#########################################################################################"
snowpack -c /scratch/ms/nl/rucs/SP_Data/Run_temp/cfgfiles/config_${point_current}.ini -e 2014-12-31T21:00

# Save final files
mv /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/output/temp.met ${p2output}${filename_part1}_${point_current}.met
mv /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/output/temp.pro ${p2output}${filename_part1}_${point_current}.pro
mv /scratch/ms/nl/rucs/SP_Data/Run_temp/point_${point_current}/output/temp.sno ${p2output}${filename_part1}_${point_current}.sno

# Move back to start directory
cd /home/ms/nl/rucs/Snowpack/Run_dir

# Remove previous point from mark_place file and run again if there are points remaining
if [ $rem_numb -gt 1 ] ; then
  python /home/ms/nl/rucs/Scripts/Python/SP_update_CPU_work_files.py $cpu_numb $host_name
 ./SP_multiple_load.sh $cpu_numb $host_name
fi

EOS_L
 # Put the script in the queue
 qsub $loadscript
