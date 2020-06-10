#!/bin/bash

#SBATCH --nodes=1                         # Number of requested nodes
#SBATCH --account=ucb164_summit1 
#SBATCH --time=04:00:00                   # Max wall time
#SBATCH --qos=normal                      # Specify testing QOS
#SBATCH --partition=shas                  # Specify Summit haswell nodes
#SBATCH --ntasks=1                       # Number of tasks per job
#SBATCH --array=1-40
##SBATCH --array=1-1

#SBATCH --job-name=smet_creation         # Job submission name
#SBATCH --output=./sbatch_out_files/%x.%j.out     # Output file name with Job ID
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.keenan@colorado.edu

# Module loading and basic setup
mkdir -p log
export LD_LIBRARY_PATH=$(pwd)/usr/lib/:${LD_LIBRARY_PATH}
ml intel; ml impi; ml proj; ml cmake; ml netcdf

ldd ./data_converter

# Define variable year and create output directory 
yr=$(echo ${SLURM_ARRAY_TASK_ID} | awk '{print 1980+($1-1)}')
rm -r output/output_${yr}
mkdir -p output/output_${yr}
cd io_files

# Create new ini file for each year
sed 's/\.\.\/\.\.\/MERRA2\/\//\.\.\/\.\.\/MERRA2\/'${yr}'\//' io.ini > io_${yr}.ini
sed -i 's/output\//output\/output_'${yr}'\//' io_${yr}.ini

# Create .smet files for each year
../data_converter ${yr}-01-01T00:30:00 ${yr}-12-31T23:30:00 60 io_${yr}.ini > ../log/${yr}.log 2>&1
