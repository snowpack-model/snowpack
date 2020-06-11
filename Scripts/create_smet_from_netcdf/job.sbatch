#!/bin/bash

#SBATCH --nodes=1                         # Number of requested nodes
#SBATCH --account=ucb164_summit1 
#SBATCH --time=00:10:00                   # Max wall time
#SBATCH --qos=normal                      # Specify testing QOS
#SBATCH --partition=shas                  # Specify Summit haswell nodes
#SBATCH --ntasks=1                       # Number of tasks per job
##SBATCH --array=1-40
#SBATCH --array=1-1

#SBATCH --job-name=smet_creation         # Job submission name
#SBATCH --output=./sbatch_out_files/%x.%j.out     # Output file name with Job ID
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.keenan@colorado.edu

command1=$(sed -n ${SLURM_ARRAY_TASK_ID}p to_exec.lst)
eval ${command1}
