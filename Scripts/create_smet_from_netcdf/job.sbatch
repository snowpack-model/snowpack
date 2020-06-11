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

