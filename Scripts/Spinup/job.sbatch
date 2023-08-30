#!/bin/bash

#
# For SLURM
#
#SBATCH --nodes=1				#	Number of requested nodes per task
#SBATCH --account=xxxxxxxx			#	Account
#SBATCH --time=00:10:00			#	Max wall time
#SBATCH --qos=normal				#	Specify QOS
#SBATCH --partition=shas			#	Specify Summit haswell nodes
#SBATCH --ntasks=1				#	Number of tasks per job
#SBATCH --job-name=spinup			#	Job submission name
#SBATCH --output=./sbatch_out_files/%x.%j.out	#	Output file name with Job ID
#SBATCH --mail-type=ALL
#SBATCH --mail-user=someuser@someaddress.edu
#SBATCH --array=1-1

#
# For PBS
#
#PBS -S /bin/bash				#	Set bash as shell
#PBS -N spinup					#	Job name
#PBS -A PXXXXXXXX				#	Account
#PBS -l walltime=0:10:00			#	Max wall time
#PBS -q regular				#	Specify compute nodes type
#PBS -j oe					#	stdout/stderr
#PBS -m abe					#	email settings
#PBS -M someuser@someaddress.edu
#PBS -l select=1:ncpus=1:mpiprocs=0		#	Number of request nodes per task
#PBS -J 1-1

#
# Check for job arrays:
#
if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
	command1=$(sed -n ${SLURM_ARRAY_TASK_ID}p to_exec.lst)
elif [ -n "${PBS_ARRAY_INDEX}" ]; then
	command1=$(sed -n ${PBS_ARRAY_INDEX}p to_exec.lst)
else
	command1=""
fi

#
# Execute job
#
if [ -n "${command1}" ]; then
	# Case using job arrays with a job manager
	eval ${command1}
	echo ${command1} >> finished.lst
else
	# Other cases
	bash to_exec.lst
fi
