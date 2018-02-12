#!/bin/bash
#This is the launching script for Alpine3D. Please make sure the user section matches your needs!
########################### Slurm directives: edit for the Slurm job manager
#SBATCH --job-name=urumqi
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
########################### End of Slurm directives
##either N or MPI, MPI_SGE or OPENMP
PARALLEL=N

##When using MPI (NOT MPI_SGE!), specify the machinefile you want to use. The machinefile format allows one
##hostname per line followed optionally by a colon and the number of processes to run on that
##host, e. g. 192.168.0.11:2 means: use 2 cores of the machine at address 192.168.0.11
##use as many hosts as you deem reasonable, not defining a machinefile results in a local run.
##Ultimately only NPROC many processes are used.
MACHINEFILE=
NPROC=8

##when using openmp but NOT relying on a job manager, how many cores to use (otherwise it will be ignored)
NCORES=2

BEGIN="2001-10-02T01:00"
END="2013-10-01T01:00"
PROG_ROOTDIR=../bin	#if you do not want to use Alpine3D from PATH, then point to where it is
#should screen messages be redirected to stdouterr.log or printed directly on the screen?
REDIRECT_LOGS=Y
########################## END OF USER CONFIGURATION

export DYLD_FALLBACK_LIBRARY_PATH=${PROG_ROOTDIR}:${DYLD_FALLBACK_LIBRARY_PATH}	#for osX
export LD_LIBRARY_PATH=${PROG_ROOTDIR}:${LD_LIBRARY_PATH}	#for Linux
EXE="${PROG_ROOTDIR}/alpine3d"
if [ ! -f ${EXE} ]; then
	EXE=`which alpine3d`
fi
N_EB=1
N_SN=1

#to combine OPENMP and MPI, run as MPI but with N_EB & N_SN > 1
if [ "${PARALLEL}" == "MPI" ]; then
	echo "Running with MPI"
	MPIEXEC=${MPIEXEC:="mpiexec"}
	MFILE=${MACHINEFILE:+"-machinefile ${MACHINEFILE}"}
	EXE="${MPIEXEC} -n ${NPROC} ${MFILE} ${EXE}"
	N_EB=1
	N_SN=1
elif [ "${PARALLEL}" == "MPI_SGE" ]; then
	echo "Running with MPI under SGE"
        MPIEXEC=${MPIEXEC:="mpiexec"}
        if [ $SLURM_TASKS_PER_NODE ]; then
        	export NSLOTS=$SLURM_TASKS_PER_NODE
        fi
        EXE="${MPIEXEC} -np ${NSLOTS} ${EXE}"
        N_EB=1
        N_SN=1
elif [ "${PARALLEL}" == "OPENMP" ]; then
	echo "Running with OPENMP"
	if [ $SLURM_CPUS_PER_TASK ]; then
		export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
	else
		export OMP_NUM_THREADS=${NSLOTS:=$NCORES}
	fi
	N_EB=$OMP_NUM_THREADS
	N_SN=$OMP_NUM_THREADS
else
	unset OMP_NUM_THREADS
	echo "Running sequentially"
fi

##Below, always use a double pound sign for comments, otherwise SGE believes this is a job manager command and refuses to run
##TOOL="valgrind --leak-check=full --show-reachable=yes --leak-resolution=high --undef-value-errors=yes --track-origins=yes --log-file=valgrind.log"

A3D_CMD="${TOOL} ${EXE} \
--iofile=./io.ini \
--enable-eb  \
--np-ebalance=${N_EB} \
--np-snowpack=${N_SN} \
--startdate=${BEGIN} --enddate=${END}"

date
if [[ ("${REDIRECT_LOGS}" == "Y") ||  ("${REDIRECT_LOGS}" == "y") ]]; then
	${A3D_CMD} > stdouterr.log 2>&1 $*
else
	${A3D_CMD} 2>&1 $*
fi
ret=$?

echo "Done Alpine3D Simulation. Return code=$ret"
date
echo
exit $ret
