
PARALLEL=OPENMP
NCORES=16

BEGIN="2017-09-01T02:00"
END="2017-08-31T23:00"
# "2018-08-31T23:00"
PROG_ROOTDIR=~/usr/bin	#if you do not want to use Alpine3D from PATH, then point to where it is

EXE="${PROG_ROOTDIR}/alpine3d"
if [ ! -f ${EXE} ]; then
	EXE=`which alpine3d`
fi

A3D_CMD="${EXE} \
--iofile=./io.ini \
--enable-eb  \
--np-ebalance=${NCORES} \
--np-snowpack=${NCORES} \
--startdate=${BEGIN} --enddate=${END}"

${A3D_CMD} 

echo "Done Alpine3D Simulation. Return code=$ret"
exit $ret
