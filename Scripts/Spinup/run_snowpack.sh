shopt -s expand_aliases		# Make sure aliases work in non-interactive shells

# Load settings
source ./spinup.rc

# Check if mawk exist, otherwise create alias
if ! command -v mawk &> /dev/null
then
	alias mawk="awk"
fi


# Check if required directories exist, otherwise create them
if [ ! -d "./snow_init/" ]; then
	mkdir snow_init
fi
if [ ! -d "./cfgfiles/" ]; then
	mkdir cfgfiles
fi


for f in ./smet/*smet
do
	sed -i.bakr -e '/^easting.*/d' -e '/^northing.*/d' -e '/^epsg/d' ${f} && rm ${f}.bakr	# Delete the coordinates, since there seems to be a slight imprecision in the easting/northing.
	stn=$(basename ${f} .smet)
	snofile="./snow_init/${stn}.sno"
	cfgfile="./cfgfiles/${stn}.ini"

	awk -v ignore_slope=${ignore_slope} '{if(/station_id/) {stn_id=$NF} else if(/station_name/) {stn=$NF} else if(/latitude/) {lat=$NF} else if(/longitude/) {lon=$NF} else if(/altitude/) {alt=$NF} else if(/slope_angle/) {sl=$NF} else if(/slope_azi/) {azi=$NF} else if(/\[DATA\]/) {getline; getline; start=$1; exit}} END {
		print "SMET 1.1 ASCII";
		print "[HEADER]";
		print "station_id       = ", stn_id;
		print "station_name     = ", stn;
		print "latitude         = ", lat;
		print "longitude        = ", lon;
		print "altitude         = ", alt;
		print "nodata           = -999";
		print "ProfileDate      = ", start;
		print "HS_Last          = 0.00";
		print "SlopeAngle       = ", (1.-ignore_slope)*sl;
		print "SlopeAzi         = ", (1.-ignore_slope)*azi;
		print "nSoilLayerData   = 0";
		print "nSnowLayerData   = 0";
		print "SoilAlbedo       = 0.09";
		print "BareSoil_z0      = 0.020";
		print "CanopyHeight     = 0.00";
		print "CanopyLeafAreaIndex = 0.00";
		print "CanopyDirectThroughfall = 1.00";
		print "WindScalingFactor = 1.00";
		print "ErosionLevel     = 0";
		print "TimeCountDeltaHS = 0.000000";
		print "fields           = timestamp Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_WP  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo";
		print "[DATA]";
	}' ${f} > ${snofile}
	echo "[GENERAL]" > ${cfgfile}
	echo "IMPORT_BEFORE		=	../base.ini" >> ${cfgfile}
	echo "IMPORT_AFTER		=	../spinup.ini" >> ${cfgfile}
	echo "[INPUT]" >> ${cfgfile}
	echo "STATION1		=	${stn}.smet" >> ${cfgfile}
	echo "SNOWFILE1		=	${stn}.sno" >> ${cfgfile}
	echo "[OUTPUT]" >> ${cfgfile}
	echo "PROF_WRITE		=       TRUE" >> ${cfgfile}
	echo "TS_WRITE		=       TRUE" >> ${cfgfile}

	exp=$(fgrep EXPERIMENT base.ini | mawk '{print $NF}')
	if [ -z "${exp}" ]; then exp="NOEXP"; fi
	echo "EXPERIMENT		=	${exp}" >> ${cfgfile}
done


> to_exec.lst
i=0
for f in ./smet/*smet
do
	echo ${stn}
	stn=$(basename ${f} .smet)
	lat_lon=$(echo ${stn} | mawk -F_ '{print $2 "_" $3}')
	if [ -d "./profiles/" ]; then
		enddate=$(find ./profiles/* | fgrep -- ${lat_lon} | mawk -F\/ '{print substr($NF,1,4) "-" substr($NF,5,2) "-" ((substr($NF,7,2)==32)?(31):(substr($NF,7,2)))}')
	else
		enddate=""
	fi
	if [ -z "${enddate}" ]; then
		if [ ! -z "${final_end}" ]; then
			enddate=${final_end}
		else
			# Determine from last time stamp in *.smet file
			enddate=$(tail -1 ${f} | mawk '{print $1}')
		fi
	fi
	let i=${i}+1
	echo "bash spinup.sh \"${which_snowpack} -c cfgfiles/${stn}.ini -e ${enddate} > log/${stn}_${exp}.log 2>&1\"" >> to_exec.lst
done


# Modify sbatch script
# First, modify the SLURM part:
sed -i "s/--array=.*/--array=1-${i}/" job.sbatch
# Second, modify the PBS part:
sed -i "s/-J.*/-J 1-${i}/" job.sbatch
