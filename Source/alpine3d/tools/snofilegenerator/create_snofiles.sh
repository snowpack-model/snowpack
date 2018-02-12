#!/bin/bash
#
# This script creates sno-files, from PREVAH land use codes.
#
# It is constructed as follows:
# 1) Set default values
# 2) Read command line parameters to overwrite default values when provided
# 3) Reads file or stdin, to get list of land use codes where sno-files should be generated for.
# 4) Cycle through all the land use codes and create sno-files
#    a- Prepare values for output to sno-files
#    b- Write actual sno-file.
#
# Examples (start with no arguments to get more help):
# 1) bash create_snofiles.sh filename=lus.asc overwrite=Y Altitude=1256 BareSoil_z0=0.03 date=20070521 date_soillayers=20000101 time=1144 soillayermapfile=a.txt
# 2) echo 12400 | bash create_snofiles.sh overwrite=Y Altitude=2000 BareSoil_z0=0.03 date=20070521 date_soillayers=20000101
#
# Author: Nander Wever (wever@slf.ch)
#


#Default values for script:
filename="lus.asc"
list_with_codes=""
snofile="output_"
stationinfofile=""
canopymapfile="canopy_map.txt"
soillayermapfile="soillayer_map.txt"
overwrite="N"


#Default values for sno file:
date=19991001
time=0000
date_soillayers=${date}
time_soillayers=${time}
HS_Last=0.0
Latitude=36.0
Longitude=16.0
Altitude=1500
SlopeAngle=0
SlopeAzi=0
nSoilLayerData=19
nSnowLayerData=0
SoilAlbedo=0.2
BareSoil_z0=0.02
CanopyHeight=0.2
CanopyLeafAreaIndex=0.2
CanopyDirectThroughfall=0
WindScalingFactor=1.00
ErosionLevel=0
TimeCountDeltaHS=0.00

SurfaceHoarIndex=0.0
DriftIndex=0.0
ThreeHourNewSnow=0.0
TwentyFourHourNewSnow=0.0

Layer_Thick=0.5
T=278.15
Vol_Frac_I=0.00
Vol_Frac_W=0.02
Vol_Frac_V=0.01
Vol_Frac_S=0.97
Rho_S=2400.0
Conduc_S=0.3
HeatCapac_S=900.0
rg=10000
rb=0.0
dd=0.0
sp=0.0
mk=0
mass_hoar=0.0
ne=1
CDot=0.0
metamo=0.0


#Internal variables used by this script
flagformissingstationinfofile=0		#Flag for whether we warned that we cannot find the stationinfofile, so we only give warning once.
flagformissingcanopymapfile=0		#Flag for whether we warned that we cannot find the canopymapfile, so we only give warning once.
flagformissingsoillayermapfile=0	#Flag for whether we warned that we cannot find the soillayermapfile, so we only give warning once.


#Check for command line parameters
#if [ "$#" -eq "0" ]; then
#	#No command line parameters provided.
#	#Put code to deal with this situation here.
#fi


#Command line parameters overwrite default values:
for var in $(echo $* |  sed 's/\s/\n/g' |grep -e "="); do
 	eval $var
done


#Make backup of some variables, so we can overwrite them and restore them later on
defLatitude=${Latitude}
defLongitude=${Longitude}
defAltitude=${Altitude}
defSlopeAngle=${SlopeAngle}
defSlopeAzi=${SlopeAzi}

defSoilAlbedo=${SoilAlbedo}
defCanopyHeight=${CanopyHeight}
defCanopyLeafAreaIndex=${CanopyLeafAreaIndex}
defCanopyDirectThroughfall=${CanopyDirectThroughfall}

defnSoilLayerData=${nSoilLayerData}


#Create list of codes from $filename or stdin
if [ ! -z "${filename}" ]; then		#If filename is provided, then read file ...
	list_with_codes=`cat ${filename} | grep ^[0-9] | tr ' \t' '\n' | sort -nu | tr '\n' ' '`
	#                ^ put file in stream
	#				   ^ only select lines with numbers (header lines, etc, can be ignored this way)
	#						 ^ Make a vertical list
	#							       ^ sort by number
	#									 ^ take unique values (every code only once)
	#										^ make horizontal list
else					#... else read stdin.
	stdin=''			#Variable to store stdin
	if [ -t "0" ]; then		#Stdin present? This test whether stdin is terminal. If so, quit script (we don't have stdin AND we don't have filename...)
		#echo "No filename given and no stdin detected."
		echo "Use:"
		echo "bash create_snofiles.sh filename= [option=...]"
		echo "   where filename is GRID-file, or list with PREVAH-codes"
		echo "-or-"
		echo "bash create_snofiles.sh [option=...]"
		echo "   and send PREVAH-codes by stdin."
		echo "note: this works slower than providing filename."
		echo ""
		echo "Options:"
		echo "  filename=: give filename with land use codes."
		echo "  snofile=: prefix for sno-files (example: snofile=RUN_, creates snofiles like: RUN_{land_use_code}.sno)"
		echo "  stationinfofile=: file with station info, in format: [station_code Latitude Longitute Altitude SlopeAngle SlopeAzi]. Default: no station map."
		echo "  canopymapfile=: file with soil albedo map file, in format: [lus-code albedo CanopyHeight LeafAreaIndex CanopyDirectThroughfall]. Default: canopy_map.txt"
		echo "  soillayermapfile=: file with soil layers per land use class. [lus-code <soil layer info in format as in sno. files>. Default: soillayer_map.txt"
		echo "  overwrite=: Y: overwrite existing snofiles without asking, N: never overwrite snofiles (skips over land use code)."
		echo "  date=: profile date, in format YYYYMMDD."
		echo "  time=: profile time, in format HHMM."
		echo "  date_soillayers=: soil layer date, in format YYYYMMDD."
		echo "  time_soillayers=: soil layer time, in format HHMM."
		echo ""
		echo "Options to override defaults of sno file (defaults are given as example):"
		echo "  HS_Last=${HS_Last}"
		echo "  Latitude=${Latitude}"
		echo "  Longitude=${Longitude}"
		echo "  Altitude=${Altitude}"
		echo "  SlopeAngle=${SlopeAngle}"
		echo "  SlopeAzi=${SlopeAzi}"
		echo "  nSoilLayerData=${nSoilLayerData}"
		echo "  nSnowLayerData=${nSnowLayerData}"
		echo "  BareSoil_z0=${BareSoil_z0}"
		echo "  CanopyHeight=${CanopyHeight}"
		echo "  CanopyLeafAreaIndex=${CanopyLeafAreaIndex}"
		echo "  CanopyDirectThroughfall=${CanopyDirectThroughfall}"
		echo "  WindScalingFactor=${WindScalingFactor}"
		echo "  ErosionLevel=${ErosionLevel}"
		echo "  TimeCountDeltaHS=${TimeCountDeltaHS}"
		echo "  SurfaceHoarIndex=${SurfaceHoarIndex} (repeated 48 times)"
		echo "  DriftIndex=${DriftIndex} (repeated 48 times)"
		echo "  ThreeHourNewSnow=${ThreeHourNewSnow} (repeated 144 times)"
		echo "  TwentyFourHourNewSnow=${TwentyFourHourNewSnow} (repeated 144 times)"
		echo ""
		echo "  T=${T}"
		echo "  Vol_Frac_I=${Vol_Frac_I}"
		echo "  Vol_Frac_W=${Vol_Frac_W}"
		echo "  Vol_Frac_V=${Vol_Frac_V}"
		echo "  Vol_Frac_S=${Vol_Frac_S}"
		echo "  Rho_S=${Rho_S}"
		echo "  Conduc_S=${Conduc_S}"
		echo "  HeatCapac_S=${HeatCapac_S}"
		echo "  rg=${rg}"
		echo "  rb=${rb}"
		echo "  dd=${dd}"
		echo "  sp=${sp}"
		echo "  mk=${mk}"
		echo "  mass_hoar=${mass_hoar}"
		echo "  ne=${ne}"
		echo "  CDot=${CDot}"
		echo "  metamo=${metamo}"
		echo ""
		echo "Usage example:"
		echo "  bash create_snofiles.sh filename=lus.asc overwrite=Y Altitude=1600 BareSoil_z0=0.03 date=20070521 time=1144"
		echo "-or-"
		echo "  bash create_snofiles.sh filename=lus.asc overwrite=Y Altitude=0 BareSoil_z0=0.03 date=20070521 date_soillayers=20000101 time=1144 soillayermapfile=a.txt"
		echo "-or-"
		echo "  echo 12400 | bash create_snofiles.sh overwrite=Y Altitude=0 BareSoil_z0=0.03 date=20070521 date_soillayers=20000101"
		echo ""
		echo "Note 1: the script priority is always to use the information from stationinfofile, canopymapfile"
		echo "and soillayermapfile. When information is not available there, it will write out default values. By"
		echo "using for example stationinfofile=\"\", the script will always use default values for station info."
		echo "Note 2: the script priority is always to look for an input file with land use codes, before reading"
		echo "stdin."

		exit
	else	#stdin is present, so read it.
		while read line
		do
			stdin=`echo ${stdin} ${line} | grep ^[0-9] | tr ' ' '\n' | sort -n | uniq | tr '\n' ' '`
		done
		list_with_codes=`echo ${stdin} | grep ^[0-9] | tr ' ' '\n' | sort -n | uniq | tr '\n' ' '`
		#                ^ put file in stream
		#				 ^ only select lines with numbers (header lines, etc, can be ignored this way)
		#					       ^ Make a vertical list
		#							     ^ sort by number
		#								       ^ take unique values (every code only once)
		#									       ^ make horizontal list
	fi
fi


#Parse codes
if [ -z "${list_with_codes}" ]; then										#If list with codes is empty
	echo "No codes found..."
	echo "No sno files created, exiting."
else
	for code in `echo ${list_with_codes} | tr ' ' '\n' | awk '{printf "%d\n", $1}' | tr '\n' ' '`		#Cycle through all the codes in the list
	#            ^ put list in stream      ^ make list vertical        ^ reformat    ^ make list horizontal again
	do
		output_snofile=`echo ${snofile}${code}.sno`			#Create output snofile name.
		if [ -e "${output_snofile}" ] && [ ${overwrite} != "Y" ]; then	#If file exists, and we are not allowed to overwrite sno files, give error message and skip code.
			echo "Code [${code}]: ${output_snofile} already exists, doing nothing..."
		else
			#Check possibility to create sno-file:
			touch ${output_snofile}			#Create file, when it doesn't exist.
			if [ ! -w "${output_snofile}" ]; then	#Check if file is writable now.
				echo "Cannot write to file ${output_snofile}. Check permissions and disk space."
				exit
			fi
			#Read out code
			lus_prefix=`echo ${code} | awk '{print int($1/10000)}'`
			lus_code=`echo ${code} | awk '{print int(($1/100)%100)}'`
			soil_depth=`echo ${code} | awk '{print int(($1/10)%10)}'`
			soil_capacity=`echo ${code} | awk '{print int($1%10)}'`
			if [ ${lus_prefix} != "1" ]; then			#If prefix is not "1", it is not a PREVAH code, so we better skip it...
				echo "Code [${code}] is not in correct format. Skipping..."
			else
				# # # # # # # # # # # # # # # # # # # # # #
				# Prepare data to go into sno file:       #
				# # # # # # # # # # # # # # # # # # # # # #
				#Prepare datum
				YYYY=`echo ${date} | awk '{printf "%04d\n", int($1/10000)}'`
				MM=`echo ${date} | awk '{printf "%02d\n", int(($1%10000)/100)}'`
				DD=`echo ${date} | awk '{printf "%02d\n", int($1%100)}'`
				HH=`echo ${time} | awk '{printf "%02d\n", int($1/100)}'`
				mm=`echo ${time} | awk '{printf "%02d\n", int($1%100)}'`

				#Copy back defaults: 	(This is because for one point, defaults might be overwritten with ones from outside (mapping files). Therefore, in the next loop we have to reset these values).
				Latitude=${defLatitude}
				Longitude=${defLongitude}
				Altitude=${defAltitude}
				SlopeAngle=${defSlopeAngle}
				SlopeAzi=${defSlopeAzi}

				SoilAlbedo=${defSoilAlbedo}
				CanopyHeight=${defCanopyHeight}
				CanopyLeafAreaIndex=${defCanopyLeafAreaIndex}
				CanopyDirectThroughfall=${defCanopyDirectThroughfall}

				#Query stationinfofile for station details
				if [ -e "${stationinfofile}" ]; then		#stationinfofile exists?
					stationinfo=`cat ${stationinfofile} | sed 's/^[ \t]*//' | grep -v ^# | grep ^${lus_code}`		#Get line from stationinfofile
					if [ `echo ${stationinfo} | awk 'END {print NF}'` == "6" ]; then				#Are there six column in it?
						#Then override defaults
						Latitude=`echo ${stationinfo} | awk '(NR==1) {print $2}'`
						Longitude=`echo ${stationinfo} | awk '(NR==1) {print $3}'`
						Altitude=`echo ${stationinfo} | awk '(NR==1) {print $4}'`
						SlopeAngle=`echo ${stationinfo} | awk '(NR==1) {print $5}'`
						SlopeAzi=`echo ${stationinfo} | awk '(NR==1) {print $6}'`
					else												#Give warning
						echo "Station ${lus_code} not in station info file ${stationinfofile}, or station info is incomplete. Using defaults."
					fi
				else						#stationinfofile does not exist
					if [ "${flagformissingstationinfofile}" == "0" ] && [ ! -z "${stationinfofile}" ]; then		# Only do this when we have not given a warning yet AND there has been a stationinfofile specified
						echo "No station info file (${stationinfofile}) found, using default station values."
						flagformissingstationinfofile=1			# Set flag to 1, so we only warn once.
					fi
				fi


				#Determine canopy values
				if [ -e "${canopymapfile}" ]; then
					canopyvalues=`cat ${canopymapfile} | sed 's/^[ \t]*//' | grep -v ^# | sort -nk1 | grep ^${lus_code} | awk '{print $1, $2, $3, $4, $5}' | head -1`
					if [ `echo ${canopyvalues} | awk 'END {print NF}'` == "5" ]; then				#Are there enough column in it?
						SoilAlbedo=`echo ${canopyvalues} | awk '{print $2}'`
						CanopyHeight=`echo ${canopyvalues} | awk '{print $3}'`
						CanopyLeafAreaIndex=`echo ${canopyvalues} | awk '{print $4}'`
						CanopyDirectThroughfall=`echo ${canopyvalues} | awk '{print $5}'`
					else
						echo "No canopy values found for code: ${lus_code}. Using defaults."
					fi
				else
					if [ "${flagformissingcanopymapfile}" == "0" ]; then		# Only do this when we have not given a warning yet
						echo "No canopy map file (${canopymapfile}) found, using default station values."
						flagformissingcanopymapfile=1				# Set flag to 1, so we only warn once.
					fi
				fi


				#Determine soil layer values
				nSoilLayerData=${defnSoilLayerData}				#Reset to default value
				if [ -e "${soillayermapfile}" ]; then				#Check if soil layer map file exists
					nSoilLayersinFile=`cat ${soillayermapfile} | grep -v ^# | awk '($1=='${lus_code}') {print $0}' | wc -l`
					if [ "${nSoilLayersinFile}" -eq 0 ]; then
						nSoilLayersinFile=""				#Make string empty, so it can be used in the check below.
					fi
				else
					if [ "${flagformissingsoillayermapfile}" == "0" ]; then		# Only do this when we have not given a warning yet
						echo "No soil layer map file (${soillayermapfile}) found, using default station values."
						flagformissingsoillayermapfile=1			# Set flag to 1, so we only warn once.
					fi
					nSoilLayersinFile=""					#Make string empty, so it can be used in the check below.
				fi
				if [ -z "${nSoilLayersinFile}" ]; then				#If empty, we have no soil layer values in map file, or the soil layer map file does not exist.
					if [ "${flagformissingsoillayermapfile}" -ne 1 ]; then	#If this is not 1, the reason nSoilLayersinFile is empty is because the lus-code is not in the file, not because the file does not exist.
						echo "No soil layer values found for code: ${lus_code}. Using defaults."
					fi
					#Determine default values for soil layers
					for nSoilLayer in `seq 1 ${nSoilLayerData}`
					do
						#Note: the construction with eval creates the variables with index nSoilLayers (e.g. Vol_Frac_W14). This is an equal approach to having an array with values.
						#Make layer birth date and time equal to profile time
						eval YYYY_soil${nSoilLayer}=`echo ${date_soillayers} | awk '{printf "%04d\n", int($1/10000)}'`
						eval MM_soil${nSoilLayer}=`echo ${date_soillayers} | awk '{printf "%02d\n", int(($1%10000)/100)}'`
						eval DD_soil${nSoilLayer}=`echo ${date_soillayers} | awk '{printf "%02d\n", int($1%100)}'`
						eval HH_soil${nSoilLayer}=`echo ${time_soillayers} | awk '{printf "%02d\n", int($1/100)}'`
						eval mm_soil${nSoilLayer}=`echo ${time_soillayers} | awk '{printf "%02d\n", int($1%100)}'`
						#Set other soil properties
						eval Layer_Thick${nSoilLayer}=${Layer_Thick}
						eval T${nSoilLayer}=${T}
						eval Vol_Frac_I${nSoilLayer}=${Vol_Frac_I}
						eval Vol_Frac_W${nSoilLayer}=${Vol_Frac_W}
						eval Vol_Frac_V${nSoilLayer}=${Vol_Frac_V}
						eval Vol_Frac_S${nSoilLayer}=${Vol_Frac_S}
						eval Rho_S${nSoilLayer}=${Rho_S}
						eval Conduc_S${nSoilLayer}=${Conduc_S}
						eval HeatCapac_S${nSoilLayer}=${HeatCapac_S}
						eval rg${nSoilLayer}=${rg}
						eval rb${nSoilLayer}=${rb}
						eval dd${nSoilLayer}=${dd}
						eval sp${nSoilLayer}=${sp}
						eval mk${nSoilLayer}=${mk}
						eval mass_hoar${nSoilLayer}=${mass_hoar}
						eval ne${nSoilLayer}=${ne}
						eval CDot${nSoilLayer}=${CDot}
						eval metamo${nSoilLayer}=${metamo}
					done
				else
					nSoilLayerData=${nSoilLayersinFile}
				fi




				# # # # # # # # # # # # # # # # # # # # # #
				# Write out to sno file:                  #
				# # # # # # # # # # # # # # # # # # # # # #
				echo "[SNOWPACK_INITIALIZATION]" > ${output_snofile}
				echo "StationName=${lus_code}" >> ${output_snofile}
				echo "ProfileDate=${YYYY} ${MM} ${DD} ${HH} ${mm}" >> ${output_snofile}
				echo "HS_Last=${HS_Last}" >> ${output_snofile}
				echo "Latitude=${Latitude}" >> ${output_snofile}
				echo "Longitude=${Longitude}" >> ${output_snofile}
				echo "Altitude=${Altitude}" >> ${output_snofile}
				echo "SlopeAngle=${SlopeAngle}" >> ${output_snofile}
				echo "SlopeAzi=${SlopeAzi}" >> ${output_snofile}
				echo "nSoilLayerData=${nSoilLayerData}" >> ${output_snofile}
				echo "nSnowLayerData=${nSnowLayerData}" >> ${output_snofile}
				echo "SoilAlbedo=${SoilAlbedo}" >> ${output_snofile}
				echo "BareSoil_z0=${BareSoil_z0}" >> ${output_snofile}
				echo "CanopyHeight=${CanopyHeight}" >> ${output_snofile}
				echo "CanopyLeafAreaIndex=${CanopyLeafAreaIndex}" >> ${output_snofile}
				echo "CanopyDirectThroughfall=${CanopyDirectThroughfall}" >> ${output_snofile}
				echo "WindScalingFactor=${WindScalingFactor}" >> ${output_snofile}
				echo "ErosionLevel=${ErosionLevel}" >> ${output_snofile}
				echo "TimeCountDeltaHS=${TimeCountDeltaHS}" >> ${output_snofile}
				echo "YYYY MM DD HH MI Layer_Thick  T  Vol_Frac_I  Vol_Frac_W  Vol_Frac_V  Vol_Frac_S Rho_S Conduc_S HeatCapac_S  rg  rb  dd  sp  mk mass_hoar ne CDot metamo" >> ${output_snofile}
				#Use this, in combination with: cat *.sno | grep hhhhhh | sed 's/hhhhhhh //' > soillayer_map.txt to generate standard soillayer map.
				#for nSoilLayer in `seq 1 ${nSoilLayerData}`
				#do
					#eval "echo hhhhhhh \${lus_code} \$YYYY_soil${nSoilLayer} \$MM_soil${nSoilLayer} \$DD_soil${nSoilLayer} \$HH_soil${nSoilLayer} \$mm_soil${nSoilLayer} \$Layer_Thick${nSoilLayer} \$T${nSoilLayer} \$Vol_Frac_I${nSoilLayer} \$Vol_Frac_W${nSoilLayer} \$Vol_Frac_V${nSoilLayer} \$Vol_Frac_S${nSoilLayer} \$Rho_S${nSoilLayer} \$Conduc_S${nSoilLayer} \$HeatCapac_S${nSoilLayer} \$rg${nSoilLayer} \$rb${nSoilLayer} \$dd${nSoilLayer} \$sp${nSoilLayer} \$mk${nSoilLayer} \$mass_hoar${nSoilLayer} \$ne${nSoilLayer} \$CDot${nSoilLayer} \$metamo${nSoilLayer}" >> ${output_snofile}
				#done


				#Write out soil layers
				if [ ! -z "${nSoilLayersinFile}" ]; then			#When we can read soil layers from file
					#Check that total volumetric content == 1:
					cat ${soillayermapfile} | grep -v ^# | awk '($1=='${lus_code}') {if ($9+$10+$11+$12!=1.0) {printf "Code ['${code}']: Warning: total volumetric content not 1 at line %d in file %s.\n", NR, "'${soillayermapfile}'"}}'
					#Then write out soil layer info:
					cat ${soillayermapfile} | grep -v ^# | awk '($1=='${lus_code}') {print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' >> ${output_snofile}
				else								#If not, write out default values
					for nSoilLayer in `seq 1 ${nSoilLayerData}`		#Write out default values
					do
						#Check if total volumetric content is 1:
						if [ $(eval "echo \$Vol_Frac_I${nSoilLayer} \$Vol_Frac_W${nSoilLayer} \$Vol_Frac_V${nSoilLayer} \$Vol_Frac_S${nSoilLayer}" | awk '{print $1+$2+$3+$4==1}') != "1" ]; then
							echo "Code [${code}]: Warning: total volumetric content not 1 for nSoilLayer==${nSoilLayer}."
						fi
						#Write out soil layer info:
						eval "echo \$YYYY_soil${nSoilLayer} \$MM_soil${nSoilLayer} \$DD_soil${nSoilLayer} \$HH_soil${nSoilLayer} \$mm_soil${nSoilLayer} \$Layer_Thick${nSoilLayer} \$T${nSoilLayer} \$Vol_Frac_I${nSoilLayer} \$Vol_Frac_W${nSoilLayer} \$Vol_Frac_V${nSoilLayer} \$Vol_Frac_S${nSoilLayer} \$Rho_S${nSoilLayer} \$Conduc_S${nSoilLayer} \$HeatCapac_S${nSoilLayer} \$rg${nSoilLayer} \$rb${nSoilLayer} \$dd${nSoilLayer} \$sp${nSoilLayer} \$mk${nSoilLayer} \$mass_hoar${nSoilLayer} \$ne${nSoilLayer} \$CDot${nSoilLayer} \$metamo${nSoilLayer}" >> ${output_snofile}
					done
				fi

				echo "SurfaceHoarIndex" >> ${output_snofile}
				echo ${SurfaceHoarIndex} | awk '{ for(i=1;i<=48;i++) printf("%s ", $1) } {printf("\n")}' >> ${output_snofile}
				echo "DriftIndex" >> ${output_snofile}
				echo ${DriftIndex} | awk '{ for(i=1;i<=48;i++) printf("%s ", $1) } {printf("\n")}' >> ${output_snofile}
				echo "ThreeHourNewSnow" >> ${output_snofile}
				echo ${ThreeHourNewSnow} | awk '{ for(i=1;i<=144;i++) printf("%s ", $1) } {printf("\n")}' >> ${output_snofile}
				echo "TwentyFourHourNewSnow" >> ${output_snofile}
				echo ${TwentyFourHourNewSnow} | awk '{ for(i=1;i<=144;i++) printf("%s ", $1) } {printf("\n")}' >> ${output_snofile}
				echo "End" >> ${output_snofile}

				#Inform user
				echo "Code [${code}]: sno-file created"
				#echo "  luscode = ${lus_code}"
				#echo "  soil_depth = ${soil_depth}"
				#echo "  soil_capacity = ${soil_capacity}"
			fi
		fi
	done
fi
