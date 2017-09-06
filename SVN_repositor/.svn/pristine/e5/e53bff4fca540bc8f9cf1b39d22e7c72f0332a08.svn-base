#!/bin/bash
# This script converts downloaded sensorscope data into SMET files for use with SNOWPACK.
# Data should be downloaded from www.climaps.com, with the settings as shown in sensorscope_downloadsettings.png. Then unpack the zip file downloaded from www.climaps.com
# The script works as follows
#    1 - First it scans the given directory to check which stations are in (station_XXXX.csv files)
#    2 - It converts meteo data to SMET format
#    3 - It converts soil moisture data to a SMET like format.
#    4 - If files already exist in the working directory, the files are merged with the new data.
#        In case of inconsistencies, the new data is used.
#
# What to set in the settings:
#   - You should specify the sensor mapping. This tells which sensor from sensorscope should be matched to which SMET variable. The variable names are called according to SMET variables.
#   - You can specify settings for the SMET header, like station_name, latitude, longitude, etc. Note that such a filename is then used for all the stations that are in the input file.
#   - You can specify a nodata value by adding nodata=<my own value> to the command line parameters.
#     Note that this script only works correctly with numerical nodata values.
#
# Example usage: 
#   bash converttosmet.sh exampledata.csv.gz sensorscope station_name="sensorscope" nodata=-99999.9
#   bash converttosmet.sh exampledata.csv.gz "" station_name="" nodata=-999 stationinfo_file="sensorscope_stn.info"
# Please be careful, this script is not so robust. For example: if the nodata values is changed between files, merging will result in files with both nodata values mixed.
# Author: Nander Wever
#
#


#Basic settings
#Sensor mapping: this maps the SMET variables to the SensorScope variables
nametimestamp="Epoch"
nameTA="SHT75 / Temperature"
nameRH="SHT75 / Humidity"
nameVW="Davis / Wind Speed"
nameDW="Davis / Wind Direction"
nameOSWR="SP-212 / Radiation"
nameISWR="SP-212 / Radiation [2]"
namePSUM="Davis / Rain"
nameTSS="TNX / Surface Temperature"

#Files:
stationinfo_file=""

#Standard settings:
nodata=-999
station_name="unknown"
latitude=${nodata}
longitude=${nodata}
altitude=${nodata}
easting=${nodata}
northing=${nodata}
epsg=${nodata}
tz=1




#
#  E  N  D    O  F    S  E  T  T  I  N  G  S 
#
##############################################################################
# Below this line, edits are not really necessary, unless the program itself
# should be changed. All the possible settings are above.



#Read command line parameters
directory=$1
smetfileprefix=$2
#Command line parameters overwrite default values:
for var in $(echo $* |  sed 's/\s/\n/g' |grep -e "="); do 
        eval $var
done

#Check for valid command line parameters
if [ -z ${directory} ]; then
	echo "ERROR: no smet output file name specified..."
	echo "use: bash converttosmet <directory> <smetfile_prefix> <additional_option=value>"
	echo "     example: bash converttosmet.sh ~/Downloads/ sensorscope station_name=\"sensorscope\" nodata=-999"
	echo "              bash converttosmet.sh ~/Downloads/ \"\" station_name=\"\" nodata=-999 stationinfo_file=\"sensorscope_stn.info\""
	echo "     The <directory> should contain station_XXXX.csv files. Use ./ to use the current working directory. Always end the directory specification with a \"/\"!"
	echo "  additional_options are for SMET file. Defaults:"
	echo "    station_name=${stations_name}    (if left empty, station_name is taken equal to station_id)"
	echo "    nodata=${nodata}"
	echo "    latitude=${latitude}"
	echo "    longitude=${longitude}"
	echo "    altitude=${altitude}"
	echo "    easting=${easting}"
	echo "    northing=${northing}"
	echo "    epsg=${epsg}"
	echo "    tz=${tz}"
	echo "    stationinfo_file=${stationinfo_file}"
	echo "         The station info file should contain: station_id station_name swiss_x swiss_y latitude longitude altitude."
	exit
fi



#Inform user
echo "Analysing: $1..."

#Look for stations
list_of_stns=`find ${directory}station*csv -printf %f\\\\n 2> /dev/null | sed 's/station_//' | sed 's/.csv//'`
if [ -z "${list_of_stns}" ]; then
	echo "ERROR: no stations found in directory: ${directory}"
	echo "   Did you put a \"/\" at the end of the directory specification?"
	exit
fi
echo "-- stations found: `echo ${list_of_stns} | tr '\n' ' '`"

#Create SMET file per stn
for stn in ${list_of_stns}
do
  echo "Parsing $stn..."


  #Determine input data file name
  datafilename=`echo ${directory}station_${stn}.csv`


  #Identify in which column the time stamps are:
  colnr_for_EPOCHtime=`cat ${datafilename} | head -1 | tr ',' '\n' | grep -n ^Epoch | awk -F\: '{print $1}'`
  if [ -z ${colnr_for_EPOCHtime} ]; then
	echo "ERROR: No field <Time> found in file..."
	exit
  fi


  #Identify in which column the GMT time stamps are (this we need to correctly set TZ in SMET file):
  #NOTE: in the new format, this is not possible anymore. So, we assume: tz="+1", without DST.
  #colnr_for_GMTtime=`cat ${datafilename} | head -1 | tr ',' '\n' | grep -n GMT\ Time | awk -F\: '{print $1}'`
  #if [ -z ${colnr_for_GMTtime} ]; then
  #	echo "ERROR: No field <GMT Time> found in file..."
  #	exit
  #fi


  #Get header, to extract variables in file
  header=`cat ${datafilename} | head -1 | tr ' ' '_' | tr ',' ' '`


  #Determine SMET file name
  if [ ! -z "${smetfileprefix}" ]; then
	smetfilename=`echo ${smetfileprefix}_${stn}.smet`
  	soilmoisturefilename=`echo ${smetfileprefix}_${stn}.vwc`
  else
	smetfilename=`echo ${stn}.smet`
  	soilmoisturefilename=`echo ${stn}.vwc`
  fi

  if [ -e "${smetfilename}" ]; then
	mergedsmetfilename="${smetfilename}"
	smetfilename="${smetfilename}.tmp"
	mergesmet=1
  else
	mergesmet=0
  fi

  #Determine VWC file name
  if [ -e "${soilmoisturefilename}" ]; then
	mergedsoilmoisturefilename="${soilmoisturefilename}"
	soilmoisturefilename="${soilmoisturefilename}.tmp"
	mergesoilmoisture=1
  else
	mergesoilmoisture=0
  fi


  #Determine TZ
  #tz=`cat ${datafilename} | awk -FGMT '{print $2}' | awk '{print $1}' | sort -u | head -1`
  tz="+1"
  export TZ="UTC"					#Set the environment variable to UTC time. This is for a correct working of awk's strftime.
  tz_shift=`echo ${tz} | awk '{print $1*60*60}'`	#This shifts the timestamp by an interval tz (given in hours, so +1), but then converted to seconds.


  #Construction of SMET-file structure
  fields="fields       = "		#Output fields line
  soilmoisturefields="fields       = "	#Output fields line
  i=0					#Counts columns in input
  columns=0				#Total number of columns in output
  nsoilsensors=0			#Total number of soil sensors

  #Here, each sensor found in ${header} is compared to the names given in the beginning of the script. If a matching SMET variable name is found, the following is done:
	# 1. The number of output columns in the SMET file is increased.
	# 2. The fields line for output to the SMET file is updated with the found variable.
	# 3. Then it is stored in which column in the SMET file, the variable is. This variable is currently not used.
	# 4. Then the mapping of the columns is done, so which column from the input should be used for which column in the output.
	# 5. Then an operation to the variable can be determined, like conversion from K to C, or something like that.
	#    Examples:
	#    Do nothing: todo_with_var[columns]=""
	#    Conversion C to K: todo_with_var[columns]="+273.15"
	#    Conversion RH: todo_with_var[columns]="/100.0"
	#    The idea is that when constructing the transfer command, the awk statement is made to print variable columns, and add the todo_with var: "echo variable | awk '{print $1${todo_with_var}}'"
  

  #First pick out the time stamp, as that one we would like to have in the first column in the SMET file.
  for sensor in ${header}				#Cycle through all header elements.
  do
        sensor=`echo ${sensor} | sed 's/_/ /g'`		#Note that spaces were replaced by _, so now convert them back to match against original sensor names.
	let i=${i}+1					#Increase number of input columns

	if [ "${sensor}" == "${nametimestamp}" ]; then	#See description above.
		let columns=${columns}+1
		fields=`echo ${fields} timestamp`
		output_colnr_for_timestamp=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]=""
		#Determine file resolution (this is used to correct precipitation amounts. It determines resolution by picking out time stamps, doubling them (except first row), putting them next together and taking the differences. Then look for the difference that occurs most often (which will be the measurement resolution).
		inputfile_resolution=`cat ${datafilename} | grep ^[0-9] | awk -F, '{print $'${input_output_colnr[columns]}'}' | awk '{if (NR==1) {printf "%s\n", $1} else {printf "%s\n%s\n", $1, $1}}' | sed '$!N;s/\n/ /' | awk '{print $2-$1}' | sort -n | uniq -c | sort -nrk1 | awk '(NR==1) {print $2}'`
		echo Resolution: ${inputfile_resolution} s.
	fi
  done
  #Then the sensors
  i=0
  for sensor in ${header}				#Cycle through all header elements.
  do
        sensor=`echo ${sensor} | sed 's/_/ /g'`		#Note that spaces were replaced by _, so now convert them back to match against original sensor names.
	let i=${i}+1					#Increase number of input columns

	#Time stamp is dealt with separately (see above)
	#if [ "${sensor}" == "${nametimestamp}" ]; then	#See description above.
		#let columns=${columns}+1
		#fields=`echo ${fields} timestamp`
		#output_colnr_for_timestamp=${columns}
		#input_output_colnr[columns]=${i}
		#todo_with_var[columns]=""
	#fi

	if [ "${sensor}" == "${nameTA}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} TA`
		output_colnr_for_TA=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]="+273.15"
	fi

	if [ "${sensor}" == "${nameRH}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} RH`
		output_colnr_for_RH=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]="/100.0"
	fi

	if [ "${sensor}" == "${nameVW}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} VW`
		output_colnr_for_VW=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]=""
	fi

	if [ "${sensor}" == "${nameDW}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} DW`
		output_colnr_for_DW=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]=""
	fi

	if [ "${sensor}" == "${nameOSWR}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} OSWR`
		output_colnr_for_OWSR=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]=""
	fi

	if [ "${sensor}" == "${nameISWR}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} ISWR`
		output_colnr_for_ISWR=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]=""
	fi

	if [ "${sensor}" == "${namePSUM}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} PSUM`
		output_colnr_for_PSUM=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]="*(3600.0/${inputfile_resolution}.0)"		#This translate [mm/time step] into [mm/hour]
	fi

	if [ "${sensor}" == "${nameTSS}" ]; then
		let columns=${columns}+1
		fields=`echo ${fields} TSS`
		output_colnr_for_TSS=${columns}
		input_output_colnr[columns]=${i}
		todo_with_var[columns]="+273.15"
	fi

	if [ "`echo ${sensor} | awk '{print $3}'`" == "VWC" ]; then
		let nsoilsensors=${nsoilsensors}+1
		soilsensor[nsoilsensors]=`echo ${sensor} | sed 's/ /_/g'`
		soilmoisturefields="${soilmoisturefields} ${soilsensor[nsoilsensors]}"
		input_output_colnr_soilmoisture[nsoilsensors]=${i}
		todo_with_var_soilmoisture[nsoilsensors]=""
	fi
  done

  #Check if stationinfo-file is present
  if [ -e "${stationinfo_file}" ]; then
	#Read station info
	station_name=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $2}'`
	easting=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $3}'`
	northing=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $4}'`
	latitude=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $5}'`
	longitude=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $6}'`
	altitude=`cat ${stationinfo_file} | sed 's/^ *//' | grep ^${stn} | awk '{print $7}'`
	if [ -z "${easting}" ]; then easting=${nodata}; fi
	if [ -z "${northing}" ]; then northing=${nodata}; fi
	if [ -z "${latitude}" ]; then latitude=${nodata}; fi
	if [ -z "${longitude}" ]; then longitude=${nodata}; fi
	if [ -z "${altitude}" ]; then altitude=${nodata}; fi
  fi

  #Create header of SMET
  echo "SMET 1.1 ASCII" > ${smetfilename}
  echo "[HEADER]" >> ${smetfilename}
  echo "station_id   = ${stn}" >> ${smetfilename}
  if [ ! -z "${station_name}" ]; then
	echo "station_name = ${station_name}" >> ${smetfilename}
  else
	echo "station_name = ${stn}" >> ${smetfilename}
  fi
  echo "latitude     = ${latitude}" >> ${smetfilename}
  echo "longitude    = ${longitude}" >> ${smetfilename}
  echo "altitude     = ${altitude}" >> ${smetfilename}
  echo "easting      = ${easting}" >> ${smetfilename}
  echo "northing     = ${northing}" >> ${smetfilename}
  echo "epsg         = ${epsg}" >> ${smetfilename}
  echo "nodata       = ${nodata}" >> ${smetfilename}
  echo "tz           = ${tz}" >> ${smetfilename}
  echo "${fields}" >> ${smetfilename}
  echo "[DATA]" >> ${smetfilename}


  #Create header of SMET
  echo "SMET 1.1 ASCII" > ${soilmoisturefilename}
  echo "[HEADER]" >> ${soilmoisturefilename}
  echo "station_id   = ${stn}" >> ${soilmoisturefilename}
  if [ ! -z "${station_name}" ]; then
	echo "station_name = ${station_name}" >> ${soilmoisturefilename}
  else
	echo "station_name = ${stn}" >> ${soilmoisturefilename}
  fi
  echo "latitude     = ${latitude}" >> ${soilmoisturefilename}
  echo "longitude    = ${longitude}" >> ${soilmoisturefilename}
  echo "altitude     = ${altitude}" >> ${soilmoisturefilename}
  echo "easting      = ${easting}" >> ${soilmoisturefilename}
  echo "northing     = ${northing}" >> ${soilmoisturefilename}
  echo "epsg         = ${epsg}" >> ${soilmoisturefilename}
  echo "nodata       = ${nodata}" >> ${soilmoisturefilename}
  echo "tz           = ${tz}" >> ${soilmoisturefilename}
  echo "${soilmoisturefields}" >> ${soilmoisturefilename}
  echo "[DATA]" >> ${soilmoisturefilename}


  #Now read meteo data from file
  #First create executecommand. We will construct an executecommand, which will be evaluated and does everything necessary.
  executecommand="cat ${datafilename} | grep ^[0-9] | sed 's/nan/${nodata}/g' | awk -F, '{print strftime(\"%Y %m %d %H %M\", \$${input_output_colnr[1]}+${tz_shift})"
  #               ^^ open file		              ^^ change nan to nodata      ^^ brake up date, into YYYY MM DD HH mm. We force this to be the first field above, so we can use index [1] here.
  #Now we already adressed the 1st column in the output file (the time stamp), now cycle through all remaining columns and add to the 
  for i in `seq 2 ${columns}`
  do
	executecommand="${executecommand}, (\$${input_output_colnr[i]}==${nodata})?${nodata}:\$${input_output_colnr[i]}${todo_with_var[i]}"
	#       ^^^ add to executecommand         ^^^ to don't mess up todo_with_var with nodata, make separation whether column is no data or not.
  done
  executecommand="${executecommand}}' | sed -e 's/ /-/' -e 's/ /-/' -e 's/ /T/' -e 's/ /\:/' | awk '("
  #   ^^^ add to executecommand         ^^^ This set-statement translates the time stamp to SMET format     ^^^ This opens an awk statement to prevent lines with only nodata values to appear in the SMET-file.

  for i in `seq 2 ${columns}`
  do
	#This part constructs the awk statement like: awk '($2!=nodata || $3!=nodata || $4 != nodata)'. This only needs to be done for columns 2 and higher, as the timestamp will never be nodata.
	if (( i==2 )); then
		executecommand="${executecommand}\$${i}!=${nodata}"
	else
		executecommand="${executecommand} || \$${i}!=${nodata}"
	fi
  done
  executecommand="${executecommand})' | awk '{print 1, \$'0'}'"
  #                 ^^^ finish awk-statement     ^^^ add a 1 as first column, to identify orginal data.
  eval ${executecommand} >> ${smetfilename}.tmp1
  # ^^^ we now constructed an command which does the translation. With eval it is executed. use echo ${executecommand} to view what it is actually doing.


  #Because sometimes, one sensor give the data later, but within the same minute merge this data to one time stamp:
  cat ${smetfilename}.tmp1 | cut -d\  -f2- | uniq -w16 -dD | tr '\n' ' ' | awk '{ for(j=1; j<=NF; j+='${columns}') if($j==$(j+'${columns}')) {print $j; for(i=j+1;i<j+'${columns}'; i++){print (($i==-999 && $(i+'${columns}')==-999) || ($i!=-999 && $(i+'${columns}') !=-999))?($i+$(i+'${columns}'))/2.0:($i+$(i+'${columns}')+999)}}}' | tr '\n' ' ' | awk '{for (i=1; i<=NF; i++) {printf "%s ", $i; if (i%'${columns}'==0) printf "\n"}}' | awk '{print 0, $0}' | sed 's/ $//' > ${smetfilename}.tmp2
  # ^^^ start pipe           ^^^ remove first column ^^^ select only timestamps which occur multiple times
  #                                                           ^^^ put everything in a row   ^^^ cycle through all blocks of data (per timestamp)
  #															^^^ if two succeeding blocks match ...
  #																		^^^ print time stamp and  ...                  ^^^^  when not only one of them is nodata                                          ^^^ print average (is nodata when both are nodata, or mean value if both are valid values.
 #																																					   ^^^ or print a single value. This construction is because we don't know which value (first or second) is valid and which one not.
  #																																									    ^^^ make a single row		^^^^ make one timestamp per row, making use of the number of columns in the file
  #																																																						  ^^^ print 0 in first column, to identify modified data.   ^^^ The final sed removes a white space at the end of the line. It is there, because of the print-loops in awk.
  cat ${smetfilename}.tmp1 ${smetfilename}.tmp2 | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 >> ${smetfilename}
  #Put now both files in the pipe		  ^^^ sort first on date, and then on first column, which contains the identifier for being original or modified data
  #							           ^^^ remove identifier in first column   and the uniq now selects single timestamps. When multiple timestamps are found, only the first occurrence is written out, which is modified data when available (because of the sorting on the first column)

  #Remove temporary files
  rm ${smetfilename}.tmp1 ${smetfilename}.tmp2




  #Now read soil moisture data from file
  if (( ${nsoilsensors}>0 )); then
	#First create executecommand. We will construct an executecommand, which will be evaluated and does everything necessary.
  	executecommand="cat ${datafilename} | grep ^[0-9] | sed 's/nan/${nodata}/g' | awk -F, '{print strftime(\"%Y %m %d %H %M\", \$${input_output_colnr[1]}+${tz_shift})"
  	#               ^^ open file		              ^^ change nan to nodata      ^^ brake up date, into YYYY MM DD HH mm. We force this to be the first field above, so we can use index [1] here.

  	#Now we already adressed the 1st column in the output file (the time stamp), now cycle through all remaining columns and add to the 
  	for i in `seq 1 ${nsoilsensors}`
  	do
		executecommand="${executecommand}, (\$${input_output_colnr_soilmoisture[i]}==${nodata})?${nodata}:\$${input_output_colnr_soilmoisture[i]}${todo_with_var_soilmoisture[i]}"
		#       ^^^ add to executecommand         ^^^ to don't mess up todo_with_var with nodata (like conversion to kelvin of nodata: -999+273.15), make separation whether column is no data or not.
  	done
  	executecommand="${executecommand}}' | sed -e 's/ /-/' -e 's/ /-/' -e 's/ /T/' -e 's/ /\:/' | awk '("
  	#   ^^^ add to executecommand         ^^^ This set-statement translates the time stamp to SMET format     ^^^ This opens an awk statement to prevent lines with only nodata values to appear in the SMET-file.

  	for i in `seq 1 ${nsoilsensors}`
  	do
		#This part constructs the awk statement like: awk '($2!=nodata || $3!=nodata || $4 != nodata)'. This only needs to be done for columns 2 and higher, as the timestamp will never be nodata.
		if (( i==1 )); then
			executecommand="${executecommand}\$(${i}+1)!=${nodata}"
		else
			executecommand="${executecommand} || \$(${i}+1)!=${nodata}"
		fi
  	done
  	executecommand="${executecommand})' | awk '{print 1, \$'0'}'"
  	#                 ^^^ finish awk-statement     ^^^ add a 1 as first column, to identify orginal data.
  	eval ${executecommand} >> ${soilmoisturefilename}.tmp1
  	# ^^^ we now constructed an command which does the translation. With eval it is executed. use echo ${executecommand} to view what it is actually doing.
  

  	#Because sometimes, one sensor give the data later, but within the same minute merge this data to one time stamp:
  	#  note: in this part, nsoilsensors+1 is used to determine the number of columns. This is because $columns for meteo data contains time stamp, but $nsoilsensors not.
  	cat ${soilmoisturefilename}.tmp1 | cut -d\  -f2- | uniq -w16 -dD | tr '\n' ' ' | awk '{ for(j=1; j<=NF; j+='${nsoilsensors}'+1) if($j==$(j+'${nsoilsensors}'+1)) {print $j; for(i=j+1;i<j+'${nsoilsensors}'+1; i++){print (($i==-999 && $(i+'${nsoilsensors}'+1)==-999) || ($i!=-999 && $(i+'${nsoilsensors}'+1) !=-999))?($i+$(i+'${nsoilsensors}'+1))/2.0:($i+$(i+'${nsoilsensors}'+1)+999)}}}' | tr '\n' ' ' | awk '{for (i=1; i<=NF; i++) {printf "%s ", $i; if (i%('${nsoilsensors}'+1)==0) printf "\n"}}' | awk '{print 0, $0}' | sed 's/ $//' > ${soilmoisturefilename}.tmp2
  	# ^^^ start pipe           ^^^ remove first column ^^^ select only timestamps which occur multiple times
  	#                                                           ^^^ put everything in a row   ^^^ cycle through all blocks of data (per timestamp)
  	#															^^^ if two succeeding blocks match ...
  	#																		^^^ print time stamp and  ...                  ^^^^  when not only one of them is nodata                                          ^^^ print average (is nodata when both are nodata, or mean value if both are valid values.
	#																																					   ^^^ or print a single value. This construction is because we don't know which value (first or second) is valid and which one not.
  	#																																									    ^^^ make a single row		^^^^ make one timestamp per row, making use of the number of columns in the file
  	#																																																						  ^^^ print 0 in first column, to identify modified data.   ^^^ The final sed removes a white space at the end of the line. It is there, because of the print-loops in awk.
	cat ${soilmoisturefilename}.tmp1 ${soilmoisturefilename}.tmp2 | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 >> ${soilmoisturefilename}
  	#Put now both files in the pipe		  ^^^ sort first on date, and then on first column, which contains the identifier for being original or modified data
  	#							           ^^^ remove identifier in first column   and the uniq now selects single timestamps. When multiple timestamps are found, only the first occurrence is written out, which is modified data when available (because of the sorting on the first column). The output is added to the file, which already contains the header.

	#Remove temporary files
	rm ${soilmoisturefilename}.tmp1 ${soilmoisturefilename}.tmp2
  fi





  #Merge files if necessary
  #Note: merging is done based on time stamp. Fields should be equal between files. When newer data differes for the same time stamp, the new data is used.
  if (( ${mergesmet} == 1 )); then
	#Inform user
	echo "Merging smet files..."

	#First check if structure is the same
  	if (( `cat ${smetfilename} ${mergedsmetfilename} | sort | grep ^fields | uniq | wc -l` != 1 )); then
		echo "ERROR: cannot merge, because fields are not the same!"
	else

		#Then write out header (take the newest one)
		cat ${smetfilename} | grep -v ^[0-9] >> ${smetfilename}.tmp1

		#Then write out data from both files, including a time stamp, and a 0 to denote the new data and a 1 to denote old data (this is used to select the newer data in case there is different data for the same time stamp). The sed transforms YYYY-MM-DDTHH:mm to YYYYMMDDHHmm, which makes direct comparison/sorting possible. 
		cat ${smetfilename} | grep ^[0-9] | awk '{print $1, 0, $0}' | sed 's/\([0-9][0-9][0-9][0-9]\)-\([0-9][0-9]\)-\([0-9][0-9]\)T\([0-9][0-9]\)\:\([0-9][0-9]\)/\1\2\3\4\5/1' > ${smetfilename}.tmp2
		cat ${mergedsmetfilename} | grep ^[0-9] | awk '{print $1, 1, $0}' | sed 's/\([0-9][0-9][0-9][0-9]\)-\([0-9][0-9]\)-\([0-9][0-9]\)T\([0-9][0-9]\)\:\([0-9][0-9]\)/\1\2\3\4\5/1' > ${smetfilename}.tmp3
  
		#Then merge files and sort:
		cat ${smetfilename}.tmp2 ${smetfilename}.tmp3 | sort -n -k 1 -k 2 > ${smetfilename}.tmp4
		#								^^^ this first select time stamps, and in case of multiple occurrences of equal timestamps, select the one from the new dataset, denoted by a 0 in column 2.

		#Now check for equal time stamps, with different data. If so, give warning.
		#The first uniq removes all really unique lines, the second one checks if still some time stamps appear more than once, which means the data is different.
		#Note, the awk is to get rid of the first 2 columns, because the second column is always different between the two files, as it indicates from which file the line is coming.
		#That's why the uniq has to check 16 characters, as the first column is now the date in the YYYY-MM-DDTHH:mm format.
		if (( `cat ${smetfilename}.tmp4 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' | uniq | uniq -w16 -dD | wc -l` > 0 )); then
			echo "WARNING: duplicate time stamps with different data encountered! Using data from the new file."
			echo "   Inconsistencies are saved in: ${smetfilename}.inconsistencies."
			cat ${smetfilename}.tmp4 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' | uniq | uniq -w16 -dD > ${smetfilename}.inconsistencies
		fi

		#Take header
		mv ${smetfilename}.tmp1 ${mergedsmetfilename}

		#Add data. First sort by time stamp, then sort by second column, which contains 0 for new data and 1 for old data. Then only compare the first 12 characters, which is the time stamp. This
                #gives only the first occurrence of a time stamp in the output (which is the new data, if data differs and also the new data when the data is equal). Then awk is used to remove
                #the first two columns, which were used for the sorting. The final sed removes a space which is added to the end of the line, due to the printf statement.
		cat ${smetfilename}.tmp4 | uniq -w12 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' >> ${mergedsmetfilename}

		#Remove temporary files
		rm ${smetfilename}.tmp2 ${smetfilename}.tmp3 ${smetfilename}.tmp4 ${smetfilename}
	fi	
  fi


  if (( ${mergesoilmoisture} == 1 )); then
	#Inform user
	echo "Merging vwc files..."

	#First check if structure is the same (note that we here can anticipate to the use of # for grepping the header)
  	if (( `cat ${soilmoisturefilename} ${mergedsoilmoisturefilename} | sort | grep fields | uniq | wc -l` != 1 )); then
		echo "ERROR: cannot merge, because fields are not the same!"
	else

		#Then write out header (take the newest one)
		cat ${soilmoisturefilename} | grep -v ^[0-9] >> ${soilmoisturefilename}.tmp1

		#Then write out data from both files, including a time stamp, and a 0 to denote the new data and a 1 to denote old data (this is used to select the newer data in case there is different data for the same time stamp). The sed transforms YYYY-MM-DDTHH:mm to YYYYMMDDHHmm, which makes direct comparison/sorting possible. 
		cat ${soilmoisturefilename} | grep ^[0-9] | awk '{print $1, 0, $0}' | sed 's/\([0-9][0-9][0-9][0-9]\)-\([0-9][0-9]\)-\([0-9][0-9]\)T\([0-9][0-9]\)\:\([0-9][0-9]\)/\1\2\3\4\5/1' > ${soilmoisturefilename}.tmp2
		cat ${mergedsoilmoisturefilename} | grep ^[0-9] | awk '{print $1, 1, $0}' | sed 's/\([0-9][0-9][0-9][0-9]\)-\([0-9][0-9]\)-\([0-9][0-9]\)T\([0-9][0-9]\)\:\([0-9][0-9]\)/\1\2\3\4\5/1' > ${soilmoisturefilename}.tmp3
  
		#Then merge files and sort:
		cat ${soilmoisturefilename}.tmp2 ${soilmoisturefilename}.tmp3 | sort -n -k 1 -k 2 > ${soilmoisturefilename}.tmp4
		#								^^^ this first select time stamps, and in case of multiple occurrences of equal timestamps, select the one from the new dataset, denoted by a 0 in column 2.

		#Now check for equal time stamps, with different data. If so, give warning.
		#The first uniq removes all really unique lines, the second one checks if still some time stamps appear more than once, which means the data is different.
		#Note, the awk is to get rid of the first 2 columns, because the second column is always different between the two files, as it indicates from which file the line is coming.
		#That's why the uniq has to check 16 characters, as the first column is now the date in the YYYY-MM-DDTHH:mm format.
		if (( `cat ${soilmoisturefilename}.tmp4 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' | uniq | uniq -w16 -dD | wc -l` > 0 )); then
			echo "WARNING: duplicate time stamps with different data encountered! Using data from the new file."
			echo "   Inconsistencies are saved in: ${soilmoisturefilename}.inconsistencies."
			cat ${soilmoisturefilename}.tmp4 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' | uniq | uniq -w16 -dD > ${soilmoisturefilename}.inconsistencies
		fi

		#Take header
		mv ${soilmoisturefilename}.tmp1 ${mergedsoilmoisturefilename}

		#Add data. First sort by time stamp, then sort by second column, which contains 0 for new data and 1 for old data. Then only compare the first 12 characters, which is the time stamp. This
                #gives only the first occurrence of a time stamp in the output (which is the new data, if data differs and also the new data when the data is equal). Then awk is used to remove
                #the first two columns, which were used for the sorting. The final sed removes a space which is added to the end of the line, due to the printf statement.
		cat ${soilmoisturefilename}.tmp4 | uniq -w12 | awk '{for (i=3; i<=NF; i++) printf "%s ", $i} {printf "\n"}' | sed 's/ $//' >> ${mergedsoilmoisturefilename}

		#Remove temporary files
		rm ${soilmoisturefilename}.tmp2 ${soilmoisturefilename}.tmp3 ${soilmoisturefilename}.tmp4 ${soilmoisturefilename}
	fi	
  fi

  unset TZ	#Unset the time zone, which was set to UTC.
done
