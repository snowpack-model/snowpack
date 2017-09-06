#!/bin/bash
#This script converts a given smet file into an XmGrace format
#to automatically, on-the-fly convert smet to agr directly in xmgrace, add the following line to your gracerc:
# DEFINE IFILTER "smet2agr.sh %s" PATTERN "*.smet"

INPUT=$1

#create general header section
printf "# Grace project file\n#\n"
printf "@version 50121\n"
printf "@page size 792, 612\n"
datum=$(date --rfc-3339=seconds)
printf "@timestamp def \"${datum}\"\n"

#create g0 graph
stat_id=$(head -20 ${INPUT} | grep "station_id" | tr -s '\t' ' ' | cut -d' ' -f 3-)
stat_name=$(head -20 ${INPUT} | grep "station_name" | tr -s '\t' ' ' | cut -d' ' -f 3-)
lat=$(head -20 ${INPUT} | grep "latitude" | tr -s '\t' ' ' | cut -d' ' -f 3-)
lon=$(head -20 ${INPUT} | grep "longitude" | tr -s '\t' ' ' | cut -d' ' -f 3-)
alt=$(head -20 ${INPUT} | grep "altitude" | tr -s '\t' ' ' | cut -d' ' -f 3-)
printf "@g0 on\n"
printf "@g0 hidden false\n"
printf "@g0 type XY\n"
printf "@g0 stacked false\n"
printf "@g0 bar hgap 0.000000\n"
printf "@g0 fixedpoint off\n"
printf "@g0 fixedpoint type 0\n"
printf "@g0 fixedpoint xy 0.000000, 0.000000\n"
printf "@g0 fixedpoint format yymmddhms general\n"
printf "@g0 fixedpoint prec 6, 6\n"
printf "@with g0\n"
printf "@    world 2454800, -1000, 2454820, 400\n"
printf "@    stack world 0, 0, 0, 0\n"
printf "@    znorm 1\n"
printf "@    view 0.150000, 0.150000, 1.150000, 0.850000\n"
printf "@    title \"${INPUT}\"\n"
printf "@    subtitle \"${stat_id} - ${stat_name} (${lat}, ${lon}, ${alt})\"\n"
printf "@    legend 0.2, 0.8\n"
printf "@    legend char size 0.500000\n"
printf "@    xaxis  label \"date\"\n"
printf "@    xaxis  ticklabel on\n"
printf "@    xaxis  ticklabel format yymmdd\n"

#create data sets metadata
columns=$(head -20 ${INPUT} | grep "fields")
echo "${columns}" | awk '
	/fields/ {
		for(i=4; i<=NF; i++) {
			f=i-4
			printf("@    s%d hidden false\n", f)
			printf("@    s%d type xy\n", f)
			printf("@    s%d comment \"%s\"\n", f, $(i))
			printf("@    s%d legend  \"%s\"\n", f, $(i))
		}
	}
'

#create data sets data
nb_sets=$(echo "${columns}" | wc -w)
#for i in {4..${nb_sets}}; do
for (( i=4 ; i<=${nb_sets} ; i++ )); do
	f=$(( i-4 ))
	printf "@target G0.S${f}\n@type xy\n"
	awk '
	BEGIN {
		field='${i}'-2
	}
	/^[[:space:]]*[0-9\-]+/ {
		printf("%s %s\n",$1, $(field))
		}' ${INPUT}
	printf "&\n"
done
