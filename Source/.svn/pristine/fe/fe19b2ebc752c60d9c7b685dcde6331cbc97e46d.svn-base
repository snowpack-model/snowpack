#!/bin/sh
#extract a given column from a smet file

INPUT=$1
FIELD=$2

if [ $# -eq 1 ]; then
	head -20 ${INPUT} | grep "fields" | cut -d'=' -f2 | tr -s ' \t' '\t' | xargs -i echo "Available fields: {}"
fi

#get generic info
stat_id=`head -20 ${INPUT} | grep "station_id" | tr -s '\t' ' ' | cut -d' ' -f 3-`
stat_name=`head -20 ${INPUT} | grep "station_name" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lat=`head -20 ${INPUT} | grep "latitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lon=`head -20 ${INPUT} | grep "longitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
alt=`head -20 ${INPUT} | grep "altitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`

#create data sets metadata
field_nr=$(head -20 ${INPUT} | grep "fields" | awk '
	/fields/ {
		found=3
		for(i=1; i<=NF; i++) {
			if($(i)=="'${FIELD}'") found=i
		}
		printf("%d\n", found-2)
	}
')

if [ ${field_nr} -eq 1 ]; then
	exit
fi

#out_name="${stat_id}_${FIELD}.dat"
out_name="${stat_id}_${alt}.dat"
printf "#${stat_id} - ${stat_name}\n#lat=${lat} - lon=${lon} - alt=${alt}\n" > ${out_name}
grep -E "^[0-9]{4}-[0-9]{2}-[0-9]{2}" ${INPUT} | tr -s ' ' | cut -d' ' -f1,${field_nr} >> ${out_name}

