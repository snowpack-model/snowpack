#!/bin/sh
#extract a given column from a smet file

if [ $# -lt 1 ]; then
	me=`basename $0`
	printf "Usage: \n"
	printf "\t$me .\n\t\t to show which fields are present for each files in the current directory\n"
	printf "\t$me {smet_file} {parameter}\n\t\t to extract the given parameter out of the given file\n"
	printf "\t$me {smet_file} {parameter} {aggregation}\n\t\t to extract the monthly aggregated given parameter out of the given file\n"
	printf "\t\t where {aggregation} is any of (AVG, MIN, MAX)\n"
	exit 0
fi

INPUT=$1

if [ $# -eq 1 ]; then
	if [ -d ${1} ]; then
		head -50 ${1}/*.smet | awk '
			BEGIN {
				name_length=17
			}
			/station_id/ {
				gsub(/\r/, "")
				station=$3
			}
			/station_name/ {
				gsub(/\r/, "")
				station_name=$3
			}
			/altitude/ {
				gsub(/\r/, "")
				altitude=$3
			}
			/fields/ {
				gsub(/\r/, "")
				station=sprintf("%s-%s", station, station_name)
				if (length(station)>name_length) station=sprintf("%s…", substr(station, 1, name_length-1))
				altitudes[station]=altitude
				for(ii=3; ii<=NF; ii++) {
					if ($(ii)=="timestamp") continue
					if ($(ii)=="julian") continue
					field[station][$(ii)] = 1
					all_fields[$(ii)]++
				}
				next
			}
			END {
				nr_fields=asorti(all_fields, fields_idx)

				format=sprintf("%%4.0f - %%-%ds\t", name_length)
				for(station in altitudes) {
					printf(format, altitudes[station], station)
					for(idx=1; idx<=nr_fields; idx++) {
						field_name=fields_idx[idx]
						field_format=sprintf("%%%ds  ", length(field_name))
						if (field[station][field_name]>0)
							printf(field_format, field_name)
						else
							printf(field_format, " ")
					}
					printf("\n")
				}
			}
		' | sort -n -k1
		exit 0
	fi

	head -50 ${INPUT} | grep "fields" | cut -d'=' -f2 | tr -s '  \t' ' ' | xargs -i echo "Available fields: {}"
	exit 0
fi

FIELD=$2
if [ $# -eq 3 ]; then
	AGG_TYPE=$3
fi

#create data sets metadata
field_nr=$(head -50 ${INPUT} | grep "fields" | awk '
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

#get generic info
stat_id=`head -50 ${INPUT} | grep "station_id" | tr -s '\t' ' ' | cut -d' ' -f 3-`
stat_name=`head -50 ${INPUT} | grep "station_name" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lat=`head -50 ${INPUT} | grep "latitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lon=`head -50 ${INPUT} | grep "longitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
alt=`head -50 ${INPUT} | grep "altitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
JULIAN=`head -50 "${INPUT}" | grep fields | grep julian`

awk '
	BEGIN {
		field="'"${field_nr}"'"
		agg_type="'"${AGG_TYPE}"'"
		if ("'"${JULIAN}"'"!="") isJulian=1
		printf("#%s - %s\n", "'"${stat_id}"'", "'"${stat_name}"'")
		printf("#lat=%s - lon=%s - alt=%s\n", "'"${lat}"'", "'"${lon}"'", "'"${alt}"'")
	}
	function getISO(ts){
		nr_secs=(ts-2440587.5)*24.*3600.
		return sprintf("%s", strftime("%FT%H:%M:%S", int(nr_secs+0.5))) #rounding to nearest second
	}
	/^[0-9][0-9][0-9][0-9]/ {
		if (agg_type=="") {
			if (isJulian==0) datum=$1
			else datum=getISO($1)
			printf("%s %s\n", datum, $(field))
		} else {
			if (isJulian==0) datum=$1
			else datum=getISO($1)
			gsub(/\-|\:|T/," ", datum); split(datum,d," ");
			key=sprintf("%s-%s-01", d[1], d[2])
			
			if (agg_type=="AVG") {
				agg[key] += $(field)
				count[key]++
			} else if (agg_type=="MIN") {
				if ($(field)<agg[key] || count[key]==0) agg[key]=$(field)
				count[key]++
			} else if (agg_type=="MAX") {
				if ($(field)>agg[key] || count[key]==0) agg[key]=$(field)
				count[key]++
			}
		}
	}
	/^units_offset/ {
		offset=$(field+2)
	}
	/^units_multiplier/ {
		multiplier=$(field+2)
	}
	/^plot_unit/ {
		unit=$(field+2)
	}
	/^plot_description/ {
		description=$(field+2)
	}
	/^\[DATA\]/ {
		printf("#offset=%g - multiplier=%g - unit=\"%s\" - description=\"%s\"\n", offset, multiplier, unit, description)
	}
	END {
		if (agg_type=="AVG") {
			n = asorti(agg, data)
			for(i=1; i<=n; i++) {
				idx=data[i]
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx]/count[idx])
				else printf("%s -999\n", idx)
			}
			#for(idx in agg) {
			#	if (count[idx]>0) printf("%s %f\n", idx, agg[idx]/count[idx])
			#	else printf("%s -999\n", idx)
			#}
		} else if (agg_type=="MIN") {
			for(idx in agg) {
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx])
				else printf("%s -999\n", idx)
			}
		} else if (agg_type=="MAX") {
			for(idx in agg) {
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx])
				else printf("%s -999\n", idx)
			}
		}
	}
' ${INPUT}


