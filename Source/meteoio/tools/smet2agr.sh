#!/bin/sh
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

#print the color table and generate the table of indices
plot_colors=`head -100 ${INPUT} | grep "plot_color" | tr -s '\t' ' ' | cut -d'=' -f2 | tr ' ' '\n'`
if [ "x$plot_colors" != "x" ]; then
	sorted_plot_colors=`echo "${plot_colors}" | sort -u | grep -v "0x000000" | grep -v "-" | grep -v "0xFFFFFF" | tr '\n' ' ' | awk '{printf("white black grey %s", $0)}' `
	color_table=`echo "${plot_colors}" | awk '
		BEGIN {
			split("'"${sorted_plot_colors}"'", sorted, " ")
		}
		!/^$/{
			for (color_idx in sorted) {
				if ($1=="0xFFFFFF") { #white is always at position 0
					printf("0 ")
					next
				}
				if ($1=="0x000000") { #black is always at position 1
					printf("1 ")
					next
				}
				if ($1=="-") { #undef is always at position 2
					printf("2 ")
					next
				}
				
				if ($1==sorted[color_idx]) {
					printf("%d ", color_idx-1)	#we must start with white=0
					next
				}
			}
		} '
	`

	echo "${sorted_plot_colors}" | tr ' ' '\n' | awk --non-decimal-data '
		BEGIN {
			printf("@map color 0 to (255, 255, 255), \"white\"\n")
			printf("@map color 1 to (0, 0, 0), \"black\"\n")
			printf("@map color 2 to (160, 160, 160), \"grey\"\n")
			col_count=2
		}
		!/^$/{
			if (substr($1,1,2)!="0x") next	#color names are handled differently
			r=int( sprintf("%f", "0x" substr($1,3,2)) )
			g=int( sprintf("%f", "0x" substr($1,5,2)) )
			b=int( sprintf("%f", "0x" substr($1,7,2)) )
			color_name=$1
			if (r==0 && g==0 && b==0) next

			col_count++
			printf("@map color %d to (%d, %d, %d), \"%s\"\n", col_count, r,g,b, color_name)
		}
	'
fi

#create g0 graph
stat_id=`head -100 ${INPUT} | grep "station_id" | tr -s '\t' ' ' | cut -d' ' -f 3-`
stat_name=`head -100 ${INPUT} | grep "station_name" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lat=`head -100 ${INPUT} | grep "latitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lon=`head -100 ${INPUT} | grep "longitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
alt=`head -100 ${INPUT} | grep "altitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
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
columns=$(head -100 ${INPUT} | grep "fields")
echo "${columns}" | awk '
	BEGIN {
		n=split("'"${color_table}"'", color_table, " ")
	}
	/fields/ {
		for(i=4; i<=NF; i++) {
			f=i-4
			printf("@    s%d hidden false\n", f)
			printf("@    s%d type xy\n", f)
			printf("@    s%d comment \"%s\"\n", f, $(i))
			printf("@    s%d legend  \"%s\"\n", f, $(i))
			if (n>0) printf("@    s%d line color %d\n", f, color_table[f+2])
		}
	}
'

#create data sets data
nb_sets=$(echo "${columns}" | wc -w)
for i in $(seq 4 ${nb_sets}); do
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
