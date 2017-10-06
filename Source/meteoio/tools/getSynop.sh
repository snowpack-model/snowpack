#!/bin/bash
#Please look for the station code on http://www.ogimet.com/indicativos.phtml.en

if [ $# -ne 3 ]; then
	printf "$0\tdownload and convert the synop data for a given station\n"
	printf "\tUsage: $0 {station_index} {start_date} {end_date}\n\n"
	printf "\tYou can get the station index on http://www.ogimet.com/indicativos.phtml.en\n"
	exit
fi

STATION=$1
START=`echo $2 | tr -d '-' | tr -d 'T' | tr -d ':' `
END=`echo $3 | tr -d '-' | tr -d 'T' | tr -d ':' `

if [ ${#START} -eq 8  ]; then
	START="${START}0000"
fi
if [ ${#START} -eq 10  ]; then
	START="${START}00"
fi
if [ ${#END} -eq 8  ]; then
	END="${END}0000"
fi
if [ ${#END} -eq 10  ]; then
	END="${END}00"
fi


CMD="http://www.ogimet.com/cgi-bin/getsynop?block=${STATION}&begin=${START}&end=${END}"
#printf "Get SYNOP data for $STATION between $START and $END\n\t-> ${CMD}\n"
OUTPUT="${STATION}.csv"
curl ${CMD} -o ${OUTPUT}

cat ${OUTPUT} | tr ',' ' ' | awk '
	BEGIN {
		printf("SMET 1.1 ASCII\n[HEADER]\n")
		printf("station_id = %s\nstation_name = %s\n", "'"${STATION}"'", "'"${STATION}"'")
		printf("latitude = \nlongitude = \naltitude = \nnodata = -999\ntz = 0\n")
		printf("fields = timestamp TA TA_max TA_min TD P PSUM PSUM_dur PSUM24 HS\n")
		printf("[DATA]\n")
	}
	function reset()
	{
		special=0
		header=1
		TA=-999
		TA_max=-999
		TA_min=-999
		TD=-999
		P=-999
		PSUM=-999
		PSUM_dur=-999
		HS=-999
	}
	function parseTA(str)
	{
		sg=1
		sg_str=substr(str, 2, 1)
		if (sg_str=="1") sg=-1
		val=substr(str, 3, 3)
		temp=sg*val/10
		if (temp>60 || temp<-80) return -999
		return temp
	}
	function parseTD(str)
	{
		sg=1
		sg_str=substr(str, 2, 1)
		if (sg_str=="1") sg=-1
		if (sg_str=="9") return -999; #this is RH
		val=substr(str, 3, 3)
		temp=sg*val/10.
		if (temp>60 || temp<-80) return -999
		return temp
	}
	function parseP(str)
	{
		val=substr(str, 2, 4)
		press=val*10; #from .1mb to Pa
		return press
	}
	function parsePSUM24(str)
	{
		val=substr(str, 2, 4)
		prec=val/10.;
		return prec
	}
	function parseHS(str)
	{
		val=substr(str, 3, 3)
		hs=val;
		return hs
	}
	function parsePSUM(str)
	{
		dur=substr(str, 5, 1)
		if (dur=="1") PSUM_dur=1*6
		if (dur=="2") PSUM_dur=2*6
		if (dur=="3") PSUM_dur=3*6
		if (dur=="4") PSUM_dur=4*6
		if (dur=="5") PSUM_dur=1
		if (dur=="6") PSUM_dur=2
		if (dur=="7") PSUM_dur=3
		if (dur=="8") PSUM_dur=9
		if (dur=="9") PSUM_dur=15
		if (dur=="/") PSUM_dur=24
		val=substr(str, 2, 3)
		if (val>=990) {
			last_dg=substr(val, 3, 1)
			prec=0.1*last_dg
		} else
			prec=val
		return prec
	}
	/'"${STATION}"'/ {
		reset()
		datum=sprintf("%s-%s-%sT%s:%s", $2, $3, $4, $5, $6)
		for (ii=12; ii<=NF; ii++) {
# 			if ($(ii)=='"${STATION}"') {
# 				header=0
# 				continue
# 			}
# 			if (header==1) continue
			if ($(ii)=="333") {
				special=1
				continue
			}
			if ($(ii)=="444") break
			if ($(ii)=="555") break
			
			if (special==1 && match($(ii),"1(1|0)[0-9][0-9][0-9]")>0) {
				TA_max=parseTA($(ii))
			}
			if (special==1 && match($(ii),"2(1|0)[0-9][0-9][0-9]")>0) {
				TA_min=parseTA($(ii))
			}
			if (special==1 && match($(ii),"7[0-9][0-9][0-9][0-9]")>0) {
				PSUM24=parsePSUM24($(ii))
			}
			if (special==1 && match($(ii),"4[0-9][0-9][0-9][0-9]")>0) {
				HS=parseHS($(ii))
			}
			
			if (special==1) continue #skip all remaining special fields
			
			if (match($(ii),"1(1|0)[0-9][0-9][0-9]")>0) {
				TA=parseTA($(ii))
			}
			if (match($(ii),"2(1|0)[0-9][0-9][0-9]")>0) {
				TD=parseTD($(ii))
			}
			if (match($(ii),"3[0-9][0-9][0-9][0-9]")>0) {
				P=parseP($(ii))
			}
			if (match($(ii),"6[0-9][0-9][0-9][0-9]")>0) {
				PSUM=parsePSUM($(ii))
			}
		}
		printf("%s %7.2f %7.2f %7.2f %7.2f %7d %7.1f %4d %7.1f %7.1f\n", datum, TA, TA_max, TA_min, TD, P, PSUM, PSUM_dur, PSUM24, HS)
		reset()
	}
' > "${STATION}".smet

