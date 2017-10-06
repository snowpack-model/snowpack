#!/bin/bash
#from a list of SMET files, this generates a KML file that can be read in google earth or map.geo.admin
#to represent on a map where the stations are
#see https://developers.google.com/kml/documentation/kml_tut for more on kml

if [ $# -lt 1 ]; then
	INPUT_DIR="."
else
	INPUT_DIR=$1
fi

cs2cs=`which cs2cs`

ls ${INPUT_DIR}/*.smet | xargs -i head -50 {} | awk '
	BEGIN {
		printf("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
		printf("<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n")
		printf("<Folder>\n<name>My Simulation</name>\n")
		printf("<Style id=\"sty0\">\n<LabelStyle>\n")
		printf("<color>ff0000ff</color>\n<scale>1</scale>\n")
		printf("</LabelStyle>\n</Style>\n")
		latitude=-99
		easting=-99
		if (length("'"${cs2cs}"'")>0) has_cs2cs=1
	}
	/\[DATA\]/ {
		if (has_cs2cs==1 && latitude==-99) { #we need to convert east/north to lat/lon
			cmd=sprintf("echo \"%f %f\" | cs2cs +init=epsg:%d +to +init=epsg:4326 -f \"%%.10f\"", easting, northing, epsg)
			cmd | getline
			longitude = $1
			latitude = $2
		}

		printf("<Placemark>\n")
		printf("<name>%s (%d)</name>\n", station_id, altitude)
		printf("<styleUrl>#sty0</styleUrl>\n")
		printf("<description>%s (%d)</description>\n", station_name, altitude)
		printf("<ExtendedData><Data name=\"fields\"><value>%s</value></Data></ExtendedData>\n", fields)
		printf("<Point><coordinates>%f, %f, %d</coordinates></Point>\n", longitude, latitude, altitude)
		printf("</Placemark>\n")
		latitude=-99
		easting=-99
		#nextfile
	}
	/station_id/ {
		gsub(/\r/, "")
		station_id = $3
	}
	/station_name/ {
		gsub(/\r/, "")
		station_name = $3
	}
	/latitude/ {
		latitude = $3
	}
	/longitude/ {
		longitude = $3
	}
	/altitude/ {
		altitude = $3
	}
	/easting/ {
		easting = $3
	}
	/northing/ {
		northing = $3
	}
	/epsg/ {
		epsg = $3
	}
	/fields/ {
		fields=""
		for (ii=3; ii<=NF;ii++) {
			if ($(ii)=="timestamp" || $(ii)=="julian") continue;
			fields=sprintf("%s%s ", fields, $(ii))
		}
		fields=substr(fields, 1, length(fields)-1) #remove the last space
	}
	END {
		printf("</Folder>\n</kml>\n")
	}
' 

