ignore_slope=1		# Ignore slope information from smet file
for f in ./smet/*smet
do
	stn=$(basename ${f} .smet)
	snofile="./snow_init/${stn}.sno"

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
done
