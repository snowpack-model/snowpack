mkdir -p COSMO-2_combi/smet

# ls COSMO-2_[0-9][0-9][0-9][0-9]

for f in COSMO-2_1996/smet/*
do
	g=$(basename ${f})
	echo "Processing: ${g}"
	
	cat COSMO-2_[0-9][0-9][0-9][0-9]/smet/${g} | awk 'BEGIN {header=1} {if(/SMET 1.1 ASCII/) {data=0}; if(data || header) {print}; if(/\[DATA\]/) {header=0; data=1}}' > COSMO-2_combi/smet/${g}
	# echo ${header}
	# echo ${data}
done
