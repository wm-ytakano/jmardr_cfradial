#!/bin/bash
#
# convert jma polametric radar (2hempa)
# usage:
# sh convert_nc.sh Z__CRJTD_YYYYMMDDHHSS_RDR_JMA_GPV_..._ANAL_girb2.bin.tgz
# output is Z__CRJTD_YYYYMMDDHHSS_RDR_JMA_GPV_...._ANAL.nc
# 
#
# set python script for convert wgrib2 to cfradial 
wgrib2ncrad=../jmardr_cfradial/convert_2hempa.py
marge_ncrad=../jmardr_cfradial/marge_netcdf_v1.py
#
#
tgz=$1
nc=${tgz%_*}.nc
files=`tar zxvf  $tgz`
declare -a v=()  
for i in $files ; do
	if [ "`echo ${i%_*} | grep 'N'`" ]; then # for remove "./" in targz file
		v+=("${i%_*}.nc")
		python3 $wgrib2ncrad  ${i%_*}_grib2.bin ${i%_*}.nc &
	fi
done
wait 
# marging
# sort uncompress files by filename
# then that may sorted by ray number or observation time schedule
v_sorted=$( printf "%s\n" "${v[@]}" | sort )
python3 $marge_ncrad ${v_sorted[@]} -o $nc
# clear file
rm  ${v_sorted[@]} *.bin
# end 