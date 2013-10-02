#!/bin/bash
#This script is to update the Macro of DV to show latest NDR read out.
#-------------------------------indiajoe

dir=$(pwd)
while :
do
    xmessage -buttons Update:10,ChangeDir:20 -default Update "Latest in : $dir"
    case $? in
	10 ) latestfile=$( find $dir -type f | xargs stat --format '%Y : %y %n ' | sort -nr | awk 'NR==2 {print $NF ; exit}' )
	echo "Read $latestfile C" > /home/tirspec/DV/macros/lastNDR.mac ;;

	20 ) read -p "Enter New Data Directory: " dir ;;
    esac
done