#!/bin/bash
#This script is just to start the 4 terminals and dv for running tirspec.
#-----------------------------------------indiajoe

binDIR=/home/tirspec/bin.sw/  #Location of binaries
logDIR=/home/tirspec/LOGS   #Location to write logs
scriptDIR=/home/tirspec/   #Location to shell scripts

xmessage -buttons Continue:10,Cancel:20 -default Continue "Before Continuing, make sure no other instances of tirspec software is running. " 
if [ $? -eq 20 ] ; then
    exit 1
fi

xmessage -buttons Continue:10,Cancel:20 -default Continue "Before Continuing, make sure the Temperature of detecter is less than 80K."
if [ $? -eq 20 ] ; then
    exit 1
fi


NOW=$(date +"%Y%m%d-%H%M%S")
#gnome-terminal --working-directory=$binDIR -e " bash -c \"./mkmac_h1_as | tee $logDIR/mkmac-$NOW.log ; exec bash \" "
gnome-terminal --working-directory=$binDIR -e " bash -c \"./run_mkmac_h1_as.sh ; exec bash \" " &

sleep 4
gnome-terminal --working-directory=$binDIR -e " bash -c \"./run_motor_server.sh ; exec bash \" " &

sleep 3
gnome-terminal --working-directory=$binDIR -e " bash -c \"./run_icserver.sh ; exec bash \" " &

sleep 10
gnome-terminal --working-directory=$binDIR -e " bash -c \"./client.sh ; exec bash \" " &

sleep 1
gnome-terminal --working-directory=$scriptDIR -e " bash -c \"./run_dv_LampAlert.sh ; exec bash \" "  &
gnome-terminal --working-directory=$scriptDIR -e " bash -c \"./lastNDR.sh ; exec bash \" "  &






