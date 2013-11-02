#!/bin/bash
#This script is to be run after the nigth's observations to generate the slope images and also sent the final images to CREST, Hoskote.
#-------------------------------------indiajoe

SlopeScript='/home/tirspec/bin/Generate4mRamp.py'
PhotometryBadpixelfile='/home/tirspec/SlopeTIRSPECphoto-BP.pl'
SpectroscopyBadPixelfile='/home/tirspec/SlopeTIRSPECspec-BP.pl'

echo "-----------------------------------"
echo "Welcome, TIRSPEC Observer...."
echo "We hope you had a great observation night."
echo "-----------------------------------"
echo
echo "Enter below the night directories you want to process and copy."
echo "You can give more than one directory by space seperated input"
echo "Eg: /data/20130125 /data/20130126 "
echo
echo -n "Enter directory (s) :"
while IFS=' ' read -ra inpnightsarray 
do
    nightarray=()
    for dir in "${inpnightsarray[@]}"; do
	if [[ -d "$dir" ]] ; then
	    nightarray+=("$dir")  #Adding to the array of nights
	else
	    echo "Warning: Cannot find the directory : $dir"
	fi
    done
    if [[ ${#nightarray[@]} > 0 ]]; then #Length of array is > 0
	break
    else
	echo "None of the directries you entered cannot be not found. Pls try again "
	echo -n "Enter directory (s) :"
    fi
done

echo "Proceeding with the processing of following directories :  ${nightarray[@]}"


echo "Which Machine at CREST do you want to download the data into ?"
echo "---------------------------------------"
echo " 0 : Do not download data to anywhere "
echo " 1 : Download to 192.168.1.23 "
echo " 2 : Download to 192.168.1.14 "
echo "---------------------------------------"
echo -n "Enter Serial number of your choice :"
while read -n 1 choice
do
    case $choice in
	0) 
	    MachineIP='NIL' 
	    echo
	    echo "No downloading of data wil be done." 
	    break ;;
	
	1) 
	    MachineIP='192.168.1.23' 
	    echo
	    echo "Final Data will be downloaded to $MachineIP"
	    break ;;	
	2) 
	    MachineIP='192.168.1.14' 
	    echo
	    echo "Final Data will be downloaded to $MachineIP"
	    echo "Sorry: Not yet setup. cannot be done..."
	    break ;;
	*)
	    echo
	    echo "Pls enter a valid choice: 0,1, or 2" 
	    ;;
    esac
done

for dir in "${nightarray[@]}"; do
    if [[ -e "$dir/SlopeimagesLog.txt" ]] ; then
	echo "$dir has previously generated slope images."
    else
	echo "Generating slope images in this directory : $dir"
	echo "This will take time, you can go to sleep..."
	echo
	date | tee -a $dir/SlopeGenerationScript.log
	#generating the slope images and also write a log
	$SlopeScript $dir | tee -a $dir/SlopeGenerationScript.log
	date | tee -a $dir/SlopeGenerationScript.log
	echo "Generation of slope images over. Now Copying the bad pixel masks to $dir"
	#Copying the Badpixel maps to the directory
	cp $PhotometryBadpixelfile $dir
	cp $SpectroscopyBadPixelfile $dir
    fi
done

#Proceeding to upload the files to CREST..

if [[ $MachineIP != 'NIL' ]] ; then
    for dir in "${nightarray[@]}"; do
	echo "Uploading $dir files to ~/TIRSPECdata/ in $MachineIP..."
	date
	rsync -e ssh -avzR $dir/Slope* observer@$MachineIP:~/TIRSPECdata/
	date
    done
    echo "All uploads to CREST over..."
fi
echo "Thanks for using TIRSPEC."
	
