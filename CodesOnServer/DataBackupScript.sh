#!/bin/bash
# This script is to semi-automate the data backup to externally connected hardisk
# It relies on the files created by the DirListUpdaterCronJob.sh cronjob
# Run this script after connecteing an external hardisk
# ------------------------------------------------indiajoe

#DataPath='/data/'
LogFile='/home/tirspec/DataBackup.log'
ToBeCopiedFile='/data/DirsToBeCopied.txt'
AlreadyCopiedFile='/data/DirsAlreadyCopiedButInTransit.txt'
DirAndChecksumFile='DirAndChecksum.txt'

# First we have to make sure only one instance of this script is running using flock
ScriptLock='/tmp/DataBackupScript.LOCK'
exec 8> $ScriptLock  #I picked 8 arbitrarily as my file handle
flock -n -e 8 || exit 1   # Silently exit with 1 if not able to obtain lock

###### Check and select one mounted external disks #####
# Create an array of mounted external storage directories and also free space in each
OLDIFS=$IFS
IFS=$'\n'  #Since the direcotry name may have spaces
External=($(df | awk '/^\/dev\/sd[^a]/{$1=$2=$3=$4=$5=""; print $0}'))
IFS=$OLDIFS

FreeSpace=($(df | awk '/^\/dev\/sd[^a]/{print $4}'))
PercentUse=($(df | awk '/^\/dev\/sd[^a]/{print $5}'))


Lastsd=$(( ${#External[@]} -1))  # Index of Last external hdisk 

if [[ ${#External[@]} == 0 ]] ; then #No external hardisk detected
    echo "No external hardisk mounted to start backup."
    echo "Mount a new external hardisk and call this backup script again."
    echo "Press Enter to close..."
    read junk	
    exit 2

elif [[ ${#External[@]} > 1 ]] ; then  # More than one external storage device mounted
    echo "More than one external hardisk detected. Choose the hardisk to backup data."
    echo "---------------------------"
    for i in $(seq 0 $Lastsd) ; do
	freespace=${FreeSpace[$i]} 
	echo "$i : ${External[$i]} with freespace= ${freespace:0:$((${#freespace}-6))} GB, and Use% = ${PercentUse[$i]}"
    done
    
    echo "---------------------------------------"
    echo -n "Enter Serial number of your choice:"
    while read -n 1 choice
    do
	if [[ $choice < ${#External[@]} && $choice > -1 ]] ; then
	    ExternalHDISK=${External[$choice]}
	    ExternalHDISKFreeSpace=${FreeSpace[$choice]}
	    break
	else
	    echo "Pls enter a valid serial number." 
	fi
    done
else   # Only one external hardisk connected
    ExternalHDISK=${External[0]}
    ExternalHDISKFreeSpace=${FreeSpace[0]}
fi
    

# Strip trailing and leading white spaces in mounted partition name
ExternalHDISK=$(echo $ExternalHDISK | sed -e 's/^ *//' -e 's/ *$//')

#### Check write permission in the hardisk ####
if date > "$ExternalHDISK"/DateofCopying.txt ; then
else
    echo "Not able to write to external hardisk $ExternalHDISK"
    echo "Mount a new external hardisk and call this backup script again."
    echo "Press Enter to close..."
    read junk	
    exit 3
fi


echo "$(date) : Starting backup to $ExternalHDISK with freespace $ExternalHDISKFreeSpace" | tee -a "$LogFile"

#### Find directories which can be copied to this hardisk ####
FilledSize=0
DirsToCopyArray=()
{ while read line
    do
    dirsize=$(echo $line | gawk '{print $NF}')
    if [[ $(($FilledSize + $dirsize)) < $(($ExternalHDISKFreeSpace - 3000000)) ]] ; then  # We leave 3 GB free
	dir=$(echo $line | gawk '{print $1}')
	DirsToCopyArray+=($dir)
    fi
    done } < "$ToBeCopiedFile"


#### Start copying the selected directories to external hardisk ####
echo "$(date) : Directories being copied to $ExternalHDISK :: ${DirsToCopyArray[@]} " | tee -a "$LogFile"

echo > "$ExternalHDISK"/ListofDirs.txt
for dir in ${DirsToCopyArray[@]} ; do
    echo "$dir" >> "$ExternalHDISK"/ListofDirs.txt
done

# Actual copying of file via rsync occcurs here. 
# We shall run it in Idle mode so that it runs only when no other program is asking for disk access. 
# For this we set ionice -c 3

# It is a shame, we cannot use ionice -c3 as tirspec user since kernal (2.6.18) is older than 2.6.25 
# Hence I am using "nice -n18 ionice -c2 -n7" instead
#ionice -c3 rsync -avr --files-from="$ExternalHDISK"/ListofDirs.txt "$ExternalHDISK/"
nice -n18 ionice -c2 -n7 rsync -avr --files-from="$ExternalHDISK"/ListofDirs.txt "$ExternalHDISK/"

if [ $? ] ; then
    echo "$(date) : Succesfully copied all data  " | tee -a "$LogFile"
else
    echo "$(date) : ERROR: Rsync exitied with $? . Something went wrong.  " | tee -a "$LogFile"
fi

#### Calculate md5sum of the directories and stor them in external hardisk ####

for dir in $(cat "$ExternalHDISK"/ListofDirs.txt) ; do
    echo "Calculating md5sum of "$dir
#    md5checksum=$(ionice -c3 find $dir -type f -exec md5sum {} + | awk '{print $1}' | sort | md5sum | cut -d' ' -f1)
    md5checksum=$(nice -n18 ionice -c2 -n7 find $dir -type f -exec md5sum {} + | awk '{print $1}' | sort | md5sum | cut -d' ' -f1)
    echo $dir $md5checksum | tee -a "$ExternalHDISK/$DirAndChecksumFile"
done

cat "$ExternalHDISK/$DirAndChecksumFile" >> "$LogFile"

#### Update the ToBeCopied as well as AlreadyCopied files ####

for dir in $(cat "$ExternalHDISK"/ListofDirs.txt) ; do
    echo $dir >> $AlreadyCopiedFile
    sed -i ":$dir:d" $ToBeCopiedFile
done

#### Ask Operator to unmount he external hardisk and call Hanle to replace Hardisk with a new one ####

for i in $(seq 20); do
    echo -en '\a' 
    sleep 0.7 
done

# Check for free space in $ExternalHDISK/
freespace=$(df | awk '$NF=="$ExternalHDISK"{print $4}')

echo "$ExternalHDISK has now only Free = ${freespace:0:$((${#freespace}-6))} GB" | tee -a "$LogFile"
echo "If it is almost full, Please Unmount the hardisk."
echo "Wait for the hardisk to get safely unmounted."
echo "Afterwards, call Hanle staff and ask them to replace with a new hardisk."
echo "Press Enter to close..."
read junk	

#---------------------------------
# Meaning of different exit codes
# 0 : Everything went well
# 1 : Not ablle to obtain lock, some other instance of this script is already running.
# 2 : No extrnal hardisk detected to be mounted
# 3 : Not able to write to the chosen external hardisk
