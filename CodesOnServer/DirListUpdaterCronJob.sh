#!/bin/bash
# This script is run by cronjob of tirspec user at 14:00hrs IST everyday.
# The purpose is to update various directory list for the Data Backup Script.
# ------------------------------------------------------indiajoe

DataPath='/data/'
Age=1  # In Units of Days, how old directories should be to conisdered to be copied

NewDirsFile='/data/NewDirsFormed.txt'
OldDirsFile='/data/OldDirsList.txt'
LogFile='/home/tirspec/DataBackup.log'
ToBeCopiedFile='/data/DirsToBeCopied.txt'

if [ ! -f "$OldDirsFile" ]; then  # Atleast an empty file should be there
    touch "$OldDirsFile"
fi

#Create a list of new directories (older that $Age) which are not in $OldDirsFile file
find "$DataPath" -maxdepth 1 -mindepth 1 -type d -mtime +$Age -print | grep -v -f  <( cut -d' ' -f 1 "$OldDirsFile")  > "$NewDirsFile"
#ALERT: Due to the previous line, directories are not allowed to have any space in their name

for dir in $(cat "$NewDirsFile"); do
    dirSize=$(du -s "$dir" | awk '{print $1}')
    echo "$dir $dirSize" >> "$ToBeCopiedFile"
    echo "$dir $dirSize" >> "$OldDirsFile"
    echo "$(date) : Added $dir :$dirSize to $ToBeCopiedFile" >> "$LogFile"
done
