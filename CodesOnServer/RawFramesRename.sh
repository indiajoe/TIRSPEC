#!/usr/bin/env bash
# This script is to rename all the raw frames an exposure of tirspec.
#-------------------------------------indiajoe

echo "Renaming files in present working directory:" $(pwd)
if [ "$#" -le 2 ] ; then  #User input check.
  echo "To rename OLD_FILENAME to NEW_FILENAME with FILENUMBER1">&2
  echo "Usage: $0 OLD_FILENAME NEW_FILENAME FILENUMBER1 FILENUMBER2 ..." >&2
  exit 1
fi

oldfilename="$1"
newfilename="$2"
for fnumber in ${@:3} ; do
    echo "**************************"
    echo "Filenumber: $fnumber"
    for file in $oldfilename-$fnumber[-.]*fits ; do
	echo "$file --> ${file/$oldfilename/$newfilename}" 
    done
    echo "----"
    echo "Do you want to rename all the files as shown above?"
    choice='NO'
    echo -n "Enter YES to proceed with renaming:"
    read choice
    if [ "$choice" == "YES" ] ; then
	echo "Renaming the files of file number: $fnumber"
	for file in $oldfilename-$fnumber[-.]*fits ; do
	    mv "$file" "${file/$oldfilename/$newfilename}" 
	done
    else
	echo "Aborting the Renameing process for file number: $fnumber"
    fi
done
echo "Finished Renaming Task."

	
