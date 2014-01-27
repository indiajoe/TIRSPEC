#!/bin/bash
#This script will help user to rename files in bulk
# Be very carefull. Backup First
#-----------------------------indiajoe

TEXTEDITOR="emacs -nw"

echo "Very IMPORTANT:  Backup your files first before proceeding."
echo "Edit only the file names you want to change in the text file. Absolutly Nothing Else!!"
echo "No spaces allowed in filename."
echo "Don't give a filename which already exists."
echo "We recommend not to change the filenumber in the filename."
echo "File to edit will open in $TEXTEDITOR . Save and exit after changes."
read -p "Press enter to continue :" junk
cp SlopeimagesLog.txt SlopeimagesLog-NEW.txt
$TEXTEDITOR SlopeimagesLog-NEW.txt

gawk '{print $1}' SlopeimagesLog.txt > temp1rename.txt
gawk '{print $1}' SlopeimagesLog-NEW.txt > temp2rename.txt
diff --side-by-side --suppress-common-lines temp1rename.txt temp2rename.txt > diff2rename.txt

while read -r line
do
    linearray=( $line )
    echo "Moving ${linearray[0]} to ${linearray[2]}"
    mv -n ${linearray[0]} ${linearray[2]}
done  < diff2rename.txt

echo "Backuping and updating SlopeimagesLog.txt"
mv SlopeimagesLog.txt SlopeimagesLog-BACKUP.txt
mv SlopeimagesLog-NEW.txt SlopeimagesLog.txt

