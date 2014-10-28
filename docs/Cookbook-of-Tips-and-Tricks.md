# Tips and Tricks

###0) Save all the steps done in a reduction attempt###
Before you start reduction of a source. Run the script command as shown in the example below

`script LogOfMyReductionOfM31.txt`

This command will will save everything that is typed or printed in that terminal to the file  _LogOfMyReductionOfM31.txt_ while you are reducing data. At the end of reduction, press **Cntrl+D** to stop and exit the recording.

###1)  List all the observations done so far on a particular object###
Say, you are looking for a source with keyword M31 in it. 

Running all the following commands standing in directory containing all the night's raw data.

Following command will list details of all observations which has the keyword M31.

`grep M31 201*/SlopeimagesLog.txt`

If you suspect, some nights the name used might be Andromeda instead of M31, you can combine all the keywords in grep as shown below.

`grep 'M31\|Andromeda' 201*/SlopeimagesLog.txt`

Once, you have all the list related to observation of your source, one can use the power of pipes in gnu/linux to extract out the details just you want.

Some examples are given below.

* List directories containing M31 observations: `grep M31 201*/SlopeimagesLog.txt | cut -b -8 | uniq`

* List directories containing M31 observations along with number of image frames in that night: `grep M31 201*/SlopeimagesLog.txt | cut -b -8 | uniq -c`

* List only spectroscopy observations of M31: `grep M31 201*/SlopeimagesLog.txt | grep ' G [LS]'`

* List only Imaging observations of M31: `grep M31 201*/SlopeimagesLog.txt | grep -v ' G [LS]'`


###2) Use flats from a different night###
There are three ways to do it. All are equally good ways to do it.

1. Copy those files as well as the lines corresponding to those images from SlopeimagesLog.txt of that night to new night directory. It is very important to append the SlopeimagesLog.txt with the new additions in directory.

2. To save hardisk space, instead of actually copying the data, you can instead create a symbolic link to those images from one directory to other. But do append the SlopeimagesLog.txt with the new additions.

3. But if your hardisk is ntfs you may not be able to create symbolic links. You can use flats form other nights still by appending those images in SlopeimagesLog.txt, but replace the filename in the first column with the relative path. Eg: ../20140329/Slope-Flat_J-23.fits.
That way, you can avoid copying those fits file to each night's directory.

Note: Do not change anything in other columns of the SlopeimagesLog.txt they are all vital to run the code. When you are copy pasting the image lines from one directory to other, care should be take to preserve all the columns while pasting.

###3) Continue the pipeline run, only on a subset of directories initially chosen###
Due to some evil reason, your script might have crashed while reducing some particular night's data.
You don't want to run all the successful runs upto that directory again, when you rerun the piepeline.

All you have to do is to backup the OUTPUTDIR/directories text file. And then remove the name of directories you don't want to run the pipeline.
After finishing the run on remaining directories, you can choose to restore the directories list.

Note: In Stage 9 of Photometry pipeline (the last actual photometry part). You may want to choose at image level instead of directories. The file you have to backup and edit is OUTPUTDIR/Images4Photo.in file. After Step 6 in Photometry pipeline, OUTPUTDIR/directories does not have any role in deciding the pipelines action.

###4) Some general tips which will help users ###

* If you know which frames are totally useless from your observation book written during night. After you get the data, simply prefix the line corresponding to those images in SlopeimageLog.txt with # to tell data pipeline to completely ignore those files.
* Set focus-follows-mouse  : This will make life a lot easier and will help to avoid many tiresome clicks while switching between windows. (this is especially difficult in OS X)
* Load the FoV_of_tirspec.vot in Aladdin and Flip X axis in ds9 to match TIRSPEC's FoV with 2MASS in Aladdin.
* Reduce all the data of one object spanning over many nights at a time. Pipeline is designed for that.
* Always download the latest nightly tarballs from github page. Updates are very frequent.
