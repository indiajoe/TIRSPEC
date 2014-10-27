#Getting Started#

I presume you have already finished installing IRAF and required python modules. If not visit [Installation page](Installing-Pipeline-Requirements.md).

**Step one:** Download the latest version of TIRSPEC data reduction pipeline from [nightly tarball](https://github.com/indiajoe/TIRSPEC/tarball/master).

Extract the downloaded file to some location.
You will find two directories in it. CodesOnServer and CodesForUser. As far as TIRSPEC data reduction is concerned, you can completely ignore codes in CodesOnServer directory. 

**Step two:** Read the CodesForUser subsection in Contents section of the [TIRSPEC's github cover page](https://github.com/indiajoe/TIRSPEC).
This section summarises the purpose of each script in CodesForUser directory.

**Step three:** Put all the TIRSPEC data of different nights under some directory. You can name that directory whatever you want except same name as any of your night's data directory. For example let us call this directory _MotherDirectory/_. Before starting reduction you should have all the different night's data directories (Eg: _20140825/, 20140930/, 20141220/_) inside this _MotherDirectory/_.

**Step four (optional):** If you want to rename some wrongly named files. You can cd to that night's directory and run RenameImgs.sh . 

Note: Never change the _Slope-_ prefix as well and _filenumber_ in the filename when you rename. You are allowed to change on the central name part of the filename.

**Step five (optional):** If you want to discard some frames of a particular night open the SlopeimagesLog.txt file inside that night's directory and prefix # to the line of that particular image. Don't worry if you don't know which frames to discard. You can do that while running the data reduction pipeline also.

**Step six :** Copy the TIRSPECscript.conf file to some place. say inside _MotherDirectory/_ . You can rename that file to something which hint which source you are reducting. For example TIRSPECscrip_M31.conf.

**Step seven :** Open the .conf file and edit the parameters in it. Parameters in it are self explanatory. If not, see the TIRSPECscript.conf section in the page [here](Understanding-various-I-O-files-of-pipeline.md#tirspecscriptconf).

**Step eight (optional):** If you are doing spectroscopy. Extract the contents in data.tar.gz into the _MotherDirectory/_. They contain wavelength calibrated Argon lamps and Vega template.

**Step nine :** cd to _MotherDirectory/_ in terminal. And run the TIRSPECdataReduction.py script with your edited .conf file as 1st argument.


That is all. See the [Cookbook of Tips and Tricks](Cookbook-of-Tips-and-Tricks.md) page from some useful tips and tricks.

To obtain enlightenment and become a Zen Master of TIRSPEC data reduction pipeline. Learn about the Input and Output files of each stage in pipeline from [here](Understanding-various-I-O-files-of-pipeline.md).