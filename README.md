TIRSPEC
=================

Contains the scripts which has to be used for the basic reduction of TIRSPEC data.
And also some add-on scripts to assist during observation.
Codes are mostly in Python.

Feel free to modify as needed and use it while observing with TIRSPEC.
Make sure, you keep the modified scripts under a seperate directory of your name, 
to prevent any conflict with default scripts for other users.

You may also find some helpful scripts in the following link: http://indiajoe.github.io/HandyTools4Astronomers/

Webpage of this repository: http://indiajoe.github.io/TIRSPEC/

Contents
================
### Major Scripts: ###
*    **Generate4mRamp.py**  *:* This script generates the slope calculation from up-the-ramp readout raw data. It produces the dark subtracted final slope images. A log file is also generated. You can load this also as a module in IPython, and use the collection of functions in it for an interative analysis of tirspec's raw data.
*    **TIRSPECdataReduction.py** (needs **TIRSPECscript.conf**) *:* This script is to help the astronomer reduce tirspec data (output of Generate4mRamp.py) . It will semi-automate and guide the user through : image selection, flat correction, aligning and combining of NIR dithered frames for final photometry/spectroscopy.
*    **AlignCombineImagesinDir.py** *:* This is a stand alone script to help astronomer align and combine his fits images.

### Minor Scripts: ### 
Intented for use on tirspec machine while observing.
*    **StarPlot.py** (needs **StarPlot.gnuplot**) *:* This script is to quickly view the surface plot and contour of a star profile in fits file. 
*    **DitherAssistant.py** *:* To help in visualising Dither pattern and generating the commands to be given to Observatory Server.
*    **TelescopeShift-4OS.py** *:* To get the input to be given to Observatory Server to move star to slit or any other location in image.
*    **TelescopeShift.py** *:* To get the input to be given to Keystone Server to move star to slit or any other location in image.
*    **lastNDR.sh** *:* To quickly load the last NDR readout frame of the directory to DV.
*    **Start_TIRSPEC.sh** (needs **run_mkmac_h1_as.sh,run_dv_LampAlert.sh**) *:* This script starts TIRSPEC software for starting observation in the night.
*    **.dv-init** *:* Some DV start up script to make life easier while observing.
      
Module Dependencies
-------------------
First few lines of the script will tell you the required python modules for each script.
Generaly speaking you will need: numpy, pyfits and pyraf
      

License
---------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.