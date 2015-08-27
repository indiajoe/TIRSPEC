""" This is the configuration file for ReduceAllImages.py """

# Text file containing list of directories to process
DIRECTORYLISTFILE = '/home/tirspec/JoeReduction/AllDataTest/directories'

# Output directory to write into
OUTDIR = '/home/tirspec/JoeReduction/AllDataTest/PipelineOutput'

# Output log file to write final fits image list of each night
OUTFITSFILELIST = 'ReducedImagesLog.txt'

# File to write photometry output
OUTPUTFILE = 'PhotometryOutput.txt'

# Name of the log file to save logs
LOG_FILENAME = 'Log_ReduceAllImagesScript.log'

# Wheter to do Raw Image Redcution
RAWIMAGEREDUCTION = True

# Wheter to do Photometry of reduced Images
DOPHOTOMETRY = False

# Observation Log filename
NIGHTLOGFILE = 'SlopeimagesLog.txt'

# Bad pixle mask name
BadPixelMaskName = 'SlopeTIRSPECphoto-BP.pl'

# Good section of array to use for stat calculation for image
STATSECTION = '[150:900,150:900]'

# Initialisation of the Pipeline 
# The .rc part!!

import os
#Create the OUTPUT directory if not already present.
try:
    os.makedirs(OUTDIR)
except OSError:
    if os.path.isdir(OUTDIR) :
        print("Output directory "+os.path.join(OUTDIR)+" already exists.\n Everything inside it will be overwritten.")
    else:
        raise
else:
    print('Created directory :'+OUTDIR)

