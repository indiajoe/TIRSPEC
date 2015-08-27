#!/usr/bin/env python
# This script is to fully automate a quick photometry reduction 
# as well as astrometry calibration of all tirspec data.
#-----------------------------------indiajoe@gmail.com
import sqlite3
import os
#import glob
#import sys, traceback 
import numpy as np
#import re
import shlex
import shutil
import subprocess
#import time
#import datetime
from astropy.io import ascii
import astropy.table as table 
from astropy.table import Column
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, ICRS
from astropy.convolution import Gaussian2DKernel

from astropy.vo.client import conesearch
from astropy.stats import gaussian_fwhm_to_sigma

import photutils as photu

from scipy import stats

# Required for certain functions... Remove it if not needed for your module
from pyraf import iraf
iraf.set(stdimage="imt2048") #Setting the image size to 2048. Not important if people use ds9 7.1

import ConfigureReduceAllImages as PC

import logging
logging.basicConfig(filename=PC.LOG_FILENAME,
                    level=logging.DEBUG,
                    )
logger = logging.getLogger(__name__)


def main():
    """ The  main funtion to run while run from terminal as a standalone pipeline """
    ##Load the list of directories to process in this run
    directories = LoadDirectories(PC)
    logger.info('Directories loaded : {0}'.format(directories))
    # Loop through each directory
    for night in directories:
        logger.info("Working on night : {0}".format(night))

        if PC.RAWIMAGEREDUCTION:
            logger.info("Doing basic Image Reduction..**")
            RawImageReduction(PC,night)

        if PC.DOPHOTOMETRY:
            logger.info("Doing Photometry..**")
            DoPhotometryofNight(PC,night)
            

def LoadDirectories(PC):
    """ Returns list of directories to process """
    with open(PC.DIRECTORYLISTFILE,'r') as dirsfile:
        directories = [dirpath.rstrip() for dirpath in dirsfile
                       if ((dirpath.strip() is not '') and (dirpath[0] !='#'))]

    for dirs in directories:
        #Create a corresponding night directory in OUTPUT directory also if not already present.
        try:
            os.makedirs(os.path.join(PC.OUTDIR,dirs))
        except OSError:
            if os.path.isdir(os.path.join(PC.OUTDIR,dirs)) :
                logger.warning("Output directory {0} already exists.\n"
                               "Everything inside it will be overwritten.".format(os.path.join(PC.OUTDIR,dirs)))
            else:
                raise
        else:
            logger.info('Created directory :'+os.path.join(PC.OUTDIR,dirs))
        
        #Also Alert user if some directory has SlopeimagesLog.txt missing in it.            
        if not os.path.isfile(os.path.join(dirs,PC.NIGHTLOGFILE)):
            logger.error('ERROR: {0} file missing.'.format(os.path.join(dirs,PC.NIGHTLOGFILE)))
            logger.error('Copy the {0} log file to the directory before procceding'.format(PC.NIGHTLOGFILE))

    return directories

def RawImageReduction(PC,night):
    """ Does the initial reduction and astrometry calibration of all images of the night """

    FiltrImgDic = dict() # Dic to save img list to create flats
    with open(os.path.join(night,PC.NIGHTLOGFILE),'r') as imglogFILE :
        # Skip blank lines and Commented out lines with #
        imglogFILElines = [imageLINE.rstrip() for imageLINE in imglogFILE 
                           if ((imageLINE.strip() is not '') and (imageLINE[0] !='#'))]
    PrevImg = None
    PrevFilt = None
    ListOfDitherSets = [[]]
    for imgline in imglogFILElines:
        # Skip if lower wheel was not open for broad band images
        if shlex.split(imgline)[6].upper() not in ['OPEN','UNKNOWN']: 
            continue
        # Skip if slit wheel was not open for broad band images
        if shlex.split(imgline)[7].upper() not in ['OPEN','UNKNOWN']: 
            continue
            # Skip if cal mirror was not out for broad band images
        if shlex.split(imgline)[8].upper() not in ['OUT','UNKNOWN']: 
            continue

        imgfile = os.path.join(night,shlex.split(imgline)[0])
        Filt = shlex.split(imgline)[5]
        FiltrImgDic.setdefault(Filt,[]).append(imgfile)

        ##find sets of same dither images
        if Filt != PrevFilt:
            SameSettings = False
            PrevFilt = Filt
        else:
            SameSettings = True
        
        # Check the new image and previous images are same
        SameImg = CheckImagesAreAligned(imgfile,PrevImg)
        PrevImg = imgfile
        if SameImg and SameSettings:
            ListOfDitherSets[-1].append(imgfile)
        else: # Start a new sub list for the new dither position
            ListOfDitherSets.append([imgfile])

    logger.info('Finished grouping images. Going forward to Master flat creation')
    ImgFiltrDic = dict((img,filt) for filt,imglist in FiltrImgDic.items() for img in imglist)
    ##Create Master flats for this night
    FiltFlatDic = {}
    for filt in FiltrImgDic:
        OutputFlatName = os.path.join(PC.OUTDIR,night,'Flat_{0}.fits'.format(filt))
        CombineImages(FiltrImgDic[filt],output=OutputFlatName,
                      method='median',scale='median',statsec=PC.STATSECTION,norm=True)
        FiltFlatDic[filt] = OutputFlatName
        logger.debug('Flat for night:{0}, Filter:{1}, Nflatfile:{2} by combining'
                     ' {3} imgs'.format(night,filt,OutputFlatName,len(FiltrImgDic[filt])))

    # Loop through each dither image set
    for imgsetlist in ListOfDitherSets:
        if len(imgsetlist) == 0:
            continue
        Filtr = ImgFiltrDic[imgsetlist[0]]
        ##Combine images of this dither set
        OutCombimg = os.path.join(PC.OUTDIR,night,os.path.basename(imgsetlist[0])[6:])   #Removing the Slope- prefix form 1st image to use as output image
        CombineImages(imgsetlist,output=OutCombimg, method='median',zero='median',
                      statsec=PC.STATSECTION)
        ## Do flat correction
        OutFinalImage = os.path.splitext(OutCombimg)[0]+'_FC.fits'
        iraf.imarith(operand1=OutCombimg,op='/',operand2=FiltFlatDic[Filtr],result=OutFinalImage)
        ## Apply Bad pixel mask correction
        FixBadPixels(PC,OutFinalImage,night)

        ## Do astrometry.net WCS calibration
        astrometrySuccess = RunAstrometryCalibration(OutFinalImage,Outputfile=OutFinalImage)
        if not astrometrySuccess:
            logger.error('ASTROMETRY_FAILED: {0}'.format(OutFinalImage))
            logger.warning('Astrometry calibration failed: {0}'.format(OutFinalImage))
            logger.info('Since astrometry failed Skipping image: {0}'.format(OutFinalImage))
            continue
        else:
            with open(os.path.join(PC.OUTDIR,night,PC.OUTFITSFILELIST),'a') as outfitsfilelist:
                outfitsfilelist.write('{0} {1}\n'.format(OutFinalImage,Filtr))


def CheckImagesAreAligned(NewImg,PrevImg,Pcoeffcut=0.3):
    """ Returns True if stars in images are aligned """
    if PrevImg is None:
        return False

    # First we have to remove uneven background from both images for looging at star coorilations
    PImgData = fits.getdata(PrevImg)
    PImgbkg = photu.Background(PImgData, (50, 50), filter_shape=(3, 3), method='median')
    PImgData -= PImgbkg.background

    NImgData = fits.getdata(NewImg)
    NImgbkg = photu.Background(NImgData, (50, 50), filter_shape=(3, 3), method='median')
    NImgData -= NImgbkg.background

    PImgthreshold =  3. * PImgbkg.background_rms

    kernel = Gaussian2DKernel(2.0 * gaussian_fwhm_to_sigma, x_size=3, y_size=3) # 2 pixel FWHM smoothing kernal
    segm = photu.detect_sources(PImgData, PImgthreshold, npixels=5, filter_kernel=kernel)
    mask = segm.astype(np.bool)
    # Remove any false positives due to bad pixels from the image edges
    mask[:50,:] = False
    mask[-50:,:] = False
    mask[:,:50] = False
    mask[:,-50:] = False
    # Also remove any negative counts area
    mask[PImgData<0] = False

    # Counts array in PrevImag at the detected sources
    PImgStarCounts = PImgData[mask]
    NImgStarCounts = NImgData[mask]

    pearsonCoeff = stats.pearsonr(PImgStarCounts,NImgStarCounts)[0]

    if pearsonCoeff > Pcoeffcut:
        Aligned = True
        logger.info('Pearson={0}: {1} is aligned with {2}.'.format(pearsonCoeff,NewImg,PrevImg))
    else:
        Aligned = False
        logger.info('Pearson={0}: {1} is not aligned with {2}.'.format(pearsonCoeff,NewImg,PrevImg))
    return Aligned

def DoPhotometryofNight(PC,night):
    """ Does photometry of all images in the night as well as create the final magnitude catlog table """
    with open(os.path.join(PC.OUTDIR,night,PC.OUTFITSFILELIST),'r') as outfitsfilelist:
        # Skip blank lines and Commented out lines with #
        imgfilterlist = [tuple(imageLINE.rstrip().split()) for imageLINE in outfitsfilelist 
                         if ((imageLINE.strip() is not '') and (imageLINE[0] !='#'))]

    for OutFinalImage,Filtr in imgfilterlist:
        MagTable = RunPhotometryofImage(OutFinalImage,Filtr)

        ## Append differential 2MASS magnitude of all sources to table
        TableHeaders = ['ra','dec','mag','magerror','Qflag']#filter,epoch,ImgFile
        filter_col = Column(name='Filter', data=[Filtr]*len(MagTable))
        epoch_col = Column(name='Epoch', data=[epoch]*len(MagTable))
        imgfile_col = Column(name='ImgFile', data=[OutFinalImage]*len(MagTable))
        TableToOutput = MagTable[TableHeaders]
        TableToOutput.add_columns([filter_col, epoch_col, imgfile_col])

        # Now append the Full output table also to an ascii file.
        OutputTableFilename = os.path.join(PC.OUTDIR,PC.OUTPUTFILE)
        try :
            PreviousFullTable = ascii.read(OutputTableFilename,delimiter=',',
                                           format='commented_header',fill_values=('--','0'))
        except IOError :
            logger.info('No previous photometry output found, hence we will be creating'
                        ' a new file {0}.'.format(OutputTableFilename))
            OutputTableToWrite = TableToOutput
        else :
            OutputTableToWrite = table.vstack([PreviousFullTable,TableToOutput], join_type='outer')
        # Writing the final appended table
        ascii.write(OutputTableToWrite, OutputTableFilename,delimiter=',',format='commented_header')

        logger.info("Photometry of {0} over.".format(OutFinalImage))
        #END of the photometry of this image..


def RunAstrometryCalibration(InputImage,Outputfile=None):
    """ Does Astrometry calibration of input fits image InputImage.
        Returns True if calibration succeded."""
    # Following lines are TIRSPEC specific for formating the time from fits header
    obs_timehdr = fits.getval(InputImage,'OBSTIME')
    obs_datehdr = fits.getval(InputImage,'OBSDATE')
    time_str = '{0} {1}'.format('-'.join(obs_datehdr.split('.')),obs_timehdr)

    time = Time(time_str,format='iso',scale='utc')
    hct_hanle = EarthLocation(lat=32.77944*u.deg, lon=78.9641*u.deg, height=4500*u.m)

    zenith = AltAz(location=hct_hanle, obstime=time, az=0*u.deg, alt=90*u.deg)
    ZenithRaDec = zenith.transform_to(ICRS)

    # Finally call astrometry.net software to calibrate fits image
    Z_ra = '{0.ra}'.format(ZenithRaDec).split()[0]
    Z_dec = '{0.dec}'.format(ZenithRaDec).split()[0]
    ret = subprocess.call(['solve-field','--no-plots', '--ra',Z_ra, '--dec',Z_dec, '--radius','85', 
                           '--scale-units','arcsecperpix', '--scale-low','0.28', '--scale-high','0.32',
                           InputImage])
    print(ret)
    if ret == 0:
        if Outputfile is not None:
            print('Copying {0} to {1}'.format(os.path.splitext(InputImage)[0]+'.new',Outputfile))
            shutil.copy(os.path.splitext(InputImage)[0]+'.new',Outputfile)
        return True
    else:
        return False

def RunPhotometryofImage(image,filt):
    """ Does photometry of input WCS calibrated image and return table of mganitudes """
    imgWCS = WCS(image)
    ## Find 2MASS bright star coordinates
    
    Catalog_2MASS = VOCatloger.RetrieveCatlog(imgwcs=imgWCS,catlog='2MASS')
    # Sort by magnitude
    Catalog_2MASS.sort('k_m')

    # Take N bright source X, Y coordinates
    N = 10
    NbrightXY = Catalog_2MASS[:N]['ra','dec']#['Xpix','Ypix']
    Gcoofilename = os.path.splitext(image)[0]+'.Gcoo'
    NbrightXY.write(Gcoofilename,format='ascii.no_header')

    ## Find median fwhn,e,PA etc from them
    fwhm,ratio,pa,sky_med,sky_std = ImageStarProperties(image,gstarcoo=NbrightXY)#Gcoofilename)

    ## Run daofind for full source list
    sources = daofind(imgdata-sky_med, fwhm=fwhm, threshold=5.*sky_std,ratio=ratio,theta=pa)

    # Match with 2MASS sources
    radec_daofind = imgWCS.wcs_pix2world(zip(sources['xcentroid'],sources['ycentroid']),1)
    Catdaofind = SkyCoord(radec_daofind[:,0],radec_daofind[:,1],frame='icrs')
    Cat2MASS = SkyCoord(Catalog_2MASS['ra'],Catalog_2MASS['dec'],frame='icrs') 
    index,dist2d,dist3d = Cat2MASS.match_to_catalog_sky(Catdaofind)

    ## Run phot for magnitude list
    
    ## Append to the Raw Instrument mag table
    
    ## Find delta mag
    

def ImageStarProperties(image,gstarcoo=None):
    """ Returns fwhm, ratio, pa, bkg, and bkg_std of input image """
    data = fits.getdata(image)
    imgWCS = WCS(image)

    bkg = photu.Background(data, (50, 50), filter_shape=(3, 3), method='median')
    threshold = bkg.background + (3. * bkg.background_rms)
    kernel = Gaussian2DKernel(2.0 * gaussian_fwhm_to_sigma, x_size=3, y_size=3) # 2 pixel FWHM smoothing kernal
    segm = photu.detect_sources(data, threshold, npixels=5, filter_kernel=kernel)
    props = segment_properties(data, segm,wcs=imgWCS)
    tbl = properties_table(props)
    catimg = SkyCoord(tbl['ra_icrs_centroid'],tbl['dec_icrs_centroid'],frame='icrs')
    catGstars = SkyCoord(gstarcoo['ra'],gstarcoo['dec'],frame='icrs') 
    # Cross match the catlogs to find same stars
    index,dist2d,dist3d = catGstars.match_to_catalog_sky(catimg)

    fwhm = np.sqrt(np.median(tbl[index]['covar_sigx2']))/gaussian_fwhm_to_sigma #WRONG FIND another way!
    ratio = np.median(tbl[index]['elongation'])
    pa = np.median(tbl[index]['orientation'])
    
    return (fwhm, ratio, pa, bkg.background, bkg.background_rms)


def CombineImages(imglist,output,method='median',zero='none' ,scale='none',norm=False,
                  reject="avsigclip", statsec='[150:900,150:900]'):
    """ Combined the input list of images and write to output fits file. """
    iraf.imcombine.unlearn()
    imglistfname = os.path.splitext(output)[0]+'.comblist'
    with open(imglistfname,'w') as imgs2combinefile:
        imgs2combinefile.write('\n'.join(imglist)+'\n')
    # Now call iraf imcombine with zero scaling
    combineoutputfile = os.path.splitext(output)[0]+'_un.fits' if norm else output

    iraf.imcombine(input='@'+imglistfname, output=combineoutputfile, combine=method, 
                   reject=reject, statsec=statsec, scale=scale, zero=zero)
    if norm:
        mediancounts = np.median(fits.getdata(combineoutputfile))
        iraf.imarith(operand1=combineoutputfile,op='/',operand2=mediancounts,result=output)



def ImgCombineWithZeroFloating(imglistfname,outputfile,cmethod="median",czero="median",
                               creject="avsigclip",cstatsection='[150:900,150:900]'):
    """ Returns the combined image with actuall average median flux, 
    It does zero scaling only for sigma rejection of stars. This is needed to remove faint 
    stars in rejection algorithm when the background sky itself is varying from frame to frame. """
    iraf.imcombine.unlearn()
    Xmin=float(cstatsection[1:-1].split(',')[0].split(':')[0])  #Everything now in fits coordinates
    Xmax=float(cstatsection[1:-1].split(',')[0].split(':')[1])
    Ymin=float(cstatsection[1:-1].split(',')[1].split(':')[0])
    Ymax=float(cstatsection[1:-1].split(',')[1].split(':')[1])

    if czero == "median" : statfunction = np.median
    elif czero == "average" : statfunction = np.mean
    else : 
        print('Error: czero should be median or average. Unknown option {0}'.format(czero))
        raise

    with open(imglistfname,'r') as imgfile:
        statlist=[]
        for img in imgfile:
            img = img.rstrip()
            statlist.append(statfunction(fits.getdata(img)[Ymin-1:Ymax,Xmin-1:Xmax]))
    print('{0} of images: {1}'.format(czero,str(statlist)))
    statAvg=np.mean(statlist)
    Zeroshifts= statAvg - np.array(statlist)
    print('Zeroshifts of images: {0} :: ImgAvg ={1}'.format(str(Zeroshifts),statAvg))
    with open(outputfile+'_zeroshifts.txt','w') as zeroshiftFILE:
        for shift in Zeroshifts: 
            zeroshiftFILE.write('{0} \n'.format(shift))
    # Now call iraf imcombine with zero scaling
    iraf.imcombine(input='@'+imglistfname, output=outputfile, combine=cmethod, 
                   reject=creject, statsec=cstatsection, zero='@'+outputfile+'_zeroshifts.txt')





def FixBadPixels(PC,images,nightdir):
    """ This will run iraf task proto.fixpix to interpolate badpixels """
    PixelMask = os.path.join(nightdir,PC.BadPixelMaskName)
    if not os.path.isfile(PixelMask):
        print("No Bad Pixel Mask file found by the name "+ PixelMask)
        print("Hence skipping Bad pixel interpolation")
        return

    iraf.proto(_doprint=0)
    iraf.fixpix.unlearn()
    iraf.fixpix(images=images,masks=PixelMask)


class VOManager(object):
    """ Class to manage and cache all calls which access VO catalogs """
    def __init__(self,search_rad=0.1*u.degree):
        self.search_rad = search_rad
        self.catnamedic = {'2MASS':'Two Micron All Sky Survey (2MASS) 1'}
        self.catalog_db = self.catnamedic['2MASS'] # By default we keep 2MASS
        self.RetrievedData = [] # List to cache the catlog retried each time
        self.MaxCacheLength = 100  # Number of maximum queries to cache
        self.imgX = 1024
        self.imgY = 1024

    def _retrievefullregion(self,center,FootPrints):
        """ Returns the full catlog region of size self.search_rad containing center and image FootPrints """
        for region in reversed(self.RetrievedData):
            if np.all(FootPrints.separation(region['center']) < self.search_rad): #All footprints inside region
                return region['table']
        else: # Not available in our cached data
            NewData = {'center':center}
            print('DOWNLOADING region around: {0}'.format(center))
            NewData['table'] = conesearch.conesearch(center,self.search_rad,catalog_db=self.catalog_db).to_table()
            self.RetrievedData.append(NewData)
            if len(self.RetrievedData) > self.MaxCacheLength: # Remove data from the top of list
                self.RetrievedData = self.RetrievedData[-self.MaxCacheLength:]
            return NewData
        
    def RetrieveCatlog(self,imgwcs,catlog='2MASS'):
        """ Returns a table of sources retreived from catalog for int input image WCS object """
        self.catalog_db = self.catnamedic[catlog]
        imgC = imgwcs.wcs_pix2world([(self.imgX/2,self.imgY/2)],1)
        center = SkyCoord(imgC[0][0], imgC[0][1], unit="deg")
        FootPrints = SkyCoord(imgwcs.calc_footprint(),unit="deg")
        FullTable = self._retrievefullregion(center,FootPrints)
        # Pixels coordinates of all stars
        XYpixels = imgwcs.wcs_world2pix(zip(FullTable['ra'],FullTable['dec']),1)
        Xpixel_col = Column(name='Xpix', data=XYpixels[:,0])
        Ypixel_col = Column(name='Ypix', data=XYpixels[:,1])
        # Add pixel column to table
        FullTable.add_columns([Xpixel_col,Ypixel_col])
        return FullTable[(FullTable['Xpix'] > 0) & (FullTable['Xpix'] < self.imgX) & 
                         (FullTable['Ypix'] > 0) & (FullTable['Ypix'] < self.imgY)]
    



VOCatloger = VOManager()

if __name__ == "__main__":
    main()
