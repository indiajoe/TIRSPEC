#!/usr/bin/env python
""" This script is to match (rotate, scale etc) and then subtract it from the target frame
This can be used for subtracting galaxy for SN template, or 
Subracting the Continuum images from Narrow band images 

IMP: Keep one ds9 open before running this script.
"""

from pyraf import iraf
iraf.imexamine.unlearn()
iraf.geotran.unlearn()
iraf.geomap.unlearn()
from astropy.io import fits
import scipy.optimize
import numpy as np
import os
import sys

def CreateMatchingXYCoords(TargetImg,Template,MatchingXYcoordsFile):
    """ Creates Xref Yref X Y, a 4 coulumn matched coordinates file MatchingXYcoordsFile."""
    iraf.display(TargetImg,1)
    print('Choose _atleast_ 9 matching stars')
    print('First Select maximum number of stars from this Target image by pressing _a_ on good stars. Press q to finish.')
    imxTarget = iraf.imexam(Stdout=1)
    iraf.display(Template,1)
    print('Press q to read from '+Template[:-5]+'_TemplXY.coo'+'. OR Select exactly same stars in same order from this Template image by pressing _a_ and then press q to finish.')
    imxTemplate = iraf.imexam(Stdout=1)
    if len(imxTemplate) == 0: # Default to the predefinted coo file
        with open(Template[:-5]+'_TemplXY.coo') as defaultcoo:
            imxTemplate = [line.rstrip() for line in defaultcoo]

    with open(MatchingXYcoordsFile,'w') as foo :    #Creating XYXY file
        i=2
        while i < len(imxTemplate) :               
            foo.write(imxTarget[i].split()[0] +' '+imxTarget[i].split()[1]+' '+imxTemplate[i].split()[0]+' '+imxTemplate[i].split()[1]+'\n')
            i=i+2
    print("Xref Yref X Y Table created in "+MatchingXYcoordsFile)


def AlignImage(TargetImg,Template,AlignedImageName):
    """Align Template Image to TargetImage to create output image AlignedImageName"""
    MatchingXYcoordsFile = AlignedImageName[:-5]+'.XYXY'
    if not os.path.isfile(MatchingXYcoordsFile) :
        CreateMatchingXYCoords(TargetImg,Template,MatchingXYcoordsFile)
    else:
        print('Using Old {0} file'.format(MatchingXYcoordsFile))

    XYtranformdb = AlignedImageName[:-5]+'rtran.db'
    iraf.geomap(input=MatchingXYcoordsFile, database=XYtranformdb, xmin=1, xmax=1024, ymin=1, ymax=1024, interactive=0,fitgeometry="general")
    #Now geotran the Template image to match the TargetImg
    iraf.geotran(input=Template,output=AlignedImageName,database=XYtranformdb,transforms=MatchingXYcoordsFile)
    #return the transfgorm .db filename just incase it is needed for the user
    return XYtranformdb
    
def dummyfunction(inp):
    """ Return the input as output """
    return inp

def ExtractTiles(FitsFile,CoordsFile,Summeryfunction=dummyfunction,hsize=5):
    """ Returns a list of tiles (of size 2*hsize x 2*hsize) from FitsFile 
    If Summeryfunction is given it will apply that function to the tiles and returns its output instead
    CoordsFile is X Y coordes in FITS/iraf coordinate system"""
    FullData=fits.getdata(FitsFile)
    output=[]
    with open(CoordsFile) as coofile:
        for cooline in coofile:
            cooline = cooline.rstrip()
            Y,X = int(float(cooline.split()[0]))-1 , int(float(cooline.split()[1]))-1
            output.append(Summeryfunction(FullData[X-hsize:X+hsize,Y-hsize:Y+hsize]))
    return output

def SkySubtractImage(Img,OutputFitsFile,SkyCoordsFile):
    """ Creates a median sky subtracted fits image """
    SkyValues = ExtractTiles(Img,SkyCoordsFile,Summeryfunction=np.median)
    MedianSky = np.median(SkyValues)
    iraf.imarith(operand1=Img,op="-",operand2=MedianSky,result=OutputFitsFile)
    return MedianSky

def MatchNSubtract(TargetImg,Template,OutputImage):
    """ Creates OutputImage =  TargetImg - Template after scaling and matching Template to TargetImg """
    
    AlignedImg = TargetImg[:-5]+"_"+Template
    TransformDBfile = AlignImage(TargetImg,Template,AlignedImg)
    
    # Now get the Good sky region coordinates
    SkyCoordsFile = TargetImg[:-5]+'_BlankSky.coo'
    if not os.path.isfile(SkyCoordsFile) :
        iraf.display(TargetImg,1)
        print ('For taking coordinates of good sky. Press _x_ over blank sky areas.')
        imx=iraf.imexam(Stdout=1)
        with open(SkyCoordsFile,'w') as foo :    #Creating blank sky coords files
            for line in imx :               
                foo.write(line.split()[0] +'  '+line.split()[1]+'\n')

    # Now get the regions in the image whose brightness has to be cancelled by scaling
    FluxCoordsFile = TargetImg[:-5]+'_FluxRegions.coo'
    if not os.path.isfile(FluxCoordsFile) :
        iraf.display(TargetImg,1)
        print ('Press _x_ over areas you want to minimise the residual flux after subtraction')
        imx=iraf.imexam(Stdout=1)
        with open(FluxCoordsFile,'w') as foo :    #Creating Flux ares which we have to remove in subtraction
            for line in imx :               
                foo.write(line.split()[0] +'  '+line.split()[1]+'\n')

    #Now we first has to remove background from both the images.
    TargetSkySubtractedFile = TargetImg[:-5]+'_SkyS.fits'
    if not os.path.isfile(TargetSkySubtractedFile):
        skyvalue = SkySubtractImage(TargetImg,TargetSkySubtractedFile,SkyCoordsFile)
    else:
        print('Warning: Using old {0} file'.format(TargetSkySubtractedFile))

    AlignedSkySubtractedFile = AlignedImg[:-5]+'_SkyS.fits'
    if not os.path.isfile(AlignedSkySubtractedFile):
        skyvalue = SkySubtractImage(AlignedImg,AlignedSkySubtractedFile,SkyCoordsFile)
    else:
        print('Warning: Using old {0} file'.format(AlignedSkySubtractedFile))

    #We shall now extract the totel Flux in each tiles from both the images
    TargetFluxinTiles = ExtractTiles(TargetSkySubtractedFile,FluxCoordsFile,Summeryfunction=np.sum,hsize=7*1.5)
    TemplateFluxinTiles = ExtractTiles(AlignedSkySubtractedFile,FluxCoordsFile,Summeryfunction=np.sum,hsize=7*1.5)
    
    def DiffSquareSum(x):
        return np.sum([(targetF - x*templateF)**2 for targetF,templateF in zip(TargetFluxinTiles,TemplateFluxinTiles)])
    
    res = scipy.optimize.minimize_scalar(DiffSquareSum)
    Scale = res.x
    print('Scaling to match the fluxes is {0}'.format(Scale))
    iraf.imarith(operand1=AlignedSkySubtractedFile,op="*",operand2=Scale,result=AlignedSkySubtractedFile[:-5]+'M.fits')

    iraf.imarith(operand1=TargetSkySubtractedFile,op="-",operand2=AlignedSkySubtractedFile[:-5]+'M.fits',result=OutputImage)
    

def main():
    if len(sys.argv)<2 :
        print('-'*10)
        print('Usage : {0} InputImageTemplate.txt'.format(sys.argv[0]))
        print('where,')
        print('     InputImageTemplate.txt is a three column file containing Target, Template and subtracted output filename')
        print(' ')
        print('-'*10)
        sys.exit(1)

    try : 
        InputFilelist=open(sys.argv[1],'r')
    except IOError :
        print ("Error: Cannot open the file "+sys.argv[1])
        sys.exit(1)
        
    for line in InputFilelist:
        line= line.rstrip()
        TargetImg = line.split()[0]
        Template = line.split()[1]
        OutputImage = line.split()[2]
        MatchNSubtract(TargetImg,Template,OutputImage)

if __name__ == "__main__":
    main()

